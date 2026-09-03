#!/usr/bin/env python3
import argparse
import json
import os
import subprocess
from pathlib import Path

import pandas as pd
import polars as pl

from drugmr import paths

# look at MR results based on dataset X (which == arg)
# IVW p_FDR < 0.05 and passes egger intercept and cochran Q
# go to that gene pQTL df["protein"] - and acess corresponding directory in dat/cis_regions/{dataset}
# grab parquet files and cmd run for ./bin/coloc.R script

# Running this script comes after X -> Y and X -> M and M -> Y runs
# Parse mediators arg (if != mediators:)
# Run the function with its current architecture
# if mediators:
# Look at mediators dir for X -> M res as long as M -> Y
# if ivw_p (M -> Y) < 0.05 & for each protein in X -> M if ivw p_FDR < 0.05 and Cochran Q and Egger intercept
# if the same protein also passes this shit in X -> Y (Then run coloc pairwise - for protein X - coloc X -> Y and X -> M)
# save in work/coloc/X_M_Y_coloc.json -> protein: traits involved (i.e. the trait prefixes in terms of the traits which passed everything above for protein X) - {}
# In python function -> read that json - check how many (No of traits) - and then parse that as an argument on muli_trait_coloc() for moloc.R
# if pp.h4 > 0.7 -> execute moloc (multi_trait_coloc() function)
# runs moloc script - for each protein in work/coloc/.json (with the following args)
# - n_o of traits for protein X
# - pheno_id
# - pqtl_dataset
# - pqtl_dir (maybe)

def pairwise_coloc(
    pqtl_dataset: str,
    local_results_dir: str,
    pqtl_dir: str,
    pheno_id: str,
    n_cases: int,
    n_controls: int,
    wald_fdr_q: float = 0.05,
    ivw_fdr_q: float = 0.05,
    cochran_q_pval: float = 0.05,
    egger_intercept_pval_min: float = 0,
    min_instruments_for_ivw: int = 3,
):
    # local_results_dir is the general Tier-2 results root (e.g. runs/<run_id>/results),
    # not just the cis-MR dir - every paths.py call below derives its own subdir from it
    pqtl_dataset = pqtl_dataset.lower()
    pqtl_dir = Path(pqtl_dir)
    coloc_script = "./bin/coloc.R"
    out_dir = paths.coloc_out(pqtl_dataset, pheno_id, local_results_dir).parent
    out_dir.mkdir(parents=True, exist_ok=True)
    df = pl.read_csv(paths.mr_out(pqtl_dataset, pheno_id, local_results_dir), separator="\t")
    results = []
    results_sensitivity = []

    # filter for proteins which passed cis-MR thresholds (params/*.yaml gates.cis_mr)
    df2 = (
        df
        .filter(
            (
                (pl.col("n_instruments") >= min_instruments_for_ivw) &
                (pl.col("IVW_FDR_q") < ivw_fdr_q) &
                (pl.col("egger_intercept_pval") > egger_intercept_pval_min) &
                (pl.col("Q_pval") > cochran_q_pval)
            )
            |
            (
                (pl.col("n_instruments") == 1) &
                (pl.col("Wald_FDR_q") < wald_fdr_q)
            )
        )
        .select("protein")
        .unique()
    )

    print(f"[TRACKING] Proteins passing cis-MR filters: {df2.height}")

    for row in df2.iter_rows(named=True):
        protein = row["protein"]
        protein_dir = pqtl_dir / protein
        gwas = protein_dir / "gwas.parquet"
        pqtl = protein_dir / "pqtl.parquet"
        protein_file = out_dir / f"{pheno_id}_{protein}_coloc.tsv"
        sensitivity_file = out_dir / f"{pheno_id}_{protein}_coloc_sensitivity.tsv"
        cmd_coloc = ["Rscript", coloc_script, pqtl_dataset, protein, pheno_id, str(gwas), str(pqtl), str(n_cases), str(n_controls), str(local_results_dir)]
        print(f"[TRACKING] Running COLOC for {protein}")
        subprocess.run(cmd_coloc, check=True)
        results.append(pd.read_csv(protein_file, sep="\t"))
        protein_file.unlink()
        results_sensitivity.append(pd.read_csv(sensitivity_file, sep="\t"))
        sensitivity_file.unlink()

        # compile into 1 master file
    master = pd.concat(results, ignore_index=True)
    out_file = out_dir / f"{pqtl_dataset}_{pheno_id}_all_coloc.tsv"
    master.to_csv(out_file, sep="\t", index=False)
    print(f"[DONE] Saved master COLOC table: {out_file}")

    # same treatment for the per-protein prior sensitivity sidecar bin/coloc.R
    # writes alongside each protein_file - collect into 1 master table so it
    # rides on the same check_output(coloc_out, ...) skip gate as everything
    # else in this step, rather than leaving 1 orphaned TSV per protein behind
    master_sensitivity = pd.concat(results_sensitivity, ignore_index=True)
    sensitivity_out_file = paths.coloc_sensitivity_out(pqtl_dataset, pheno_id, local_results_dir)
    master_sensitivity.to_csv(sensitivity_out_file, sep="\t", index=False)
    print(f"[DONE] Saved master COLOC sensitivity table: {sensitivity_out_file}")




# include -> HyPrColoc here and also include filter
# -> * passes cochran Q and IVW



def coloc_with_mediators(
    pqtl_dataset: str,
    local_results_dir: str,
    pqtl_dir: str,
    pheno_id: str,
    n_cases: int,
    n_controls: int,
    mediator_manifest: str,
    ivw_fdr_q: float = 0.05,
    pp4_threshold: float = 0.7,
):
    standard_coloc = "./bin/coloc.R"
    moloc = "./bin/moloc.R"
    Ms = Path(mediator_manifest)

    # pqtl + gwas dir
    cis_regions = Path(f"./dat/cis_regions/{pqtl_dataset}")

    # out_dir
    out_dir = paths.coloc_out(pqtl_dataset, pheno_id, local_results_dir).parent
    os.makedirs(out_dir, exist_ok=True)

    moloc_json_dir = Path("./work/coloc")
    moloc_json_dir.mkdir(parents=True, exist_ok=True)
    moloc_json = {}

    Ms = pl.read_csv(mediator_manifest)
    mediators = Ms["pheno_id"].to_list()
    res = paths.network_mr_mediation_estimates_out(pqtl_dataset, pheno_id, local_results_dir)

    # ukb_ppp_AD_networkMR.tsv
    df = pl.read_csv(res, separator="\t")
    
    # results list (for pairwise coloc runs)
    results_pairwise = []
    results_mediator_pairwise = []
    results_pairwise_sensitivity = []
    results_mediator_pairwise_sensitivity = []

    # compile candidate proteins for coloc
    proteins = df["protein"].to_list()

    # PASSING PARAMETERS FOR CANDIDATE PROTEINS
    # X_M_IVW_FDR_q
    # X_Y_IVW_FDR_q
    # M_Y_IVW_pval
    for row in df.iter_rows(named=True):
        # MAYBE AT A LATER STAGE -> SAVE Cochran_Q pval and Egger Intercept pval onto mediator output CSV
        X_M_IVW_FDR_q = row["X_M_IVW_FDR_q"]
        X_Y_IVW_FDR_q = row["X_Y_IVW_FDR_q"]
        # X_Y_cochran_q = row["X_Y_IVW_FDR_q"]
        M_Y_IVW_pval = row["M_Y_IVW_pval"]

        if X_Y_IVW_FDR_q < ivw_fdr_q:
            outcome = row["pheno_id"]
            protein = row["protein"]
            if str(outcome) == pheno_id:
                print("[TRACKING] Pheno IDs match!")
            else:
                print("[CONCERN] Yowza! Something is wrong cuz outcome pheno IDs do not match!")

            protein_dir = cis_regions / protein

            # re-doing the pairwise coloc logic from the function above
            gwas = protein_dir / "gwas.parquet"
            pqtl = protein_dir / "pqtl.parquet"

            protein_file = out_dir / f"{pheno_id}_{protein}_coloc.tsv"
            protein_sensitivity_file = out_dir / f"{pheno_id}_{protein}_coloc_sensitivity.tsv"
            cmd_coloc = ["Rscript", standard_coloc, pqtl_dataset, protein, pheno_id, str(gwas), str(pqtl), str(n_cases), str(n_controls), str(local_results_dir)]
            print(f"[TRACKING] Running COLOC for {protein}")
            subprocess.run(cmd_coloc, check=True)
            pairwise_df = pl.read_csv(protein_file, separator="\t")
            results_pairwise.append(pairwise_df)
            protein_file.unlink()
            results_pairwise_sensitivity.append(pl.read_csv(protein_sensitivity_file, separator="\t"))
            protein_sensitivity_file.unlink()

            # second if (not strictly necessary -> as we can carry on with only X -> Y coloc)
            if X_M_IVW_FDR_q < ivw_fdr_q:
                mediator = row["mediator"]
                if mediator in mediators:
                    print(f"[TRACKING] Mediator {mediator} tracked!")
                else:
                    print(f"[CONCERN] Mediator {mediator} not found...")

                # cont...
                m = protein_dir / "mediators" / f"{mediator}.parquet"
                mediator_file = out_dir / f"{mediator}_{protein}_coloc.tsv"
                mediator_sensitivity_file = out_dir / f"{mediator}_{protein}_coloc_sensitivity.tsv"

                # n_cases and n_controls == n/a because mediator == quant trait
                cmd_coloc = [ "Rscript", standard_coloc, pqtl_dataset, protein, mediator, str(m), str(pqtl), str(n_cases), str(n_controls), str(local_results_dir)]
                print(f"[TRACKING] Running COLOC for {protein}")
                subprocess.run(cmd_coloc, check=True)
                mediator_df = pl.read_csv(mediator_file, separator="\t")
                results_mediator_pairwise.append(mediator_df)
                mediator_file.unlink()
                results_mediator_pairwise_sensitivity.append(pl.read_csv(mediator_sensitivity_file, separator="\t"))
                mediator_sensitivity_file.unlink()

                # open those two coloc results and make sure pp4 for the same protein == > 0.7
                # if pp4 > 0.7 on both: (PP.H4.abf)

                y_row = pairwise_df.filter(pl.col("protein_id") == protein)
                m_row = mediator_df.filter(pl.col("protein_id") == protein)

                if y_row.height > 0 and m_row.height > 0:
                    y_pp4 = y_row["PP.H4.abf"][0]
                    m_pp4 = m_row["PP.H4.abf"][0]
                    if y_pp4 > pp4_threshold and m_pp4 > pp4_threshold:
                        print(f"[TRACKING] {protein} passed pairwise COLOC for {outcome} and {mediator}")

                        if protein not in moloc_json:
                            moloc_json[protein] = [pheno_id]

                        if mediator not in moloc_json[protein]:
                            moloc_json[protein].append(mediator)

                # we need to also check whether > 1 mediator colocalises oin that same protein and save the correspondign according json file for moloc

    # compile into master files, same treatment as pairwise_coloc(): 1 row per
    # protein/mediator pair rather than 1 orphaned TSV per pair left in out_dir.
    # coloc_out() is what local.py's/hpc.py's check_output(coloc_out, "COLOC", ...)
    # gates the whole step on, so writing it here too (this branch previously
    # never wrote it at all) means a mediators run now also satisfies that gate.
    if results_pairwise:
        pl.concat(results_pairwise).write_csv(
            str(paths.coloc_out(pqtl_dataset, pheno_id, local_results_dir)), separator="\t"
        )
        print(f"[DONE] Saved master COLOC table: {paths.coloc_out(pqtl_dataset, pheno_id, local_results_dir)}")

        pl.concat(results_pairwise_sensitivity).write_csv(
            str(paths.coloc_sensitivity_out(pqtl_dataset, pheno_id, local_results_dir)), separator="\t"
        )
        print(f"[DONE] Saved master COLOC sensitivity table: {paths.coloc_sensitivity_out(pqtl_dataset, pheno_id, local_results_dir)}")
    else:
        print("[CONCERN] No proteins passed X_Y_IVW_FDR_q - no COLOC master table to save.")

    if results_mediator_pairwise:
        pl.concat(results_mediator_pairwise).write_csv(
            str(paths.coloc_mediator_out(pqtl_dataset, pheno_id, local_results_dir)), separator="\t"
        )
        print(f"[DONE] Saved master mediator COLOC table: {paths.coloc_mediator_out(pqtl_dataset, pheno_id, local_results_dir)}")

        pl.concat(results_mediator_pairwise_sensitivity).write_csv(
            str(paths.coloc_mediator_sensitivity_out(pqtl_dataset, pheno_id, local_results_dir)), separator="\t"
        )
        print(f"[DONE] Saved master mediator COLOC sensitivity table: {paths.coloc_mediator_sensitivity_out(pqtl_dataset, pheno_id, local_results_dir)}")
    else:
        print("[CONCERN] No proteins passed X_M_IVW_FDR_q - no mediator COLOC master table to save.")

    moloc_json = {
        protein: traits
        for protein, traits in moloc_json.items()
        if len(traits) >= 2
    }

    moloc_json_file = moloc_json_dir / f"{pqtl_dataset}_{pheno_id}_moloc.json"

    with open(moloc_json_file, "w") as f:
        json.dump(moloc_json, f, indent=4)

    print(f"[TRACKING] Saved MOLOC JSON: {moloc_json_file}")

    # Run moloc.R
    print(f"[TRACKING] Running MOLOC for all proteins within {moloc_json_file}!")

    cmd_moloc = f"""
set -euo pipefail
Rscript {moloc} \
    {pheno_id} \
    {pqtl_dataset} \
    {moloc_json_file}
"""
    
    subprocess.run(cmd_moloc, shell=True, check=True, executable="/bin/bash")


def main():
    # if mediators: true - do NOT run this one - if mediators = true run coloc_with_mediators():
    p = argparse.ArgumentParser()
    p.add_argument("--pqtl_dataset", required=True, choices=["ukb_ppp", "decode", "wu_csf", "wingo_brain"])
    p.add_argument("--local_results_dir", required=True)
    p.add_argument("--pqtl_dir", required=True)
    p.add_argument("--pheno_id", required=True)
    p.add_argument("--n_cases", required=True, type=int)
    p.add_argument("--n_controls", required=True, type=int)
    p.add_argument("--mediators", action="store_true")
    p.add_argument("--mediator_manifest", required=False)
    p.add_argument("--wald_fdr_q", type=float, default=0.05)
    p.add_argument("--ivw_fdr_q", type=float, default=0.05)
    p.add_argument("--cochran_q_pval", type=float, default=0.05)
    p.add_argument("--egger_intercept_pval_min", type=float, default=0)
    p.add_argument("--min_instruments_for_ivw", type=int, default=3)
    p.add_argument("--pp4_threshold", type=float, default=0.7)
    args = p.parse_args()

    # if mediators: true:
    if args.mediators:
        if args.mediator_manifest is None:
            raise ValueError("--mediator_manifest is required when --mediators is used")

        coloc_with_mediators(
            pqtl_dataset=args.pqtl_dataset,
            local_results_dir=args.local_results_dir,
            pqtl_dir=args.pqtl_dir,
            pheno_id=args.pheno_id,
            n_cases=args.n_cases,
            n_controls=args.n_controls,
            mediator_manifest=args.mediator_manifest,
            ivw_fdr_q=args.ivw_fdr_q,
            pp4_threshold=args.pp4_threshold,
        )

    # else:
    else:
        pairwise_coloc(
            pqtl_dataset=args.pqtl_dataset,
            local_results_dir=args.local_results_dir,
            pqtl_dir=args.pqtl_dir,
            pheno_id=args.pheno_id,
            n_cases=args.n_cases,
            n_controls=args.n_controls,
            wald_fdr_q=args.wald_fdr_q,
            ivw_fdr_q=args.ivw_fdr_q,
            cochran_q_pval=args.cochran_q_pval,
            egger_intercept_pval_min=args.egger_intercept_pval_min,
            min_instruments_for_ivw=args.min_instruments_for_ivw,
        )

if __name__ == "__main__":
    main()