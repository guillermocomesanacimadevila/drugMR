#!/usr/bin/env python3
import polars as pl
import pandas as pd
import subprocess
import argparse
from pathlib import Path
import os 
import json

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

def pairwise_coloc(pqtl_dataset: str, local_results_dir: str, pqtl_dir: str, pheno_id: str, n_cases: int, n_controls: int):
    pqtl_dataset = pqtl_dataset.lower()
    pqtl_dir = Path(pqtl_dir)
    local_results_dir = Path(local_results_dir)
    coloc_script = "./bin/coloc.R"
    out_dir = Path("./results/coloc") / pqtl_dataset
    out_dir.mkdir(parents=True, exist_ok=True)
    df = pl.read_csv(f"{local_results_dir}/{pqtl_dataset}_{pheno_id}_all_MR.tsv", separator="\t")
    results = []

    # filter for proteins which passed cis-MR thresholds
    df2 = (
        df
        .filter(
            (pl.col("IVW_FDR_q") < 1) & # 0.05 -> 1 for CI/CD testing ########## CHANGE THIS AFTER CI/CD TESTING
            (pl.col("egger_intercept_pval") > 0.05) & ########## CHANGE THIS AFTER CI/CD TESTING
            (pl.col("Q_pval") > 0.05) ########## CHANGE THIS AFTER CI/CD TESTING
        )
        .select("protein")
    )

    print(f"[TRACKING] Proteins passing cis-MR filters: {df2.height}")

    for row in df2.iter_rows(named=True):
        protein = row["protein"]
        protein_dir = pqtl_dir / protein
        gwas = protein_dir / "gwas.parquet"
        pqtl = protein_dir / "pqtl.parquet"
        protein_file = out_dir / f"{pheno_id}_{protein}_coloc.tsv"
        cmd_coloc = ["Rscript", coloc_script, pqtl_dataset, protein, pheno_id, str(gwas), str(pqtl), str(n_cases), str(n_controls)]
        print(f"[TRACKING] Running COLOC for {protein}")
        subprocess.run(cmd_coloc, check=True)
        results.append(pd.read_csv(protein_file, sep="\t"))
        protein_file.unlink()

        # compile into 1 master file 
    master = pd.concat(results, ignore_index=True)
    out_file = out_dir / f"{pqtl_dataset}_{pheno_id}_all_coloc.tsv"
    master.to_csv(out_file, sep="\t", index=False)
    print(f"[DONE] Saved master COLOC table: {out_file}")


def coloc_with_mediators(pqtl_dataset: str, local_results_dir: str, pqtl_dir: str, pheno_id: str, n_cases: int, n_controls: int, mediator_manifest: str):
    standard_coloc = "./bin/coloc.R"
    moloc = "./bin/moloc.R"
    Ms = Path(mediator_manifest)

    # pqtl + gwas dir
    cis_regions = Path(f"./dat/cis_regions/{pqtl_dataset}")

    # out_dir
    out_dir = Path("./results/coloc") / pqtl_dataset
    os.makedirs(out_dir, exist_ok=True)

    moloc_json_dir = Path("./work/coloc")
    moloc_json_dir.mkdir(parents=True, exist_ok=True)
    moloc_json = {}

    Ms = pl.read_csv(mediator_manifest)
    mediators = Ms["pheno_id"].to_list()
    network_mr_outputs = Path(f"./results/networkMR/mediation_estimates/{pqtl_dataset}")
    res = network_mr_outputs / f"{pqtl_dataset}_{pheno_id}_networkMR.tsv"

    # ukb_ppp_AD_networkMR.tsv
    df = pl.read_csv(res, separator="\t")
    
    # results list (for pairwise coloc runs)
    results_pairwise = []
    results_mediator_pairwise = []

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
        M_Y_IVW_pval = row["M_Y_IVW_pval"]

        if X_Y_IVW_FDR_q < 1: ####### CHANGE THESE PARAMS BACK TO 0.05 AFETR CI/CD TESTING
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
            cmd_coloc = ["Rscript", standard_coloc, pqtl_dataset, protein, pheno_id, str(gwas), str(pqtl), str(n_cases), str(n_controls)]
            print(f"[TRACKING] Running COLOC for {protein}")
            subprocess.run(cmd_coloc, check=True)
            pairwise_df = pl.read_csv(protein_file, separator="\t")
            results_pairwise.append(pairwise_df)

            # second if (not strictly necessary -> as we can carry on with only X -> Y coloc)
            if X_M_IVW_FDR_q < 1: ####### CHANGE THESE PARAMS BACK TO 0.05 AFETR CI/CD TESTING
                mediator = row["mediator"]
                if mediator in mediators:
                    print(f"[TRACKING] Mediator {mediator} tracked!")
                else:
                    print(f"[CONCERN] Mediator {mediator} not found...")

                # cont...
                m = protein_dir / "mediators" / f"{mediator}.parquet"
                mediator_file = out_dir / f"{mediator}_{protein}_coloc.tsv"

                # n_cases and n_controls == n/a because mediator == quant trait
                cmd_coloc = [ "Rscript", standard_coloc, pqtl_dataset, protein, mediator, str(m), str(pqtl), str(n_cases), str(n_controls)]
                print(f"[TRACKING] Running COLOC for {protein}")
                subprocess.run(cmd_coloc, check=True)
                mediator_df = pl.read_csv(mediator_file, separator="\t")
                results_mediator_pairwise.append(mediator_df)

                # open those two coloc results and make sure pp4 for the same protein == > 0.7
                # if pp4 > 0.7 on both: (PP.H4.abf)

                y_row = pairwise_df.filter(pl.col("protein_id") == protein)
                m_row = mediator_df.filter(pl.col("protein_id") == protein)

                if y_row.height > 0 and m_row.height > 0:
                    y_pp4 = y_row["PP.H4.abf"][0]
                    m_pp4 = m_row["PP.H4.abf"][0]
                    if y_pp4 > 0.01 and m_pp4 > 0.01: ####### CHANGE THESE PARAMS BACK TO 0.7 AFETR CI/CD TESTING
                        print(f"[TRACKING] {protein} passed pairwise COLOC for {outcome} and {mediator}")

                        if protein not in moloc_json:
                            moloc_json[protein] = [pheno_id]

                        if mediator not in moloc_json[protein]:
                            moloc_json[protein].append(mediator)

                # we need to also check whether > 1 mediator colocalises oin that same protein and save the correspondign according json file for moloc

    moloc_json = {
        protein: traits
        for protein, traits in moloc_json.items()
        if len(traits) > 2
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
    p.add_argument("--pqtl_dataset", required=True, choices=["ukb_ppp", "decode"])
    p.add_argument("--local_results_dir", required=True)
    p.add_argument("--pqtl_dir", required=True)
    p.add_argument("--pheno_id", required=True)
    p.add_argument("--n_cases", required=True, type=int)
    p.add_argument("--n_controls", required=True, type=int)
    p.add_argument("--mediators", action="store_true")
    p.add_argument("--mediator_manifest", required=False)
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
        )

if __name__ == "__main__":
    main()