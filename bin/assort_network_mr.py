#!/usr/bin/env python3
import argparse
import os
import subprocess
from pathlib import Path

import polars as pl

from drugmr import NetworkMR, paths

# what do we need to run it for CI/CD testing
# AD GWAS (pheno_id)
# ukb_ppp pQTLs (dir) -> after SLURM ETL pipeline
# cp -R SCZ GWAS (with different header names onto a mediators directory)
# adjust params to be permisive for CI/CD testing rather than standard significance thresholds

# DS NetworkMR pipeline
# dictionary in jupyter notebook {M_id: 'User/Path/...'}
# FROM NOTEBOOK -> MAKE A MEDIATOR MANIFEST
# For each protein part of dataset X
# Run cis-MR (twice) -> For each X -> M
# Also run M -> Y (whole genome) 
# dat/derived/<pheno_id>/network_mr/M_Y/....csv (Genome-wide, shared across every pqtl_dataset -
#   one CSV with MR outputs where 1 entry == univariable MR from a mediator M on Y)
# <local_results_dir>/network_mr/X_M/mediator_1/....csv (1 entry == univariable cis-MR - 1 protein vs that mediator)
# <local_results_dir>/network_mr/X_M/mediator_2/....csv (1 entry == univariable cis-MR - 1 protein vs that mediator)
# <local_results_dir>/network_mr/X_M/mediator_N/....csv (1 entry == univariable cis-MR - 1 protein vs that mediator)
# <local_results_dir>/network_mr/mediation_estimates/...csv (massive CSV with a given protein that FDR significant in X->M and X->Y and also if IVW_p < 0.05 in X->Y run NetworkMR package - here the output of NetworkMR package)

# ARGS
# pheno_id   
# pheno_gwas 
# ref_bfile  
# out_dir    
# mediator_dict (or another data type which allows for >1 value) THIS IT FOR MEDIATOR SAMPLE SIZE!!!!!
# pqtl_dataset
# pqtl_dir

# X -> Y == already done
# For each M -> run X -> M
# For each M -> run M -> Y (genome-wide)
# Inherit NMR library from drugmr/ and run NetworkMR

# for each mediatror 

def run_genomewide_mr(ref_bfile: str, pheno_id: str, pheno_gwas: str):
    genomewide_mr = "./bin/genomewide_mr.R"
    mediator_dir = paths.qc_mediator_dir()
    out_dir = paths.network_mr_m_y_dir(pheno_id)
    out_dir.mkdir(parents=True, exist_ok=True)
    all_results = []
    for file in mediator_dir.glob("*.tsv"):
        M_id = file.stem
        M = pl.read_csv(file, separator="\t", n_rows=1)
        M_N = M["N"][0]

        cmd = f"""
set -euo pipefail
Rscript {genomewide_mr} \\
    '{M_id}' \\
    '{file}' \\
    '{M_N}' \\
    '{pheno_id}' \\
    '{pheno_gwas}' \\
    '{ref_bfile}' \\
    '{out_dir}'
"""

        subprocess.run(cmd, shell=True, check=True, executable="/bin/bash")
        res_file = out_dir / f"{M_id}_{pheno_id}_genomewide_MR.tsv"

        if res_file.exists():
            result = pl.read_csv(res_file, separator="\t", null_values=["NA", "NaN", "nan", ""])
            numeric_cols = ["mediator_N", "n_instruments", "IVW_beta", "IVW_se", "IVW_pval", "Egger_beta", "Egger_se", "Egger_pval", "egger_intercept", "egger_intercept_pval", "WME_beta", "WME_se", "WME_pval", "Wald_beta", "Wald_se", "Wald_pval", "Q", "Q_df", "Q_pval"]
            result = result.with_columns([pl.col(col).cast(pl.Float64, strict=False) for col in numeric_cols if col in result.columns])
            all_results.append(result)

    if all_results:
        compiled = pl.concat(all_results, how="diagonal_relaxed")
        compiled_file = out_dir / f"{pheno_id}_mediator_genomewide_MR.tsv"
        compiled.write_csv(compiled_file, separator="\t")
        print(f"[TRACKING] Saved compiled genome-wide mediator MR results: {compiled_file}")
        
        for file in out_dir.glob(f"*_{pheno_id}_genomewide_MR.tsv"):
            if file.name != compiled_file.name:
                file.unlink()
    else:
        print("[CONCERN] No genome-wide mediator MR results generated.")


def run_cis_mr_X_M(pqtl_dataset: str, pqtl_dir: str, ref_bfile: str, local_results_dir: str = "results"):
    ref_bfile = Path(ref_bfile)
    pqtl_dir = Path(pqtl_dir)
    out_dir = paths.network_mr_x_m_dir(pqtl_dataset, local_results_dir)
    protein_dir = Path(f"./dat/cis_regions/{pqtl_dataset}")
    mediator_gwas = paths.qc_mediator_dir()
    out_dir.mkdir(parents=True, exist_ok=True)

    for protein_path in protein_dir.iterdir():
        protein = protein_path.name
        print(f"[TRACKING] Processing protein {protein}")
        pqtl_file = protein_path / "pqtl.parquet"
        gwas_file = protein_path / "gwas.parquet"
        pqtl = pl.read_parquet(pqtl_file)
        gwas = pl.read_parquet(gwas_file)
        chr = pqtl["CHR"][0]
        start = min(pqtl["BP"].to_list())
        end = max(pqtl["BP"].to_list())
        mediator_out = protein_path / "mediators"
        mediator_out.mkdir(parents=True, exist_ok=True)

        for mediator_file in mediator_gwas.iterdir():
            if mediator_file.suffix not in [".tsv", ".txt"]:
                continue

            mediator_id = mediator_file.stem
            print(f"[TRACKING] Processing mediator {mediator_id} for {protein}")
            mediator = pl.read_csv(mediator_file, separator="\t")
            mediator = mediator.filter(
                (pl.col("CHR") == chr) &
                (pl.col("BP") >= start) &
                (pl.col("BP") <= end)
            )

            if mediator.height == 0:
                print(f"[CONCERN] No mediator SNPs in cis-region for {protein} / {mediator_id}")
                continue

            overlap_snps = (
                set(pqtl["SNP"].to_list()) &
                set(gwas["SNP"].to_list()) &
                set(mediator["SNP"].to_list())
            )

            if len(overlap_snps) == 0:
                print(f"[CONCERN] No overlapping SNPs for {protein} / {mediator_id}")
                continue

            mediator_keep = mediator.filter(
                pl.col("SNP").is_in(overlap_snps)
            )

            mediator_keep.write_parquet(
                mediator_out / f"{mediator_id}.parquet"
            )

            print(
                f"[TRACKING] Saved mediator {mediator_id}.parquet for {protein} "
                f"with {len(overlap_snps)} overlapping SNPs"
            )

    # run cis-MR: Protein -> Mediator
    for mediator_file in mediator_gwas.iterdir():
        if mediator_file.suffix not in [".tsv", ".txt"]:
            continue

        mediator_id = mediator_file.stem
        print(f"[TRACKING] Running cis-MR X -> {mediator_id}")

        # NOTE: this reads via paths.qc_out(mediator_id), which mirrors the *main*
        # pheno QC path template (dat/derived/{{id}}/qc_gwas/{{id}}.tsv) rather than
        # paths.qc_mediator_out(mediator_id) (the actual flat mediator QC location
        # used everywhere else in this file). That mismatch pre-dates this migration
        # (it was "./results/QC/{{mediator_id}}/{{mediator_id}}.tsv" vs the real
        # "./results/QC/mediators/{{mediator_id}}.tsv") - preserved as-is here rather
        # than silently fixed, since fixing it changes behavior rather than just
        # relocating a path. Flag for follow-up: this file likely never resolves.
        cmd = f"""
set -euo pipefail
Rscript bin/cis_mr.R \\
    {pqtl_dataset} \\
    dat/cis_regions/{pqtl_dataset} \\
    {mediator_id} \\
    {paths.qc_out(mediator_id)} \\
    {ref_bfile}
mv ./results/cis-MR/{pqtl_dataset}_{mediator_id}_all_MR.tsv \\
   {out_dir}/{pqtl_dataset}_{mediator_id}_all_MR.tsv
"""

        subprocess.run(cmd, shell=True, check=True, executable="/bin/bash")
            

# B_XM: float  
# SE_XM: float
# B_XY: float
# SE_XY: float
# B_MY: float 
# SE_MY: float



def perform_network_mr(
    pheno_id: str,
    pqtl_dataset: str,
    local_results_dir: str = "results",
    ivw_fdr_q: float = 0.05,
    egger_intercept_pval_min: float = 0,
    cochran_q_pval: float = 0.05,
    m_y_pval_threshold: float = 0.05,
):
    mediator_dir = paths.qc_mediator_dir() # ***** PROBS NEED TO CHANGE THIS TO THE MANIFEST ITSELF FOR IT NOT TO BREAK WITH >1 RUN/S
    mediators = [file.stem for file in mediator_dir.glob("*.tsv")]
    results = []

    # out_dir for X->Ms
    X_to_M = paths.network_mr_x_m_dir(pqtl_dataset, local_results_dir)

    # read M -> Y results (Tier 1, pheno-scoped only)
    df_M_to_Y = pl.read_csv(paths.network_mr_m_y_out(pheno_id), separator="\t", null_values=["NA", "NaN", "nan", ""])

    # read X -> Y results once
    df_X_to_Y = pl.read_csv(paths.mr_out(pqtl_dataset, pheno_id, local_results_dir), separator="\t", null_values=["NA", "NaN", "nan", ""])

    for m in mediators:
        row_M_to_Y = df_M_to_Y.filter(pl.col("mediator") == m)

        if row_M_to_Y.height == 0:
            print(f"[CONCERN] No M -> Y result for {m}")
            continue

        n_instruments = int(row_M_to_Y["n_instruments"][0])

        if n_instruments == 1:
            m_y_method = "Wald ratio"
            m_y_beta = row_M_to_Y["Wald_beta"][0]
            m_y_se = row_M_to_Y["Wald_se"][0]
            m_y_pval = row_M_to_Y["Wald_pval"][0]
        else:
            m_y_method = "IVW"
            m_y_beta = row_M_to_Y["IVW_beta"][0]
            m_y_se = row_M_to_Y["IVW_se"][0]
            m_y_pval = row_M_to_Y["IVW_pval"][0]

        if m_y_beta is None or m_y_se is None or m_y_pval is None:
            print(f"[CONCERN] Missing {m_y_method} result for {m}. Skipping.")
            continue

        # declare cis-MR result output for mediator M
        cis_X_to_M = paths.network_mr_x_m_out(pqtl_dataset, m, local_results_dir)

        if not cis_X_to_M.exists():
            print(f"[CONCERN] Missing X -> M cis-MR file for {m}")
            continue

        df_X_to_M = pl.read_csv(cis_X_to_M, separator="\t", null_values=["NA", "NaN", "nan", ""])

        if m_y_pval < m_y_pval_threshold:
            print(f"[TRACKING] All good! M -> Y {m_y_method} p-value passed!")
        else:
            print(f"[CONCERN] {m} -> {pheno_id} {m_y_method} p-value failed. Skipping.")
            continue

        proteins = set(df_X_to_M["protein"].to_list()).intersection(set(df_X_to_Y["protein"].to_list()))

        # next condition for NetworkMR - for each protein in X -> AD
        # (IVW_FDR < 0.05 and Cochran Q-p > 0.05 and Egger intercept > 0.05)
        # and the same thing on X -> M for this given M
        for p in proteins:
            row_X_to_M = df_X_to_M.filter(pl.col("protein") == p)
            row_X_to_Y = df_X_to_Y.filter(pl.col("protein") == p)

            # X -> M
            cis_ivw_p = row_X_to_M["IVW_FDR_q"][0]
            egger_i_p = row_X_to_M["egger_intercept_pval"][0]
            cochan_p = row_X_to_M["Q_pval"][0]

            # X -> Y
            cis_ivw_p_x_y = row_X_to_Y["IVW_FDR_q"][0]
            egger_i_p_x_y = row_X_to_Y["egger_intercept_pval"][0]
            cochan_p_x_y = row_X_to_Y["Q_pval"][0]

            if cis_ivw_p is None or egger_i_p is None or cochan_p is None:
                print(f"[CONCERN] Missing X -> M IVW sensitivity results for {p} / {m}. Skipping.")
                continue

            if cis_ivw_p_x_y is None or egger_i_p_x_y is None or cochan_p_x_y is None:
                print(f"[CONCERN] Missing X -> Y IVW sensitivity results for {p}. Skipping.")
                continue

            if cis_ivw_p < ivw_fdr_q and egger_i_p > egger_intercept_pval_min and cochan_p > cochran_q_pval:
                print(f"[TRACKING] Protein {p} -> passed X -> M cis-MR!")

                if cis_ivw_p_x_y < ivw_fdr_q and egger_i_p_x_y > egger_intercept_pval_min and cochan_p_x_y > cochran_q_pval:
                    print(f"[TRACKING] Protein {p} -> passed X -> Y cis-MR too! Carried forward for NetworkMR!")

                    # networkMR
                    res = NetworkMR(
                        B_XM=row_X_to_M["IVW_beta"][0],
                        SE_XM=row_X_to_M["IVW_se"][0],
                        B_XY=row_X_to_Y["IVW_beta"][0],
                        SE_XY=row_X_to_Y["IVW_se"][0],
                        B_MY=m_y_beta,
                        SE_MY=m_y_se,
                    )

                    # store networkMR res
                    results.append({
                        "protein": p,
                        "mediator": m,
                        "pheno_id": pheno_id,
                        "pqtl_dataset": pqtl_dataset,
                        "X_M_IVW_beta": row_X_to_M["IVW_beta"][0],
                        "X_M_IVW_se": row_X_to_M["IVW_se"][0],
                        "X_M_IVW_FDR_q": cis_ivw_p,
                        "X_Y_IVW_beta": row_X_to_Y["IVW_beta"][0],
                        "X_Y_IVW_se": row_X_to_Y["IVW_se"][0],
                        "X_Y_IVW_FDR_q": cis_ivw_p_x_y,
                        "M_Y_method": m_y_method,
                        "M_Y_n_instruments": n_instruments,
                        "M_Y_beta": m_y_beta,
                        "M_Y_se": m_y_se,
                        "M_Y_pval": m_y_pval,
                        **res
                    }) 
    
    # saving networkMR res
    mediation_out = paths.network_mr_mediation_estimates_out(pqtl_dataset, pheno_id, local_results_dir)
    mediation_out.parent.mkdir(parents=True, exist_ok=True)
    if results:
        pl.DataFrame(results).write_csv(mediation_out, separator="\t")
    else:
        print("[TRACKING] No protein-mediator pairs passed filters for NetworkMR.")


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--pheno_id", required=True)
    p.add_argument("--pheno_gwas", required=True)
    p.add_argument("--ref_bfile", required=True)
    p.add_argument("--pqtl_dataset", required=True)
    p.add_argument("--pqtl_dir", required=True)
    p.add_argument("--run_genomewide_mr", action="store_true")
    p.add_argument("--run_cis_mr_X_M", action="store_true")
    p.add_argument("--run_network_mr", action="store_true")
    p.add_argument("--local_results_dir", default="results")
    p.add_argument("--ivw_fdr_q", type=float, default=0.05)
    p.add_argument("--egger_intercept_pval_min", type=float, default=0)
    p.add_argument("--cochran_q_pval", type=float, default=0.05)
    p.add_argument("--m_y_pval_threshold", type=float, default=0.05)
    args = p.parse_args()

    if args.run_genomewide_mr:
        run_genomewide_mr(
            ref_bfile=args.ref_bfile,
            pheno_id=args.pheno_id,
            pheno_gwas=args.pheno_gwas,
        )

    if args.run_cis_mr_X_M:
        run_cis_mr_X_M(
            pqtl_dataset=args.pqtl_dataset,
            pqtl_dir=args.pqtl_dir,
            ref_bfile=args.ref_bfile,
            local_results_dir=args.local_results_dir,
        )

    if args.run_network_mr:
        perform_network_mr(
            pheno_id=args.pheno_id,
            pqtl_dataset=args.pqtl_dataset,
            local_results_dir=args.local_results_dir,
            ivw_fdr_q=args.ivw_fdr_q,
            egger_intercept_pval_min=args.egger_intercept_pval_min,
            cochran_q_pval=args.cochran_q_pval,
            m_y_pval_threshold=args.m_y_pval_threshold,
        )


if __name__ == "__main__":
    main()
