#!/usr/bin/env python3
import polars as pl
import pandas as pd
import subprocess
import argparse
from pathlib import Path

# look at MR results based on dataset X (which == arg)
# IVW p_FDR < 0.05 and passes egger intercept and cochran Q
# go to that gene pQTL df["protein"] - and acess corresponding directory in dat/cis_regions/{dataset}
# grab parquet files and cmd run for ./bin/coloc.R script

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
            (pl.col("IVW_FDR_q") < 1) & # 0.05 -> 1 for CI/CD testing
            (pl.col("egger_intercept_pval") > 0.05) &
            (pl.col("Q_pval") > 0.05)
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

def main():
    p = argparse.ArgumentParser()
    p.add_argument("--pqtl_dataset", required=True, choices=["ukb_ppp", "decode"])
    p.add_argument("--local_results_dir", required=True)
    p.add_argument("--pqtl_dir", required=True)
    p.add_argument("--pheno_id", required=True)
    p.add_argument("--n_cases", required=True, type=int)
    p.add_argument("--n_controls", required=True, type=int)
    args = p.parse_args()
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