import argparse
import subprocess
from pathlib import Path

import pandas as pd
import polars as pl

from drugmr import paths

# look at MR results based on dataset X (which == arg)
# IVW p_FDR < 0.05 and passes egger intercept and cochran Q
# go to that gene pQTL df["protein"] - and acess corresponding directory in dat/cis_regions/{dataset}
# grab parquet files and cmd run for ./bin/coloc.R script

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
    master = master.rename(columns={"protein_id": "protein"})
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


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--pqtl_dataset", required=True, choices=["ukb_ppp", "decode", "wu_csf", "wingo_brain"])
    p.add_argument("--local_results_dir", required=True)
    p.add_argument("--pqtl_dir", required=True)
    p.add_argument("--pheno_id", required=True)
    p.add_argument("--n_cases", required=True, type=int)
    p.add_argument("--n_controls", required=True, type=int)
    p.add_argument("--wald_fdr_q", type=float, default=0.05)
    p.add_argument("--ivw_fdr_q", type=float, default=0.05)
    p.add_argument("--cochran_q_pval", type=float, default=0.05)
    p.add_argument("--egger_intercept_pval_min", type=float, default=0)
    p.add_argument("--min_instruments_for_ivw", type=int, default=3)
    args = p.parse_args()

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
