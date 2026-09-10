import argparse
from pathlib import Path

import polars as pl

from drugmr.paths import mr_out, pwcoco_out, pwcoco_raw_prefix
from drugmr.pwcoco import PWCoCo
from drugmr.utils import grab_cis_mr_hits

""" Workflow wrapper (pre-Snakemake) for PWCoCo """

# real native pwcoco tool output columns (confirmed against a real production
# pwcoco.tsv) - used only when zero targets pass the cis-MR gate, so the
# written file still has the columns every downstream reader expects, instead
# of a schema-less empty file that raises polars.exceptions.NoDataError on
# the very next read.
_PWCOCO_COLUMNS = ["SNP1", "SNP2", "nsnps", "H0", "H1", "H2", "H3", "H4", "log_abf_all", "protein", "outcome_trait"]


def resolve_maf_col(df: pl.DataFrame) -> str:
    return "FRQ" if "FRQ" in df.columns else "MAF"


def run_pwcoco(
    ref_bfile: str,
    n_cases: int,
    n_controls: int,
    pqtl_dataset: str,
    pheno_id: str,
    local_results_dir: str = "results",
    wald_fdr_q: float = 0.05,
    ivw_fdr_q: float = 0.05,
    cochran_q_pval: float = 0.05,
    egger_intercept_pval_min: float = 0,
    min_instruments_for_ivw: int = 3,
    cis_regions_dir: str | None = None,
):

    """
    PWCoCo wrapper for workflow
    """

    pwcoco = PWCoCo()
    # same canonical cis-MR pass/fail rule standard COLOC uses (bin/coloc_targets.py's
    # pairwise_coloc(), via drugmr.utils.select_cis_mr_passing_proteins) - a target
    # that fails cis-MR for COLOC can no longer silently still reach PWCoCo
    targets = grab_cis_mr_hits(
        csv_file=mr_out(pqtl_dataset, pheno_id, local_results_dir),
        wald_fdr_q=wald_fdr_q,
        ivw_fdr_q=ivw_fdr_q,
        cochran_q_pval=cochran_q_pval,
        egger_intercept_pval_min=egger_intercept_pval_min,
        min_instruments_for_ivw=min_instruments_for_ivw,
    )

    results = []
    for target in targets:
        # 1 bad target (a missing/degenerate cis-region, a harmonisation
        # failure, a native pwcoco crash) must not lose every other
        # already-computed target - same reasoning as cis_mr.R's own
        # per-protein tryCatch
        try:
            cis_regions = Path(cis_regions_dir) / target if cis_regions_dir else Path(f"./dat/cis_regions/{pqtl_dataset}/{target}")
            gwas = pl.read_parquet(cis_regions / "gwas.parquet")
            pqtl = pl.read_parquet(cis_regions / "pqtl.parquet")

            # harmonise sumstats (pre-PWCoCo)
            # GWAS
            gwas = pwcoco.harmonise_sumstats(
                df=gwas,
                snp_col="SNP",
                a1_col="A1",
                a2_col="A2",
                maf_col=resolve_maf_col(gwas),
                beta_col="BETA",
                se_col="SE",
                p_col="P",
                n_col="N"
            )

            # pQTL
            pqtl = pwcoco.harmonise_sumstats(
                df=pqtl,
                snp_col="SNP",
                a1_col="A1",
                a2_col="A2",
                maf_col=resolve_maf_col(pqtl),
                beta_col="BETA",
                se_col="SE",
                p_col="P",
                n_col="N"
            )

            out_prefix = pwcoco_raw_prefix(pqtl_dataset, pheno_id, target, local_results_dir)
            out_prefix.parent.mkdir(parents=True, exist_ok=True)

            # n_1 mirrors coloc.R's own convention (max N across the pQTL
            # cis-region, not first-row) so both paths report the same
            # sample size for the same protein
            n_pqtl = pqtl["n"].max()
            if n_pqtl is None:
                raise ValueError(f"pQTL N column is entirely null for {target}")

            # run PWCoCo
            pwcoco.pwcoco(
                ref_bfile=ref_bfile,
                sumstats_1=pqtl,
                sumstats_2=gwas,
                threads=8,
                n_1=int(n_pqtl),
                n_2=(n_cases + n_controls),
                n2_case=n_cases,
                out_dir=str(out_prefix),
            )

            # read this target's own .coloc output back and tag it - deliberately
            # stale .coloc files for this (pqtl_dataset, pheno_id) can't get mixed in
            coloc_file = Path(f"{out_prefix}.coloc")
            if coloc_file.exists():
                result = pl.read_csv(coloc_file, separator="\t").with_columns(
                    pl.lit(target).alias("protein"),
                    pl.lit(pheno_id).alias("outcome_trait"),
                )
                results.append(result)
            else:
                print(f"[CONCERN] Expected PWCoCo output not found: {coloc_file}")
        except Exception as error:
            print(f"[CONCERN] PWCoCo failed for {target} - continuing without it: {error}")

    master_file = pwcoco_out(pqtl_dataset, pheno_id, local_results_dir)
    master_file.parent.mkdir(parents=True, exist_ok=True)
    if results:
        master_df = pl.concat(results, how="diagonal_relaxed")
        master_df = master_df.drop([c for c in ("Dataset1", "Dataset2") if c in master_df.columns])
    else:
        print("[TRACKING] No targets passed cis-MR filters - writing an empty PWCoCo table")
        master_df = pl.DataFrame(schema=_PWCOCO_COLUMNS)
    master_df.write_csv(master_file, separator="\t")
    print(f"[DONE] Saved master PWCoCo table: {master_file}")


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--pqtl_dataset", required=True)
    p.add_argument("--pheno_id", required=True)
    p.add_argument("--ref_bfile", required=True)
    p.add_argument("--n_cases", required=True, type=int)
    p.add_argument("--n_controls", required=True, type=int)
    p.add_argument("--local_results_dir", default="results")
    p.add_argument("--wald_fdr_q", type=float, default=0.05)
    p.add_argument("--ivw_fdr_q", type=float, default=0.05)
    p.add_argument("--cochran_q_pval", type=float, default=0.05)
    p.add_argument("--egger_intercept_pval_min", type=float, default=0)
    p.add_argument("--min_instruments_for_ivw", type=int, default=3)
    p.add_argument("--cis_regions_dir", default=None)
    args = p.parse_args()

    run_pwcoco(
        ref_bfile=args.ref_bfile,
        n_cases=args.n_cases,
        n_controls=args.n_controls,
        pqtl_dataset=args.pqtl_dataset,
        pheno_id=args.pheno_id,
        local_results_dir=args.local_results_dir,
        wald_fdr_q=args.wald_fdr_q,
        ivw_fdr_q=args.ivw_fdr_q,
        cochran_q_pval=args.cochran_q_pval,
        egger_intercept_pval_min=args.egger_intercept_pval_min,
        min_instruments_for_ivw=args.min_instruments_for_ivw,
        cis_regions_dir=args.cis_regions_dir,
    )


if __name__ == "__main__":
    main()
