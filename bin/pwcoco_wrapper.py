import argparse
from pathlib import Path
from drugmr.pwcoco import PWCoCo
from drugmr.utils import grab_cis_mr_hits
from drugmr.paths import pwcoco_out, pwcoco_raw_prefix, mr_out
import polars as pl

""" Workflow wrapper (pre-Snakemake) for PWCoCo """


def resolve_maf_col(df: pl.DataFrame) -> str:
    return "FRQ" if "FRQ" in df.columns else "MAF"


def run_pwcoco(ref_bfile: str, n_cases: int, n_controls: int, pqtl_dataset: str, pheno_id: str, local_results_dir: str = "results", cochran_q_thresh: float = 0.05, causal_thresh: float = 0.05):

    """
    PWCoCo wrapper for workflow
    """

    pwcoco = PWCoCo()
    targets = grab_cis_mr_hits(
        csv_file=mr_out(pqtl_dataset, pheno_id, local_results_dir),
        cochran_q_thresh=cochran_q_thresh,
        causal_thresh=causal_thresh,
    )

    results = []
    for target in targets:
        cis_regions = Path(f"./dat/cis_regions/{pqtl_dataset}/{target}")
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

        # run PWCoCo
        pwcoco.pwcoco(
            ref_bfile=ref_bfile,
            sumstats_1=pqtl,
            sumstats_2=gwas,
            threads=8,
            n_1=int(pqtl["n"][0]),
            n_2=(n_cases + n_controls),
            n2_case=n_cases,
            out_dir=str(out_prefix),
        )

        # read this target's own .coloc output back and tag it - deliberately
        # stale .coloc files for this (pqtl_dataset, pheno_id) can't get mixed in
        coloc_file = Path(f"{out_prefix}.coloc")
        if coloc_file.exists():
            result = pl.read_csv(coloc_file, separator="\t").with_columns(
                pl.lit(target).alias("protein")
            )
            results.append(result)
        else:
            print(f"[CONCERN] Expected PWCoCo output not found: {coloc_file}")

    master_file = pwcoco_out(pqtl_dataset, pheno_id, local_results_dir)
    master_file.parent.mkdir(parents=True, exist_ok=True)
    pl.concat(results, how="diagonal_relaxed").write_csv(master_file, separator="\t")
    print(f"[DONE] Saved master PWCoCo table: {master_file}")


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--pqtl_dataset", required=True, choices=["ukb_ppp", "decode", "wu_csf", "wingo_brain"])
    p.add_argument("--pheno_id", required=True)
    p.add_argument("--ref_bfile", required=True)
    p.add_argument("--n_cases", required=True, type=int)
    p.add_argument("--n_controls", required=True, type=int)
    p.add_argument("--local_results_dir", default="results")
    p.add_argument("--cochran_q_pval", type=float, default=0.05)
    p.add_argument("--wald_fdr_q", type=float, default=0.05)
    args = p.parse_args()

    run_pwcoco(
        ref_bfile=args.ref_bfile,
        n_cases=args.n_cases,
        n_controls=args.n_controls,
        pqtl_dataset=args.pqtl_dataset,
        pheno_id=args.pheno_id,
        local_results_dir=args.local_results_dir,
        cochran_q_thresh=args.cochran_q_pval,
        causal_thresh=args.wald_fdr_q,
    )


if __name__ == "__main__":
    main()
