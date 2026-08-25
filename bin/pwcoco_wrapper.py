import argparse
from pathlib import Path
from drugmr.pwcoco import PWCoCo
from drugmr.utils import grab_cis_mr_hits
from drugmr.paths import pwcoco_out, pwcoco_raw_dir, pwcoco_raw_prefix, mr_out
import polars as pl 

""" Workflow wrapper (pre-Snakemake) for PWCoCo """

def run_pwcoco(ref_bile: str, n_cases: int, n_controls: int, pqtl_dataset: str, pheno_id: str, cochran_q_thresh: float = 0.05, causal_thresh: float = 0.05):

    """
    PWCoCo wrapper for workflow
    """

    pwcoco = PWCoCo()
    df = mr_out(pheno_id=pheno_id, pqtl_dataset=pqtl_dataset)
    df = pl.read_csv(df, separator="\t")
    targets = grab_cis_mr_hits(
        csv_file=df,
        cochran_q_thresh=cochran_q_thresh,
        causal_thresh=causal_thresh,
    )

    for target in targets:
        cis_regions = f"./dat/cis_regions/{pqtl_dataset}/{target}"
        cis_regions = Path(cis_regions)
        for file in cis_regions.glob("*.parquet"):
            gwas = "gwas.parquet"
            pqtl = "pqtl.parquet"

            # harmonise sumstats (pre-PWCoCo)
            # GWAS
            gwas = pwcoco.harmonise_sumstats(
                df=gwas,
                snp_col="SNP",
                a1_col="A1",
                a2_col="A2",
                maf_col="FRQ",
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
                maf_col="FRQ",
                beta_col="BETA",
                se_col="SE",
                p_col="P",
                n_col="N"
            )

            # run PWCoCo
            pwcoco.pwcoco(
                ref_bfile=ref_bile,
                sumstats_1=pqtl,
                sumstats_2=gwas,
                threads=8,
                n_2=(n_cases + n_controls),
                n2_case=n_cases,
                # out_dir=
                # n_1
            )
        