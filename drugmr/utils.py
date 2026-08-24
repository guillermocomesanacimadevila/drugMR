#!/usr/bin/env python3
# utils.py
import polars as pl
import subprocess
import numpy as np 

# generic 3+ trait SNP matcher for multi-trait coloc-style analyses (HyPrColoc,
# MOLOC, ...) - each df in `datasets` needs SNP / A1 / A2 / BETA / SE columns and
# 1 row per SNP already (dedup by whatever the caller's significance criterion is,
# e.g. lowest P). Every non-reference dataset gets its effect allele aligned onto
# `reference`'s A1/A2 (BETA flipped where A1/A2 are swapped, SNP dropped where
# neither allele pairing resolves), then all datasets are reduced to the SNPs
# shared across the lot. Returns {name: df[SNP, BETA, SE]}, row-aligned by SNP.
def extract_common_snps(datasets: dict, reference: str):

    if reference not in datasets:
        raise ValueError(f"Reference dataset '{reference}' not found in datasets: {list(datasets.keys())}")

    ref_alleles = datasets[reference].select([
        "SNP",
        pl.col("A1").alias("ref_a1"),
        pl.col("A2").alias("ref_a2")
    ])

    aligned = {}
    for name, df in datasets.items():
        if name == reference:
            aligned[name] = df.select(["SNP", "BETA", "SE"])
            continue

        df = df.join(ref_alleles, on="SNP", how="inner")
        matched = (pl.col("A1") == pl.col("ref_a1")) & (pl.col("A2") == pl.col("ref_a2"))
        swapped = (pl.col("A1") == pl.col("ref_a2")) & (pl.col("A2") == pl.col("ref_a1"))

        aligned[name] = (
            df
            .filter(matched | swapped)
            .with_columns(
                pl.when(swapped).then(-pl.col("BETA")).otherwise(pl.col("BETA")).alias("BETA")
            )
            .select(["SNP", "BETA", "SE"])
        )

    shared_snps = None

    for df in aligned.values():
        snps = df.select("SNP")
        shared_snps = snps if shared_snps is None else shared_snps.join(snps, on="SNP", how="inner")

    return {name: shared_snps.join(df, on="SNP", how="inner").sort("SNP") for name, df in aligned.items()}


def impute_ld(ref_bfile, snp_1, snp_2):
    cmd = f"""
plink \
    --bfile {ref_bfile} \
    --ld {snp_1} {snp_2}
"""
    return subprocess.run(cmd, shell=True, check=False, executable="/bin/bash", capture_output=True, text=True)


def quick_f_statistic(beta_exposure, se_exposure):
    return (beta_exposure / se_exposure)**2


def lambda_sample_overlap(
        n_overlap,
        n_exposure_total,
        n_outcome_total):
    lambda_res = (
        n_overlap / np.sqrt(n_exposure_total * n_outcome_total)
    )
    return lambda_res


def sample_overlap_relative_bias(lambda_funct, f_statistic):
    raw = lambda_funct / f_statistic
    percent = raw * 100
    return raw, percent