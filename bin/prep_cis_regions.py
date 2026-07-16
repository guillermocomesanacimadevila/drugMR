#!/usr/bin/env python3
import polars as pl
import os 
import argparse
from pathlib import Path

# grab .parquet files from pQTLs
# add fixed N based on either sample size (either deCODE or UKB-PPP)
# for each parquet file - create a directory specific to it 
# map grab exactly those same SNPs on the .parquet file and map the same SNPs on outcome GWAS
# save locus from GWAS specific to protein X onto the same dir as the protein
# make sure its harmonised to LDSC format
# these will be the ones used for MR and COLOC

decode_n = 35559
ukb_ppp_n = 54219

# pre-established args from notebook
# * pheno_id
# * pqtl_dir
# * ref_bfile
# * pqtl_dataset
# * colnames (pQTL and GWAS)

def define_loci_from_cis_regions(pqtl_dataset: str, pheno_id: str, pqtl_dir: str):
    gwas = pl.read_csv(f"./results/QC/{pheno_id}/{pheno_id}.tsv", separator="\t")
    pqtl_dir = Path(pqtl_dir)
    pqtl_dataset = pqtl_dataset.lower()
    for file in pqtl_dir.glob("*.parquet"):
        gene = file.stem.split("_")[0]
        protein = file.stem.split("_")[1]   
        out_dir = Path(f"./dat/cis_regions/{pqtl_dataset}/{gene}_{protein}")
        os.makedirs(out_dir, exist_ok=True)
        # we need to move both .parquet files (pQTL and GWAS) into that new dir
        df = pl.read_parquet(file)

        if df.height == 0:
            print(f"[SKIP] {gene}_{protein}: empty pQTL parquet")
            continue

        if pqtl_dataset == "ukb_ppp":
            df = df.with_columns(
                pl.lit(ukb_ppp_n).alias("N")
            )
        elif pqtl_dataset == "decode":
            df = df.with_columns(
                pl.lit(decode_n).alias("N")
            )
        
        # pos
        chr = df.select(pl.col("CHR").cast(pl.Int64).unique()).item()
        start = df.select(pl.col("BP").min()).item()
        end = df.select(pl.col("BP").max()).item()
        df2 = gwas.filter((pl.col("CHR").cast(pl.Int64) == chr) & (pl.col("BP").is_between(start, end)))

        # remove duplicates
        df = (df.sort("P").unique(subset=["SNP"], keep="first"))
        df2 = (df2.sort("P").unique(subset=["SNP"], keep="first"))

        # match SNPs with pQTL
        pqtl_matched = df.join(df2.select("SNP"), on="SNP", how="inner")
        gwas_matched = df2.join(df.select("SNP"), on="SNP", how="inner")

        # save onto out_dir
        pqtl_matched.write_parquet(out_dir / "pqtl.parquet")
        gwas_matched.write_parquet(out_dir / "gwas.parquet")
        print(
            f"{gene}_{protein}: "
            f"pQTL={df.height}, GWAS_region={df2.height}, matched={pqtl_matched.height}"
        )

def main():
    p = argparse.ArgumentParser()
    p.add_argument("--pqtl_dataset", required=True, choices=["ukb_ppp", "decode"])
    p.add_argument("--pheno_id", required=True)
    p.add_argument("--pqtl_dir", required=True)
    args = p.parse_args()
    define_loci_from_cis_regions(
        pqtl_dataset=args.pqtl_dataset,
        pheno_id=args.pheno_id,
        pqtl_dir=args.pqtl_dir,
    )

if __name__ == "__main__":
    main()