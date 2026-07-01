#!/usr/bin/env python3
import polars as pl
from pathlib import Path
import os 
import argparse
import subprocess
from urllib.parse import urlparse, parse_qs

# load manifest
# THIS wget command works! 
# wget 'https://download.decode.is/s3/download?token=e1530773-978a-4039-919b-4f9edcace104&file=10000_28_CRYBB2_CRBB2.txt.gz'
# for file in that manifest.csv onto dat/pqtls/deCODE
# download - gunzip
# rename
# grab cis-region
# QC cis-region
# rename cols to LDSC
# convert to .parquet and save in dat/pqtls/deCODE
# remove old file (the big file) - we only keep the cis.region .parquet
# go onto next protein
# rename 

def cmd(url, out_dir):
    out_dir = Path(out_dir)
    filename = parse_qs(urlparse(url).query)["file"][0]
    txt_file = filename.replace(".gz", "")
    cmd = f"""
set -euo pipefail
cd "{out_dir}"
wget -q '{url}' -O "{filename}"
gunzip -f "{filename}"
    """
    subprocess.run(cmd, shell=True, check=True, executable="/bin/bash")
    return out_dir / txt_file

def extract_cis_region(file: str, gene: str, chrom: str, start: int, end: int, window: int, out_dir: str):
    out_dir = Path(out_dir)
    df = pl.read_csv(file, separator="\t")
    cis_start = max(0, int(start) - window)
    cis_end = int(end) + window

    # extract cis region around gene coords +/- window
    df_cis = (
        df
        .filter(
            (pl.col("Chrom").cast(pl.Utf8) == str(chrom)) &
            (pl.col("Pos") >= cis_start) &
            (pl.col("Pos") <= cis_end)
        )
    )

    # qc locus
    df_cis = (
        df_cis
        .drop(["minus_log10_pval", "Name", "N"])
        .rename({
            "Chrom": "CHR",
            "Pos": "BP",
            "effectAllele": "A1",
            "otherAllele": "A2",
            "ImpMAF": "FRQ",
            "rsids": "SNP",
            "Beta": "BETA",
            "Pval": "P"
        })
    )

    # remove INDELs
    bases = ["A", "C", "T", "G"]
    a1 = pl.col("A1").str.to_uppercase()
    a2 = pl.col("A2").str.to_uppercase()
    ok_len = (a1.str.len_chars() == 1) & (a2.str.len_chars() == 1)
    ok_bases = a1.is_in(bases) & a2.is_in(bases)
    no_gap = ~a1.str.contains("-") & ~a2.str.contains("-")
    df_cis = df_cis.filter(ok_len & ok_bases & no_gap)

    # remove empty entries / NANs
    df_cis = df_cis.drop_nulls()
    df_cis.write_parquet(out_dir / f"{gene}.parquet")
    print(f"[DONE] Saved {gene}: {df_cis.height:,} SNPs")


def decode_preprocessing_pipeline(out_dir: str):
    out_dir = Path(out_dir) # "./dat/pqtls/deCODE"
    out_dir.mkdir(parents=True, exist_ok=True)
    manifest = "./results/deCODE_manifest/deCODE_eur_pgwas_manifest.csv"
    df = pl.read_csv(manifest)
    for row in df.iter_rows(named=True):
        url = row["download_url"]
        gene = row["gene_symbol"]
        chrom = row["chr"]
        start = row["start"]
        end = row["end"]
        file = cmd(url, out_dir) # download file into out_dir
        extract_cis_region(
            file=file,
            gene=gene,
            chrom=chrom,
            start=start,
            end=end,
            window=1_000_000,
            out_dir=out_dir
        )
        Path(file).unlink()


        

