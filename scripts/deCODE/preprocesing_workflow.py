#!/usr/bin/env python3
import polars as pl
from pathlib import Path
import os 
import argparse
import subprocess
from urllib.parse import urlparse, parse_qs

# load manifest
# THIS wget command works! 
# wget 'https://download.decode.is/s3/download?token=e1530773-978a-4039-919b-4f9edcace104&file=10000_28_CRYBB2_CRBB2.txt.gz'
# for file in that manifest.csv onto dat/pqtls/deCODE
# download compressed file
# grab cis-region
# QC cis-region
# rename cols to LDSC
# convert to .parquet and save in dat/pqtls/deCODE
# remove old file (the big file) - we only keep the cis.region .parquet
# go onto next protein
# rename 

# 4907 - total
# 4742 - autosomes


def cmd(url, out_dir):
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    filename = parse_qs(urlparse(url).query)["file"][0]
    print(f"[TRACKING] Downloading: {filename}")
    cmd = f"""
set -euo pipefail
cd "{out_dir}"
wget -q --tries=5 --timeout=60 '{url}' -O "{filename}"
    """
    subprocess.run(cmd, shell=True, check=True, executable="/bin/bash")
    return out_dir / filename

def extract_cis_region(file: str, gene: str, chrom: str, start: int, end: int, window: int, out_dir: str, aptamer_id: str):
    out_dir = Path(out_dir)
    df = pl.scan_csv(file, separator="\t").select(["Chrom", "Pos", "effectAllele", "otherAllele", "ImpMAF", "rsids", "Beta", "Pval"])
    df = df.with_columns([
        pl.col("Chrom").cast(pl.Utf8).str.replace("^chr", "").alias("Chrom"),
        pl.col("Pos").cast(pl.Int64).alias("Pos"),
    ])
    cis_start = max(0, int(start) - window)
    cis_end = int(end) + window
    print(
        f"[TRACKING] {gene} ({aptamer_id}) | "
        f"chr{chrom}:{cis_start:,}-{cis_end:,}"
    )

    # extract cis region around gene coords +/- window
    df_cis = (
        df
        .filter(
            (pl.col("Chrom") == str(chrom)) &
            (pl.col("Pos") >= cis_start) &
            (pl.col("Pos") <= cis_end)
        )
    )

    # qc locus
    df_cis = (
        df_cis
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
    df_cis = df_cis.drop_nulls().collect(engine="streaming")
    df_cis.write_parquet(out_dir / f"{gene}_{aptamer_id}.parquet", compression="zstd", compression_level=3)
    print(f"[DONE] Saved {gene}: {df_cis.height:,} SNPs")

def decode_preprocessing_pipeline(out_dir: str = "./dat/pQTL/deCODE"):
    out_dir = Path(out_dir) # "./dat/pqtl/deCODE"
    out_dir.mkdir(parents=True, exist_ok=True)
    scratch_dir = Path(os.environ.get("SLURM_TMPDIR", os.environ.get("TMPDIR", out_dir))) / f"decode_{os.environ.get('SLURM_JOB_ID', 'local')}"
    scratch_dir.mkdir(parents=True, exist_ok=True)
    manifest = "./results/deCODE_manifest/deCODE_eur_pgwas_manifest.csv"
    df = pl.read_csv(manifest)

    print(f"[TRACKING] Total Targets: {df.height:,}")
    print("[TRACKING] Dropping non-autosomal targets...")

    # drop non-autosomal targets
    df = df.filter(~pl.col("chr").cast(pl.Utf8).is_in(["X", "Y"]))
    print(f"[TRACKING] Autosomal Targets: {df.height:,}")

    # df = df.head(3) # testing only
    # print(df)

    print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
    print("[TRACKING] deCODE preprocessing pipeline")
    print(f"[TRACKING] Proteins/aptamers: {df.height:,}")
    print(f"[TRACKING] Output directory: {out_dir}")
    print(f"[TRACKING] Temporary directory: {scratch_dir}")
    print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")

    for i, row in enumerate(df.iter_rows(named=True), start=1):
        url = row["download_url"]
        gene = row["gene_symbol"]
        chrom = row["chr"]
        start = row["start"]
        end = row["end"]
        aptamer_id = row["aptamer_id"]
        parquet_file = out_dir / f"{gene}_{aptamer_id}.parquet"
        print(f"\n[{i:,}/{df.height:,}] {gene} ({aptamer_id})")

        # skip if this protein-aptamer has already been processed
        if parquet_file.exists():
            print(f"[SKIPPING] Already exists: {parquet_file.name}")
            continue

        print("[TRACKING] Downloading...")
        file = cmd(url, scratch_dir) # download compressed file into node scratch
        print("[TRACKING] Extracting cis-region...")
        extract_cis_region(
            file=file,
            gene=gene,
            chrom=chrom,
            start=start,
            end=end,
            window=1_000_000,
            out_dir=out_dir,
            aptamer_id=aptamer_id
        )
        print("[TRACKING] Removing raw file...")
        Path(file).unlink()

    print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
    print("[DONE] deCODE preprocessing completed.")
    print(f"[DONE] Processed {df.height:,} protein-aptamers.")
    print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")


if __name__ == "__main__":
    decode_preprocessing_pipeline()