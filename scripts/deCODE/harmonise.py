#!/usr/bin/env python3
import polars as pl
from pathlib import Path

# remove NaNs
# only keep biallelic SNPs

def harmonise_decode_pqtls():
    dir = Path("./dat/pQTL/deCODE")
    files = sorted(dir.glob("*.parquet"))
    print(f"[TRACKING] Found {len(files)} deCODE pQTL files")
    for i, file in enumerate(files, start=1):
        before = pl.scan_parquet(file).select(pl.len()).collect().item()
        df = (
            pl.scan_parquet(file)
            .with_columns(
                pl.col("A1").str.to_uppercase().str.strip_chars(),
                pl.col("A2").str.to_uppercase().str.strip_chars(),
                pl.col("SNP").str.strip_chars(),
            )
            .filter(
                pl.col("CHR").is_not_null()
                & pl.col("BP").is_not_null()
                & pl.col("A1").is_not_null()
                & pl.col("A2").is_not_null()
                & pl.col("FRQ").is_not_null()
                & pl.col("SNP").is_not_null()
                & pl.col("BETA").is_not_null()
                & pl.col("P").is_not_null()
                & pl.col("A1").is_in(["A", "C", "G", "T"])
                & pl.col("A2").is_in(["A", "C", "G", "T"])
                & (pl.col("A1") != pl.col("A2"))
                & (~pl.col("SNP").is_in(["NA", "NaN", "nan", ".", ""]))
                & pl.col("FRQ").is_finite()
                & pl.col("BETA").is_finite()
                & pl.col("P").is_finite()
                & (pl.col("FRQ") > 0)
                & (pl.col("FRQ") < 1)
                & (pl.col("P") >= 0)
                & (pl.col("P") <= 1)
            )
            .unique(
                subset=["CHR", "BP", "A1", "A2"],
                keep="first",
            )
            .sort(["CHR", "BP"])
            .collect(engine="streaming")
        )

        after = df.height
        df.write_parquet(file, compression="zstd", statistics=True,)
        print(f"[TRACKING] {i}/{len(files)} {file.name}: {before:,} -> {after:,} ({before-after:,} removed)")
    print("[TRACKING] deCODE pQTL harmonisation completed")

if __name__ == "__main__":
    harmonise_decode_pqtls()