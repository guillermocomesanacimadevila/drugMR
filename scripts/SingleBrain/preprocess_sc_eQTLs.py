#!/usr/bin/env python3
from zenodo_get import download
from pathlib import Path
import polars as pl
import os
import subprocess

# https://zenodo.org/records/14908182
# test
# Ast_disease_eqtl_full_assoc.tsv.gz
# f{cell_type}_disease_eqtl_full_assoc.tsv.gz
# 14908182
# Ast: astrocytes
# Ext: excitatory neurons
# IN: inhibitory neurons
# MG: microglia
# OD: oligodendrocytes
# OPC: oligodendrocyte progenitor cell
# OLD ONE -> 16051904

# for each cell-type full_asoc.tsv.gz
# download and store in dir X -> and gunzip
# harmonise to LDSC format
# overwrite the original TSV -> AND save as f{cell_type}.parquet

# SMR-ready subworkflow

def preprocess_single_brain_eqtls():
    # REMINDER FOR MYSELF -> PRINT THE SIZE OF THE FILE AT EACH STAGE
    # du -sh at - .tsv.gz -> .tsv -> .parquet (with QC)
    print("[TRACKING] Starting SingleBrain preprocessing...")
    cell_types = ["Ast", "Ext", "IN", "MG", "OD", "OPC", "End"]
    zenodo_doi = "14908182"
    out_dir = "./dat/sc-eQTL/SingleBrain"
    out_dir = Path(out_dir)
    os.makedirs(out_dir, exist_ok=True)
    bases = ["A", "C", "T", "G"]
    N = 983

    for cell_type in cell_types:
        print("\n" + "~" * 80)
        print(f"[TRACKING] Processing {cell_type}")
        print("~" * 80)
        gz_path = out_dir / f"{cell_type}_eqtl_full_assoc.tsv.gz"
        tsv_path = out_dir / f"{cell_type}_eqtl_full_assoc.tsv"
        parquet_path = out_dir / f"{cell_type}.parquet"
        print("[TRACKING] Downloading...")
        download(
            record_or_doi=f"{zenodo_doi}",
            output_dir=out_dir,
            file_glob=f"{cell_type}_eqtl_full_assoc.tsv.gz",
        )
        print(f"[TRACKING] Download complete: {gz_path}")
        print("[TRACKING] Decompressing...")
        subprocess.run(f"gunzip -f {gz_path}", shell=True, check=True, executable="/bin/bash",)
        print(f"[TRACKING] Decompressed: {tsv_path}")
        print("[TRACKING] Reading TSV lazily...")

        df = pl.scan_csv(
            tsv_path,
            separator="\t",
            null_values=["NA"],
            infer_schema_length=10000,
            schema_overrides={
                "feature": pl.Utf8,
                "variant_id": pl.Utf8,
                "chr": pl.Utf8,
                "pos": pl.Int64,
                "ref": pl.Utf8,
                "alt": pl.Utf8,
                "Allele": pl.Utf8,
                "fixed_beta": pl.Float64,
                "fixed_sd": pl.Float64,
                "fixed_z": pl.Float64,
                "Fixed_P": pl.Float64,
            },
        )

        print("[TRACKING] Harmonising columns...")

        df = (
            df
            .select([
                "feature",
                "variant_id",
                "chr",
                "pos",
                "ref",
                "alt",
                "Allele",
                "fixed_beta",
                "fixed_sd",
                "fixed_z",
                "Fixed_P",
            ])
            .rename({
                "feature": "GENE",
                "variant_id": "SNP",
                "chr": "CHR",
                "pos": "BP",
                "ref": "A1",
                "alt": "A2",
                "Allele": "EA",
                "fixed_beta": "BETA",
                "fixed_sd": "SE",
                "fixed_z": "Z",
                "Fixed_P": "P",
            })
        )

        print("[TRACKING] Casting datatypes...")

        df = (
            df
            .with_columns([
                pl.col("CHR").str.replace("^chr", "").cast(pl.UInt8),
                pl.col("BP").cast(pl.Int64),
                pl.col("BETA").cast(pl.Float64),
                pl.col("SE").cast(pl.Float64),
                pl.col("Z").cast(pl.Float64),
                pl.col("P").cast(pl.Float64),
                pl.col("A1").str.to_uppercase(),
                pl.col("A2").str.to_uppercase(),
                pl.col("EA").str.to_uppercase(),
            ])
        )

        print("[TRACKING] Removing missing values...")

        df = df.drop_nulls([
            "CHR",
            "BP",
            "A1",
            "A2",
            "EA",
            "BETA",
            "SE",
            "Z",
            "P",
        ])

        print("[TRACKING] Performing SNP QC...")

        df = (
            df
            .filter(
                (pl.col("A1").str.len_chars() == 1) &
                (pl.col("A2").str.len_chars() == 1) &
                (pl.col("A1").is_in(bases)) &
                (pl.col("A2").is_in(bases)) &
                (pl.col("SE") > 0) &
                (pl.col("P") > 0) &
                (pl.col("P") <= 1)
            )
            .with_columns(
                pl.lit(N).alias("N")
            )
        )

        print("[TRACKING] Writing parquet...")
        df.sink_parquet(parquet_path)
        print(f"[TRACKING] Saved: {parquet_path}")
        print("[TRACKING] Removing intermediate TSV...")
        tsv_path.unlink()
        print(f"[TRACKING] Finished {cell_type}")
    print("\n[TRACKING] All SingleBrain cell types processed successfully.")

if __name__ == "__main__":
    preprocess_single_brain_eqtls()