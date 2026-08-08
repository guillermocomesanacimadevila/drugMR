#!/usr/bin/env python3
import subprocess
from pathlib import Path
import polars as pl 
import tarfile
import os 

# https://yanglab.westlake.edu.cn/data/SMR/GTEx_V10_cis_eqtl_summary/Brain_Amygdala.zip
# smr \
#   --beqtl-summary BrainMeta_cis_eQTL_chr1 \
#   --query 1 \
#   --out BrainMeta_cis_eQTL_chr1_txt

def reverse_etl(file, out_name):

    cmd = f"""
smr \
    --beqtl-summary {file} \
    --query 1 \
    --out {out_name}
    """

    return subprocess.run(cmd, shell=True, check=True, executable="/bin/bash")


def etl_workflow_metabrain(path: str | Path):
    path = Path(path)
    for file in os.listdir(path):
        if file.endswith(".besd"):
            file = path / file
            file = file.with_suffix("")
            out_name = file
            reverse_etl(file=file, out_name=out_name)
            df = pl.scan_csv(f"{out_name}.txt", separator="\t")
            bases = ["A", "C", "T", "G"]
            a1 = pl.col("A1").str.to_uppercase()
            a2 = pl.col("A2").str.to_uppercase()
            ok_len = (a1.str.len_chars() == 1) & (a2.str.len_chars() == 1)
            ok_bases = a1.is_in(bases) & a2.is_in(bases)
            no_gap = ~a1.str.contains("-") & ~a2.str.contains("-")

            df = (
                df
                .drop_nulls()
                .filter(ok_len & ok_bases & no_gap)
                .with_columns(
                    pl.lit(2865).alias("N")
                )
            )

            df.sink_parquet(f"{out_name}.parquet")
            Path(f"{out_name}.txt").unlink()

    df = pl.scan_parquet(
        str(path / "BrainMeta_cis_eQTL_chr*.parquet")
    )

    df.sink_parquet(
        path / "BrainMeta_cis_eQTL.parquet"
    )


def clean_dir(dir):
    dir = Path(dir)
    for file in dir.glob("*.parquet"):
        if any(file.name.endswith(f"chr{i}.parquet") for i in range(1, 23)):
            os.remove(file)


def qc_cis_eqtls(dir):
    dir = Path(dir)
    eqtl = pl.scan_parquet(dir / "BrainMeta_cis_eQTL.parquet")
    eqtl = eqtl.rename({
        "Chr": "CHR",
        "Freq": "FRQ",
        "Probe": "PROBE",
        "Probe_Chr": "PROBE_CHR",
        "Probe_bp": "PROBE_BP",
        "Gene": "GENE",
        "Orientation": "ORIENTATION",
        "b": "BETA",
        "p": "P"
    })

    eqtl = eqtl.with_columns((pl.col("BETA") / pl.col("SE")).alias("Z"))
    return eqtl


if __name__ == "__main__":
    path = "./dat/bulk-eQTL/MetaBrain/"
    df = etl_workflow_metabrain(path) 
    clean_dir(path)
    qc_cis_eqtls(path)
