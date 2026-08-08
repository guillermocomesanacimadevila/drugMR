#!/usr/bin/env python3

from pathlib import Path
import os
import polars as pl
import subprocess


def wget_it(dir, url, file):
    cmd = f"cd {dir} && wget {url}/{file}"
    return subprocess.run(cmd, shell=True, check=True, executable="/bin/bash")


def unzip_it(file, dir):
    cmd = f"unzip -o {file} -d {dir}"
    return subprocess.run(cmd, shell=True, check=True, executable="/bin/bash")


def reverse_etl(file, out_name):
    cmd = f"smr --beqtl-summary {file} --query 1 --out {out_name}"
    return subprocess.run(cmd, shell=True, check=True, executable="/bin/bash")


def wget_files(dir, url):
    dir = Path(dir)

    regions = {
        "Brain_Amygdala.zip": 180,
        "Brain_Anterior_cingulate_cortex_BA24.zip": 233,
        "Brain_Caudate_basal_ganglia.zip": 300,
        "Brain_Cerebellar_Hemisphere.zip": 277,
        "Brain_Cerebellum.zip": 266,
        "Brain_Cortex.zip": 270,
        "Brain_Frontal_Cortex_BA9.zip": 269,
        "Brain_Hippocampus.zip": 255,
        "Brain_Hypothalamus.zip": 257,
        "Brain_Nucleus_accumbens_basal_ganglia.zip": 285,
        "Brain_Putamen_basal_ganglia.zip": 254,
        "Brain_Spinal_cord_cervical_c-1.zip": 204,
        "Brain_Substantia_nigra.zip": 183,
        "Whole_Blood.zip": 800
    }

    for file in regions:
        region = file.replace(".zip", "")
        region_dir = dir / region
        zip_file = region_dir / file
        os.makedirs(region_dir, exist_ok=True)

        if list(region_dir.glob("*.besd")):
            if zip_file.exists():
                zip_file.unlink()
            continue

        wget_it(region_dir, url, file)
        unzip_it(zip_file, region_dir)
        zip_file.unlink()

    return regions


def clean_dir(dir):
    dir = Path(dir)

    for file in dir.glob("*.parquet"):
        if "cis_qtl_pairs" in file.name:
            file.unlink()

    for file in dir.glob("*.txt"):
        file.unlink()

    for file in dir.glob("*.zip"):
        file.unlink()


def etl_workflow_gtex(path, regions):
    path = Path(path)

    for region, n in regions.items():
        region = region.replace(".zip", "")
        region_dir = path / region
        final = region_dir / f"{region}.parquet"

        if final.exists():
            continue

        for file in region_dir.glob("*.besd"):
            file = file.with_suffix("")
            txt = Path(f"{file}.txt")
            parquet = Path(f"{file}.parquet")

            if parquet.exists():
                if txt.exists():
                    txt.unlink()
                continue

            if not txt.exists():
                reverse_etl(file, file)

            df = pl.scan_csv(txt, separator="\t")
            a1 = pl.col("A1").str.to_uppercase()
            a2 = pl.col("A2").str.to_uppercase()
            bases = ["A", "C", "T", "G"]

            df = (
                df
                .drop_nulls()
                .filter(
                    (a1.str.len_chars() == 1) &
                    (a2.str.len_chars() == 1) &
                    a1.is_in(bases) &
                    a2.is_in(bases)
                )
                .with_columns(pl.lit(n).alias("N"))
            )

            df.sink_parquet(parquet)
            txt.unlink()

        df = pl.scan_parquet(
            str(region_dir / "*.cis_qtl_pairs.*.parquet")
        )

        df.sink_parquet(final)
        clean_dir(region_dir)


if __name__ == "__main__":
    path = "./dat/bulk-eQTLs/GTEx_v10"
    url = "https://yanglab.westlake.edu.cn/data/SMR/GTEx_V10_cis_eqtl_summary"
    regions = wget_files(path, url)
    etl_workflow_gtex(path, regions)