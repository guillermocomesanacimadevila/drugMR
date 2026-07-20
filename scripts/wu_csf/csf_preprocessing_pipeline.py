#!/usr/bin/env python3
import polars as pl 
import subprocess
from pathlib import Path
import os 

# load autosomal manifest
# create stable gene_uniprot_seqid ID
# download each GCST file once
# extract cis-region for each mapped gene
# QC cis-region
# rename cols
# convert to .parquet
# remove downloaded file
# go onto next aptamer


def cmd(url, out_dir, gcst):
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    filename = f"{gcst}.tsv.gz"
    print(f"[TRACKING] Downloading: {filename}")
    cmd = f"""
set -euo pipefail
cd "{out_dir}"
wget -q --tries=5 --timeout=60 '{url}' -O "{filename}"
    """
    subprocess.run(cmd, shell=True, check=True, executable="/bin/bash")
    return out_dir / filename


def extract_cis_region(file: str, gene: str, chrom: str, start: int, end: int, window: int, out_dir: str, id_target: str):
    out_dir = Path(out_dir)
    df = pl.scan_csv(file, separator="\t").select(["chromosome", "base_pair_location", "effect_allele", "other_allele", "effect_allele_frequency", "rsid", "beta", "standard_error", "p_value"])
    df = df.with_columns([
        pl.col("chromosome").cast(pl.Utf8).str.replace("^chr", "").alias("chromosome"),
        pl.col("base_pair_location").cast(pl.Int64).alias("base_pair_location")
    ])
    cis_start = max(1, int(start) - window)
    cis_end = int(end) + window
    print(
        f"[TRACKING] {gene} ({id_target}) | "
        f"chr{chrom}:{cis_start:,}-{cis_end:,}"
    )

    # extract cis region around gene coords +/- window
    df_cis = (
        df
        .filter(
            (pl.col("chromosome") == str(chrom)) &
            (pl.col("base_pair_location") >= cis_start) &
            (pl.col("base_pair_location") <= cis_end)
        )
    )

    # qc locus
    df_cis = (
        df_cis
        .rename({
            "chromosome": "CHR",
            "base_pair_location": "BP",
            "standard_error": "SE",
            "effect_allele": "A1",
            "other_allele": "A2",
            "effect_allele_frequency": "FRQ",
            "rsid": "SNP",
            "beta": "BETA",
            "p_value": "P"
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
    df_cis.write_parquet(out_dir / f"{id_target}.parquet", compression="zstd", compression_level=3)
    print(f"[DONE] Saved {gene}: {df_cis.height:,} SNPs")


def csf_preprocesing_pipeline(out_dir: str = "./dat/pQTL/wu_csf"):
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    scratch_dir = Path(os.environ.get("SLURM_TMPDIR", os.environ.get("TMPDIR", out_dir))) / f"wu_csf_{os.environ.get('SLURM_JOB_ID', 'local')}"
    scratch_dir.mkdir(parents=True, exist_ok=True)
    autosomal_manifest = "./results/WS_CSF_manifest/wu_csf_manifest_autosomal.csv"
    df = pl.read_csv(autosomal_manifest)
    # create identifier for each row
    df = df.with_columns((pl.col("Gene").cast(pl.Utf8) + "_" + pl.col("UniProt").cast(pl.Utf8) + "_" + pl.col("matched_SeqId").cast(pl.Utf8)).alias("ID"))

    print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
    print("[TRACKING] WU CSF preprocessing pipeline")
    print(f"[TRACKING] Aptamer-gene pairs: {df.height:,}")
    print(f"[TRACKING] Unique aptamers: {df['GCST'].n_unique():,}")
    print(f"[TRACKING] Unique autosomal targets: {df['Gene'].n_unique():,}")
    print(f"[TRACKING] Output directory: {out_dir}")
    print(f"[TRACKING] Temporary directory: {scratch_dir}")
    print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")

    # one download per GCST
    gcst_groups = df.partition_by("GCST", maintain_order=True)
    for i, gcst_df in enumerate(gcst_groups, start=1):
        gcst = gcst_df["GCST"][0]
        url = gcst_df["url"][0]
        parquet_files = [out_dir / f"{id_target}.parquet" for id_target in gcst_df["ID"].to_list()]
        print(f"\n[{i:,}/{len(gcst_groups):,}] {gcst}")
        # skip if all gene cis-regions already exist
        if all(file.exists() for file in parquet_files):
            print(f"[SKIPPING] Already processed: {gcst}")
            continue
        print("[TRACKING] Downloading...")
        file = cmd(url, scratch_dir, gcst)

        # extract each mapped gene from the same aptamer file
        for row in gcst_df.iter_rows(named=True):
            gene = row["Gene"]
            chrom = row["CHR"]
            start = row["START"]
            end = row["END"]
            id_target = row["ID"]
            parquet_file = out_dir / f"{id_target}.parquet"

            if parquet_file.exists():
                print(f"[SKIPPING] Already exists: {parquet_file.name}")
                continue

            print("[TRACKING] Extracting cis-region...")
            extract_cis_region(
                file=file,
                gene=gene,
                chrom=chrom,
                start=start,
                end=end,
                window=1_000_000,
                out_dir=out_dir,
                id_target=id_target
            )

        print("[TRACKING] Removing raw file...")
        Path(file).unlink()

    print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
    print("[DONE] WU CSF preprocessing completed.")
    print(f"[DONE] Processed {df['GCST'].n_unique():,} aptamer files.")
    print(f"[DONE] Generated {df.height:,} aptamer-gene cis-regions.")
    print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")

if __name__ == "__main__":
    csf_preprocesing_pipeline()