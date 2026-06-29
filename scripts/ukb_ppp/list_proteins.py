#!/usr/bin/env python3
import synapseclient
import polars as pl
import os
from pathlib import Path
import shutil
import argparse
import subprocess
import requests

# syn51364943 for UKB-PPP
synapse_id = "syn51364943"
# syn = synapseclient.Synapse()
# syn.login() # pull from base/syn...
ancestry = "European (discovery)" # arg

# cis-region function
def extract_cis_regions(protein_parquet):
    ncbi_hg38 = "./dat/ref/NCBI/NCBI_genes_grch38.tsv"
    ncbi = pl.read_csv(
        ncbi_hg38,
        separator="\t",
        schema_overrides={
            "Accession": pl.Utf8,
            "Begin": pl.Int64,
            "End": pl.Int64,
            "Chromosome": pl.Utf8,
            "Orientation": pl.Utf8,
            "Name": pl.Utf8,
            "Symbol": pl.Utf8,
            "Gene ID": pl.Utf8,
            "Gene Type": pl.Utf8,
            "Transcripts accession": pl.Utf8,
            "Protein accession": pl.Utf8,
            "Protein length": pl.Utf8,
            "Locus tag": pl.Utf8,
        },
    )
    
    protein_parquet = Path(protein_parquet)
    gene = protein_parquet.name.split("_")[0]
    gene_row = (
        ncbi
        .filter(pl.col("Symbol") == gene)
        .select([
            "Symbol",
            "Chromosome",
            "Begin",
            "End",
            "Orientation",
            "Gene Type",
        ])
    )

    if gene_row.height == 0:
        print(f"[CONCERN] No NCBI coordinates found for {gene}. Skipping cis extraction.")
        return

    gene_row = gene_row.head(1)
    chr_ = str(gene_row["Chromosome"][0])
    start = int(gene_row["Begin"][0])
    end = int(gene_row["End"][0])
    cis_start = max(0, start - 1_000_000)
    cis_end = end + 1_000_000
    print(f"[TRACKING] {gene}: chr{chr_}:{cis_start}-{cis_end}")
    df = pl.read_parquet(protein_parquet)
    df_cis = (
        df
        .filter(
            (pl.col("CHR").cast(pl.Utf8) == chr_) &
            (pl.col("BP") >= cis_start) &
            (pl.col("BP") <= cis_end)
        )
    )

    df_cis.write_parquet(protein_parquet)
    print(f"[DONE] Saved cis parquet: {protein_parquet}")

def download_into_hpc(syn):

    """
    Steps:
    1. Log into SLURM HPC (make sure synapse client file in base dir & synapseclient package update within GRCH Docker container)
    2. Ensure .sif container within ./drugMR == ok
    3. Re-make the loop above but with the ukbppp_dl package
    4. Run a re-adapted version of scripts/ukb_ppp/ukb_ppp.py for 1 single file and that rm -rf the previous one
    5. Map rsIDs
    """

    # Establish out_dir
    out_dir = os.path.expanduser("~/drugMR/dat/pQTL/ukb_ppp/") # all .parquet files here
    os.makedirs(out_dir, exist_ok=True)

    ukb_ppp_snps = Path(os.path.expanduser("~/drugMR/dat/ukbb_ppp_snps"))

    # download files
    for item in syn.getChildren(synapse_id):
        if item["name"] == "UKB-PPP pGWAS summary statistics":
            print(f"Entering: {item['name']}\n")
            for folder in syn.getChildren(item["id"]):
                if folder["name"] == ancestry:
                    print(f"Entering: {folder['name']}\n")
                    rows = []
                    for file in syn.getChildren(folder["id"]):
                        downloaded = syn.get(file["id"], downloadLocation=out_dir)
                        downloaded_path = downloaded.path
                        print(f"Downloaded: {downloaded_path}")

                        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                        #  preprocessing one protein at a time
                        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

                        protein_name = file["name"].replace(".tar", "")
                        parts = protein_name.split("_")
                        gene = parts[0]
                        uniprot = parts[1]
                        protein_id = f"{gene}_{uniprot}"
                        protein_dir = Path(out_dir) / protein_name
                        inner_dir = protein_dir / protein_name
                        out_file = Path(out_dir) / f"{protein_id}.parquet"

                        # STEPS
                        # 1. tar -xvf {file}
                        # 2. cd into protein dir
                        # 3. gunzip all chr files
                        cmd_proteins = f"""
set -euo pipefail
mkdir -p "{protein_dir}"
tar -xvf "{downloaded_path}" -C "{protein_dir}"
cd "{inner_dir}"
for gz_file in *.gz; do
    [ -e "$gz_file" ] || continue
    gunzip -f "$gz_file"
done
"""
                        subprocess.run(cmd_proteins, shell=True, check=True, executable="/bin/bash")

                        # merge into one single file across all chromosomes
                        chr_dfs = []
                        total_before_map = 0
                        total_mapped = 0
                        chr_files = sorted(inner_dir.glob("discovery_chr*"))

                        for f in chr_files:
                            chr_name = f.name.split("_")[1]
                            snp_map_file = list(ukb_ppp_snps.glob(f"olink_rsid_map_*_{chr_name}_patched_v2.tsv*"))
                            if len(snp_map_file) == 0:
                                print(f"[CONCERN] No SNP map found for {chr_name}. Skipping {f.name}")
                                continue

                            snp_map_file = snp_map_file[0]
                            print(f"[TRACKING] Mapping {protein_name} {chr_name} using {snp_map_file.name}")
                            df_chr = pl.read_csv(f, separator=" ", has_header=True)
                            rsid_map = (
                                pl.read_csv(snp_map_file, separator="\t", has_header=True)
                                .select(["ID", "rsid"])
                            )

                            n_before = df_chr.height
                            df_chr = df_chr.join(rsid_map, on="ID", how="left")
                            n_mapped = df_chr.filter(pl.col("rsid").is_not_null()).height
                            total_before_map += n_before
                            total_mapped += n_mapped
                            df_chr = df_chr.with_columns(
                                pl.when(pl.col("rsid").is_not_null())
                                .then(pl.col("rsid"))
                                .otherwise(pl.col("ID"))
                                .alias("SNP")
                            )

                            chr_dfs.append(df_chr)

                        if len(chr_dfs) == 0:
                            print(f"[CONCERN] No chromosome files mapped for {protein_name}. Skipping protein.")
                        else:
                            df = pl.concat(chr_dfs)

                            pct_mapped = (total_mapped / total_before_map) * 100 if total_before_map > 0 else 0

                            print(f"[TRACKING] {protein_name}: {total_mapped:,} / {total_before_map:,} SNPs mapped to rsID ({pct_mapped:.2f}%)")
                            print(df.head())

                            # rename cols to Bayesian COLOC / MR format
                            # KEEP ALL EXCEPT "EXTRA"
                            df = (
                                df
                                .drop(["EXTRA", "rsid"])
                                .rename({
                                    "CHROM": "CHR",
                                    "GENPOS": "BP",
                                    "ALLELE1": "A1",
                                    "ALLELE0": "A2",
                                    "A1FREQ": "FRQ",
                                })
                                .with_columns(
                                    P = 10 ** (-pl.col("LOG10P"))
                                )
                                .drop("LOG10P")
                            )

                            n_before = df.height
                            print(f"TOTAL SNPs before removing missing IDs = {n_before:,}")
                            df = df.drop_nulls(subset=["SNP"])
                            n_after = df.height
                            print(f"TOTAL SNPs after removing missing IDs = {n_after:,}")
                            print(f"SNPs removed = {n_before - n_after:,}")
                            df = df.with_columns(pl.lit(gene).alias("GENE"), pl.lit(uniprot).alias("UNIPROT"), pl.lit(protein_id).alias("PROTEIN_ID"))
                            df.write_parquet(out_file)
                            print(f"[DONE] Saved parquet: {out_file}")
                            extract_cis_regions(out_file)

                        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                        # cleanup: remove preliminary files
                        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

                        if Path(downloaded_path).exists():
                            os.remove(downloaded_path)
                            print(f"Deleted tar: {downloaded_path}")

                        if protein_dir.exists():
                            shutil.rmtree(protein_dir, ignore_errors=True)
                            print(f"Deleted extracted dir: {protein_dir}")

                       # continue

                    break
            break


# big to do's
# Rebuild singularity in SLURM
# Fix save -> out_dir (onto the cloud - not pull locally)

# Rebuild .sif file

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--synapse-username", required=True, type=str)
    parser.add_argument("--synapse-token", required=True, type=str)
    args = parser.parse_args()
    syn = synapseclient.Synapse()
    syn.login(email=args.synapse_username, authToken=args.synapse_token)
    download_into_hpc(syn)

if __name__ == "__main__":
    main()