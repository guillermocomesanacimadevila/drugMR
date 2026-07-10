#!/usr/bin/env python3
import synapseclient
import polars as pl
import os
from pathlib import Path
import shutil
import argparse
import subprocess
import requests
import time




# ONLY FOR MISSING TARGETS ON THE ORIGINAL PREPROCESSING PIPELINE 



# syn51364943 for UKB-PPP
synapse_id = "syn51364943"
# syn = synapseclient.Synapse()
# syn.login() # pull from base/syn...
ancestry = "European (discovery)" # arg

# cis-region function
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
    gene = protein_parquet.name.split("_")[0] # gene ID
    print(f"[STEP] Looking up NCBI coordinates for {gene}")

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
    cis_start = max(0, start - 1_000_000) # +/- 1Mb
    cis_end = end + 1_000_000 # +/- 1Mb
    print(f"[TRACKING] {gene}: chr{chr_}:{cis_start}-{cis_end}")

    df = pl.read_parquet(protein_parquet)
    n_before_cis = df.height

    df_cis = (
        df
        .filter(
            (pl.col("CHR").cast(pl.Utf8) == chr_) &
            (pl.col("BP") >= cis_start) &
            (pl.col("BP") <= cis_end)
        )
    )

    n_after_cis = df_cis.height
    print(f"[TRACKING] {gene}: cis SNPs kept = {n_after_cis:,} / {n_before_cis:,}")
    df_cis.write_parquet(protein_parquet)
    print(f"[DONE] Saved cis parquet: {protein_parquet}")

def download_into_hpc(syn, args):

    """
    Steps:
    1. Log into SLURM HPC (make sure synapse client file in base dir & synapseclient package update within GRCH Docker container)
    2. Ensure .sif container within ./drugMR == ok
    3. Re-make the loop above but with the ukbppp_dl package
    4. Run a re-adapted version of scripts/ukb_ppp/ukb_ppp.py for 1 single file and that rm -rf the previous one
    5. Map rsIDs
    """

    pipeline_start = time.time()

    # Establish out_dir
    out_dir = os.path.expanduser("~/drugMR/dat/pQTL/ukb_ppp/missing/") # all .parquet files here
    os.makedirs(out_dir, exist_ok=True)
    ukb_ppp_snps = Path(os.path.expanduser("~/drugMR/dat/ukbb_ppp_snps"))

    print("\n" + "~"*80)
    print("[START] UKB-PPP pQTL processing")
    print(f"[START] Output dir: {out_dir}")
    print(f"[START] SNP map dir: {ukb_ppp_snps}")
    print(f"[START] Missing targets file: {args.gene_list}")
    print("~"*80 + "\n")

    # download files
    for item in syn.getChildren(synapse_id):
        if item["name"] == "UKB-PPP pGWAS summary statistics":
            print(f"Entering: {item['name']}\n")
            for folder in syn.getChildren(item["id"]):
                if folder["name"] == ancestry:
                    print(f"Entering: {folder['name']}\n")
                    protein_files_all = list(syn.getChildren(folder["id"]))
                    total_available = len(protein_files_all)

                    wanted_genes = {
                        line.strip()
                        for line in Path(args.gene_list).read_text().splitlines()
                        if line.strip()
                    }

                    protein_files = [
                        f for f in protein_files_all
                        if f["name"].split("_")[0] in wanted_genes
                    ]

                    found_genes = {f["name"].split("_")[0] for f in protein_files}
                    not_found = sorted(wanted_genes - found_genes)

                    # log info 
                    print(f"[START] Total available proteins: {total_available:,}")
                    print(f"[START] Requested missing targets: {len(wanted_genes):,}")
                    print(f"[START] Matching Synapse files: {len(protein_files):,}")
                    print(f"[START] Targets not found in Synapse: {not_found}")

                    rows = []
                    for local_idx, file in enumerate(protein_files, start=1):
                        protein_global_idx = local_idx
                        protein_start = time.time()

                        print("\n" + "="*80)
                        print(f"[PROTEIN {protein_global_idx:,}/{len(protein_files):,}] {file['name']}")
                        print(f"[BATCH PROGRESS] {local_idx:,}/{len(protein_files):,} in this batch")
                        print("="*80)

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

                        if out_file.exists():
                            print(f"[SKIP] {out_file} already exists.")
                            elapsed = time.time() - protein_start
                            print(f"[DONE] Protein {protein_global_idx:,}/{len(protein_files):,} skipped in {elapsed/60:.2f} min")
                            continue

                        print(f"[STEP 1/7] Downloading {protein_name}")
                        downloaded = syn.get(file["id"], downloadLocation=out_dir)
                        downloaded_path = downloaded.path
                        print(f"Downloaded: {downloaded_path}")

                        # STEPS
                        # 1. tar -xvf {file}
                        # 2. cd into protein dir
                        # 3. gunzip all chr files
                        print(f"[STEP 2/7] Extracting tar and gunzipping chromosome files")

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
                        print(f"[STEP 3/7] Mapping chromosome files to rsIDs")

                        chr_dfs = []
                        total_before_map = 0
                        total_mapped = 0
                        chr_files = sorted(inner_dir.glob("discovery_chr*"))

                        print(f"[TRACKING] Found {len(chr_files):,} chromosome files for {protein_name}")

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
                            print(f"[TRACKING] {chr_name}: {n_mapped:,} / {n_before:,} SNPs mapped")

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
                            print(f"[STEP 4/7] Merging chromosomes")

                            df = pl.concat(chr_dfs)
                            pct_mapped = (total_mapped / total_before_map) * 100 if total_before_map > 0 else 0
                            print(f"[TRACKING] {protein_name}: {total_mapped:,} / {total_before_map:,} SNPs mapped to rsID ({pct_mapped:.2f}%)")
                            print(df.head())
                            print(f"[STEP 5/7] Renaming columns and calculating P values")

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

                            df = df.with_columns(
                                pl.lit(gene).alias("GENE"),
                                pl.lit(uniprot).alias("UNIPROT"),
                                pl.lit(protein_id).alias("PROTEIN_ID")
                            )

                            df.write_parquet(out_file)
                            print(f"[DONE] Saved full parquet: {out_file}")

                            print(f"[STEP 6/7] Extracting cis-region")
                            extract_cis_regions(out_file)

                        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                        # cleanup: remove preliminary files
                        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

                        print(f"[STEP 7/7] Cleaning temporary files")

                        if Path(downloaded_path).exists():
                            os.remove(downloaded_path)
                            print(f"Deleted tar: {downloaded_path}")

                        if protein_dir.exists():
                            shutil.rmtree(protein_dir, ignore_errors=True)
                            print(f"Deleted extracted dir: {protein_dir}")

                        elapsed = time.time() - protein_start
                        print(f"[DONE] Protein {protein_global_idx:,}/{len(protein_files):,} completed in {elapsed/60:.2f} min")

                       # continue

                    break
            break

    total_elapsed = time.time() - pipeline_start
    print("\n" + "="*80)
    print("[FINISHED] UKB-PPP batch finished")
    print(f"[FINISHED] Runtime: {total_elapsed/3600:.2f} hours")
    print("="*80 + "\n")


# big to do's
# Rebuild singularity in SLURM
# Fix save -> out_dir (onto the cloud - not pull locally)

# Rebuild .sif file

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--synapse-username", required=True, type=str)
    parser.add_argument("--synapse-token", required=True, type=str)
    parser.add_argument("--gene-list", required=True, type=str)
    args = parser.parse_args()
    syn = synapseclient.Synapse()
    syn.login(email=args.synapse_username, authToken=args.synapse_token)
    download_into_hpc(syn, args)

if __name__ == "__main__":
    main()