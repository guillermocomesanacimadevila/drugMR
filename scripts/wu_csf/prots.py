#!/usr/bin/env python3
import polars as pl
from pathlib import Path
import json
import requests
from time import sleep
from rapidfuzz import process, fuzz

# ETL pipeline workflow plan 
# fastexcel -> onto image
# rapidfuzz -> onto image 
# Download .xlxs from Dropbox dir/
# Extract all GC IDs from EBI GWAS catalogue alongside full names for GWAS phenotype
# Match with .xlsx to check overlap - rename the ones in EBI -> acronym
# Map each apptainer to its corresponding gene ID + Uniprot
# Preserve (ONLY cis-pQTLs)
# Map each gene to NCBI gene manifest
# Extract wg_csf manifest file with gene names, coords to hg38, etc... -> 100% overlap
# Download each file from the summary stats dir/
# Map apptainer
# ETL pipeline to -> .parquet (cis-region + biallelic + -INDELs) etc...
# Store as geneID_uniprotID.parquet (unlink() raw file) -> do the same for all proteins (cis-only)

def get_json_file():
    out_dir = "./results/WS_CSF_manifest"
    out_dir = Path(out_dir)
    start = 90421033
    end = 90428040
    proteins = {}
    for i in range(start, end + 1):
        accession = f"GCST{i}"
        r = requests.get(f"https://www.ebi.ac.uk/gwas/rest/api/studies/{accession}")
        if not r.ok:
            print(f"{accession}: FAILED ({r.status_code})")
            continue
        protein = (r.json().get("diseaseTrait", {}).get("trait", "").removesuffix(" levels"))
        lower = ((i - 1) // 1000) * 1000 + 1
        upper = lower + 999
        url = (
            "https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/"
            f"GCST{lower:08d}-GCST{upper:08d}/"
            f"{accession}/harmonised/"
            f"{accession}.h.tsv.gz"
        )
        proteins[accession] = {"protein": protein, "url": url}
        print(f"{accession} -> {protein}\n{url}\n")
        sleep(0.05)
    json_file = out_dir / "wu_csf_proteins.json"
    with open(json_file, "w") as f:
        json.dump(proteins, f, indent=4)
    print(f"\nSaved to {json_file}")

if __name__ == "__main__":
    get_json_file()
