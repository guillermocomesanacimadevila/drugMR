#!/usr/bin/env python3
import polars as pl
import argparse
from pathlib import Path

# PULL ALL FILE NAMES WITHIN pQTL DECODE PORTAL
# ONLY IN FORMAT (id1_id2_gene....)# EXTRACT GENE NAME
# BUT THIS SCRIPT IS JUST TO GRAB GENE NAMES AND OVERLAP THEM WITH NCBI REFERENCE
# % MISSINGNESS

def examine_decode_missingness():
    ncbi_hg38 = "../dat/NCBI/NCBI_genes_grch38.tsv"
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
    ncbi_genes = set(ncbi["Symbol"].unique().to_list())
    decode_genes = set()
    with open("../dat/deCODE/decode_eur_primary_somascan_urls.txt") as f:
        for line in f:
            filename = Path(line.strip()).name
            gene = filename.split("_")[2]
            decode_genes.add(gene)
    decode_missing_in_ncbi = decode_genes - ncbi_genes
    print(f"NCBI genes: {len(ncbi_genes)}")
    print(f"deCODE genes: {len(decode_genes)}")
    print(f"deCODE genes missing in NCBI: {len(decode_missing_in_ncbi)}")
    print(sorted(decode_missing_in_ncbi)[:50])


if __name__ == "__main__":
    examine_decode_missingness()


# FOR TARGET PROTEIN WITH >1 APTAMER
# EXTRACT LOCUS COORDS FOR GENE THAT ENCODES PROTEIN Y
# GRAB LEAD SNP FOR BOTH AND COMPARE IT - LOWEST PVAL APTANER FOR SAME TARGET == KEEP
# OTHER == DISCARD 

# WORKFLOW STEPS (deCODE pQTLs -> GRCh38/hg38)
# id1_id2_gene - extract gene name - overlap with NCBI - % overlap
# for file in decode portal
# gunzip file.txt.gz
# extract gene name
# download file from decode portal 
# df.remove(df.isna())
# harmonise to LDSC-like colnames
# find gene name (and coordinates - chr:start-end on NCBI ref)
# grab cis-region (+/- 1mb from start and end coords)
# save as .parquet
# remove previous txt.gz file for protein X
# carry on with the next protein