#!/usr/bin/env python3
import polars as pl
import argparse
from pathlib import Path
from urllib.parse import urlparse, parse_qs

# SwissProt Accessions for proteins

# PULL ALL FILE NAMES WITHIN pQTL DECODE PORTAL
# ONLY IN FORMAT (id1_id2_gene....)# EXTRACT GENE NAME
# BUT THIS SCRIPT IS JUST TO GRAB GENE NAMES AND OVERLAP THEM WITH NCBI REFERENCE
# % MISSINGNESS

# BAGE3 - Chr 21: 7,000,001 - 10,900,000
# ERVV1 - Chr 19: 53,013,921-53,016,123 
# GAGE2D - Chr X : 47,600,001–50,100,000

def parse_decode_gene(filename: str):
    filename = filename.replace(".txt.gz", "")
    gene = filename.split("_")[2].strip().upper()

    # manual fixes for awkward deCODE filenames
    if "ERVV_1_" in filename:
        gene = "ERVV1"
    elif "HLA_G_" in filename:
        gene = "HLA-G"
    elif "HLA_DRB3_" in filename:
        gene = "HLA-DRB3"
    elif "HLA_DQA2_" in filename:
        gene = "HLA-DQA2"
    elif "KRTAP2_4_" in filename:
        gene = "KRTAP2-4"

    return gene

def parse_aptamer_id(filename: str):
    filename = filename.replace(".txt.gz", "")
    parts = filename.split("_")
    aptamer_id = "_".join(parts[:2])
    return aptamer_id

def examine_decode_missingness():
    ncbi_hg38 = "../dat/NCBI/NCBI_genes_grch38_with_synonyms.tsv"
    urls_file = "../dat/deCODE/decode_eur_primary_somascan_urls.txt"
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
            "Synonyms": pl.Utf8,
        },
    )

    ncbi_genes = set(ncbi["Symbol"].drop_nulls().str.to_uppercase().unique().to_list())
    ncbi_synonyms = set()
    for syn in ncbi["Synonyms"].drop_nulls().to_list():
        for x in syn.split(","):
            x = x.strip().upper()
            if x:
                ncbi_synonyms.add(x)

    decode_genes = set()
    rows = []
    with open(urls_file) as f:
        for line in f:
            url = line.strip()
            filename = parse_qs(urlparse(url).query)["file"][0]
            gene = parse_decode_gene(filename)
            aptamer_id = parse_aptamer_id(filename)
            protein_id = f"{gene}_{aptamer_id}"

            decode_genes.add(gene)

            rows.append({
                "protein_id": protein_id,
                "gene_symbol": gene,
                "aptamer_id": aptamer_id,
                "filename": filename,
                "download_url": url,
            })

    matched_by_symbol = decode_genes & ncbi_genes
    missing_after_symbol = decode_genes - ncbi_genes
    rescued_by_synonym = missing_after_symbol & ncbi_synonyms
    decode_missing_in_ncbi = missing_after_symbol - ncbi_synonyms

    print(f"NCBI genes: {len(ncbi_genes)}")
    print(f"NCBI synonyms: {len(ncbi_synonyms)}")
    print(f"deCODE genes: {len(decode_genes)}")
    print()
    print(f"deCODE genes matched by NCBI Symbol: {len(matched_by_symbol)}")
    print(f"deCODE genes rescued by NCBI Synonyms: {len(rescued_by_synonym)}")
    print(f"deCODE genes missing in NCBI after Symbol + Synonyms: {len(decode_missing_in_ncbi)}")
    print()
    print(f"Symbol match %: {len(matched_by_symbol) / len(decode_genes) * 100:.2f}")
    print(f"Symbol + synonym match %: {(len(matched_by_symbol) + len(rescued_by_synonym)) / len(decode_genes) * 100:.2f}")
    print()
    print(sorted(decode_missing_in_ncbi))

    # assemble deCODE target manifest.tsv (deCODE_eur_pgwas_manifest.csv)
    # - protein_id, gene name, aptamer_id, chr, start, end, download link from decode_eur_primary_somascan_urls.txt
    # - if deCODE gene is an NCBI synonym, keep deCODE gene name but use coordinates from matched NCBI symbol

    manifest = pl.DataFrame(rows)

    symbol_lookup = (
        ncbi
        .with_columns([
            pl.col("Symbol").str.to_uppercase().alias("gene_key"),
            pl.lit("symbol").alias("match_type"),
        ])
        .select([
            "gene_key",
            "match_type",
            "Symbol",
            "Chromosome",
            "Begin",
            "End",
            "Orientation",
            "Gene Type",
            "Protein accession"
        ])
    )

    synonym_lookup = (
        ncbi
        .filter(pl.col("Synonyms").is_not_null())
        .with_columns(pl.col("Synonyms").str.split(",").alias("synonym_list"))
        .explode("synonym_list")
        .with_columns([
            pl.col("synonym_list").str.strip_chars().str.to_uppercase().alias("gene_key"),
            pl.lit("synonym").alias("match_type"),
        ])
        .select([
            "gene_key",
            "match_type",
            "Symbol",
            "Chromosome",
            "Begin",
            "End",
            "Orientation",
            "Gene Type",
            "Protein accession"
        ])
    )

    ncbi_lookup = (
        pl.concat([symbol_lookup, synonym_lookup])
        .filter(pl.col("gene_key").is_not_null() & (pl.col("gene_key") != ""))
        .unique(subset=["gene_key"], keep="first")
    )

    manifest = (
        manifest
        .with_columns(pl.col("gene_symbol").str.to_uppercase().alias("gene_key"))
        .join(ncbi_lookup, on="gene_key", how="left")
        .drop("gene_key")
        .rename({
            "Symbol": "ncbi_symbol",
            "Chromosome": "chr",
            "Begin": "start",
            "End": "end",
            "Orientation": "orientation",
            "Gene Type": "gene_type",
            "Protein accession": "ncbi_protein_accession"
        })
    )

    out_dir = Path("./results/deCODE_manifest")
    out_dir.mkdir(parents=True, exist_ok=True)
    out_file = out_dir / "deCODE_eur_pgwas_manifest.csv"
    manifest.write_csv(out_file)

    print()
    print(f"[DONE] Saved deCODE manifest: {out_file}")
    print(f"[DONE] Manifest rows: {manifest.height}")
    print(f"[DONE] Unique deCODE genes: {manifest['gene_symbol'].n_unique()}")
    print(f"[DONE] Unique deCODE protein/aptamer IDs: {manifest['protein_id'].n_unique()}")
    print(f"[DONE] Missing coordinates: {manifest.filter(pl.col('chr').is_null()).height}")
    print(f"[DONE] Symbol matches in manifest: {manifest.filter(pl.col('match_type') == 'symbol').height}")
    print(f"[DONE] Synonym matches in manifest: {manifest.filter(pl.col('match_type') == 'synonym').height}")


if __name__ == "__main__":
    examine_decode_missingness()


# FOR TARGET PROTEIN WITH >1 APTAMER
# EXTRACT LOCUS COORDS FOR GENE THAT ENCODES PROTEIN Y
# GRAB LEAD SNP FOR BOTH AND COMPARE IT - LOWEST PVAL APTAMER FOR SAME TARGET == KEEP
# OTHER == DISCARD

# WORKFLOW STEPS (deCODE pQTLs -> GRCh38/hg38)
# id1_id2_gene - extract gene name - overlap with NCBI - % overlap
# for file in decode portal
# gunzip file.txt.gz
# extract gene name
# download file from decode portal
# df.remove(df.isna())
# harmonise to LDSC-like colnames
# find gene name (and coordinates - chr:start-end on NCBI ref)
# grab cis-region (+/- 1mb from start and end coords)
# save as .parquet
# remove previous txt.gz file for protein X
# carry on with the next protein