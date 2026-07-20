#!/usr/bin/env python3
import polars as pl 
from pathlib import Path


# FLJ44635 - NHSL2 - X
# GAGE2B  

# Start: 47,600,001 bp
# End: 50,100,000 bp
# CHR X

# rev

# remove rev 
# map coords of the other two genes in chrX
# remove all chrX - only keep autosomal targets 


# candidate_genes
# candidate_uniprots

def map_target_coords():
    # current manifest (without coords)
    manifest = "./results/WS_CSF_manifest/wu_csf_gcst_safe_mapping.csv"
    manifest = pl.read_csv(manifest)
    print(f"Shape of manifest df: {manifest.shape}")
    # candidate_genes
    # candidate_uniprots
    targets = manifest["EntrezGeneSymbol"].to_list()
    print(f"Unique targets (before initial preprocessing): {len(manifest['EntrezGeneSymbol'].unique()):,}")

    # ncbi stuff
    ncbi_hg38 = "../dat/NCBI/NCBI_genes_grch38_with_synonyms.tsv"
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

    # check overlap
    # ncbi genes + synonims
    ncbi_genes = set(ncbi["Symbol"].drop_nulls().str.to_uppercase().unique().to_list())
    ncbi_synonyms = set()
    for syn in ncbi["Synonyms"].drop_nulls().to_list():
        for x in syn.split(","):
            x = x.strip().upper()
            if x:
                ncbi_synonyms.add(x)

    # manual coords for genes not found in NCBI
    manual_coords = {
        "FLJ44635": {
            "Chromosome": "X",
            "Begin": 71910845,
            "End": 72161750,
            "Orientation": None,
        },
        "GAGE2B": {
            "Chromosome": "X",
            "Begin": 49331615,
            "End": 49338933,
            "Orientation": None,
        },
    }

    # clean EntrezGeneSymbol
    manifest = manifest.with_columns(pl.col("EntrezGeneSymbol").str.to_uppercase().str.replace_all(r"\s*\|\s*", "|").alias("EntrezGeneSymbol"))

    # remove missing / non-human targets
    manifest = manifest.filter(pl.col("EntrezGeneSymbol").is_not_null() & ~pl.col("EntrezGeneSymbol").is_in(["NONE", "REV"]))

    # update targets after cleaning
    targets = manifest["EntrezGeneSymbol"].to_list()

    count = 0
    for i in targets:
        genes = [x.strip().upper() for x in i.split("|") if x.strip()]
        matched = any(gene in ncbi_genes or gene in ncbi_synonyms or gene in manual_coords for gene in genes)
        if matched:
            count += 1
        else:
            print(i)

    print(f"Percentage overlap = {(count / len(targets)) * 100:.2f}%")
    print(f"Unique targets (after initial preprocessing): {len(manifest['EntrezGeneSymbol'].unique()):,}")

    # lookup table
    lookup = {}
    for row in ncbi.iter_rows(named=True):
        if row["Symbol"] is None:
            continue
        info = {"CHR": row["Chromosome"], "START": row["Begin"], "END": row["End"], "ORIENTATION": row["Orientation"]}
        lookup[row["Symbol"].upper()] = info
        if row["Synonyms"] is not None:
            for syn in row["Synonyms"].split(","):
                syn = syn.strip().upper()
                if syn:
                    lookup[syn] = info

    # manual coords
    for gene, info in manual_coords.items():
        lookup[gene] = {"CHR": info["Chromosome"], "START": info["Begin"], "END": info["End"], "ORIENTATION": info["Orientation"]}

    # one name per file for downloads
    manifest = manifest.with_columns(pl.col("EntrezGeneSymbol").str.split("|").list.first().alias("Gene"))

    # save download manifest (one row per GCST)
    manifest.write_csv("./results/WS_CSF_manifest/wu_csf_download_manifest.csv")
    print(f"Download rows: {len(manifest):,}")
    print(f"Unique download files: {manifest['GCST'].n_unique():,}")

    # map coords (one row per gene)
    rows = []
    for row in manifest.iter_rows(named=True):
        genes = [x.strip().upper() for x in row["EntrezGeneSymbol"].split("|")]
        for gene in genes:
            if gene not in lookup:
                print(gene)
                continue
            rows.append({
                **row,
                "Gene": gene,
                "CHR": lookup[gene]["CHR"],
                "START": lookup[gene]["START"],
                "END": lookup[gene]["END"],
                "ORIENTATION": lookup[gene]["ORIENTATION"],
            })

    coords_manifest = pl.DataFrame(rows)
    print(f"Mapped rows: {len(coords_manifest):,}")
    print(f"Unique mapped targets: {coords_manifest['Gene'].n_unique():,}")
    coords_manifest.write_csv("./results/WS_CSF_manifest/wu_csf_manifest_with_coords.csv")
    # autosomal only
    autosomal = coords_manifest.filter(~pl.col("CHR").is_in(["X", "Y"]))
    print(f"Autosomal rows: {len(autosomal):,}")
    print(f"Unique autosomal targets: {autosomal['Gene'].n_unique():,}")
    autosomal.write_csv("./results/WS_CSF_manifest/wu_csf_manifest_autosomal.csv")

if __name__ == "__main__":
    map_target_coords()