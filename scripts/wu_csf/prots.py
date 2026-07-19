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


def assemble_manifest():
    apptamer_specs = "./results/WS_CSF_manifest/aptamer_info.xlsx"
    df = pl.read_excel(apptamer_specs)
    out_dir = Path("./results/WS_CSF_manifest")
    proteins = json.load(open(out_dir / "wu_csf_proteins.json"))
    protein_names = {accession: info["protein"] for accession, info in proteins.items()}
    gcst = []
    urls = []
    for row in df.iter_rows(named=True):
        name = row["TargetFullName"]
        result = process.extractOne(name, protein_names, scorer=fuzz.WRatio)
        if result is None:
            print(f"{name} -> NO MATCH")
            gcst.append(None)
            urls.append(None)
            continue
        protein, score, accession = result
        url = proteins[accession]["url"]
        print(f"{name} -> {accession} -> {protein} ({score:.1f})")
        print(url)
        # slap accession ID onto excel
        if score == 100:
            gcst.append(accession)
            urls.append(url)
        else:
            gcst.append(None)
            urls.append(None)
    df = df.with_columns(pl.Series("GCST", gcst), pl.Series("url", urls))
    df.write_excel(out_dir / "aptamer_info_gcst.xlsx")


def map_target_info_from_ncbi():
    # EntrezGeneSymbol
    apptamer_specs = "./results/WS_CSF_manifest/aptamer_info_gcst.xlsx"
    df = pl.read_excel(apptamer_specs)
    targets = df["EntrezGeneSymbol"].to_list()
    score = 0
    missing = []

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

    non_human = {
        "GFP",
        "REV",
        "NODH",
        "CYSH",
        "MAGAININS",
        "APCB",
        "APCA",
    }

    symbol_lookup = {}
    synonym_lookup = {}
    for g in ncbi.iter_rows(named=True):
        if g["Symbol"] is not None:
            symbol_lookup[g["Symbol"].upper()] = g
        if g["Synonyms"] is not None:
            for syn in g["Synonyms"].split(","):
                syn = syn.strip().upper()
                if syn and syn not in synonym_lookup:
                    synonym_lookup[syn] = g

    # check how many present in NCBI
    for t in targets:
        if t is None or t.upper() == "NONE":
            missing.append(t)
            continue
        genes = [x.strip().upper() for x in t.split("|")]
        if all(gene in symbol_lookup or gene in synonym_lookup for gene in genes):
            score += 1
        else:
            missing.append(t)

    # print score
    print(f"Found {score:,}/{len(targets):,} ({(score / len(targets)) * 100:.2f}%) targets in NCBI")
    if len(missing) > 0:
        print("\nMissing targets:")
        for t in missing:
            print(t)

    # NOW MAP THE CORRESPONDING COORDS AND GENERATE MANIFEST
    gene_key = []
    match_type = []
    symbol = []
    chromosome = []
    begin = []
    end = []
    orientation = []
    gene_type = []
    protein_accession = []
    for row in df.iter_rows(named=True):
        gene = row["EntrezGeneSymbol"]
        if gene is None or gene.upper() == "NONE":
            gene_key.append(None)
            match_type.append("missing")
            symbol.append(None)
            chromosome.append(None)
            begin.append(None)
            end.append(None)
            orientation.append(None)
            gene_type.append(None)
            protein_accession.append(None)
            continue
        genes = [x.strip().upper() for x in gene.split("|")]
        if all(gene_part in non_human for gene_part in genes):
            gene_key.append(gene)
            match_type.append("non_human")
            symbol.append(None)
            chromosome.append(None)
            begin.append(None)
            end.append(None)
            orientation.append(None)
            gene_type.append(None)
            protein_accession.append(None)
            print(f"{gene} -> NON-HUMAN/CONTROL")
            continue
        matches = []
        for gene_part in genes:
            if gene_part in symbol_lookup:
                matches.append(("symbol", symbol_lookup[gene_part]))
            elif gene_part in synonym_lookup:
                matches.append(("synonym", synonym_lookup[gene_part]))
        if len(matches) == len(genes):
            gene_key.append("|".join(genes))
            match_type.append("|".join([x[0] for x in matches]))
            symbol.append("|".join([str(x[1]["Symbol"]) for x in matches]))
            chromosome.append("|".join([str(x[1]["Chromosome"]) for x in matches]))
            begin.append("|".join([str(x[1]["Begin"]) for x in matches]))
            end.append("|".join([str(x[1]["End"]) for x in matches]))
            orientation.append("|".join([str(x[1]["Orientation"]) for x in matches]))
            gene_type.append("|".join([str(x[1]["Gene Type"]) for x in matches]))
            protein_accession.append("|".join([str(x[1]["Protein accession"]) for x in matches]))
            print(f"{gene} -> {'|'.join([str(x[1]['Symbol']) for x in matches])}")
        else:
            gene_key.append(gene)
            match_type.append("unresolved")
            symbol.append(None)
            chromosome.append(None)
            begin.append(None)
            end.append(None)
            orientation.append(None)
            gene_type.append(None)
            protein_accession.append(None)
            print(f"{gene} -> NO MATCH")

    df = df.with_columns(
        pl.Series("gene_key", gene_key),
        pl.Series("match_type", match_type),
        pl.Series("NCBI_Symbol", symbol),
        pl.Series("Chromosome", chromosome),
        pl.Series("Begin", begin),
        pl.Series("End", end),
        pl.Series("Orientation", orientation),
        pl.Series("Gene_Type", gene_type),
        pl.Series("Protein_accession", protein_accession),
    )

    out_dir = Path("./results/WS_CSF_manifest")
    manifest_file = out_dir / "wu_csf_manifest.csv"
    df.write_csv(manifest_file)
    print(f"\nSaved to {manifest_file}")

    # for each target which == prsent in NCBI and/or synonim
    # store in temp variables:
    # "gene_key",
    # "match_type",
    # "Symbol",
    # "Chromosome",
    # "Begin",
    # "End",
    # "Orientation",
    # "Gene Type",
    # "Protein accession"
    # append to excel file df
    # THEN IN ANOTHER SCRIPT
    # establish path -> path("https://ftp.ebi.ac.uk/.../harmonised/GCST90428040.h.tsv.gz")

if __name__ == "__main__":
    # get_json_file()
    # assemble_manifest()
    map_target_info_from_ncbi()