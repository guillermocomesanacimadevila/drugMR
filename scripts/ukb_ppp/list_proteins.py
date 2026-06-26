#!/usr/bin/env python3
import synapseclient
import polars as pl 
import os
import subprocess

# syn51364943 for UKB-PPP
synapse_id = "syn51364943"
syn = synapseclient.Synapse()
syn.login() # pull from base/syn...
ancestry = "European (discovery)"

for item in syn.getChildren(synapse_id):
    if item["name"] == "UKB-PPP pGWAS summary statistics":
        print(f"Entering: {item['name']}\n")
        for folder in syn.getChildren(item["id"]):
            if folder["name"] == ancestry:
                print(f"Entering: {folder['name']}\n")
                rows = []
                for file in syn.getChildren(folder["id"]):
                    name = file["name"]
                    parts = name.replace(".tar", "").split("_")
                    rows.append({
                        "synapse_id": file["id"],
                        "filename": name,
                        "gene_symbol": parts[0],
                        "uniprot_id": parts[1],
                        "olink_id": parts[2],
                        "version": parts[3],
                        "panel": "_".join(parts[4:]),
                    })
                df = pl.DataFrame(rows)
                out_dir = "../results/UKB-PPP_manifest"
                os.makedirs(out_dir)
                df.write_csv(f"{out_dir}/ukbppp_eur_pgwas_manifest.csv")
                print(df.head())
                print(f"\nTotal files: {df.height}")
                print(f"Total unique proteins = {len(df['uniprot_id'].unique()):,}")
                print(f"Total unique genes = {len(df['gene_symbol'].unique()):,}")
                break
        break

def download_into_hpc():
    return 