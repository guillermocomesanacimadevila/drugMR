#!/usr/bin/env python3
import polars as pl
import argparse
import os
from pathlib import Path
from drugmr import paths

# TO DO'S
# Coloc outdir for pqtl dataset X
# Map the original ingo fro mthe original pQTL -> define risk allele 
# Go to original GWAS + pQTL dir/ 
# map Ps - Betas - A1 - A2
# slap onto results/info 

def compile_top_cis_hits(pheno_id: str, pqtl_dataset: str, local_results_dir: str = "results"):
    coloc_res = pl.read_csv(paths.coloc_out(pqtl_dataset, pheno_id, local_results_dir), separator="\t")
    if "protein_id" in coloc_res.columns:
        coloc_res = coloc_res.rename({"protein_id": "protein"})
    top_hits = []
    for row in coloc_res.iter_rows(named=True):
        target = row["protein"]
        snp_id = str(row["top_snp"]) ########## check this for the exact colname
        # dir for the gwas and pQTL
        cis_r = Path(f"./dat/cis_regions/{pqtl_dataset}/{target}")
        gwas = cis_r / "gwas.parquet"
        pqtl = cis_r / "pqtl.parquet"
        gwas = pl.read_parquet(gwas)
        pqtl = pl.read_parquet(pqtl)
        gwas_snp = gwas.filter(pl.col("SNP").cast(pl.Utf8) == snp_id)
        pqtl_snp = pqtl.filter(pl.col("SNP").cast(pl.Utf8) == snp_id)
        gwas_row = gwas_snp.row(0, named=True)
        pqtl_row = pqtl_snp.row(0, named=True)
        # GWAS values
        a1 = str(gwas_row["A1"]).upper()
        a2 = str(gwas_row["A2"]).upper()
        a1_frq = float(gwas_row["FRQ"])
        gwas_beta = float(gwas_row["BETA"])
        gwas_p = float(gwas_row["P"])
        if gwas_p == 0:
            gwas_p = 1e-300
        print(f"[TRACKING] Assorting SNP directionality table for target {target}")

        # make A1 the GWAS risk allele
        if gwas_beta < 0:
            a1, a2 = a2, a1
            a1_frq = 1 - a1_frq
            gwas_beta = -gwas_beta

        # pQTL values
        pqtl_a1 = str(pqtl_row["A1"]).upper()
        pqtl_a2 = str(pqtl_row["A2"]).upper()
        pqtl_beta = float(pqtl_row["BETA"])
        pqtl_p = float(pqtl_row["P"])
        if pqtl_p == 0:
            pqtl_p = 1e-300

        # align pQTL beta to the GWAS risk allele
        if pqtl_a1 == a1 and pqtl_a2 == a2:
            pass

        elif pqtl_a1 == a2 and pqtl_a2 == a1:
            pqtl_beta = -pqtl_beta

        else:
            print(
                f"{target}: allele mismatch for {snp_id} | "
                f"GWAS {a1}/{a2} | "
                f"pQTL {pqtl_a1}/{pqtl_a2}"
            )
            continue

        top_hits.append({
            "protein": target,
            "SNP": snp_id,
            "A1": a1,                  
            "A2": a2,                  
            "FRQ": a1_frq,             
            "GWAS_BETA": gwas_beta,    
            "GWAS_P": gwas_p,
            "pQTL_BETA": pqtl_beta,    
            "pQTL_P": pqtl_p,
        })

    top_hits = pl.DataFrame(top_hits)
    output = paths.target_stats_out(pqtl_dataset, pheno_id, local_results_dir)
    output.parent.mkdir(parents=True, exist_ok=True)
    top_hits.write_csv(output, separator="\t")
    return top_hits

def main():
    p = argparse.ArgumentParser()
    p.add_argument("--pheno_id", required=True)
    p.add_argument("--pqtl_dataset", required=True)
    p.add_argument("--local_results_dir", default="results")
    args = p.parse_args()
    compile_top_cis_hits(
        pheno_id=args.pheno_id,
        pqtl_dataset=args.pqtl_dataset,
        local_results_dir=args.local_results_dir
    )

if __name__ == "__main__":
    main()