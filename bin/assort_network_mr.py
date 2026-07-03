#!/usr/bin/env python3
import argparse
import polars as pl
from pathlib import Path
import subprocess
from drugmr import NetworkMR
import os 

# DS NetworkMR pipeline
# dictionary in jupyter notebook {M_id: 'User/Path/...'}
# FROM NOTEBOOK -> MAKE A MEDIATOR MANIFEST
# For each protein part of dataset X
# Run cis-MR (twice) -> For each X -> M
# Also run M -> Y (whole genome) 
# results/networkMR/ 3 subdirectories
# results/networkMR/M_Y/....csv (Genome-wide - one CSV with MR outputs where 1 entry == univariable MR from a mediator M on Y)
# results/networkMR/X_M/mediator_1/....csv (1 entry == univariable cis-MR - 1 protein vs that mediator)
# results/networkMR/X_M/mediator_2/....csv (1 entry == univariable cis-MR - 1 protein vs that mediator)
# results/networkMR/X_M/mediator_N/....csv (1 entry == univariable cis-MR - 1 protein vs that mediator)
# results/networkMR/mediation_estimates/...csv (massive CSV with a given protein that FDR significant in X->M and X->Y and also if IVW_p < 0.05 in X->Y run NetworkMR package - here the output of NetworkMR package)

# ARGS
# pheno_id   
# pheno_gwas 
# ref_bfile  
# out_dir    
# mediator_dict (or another data type which allows for >1 value) THIS IT FOR MEDIATOR SAMPLE SIZE!!!!!
# pqtl_dataset
# pqtl_dir

# X -> Y == already done
# For each M -> run X -> M
# For each M -> run M -> Y (genome-wide)
# Inherit NMR library from drugmr/ and run NetworkMR

# for each mediatror 

def run_genomewide_mr(ref_bfile: str, pheno_id: str, pheno_gwas: str):
    genomewide_mr = "./bin/genomewide_mr.R"
    mediator_dir = Path("./dat/gwas/mediators")
    out_dir = Path(f"./results/networkMR/M_Y/{pheno_id}")
    out_dir.mkdir(parents=True, exist_ok=True)
    all_results = []
    for file in mediator_dir.glob("*.tsv"):
        M_id = file.stem
        M = pl.read_csv(file, separator="\t", n_rows=1)
        M_N = M["N"][0]

        cmd = f"""
set -euo pipefail
Rscript {genomewide_mr} \\
    '{M_id}' \\
    '{file}' \\
    '{M_N}' \\
    '{pheno_id}' \\
    '{pheno_gwas}' \\
    '{ref_bfile}' \\
    '{out_dir}'
"""

        subprocess.run(cmd, shell=True, check=True, executable="/bin/bash")
        res_file = out_dir / f"{M_id}_{pheno_id}_genomewide_MR.tsv"

        if res_file.exists():
            all_results.append(pl.read_csv(res_file, separator="\t"))

    if all_results:
        compiled = pl.concat(all_results, how="diagonal")
        compiled_file = out_dir / f"{pheno_id}_mediator_genomewide_MR.tsv"
        compiled.write_csv(compiled_file, separator="\t")
        print(f"[TRACKING] Saved compiled genome-wide mediator MR results: {compiled_file}")
        
        for file in out_dir.glob(f"*_{pheno_id}_genomewide_MR.tsv"):
            if file.name != compiled_file.name:
                file.unlink()
    else:
        print("[CONCERN] No genome-wide mediator MR results generated.")


def run_cis_mr_X_M(pqtl_dataset: str, pqtl_dir: str, ref_bfile: str):
    ref_bfile = Path(ref_bfile)
    pqtl_dir = Path(pqtl_dir)
    out_dir = Path(f"./results/networkMR/X_M/{pqtl_dataset}")
    protein_dir = Path(f"./dat/cis_regions/{pqtl_dataset}")
    mediator_gwas = Path("./dat/gwas/mediators")
    out_dir.mkdir(parents=True, exist_ok=True)

    for protein_path in protein_dir.iterdir():
        protein = protein_path.name
        print(f"[TRACKING] Processing protein {protein}")
        pqtl_file = protein_path / "pqtl.parquet"
        gwas_file = protein_path / "gwas.parquet"
        pqtl = pl.read_parquet(pqtl_file)
        gwas = pl.read_parquet(gwas_file)
        chr = pqtl["CHR"][0]
        start = min(pqtl["BP"].to_list())
        end = max(pqtl["BP"].to_list())
        mediator_out = protein_path / "mediators"
        mediator_out.mkdir(parents=True, exist_ok=True)

        for mediator_file in mediator_gwas.iterdir():
            if mediator_file.suffix not in [".tsv", ".txt"]:
                continue
            mediator_id = mediator_file.stem
            print(f"[TRACKING] Processing mediator {mediator_id} for {protein}")
            mediator = pl.read_csv(mediator_file, separator="\t")

            mediator = mediator.filter(
                (pl.col("CHR") == chr) &
                (pl.col("BP") >= start) &
                (pl.col("BP") <= end)
            )

            if mediator.height == 0:
                print(f"[CONCERN] No mediator SNPs in cis-region for {protein} / {mediator_id}")
                continue

            overlap_snps = (
                set(pqtl["SNP"].to_list()) &
                set(gwas["SNP"].to_list()) &
                set(mediator["SNP"].to_list())
            )

            if len(overlap_snps) == 0:
                print(f"[CONCERN] No overlapping SNPs for {protein} / {mediator_id}")
                continue

            mediator_keep = mediator.filter(pl.col("SNP").is_in(overlap_snps))
            mediator_keep.write_parquet(mediator_out / f"{mediator_id}.parquet")

            print(
                f"[TRACKING] Saved mediator {mediator_id}.parquet for {protein} "
                f"with {len(overlap_snps)} overlapping SNPs"
            )
            # ukb_ppp_AD_all_MR.tsv


            # slap this onto a different function
            # ***************
            # ***************
            # ***************

            cmd_cis_mr = f"""
set -euo pipefail
Rscript bin/cis_mr.R \\
    {pqtl_dataset} \\
    dat/cis_regions/{pqtl_dataset} \\
    {mediator_id} \\
    ./results/QC/{mediator_id}/{mediator_id}.tsv \\
    {ref_bfile}  
mkdir -p ./results/networkMR/X_M     
mv ./results/cis-MR/{pqtl_dataset}_{mediator_id}_all_MR.tsv {out_dir}           
"""
            subprocess.run(cmd_cis_mr, shell=True, check=True, executable="/bin/bash")
            
            # ***************
            # ***************
            # ***************
            


def perform_network_mr():

    # for each mediator 
    # pull info from "./results/networkMR/..."
    return