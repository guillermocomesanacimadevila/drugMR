#!/usr/bin/env python3
import argparse
import os
from pathlib import Path

import polars as pl

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

    # top_snp per protein: COLOC's own recorded value is used when available - it's
    # actually the pQTL's own lead (lowest-P) SNP for that cis-region (bin/coloc.R:
    # `top_snp <- protein$SNP[1]`, computed from the pQTL data alone, before COLOC
    # itself even runs), not something standard COLOC "produces". A protein that only
    # colocalises under PWCoCo (see project_pwcoco_wiring memory) has no coloc_out row
    # at all and would otherwise silently disappear from target_stats_out() / the
    # dashboard's harmonised SNP-allele-beta columns - so its lead SNP is recomputed
    # here directly from pqtl.parquet, the exact same "sort by P, take the top row"
    # logic coloc.R itself uses, rather than depending on either method's own output.
    target_top_snp = {
        row["protein"]: str(row["top_snp"])
        for row in coloc_res.iter_rows(named=True)
    }

    pwcoco_file = paths.pwcoco_out(pqtl_dataset, pheno_id, local_results_dir)
    if pwcoco_file.exists():
        pwcoco_res = pl.read_csv(pwcoco_file, separator="\t")
        pwcoco_only_proteins = set(pwcoco_res["protein"].unique().to_list()) - set(target_top_snp.keys())
        for pwcoco_target in pwcoco_only_proteins:
            pwcoco_pqtl = pl.read_parquet(Path(f"./dat/cis_regions/{pqtl_dataset}/{pwcoco_target}") / "pqtl.parquet")
            target_top_snp[pwcoco_target] = str(pwcoco_pqtl.sort("P").row(0, named=True)["SNP"])

    top_hits = []
    for target, snp_id in target_top_snp.items():
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
            "outcome_trait": pheno_id,
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