#!/usr/bin/env python3
import argparse
from pathlib import Path
import polars as pl 
import subprocess
import os 

# out_dir = f"./results/SMR/{eqtl_dataset}/{pheno_id}"
# final_targets_file = out_dir / f"{pqtl_dataset}_{pheno_id}_final_multi_omics_targets.tsv"

# ------------------------------------
# ------------------------------------
# PRELIMINARY COLOCALISATION FUNCTIONS 
# ------------------------------------
# ------------------------------------

# additional scripts
# GWAS -> eQTL coloc
# GWAS - pQTL - eQTL coloc

def gwas_eqtl_coloc(pheno_id: str, pqtl_dataset: str, eqtl_dataset: str, n_cases: int, n_controls: int):
    moloc_src = Path("./bin/moloc_qtl.R")
    coloc_src = Path("./bin/coloc_eqtl.R")

    # checks
    if len(pheno_id) == 0 or len(pqtl_dataset) == 0 or len(eqtl_dataset) == 0:
        print("[CONCERN] Yowza! You need to properly define your pheno_id, pqtl_dataset and eqtl_dataset")
        return

    cmd_coloc = f"""
set -euo pipefail
Rscript {coloc_src} \\
  {pheno_id} \\
  {pqtl_dataset} \\
  {eqtl_dataset} \\
  {n_cases} \\
  {n_controls}
"""

    cmd_moloc = f"""
set -euo pipefail
Rscript {moloc_src} \\
  {pheno_id} \\
  {pqtl_dataset} \\
  {eqtl_dataset}
"""

    # Multi-omics coloc + GWAS-eQTL coloc
    subprocess.run(cmd_coloc, shell=True, check=True, executable="/bin/bash")
    subprocess.run(cmd_moloc, shell=True, check=True, executable="/bin/bash")


# need to mkae prints as to % vs origunak loci when it comes to SNP overlap 

def target_multi_omics_moloc(pheno_id: str, eqtl_dataset: str, pqtl_dataset: str):
    res_dir = Path(f"./results/SMR/{eqtl_dataset}/{pheno_id}")
    targets = res_dir / f"{pqtl_dataset}_{pheno_id}_final_multi_omics_targets.tsv"
    pqtl_dir = Path(f"./dat/cis_regions/{pqtl_dataset}")
    eqtl_dir = Path(f"./dat/sc-eQTL/{eqtl_dataset}")
    multi_omics_dir = Path(f"./dat/multi_omics_targets/{pqtl_dataset}")
    t = pl.read_csv(targets, separator="\t")
    print(f"[TRACKING] {t.height} target x cell-type multi-omics results loaded")
    all_prepared_targets = []

    # q_FDR - p_HEIDI
    for row in t.iter_rows(named=True):
        q_smr = row["q_SMR"]
        p_heidi = row["p_HEIDI"]
        cell_type = row["cell_type"]
        protein = row["protein"]

        if q_smr is None or p_heidi is None:
            continue

        if q_smr >= 0.05 or p_heidi <= 0.01:
            continue

        # protein == BLNK_Q8WV28
        # gene == BLNK
        gene = protein.split("_")[0]
        probe_id = row["probeID"].split(".")[0]

        print("------------------------------------------------------------")
        print(f"[TRACKING] Preparing {protein} in {cell_type}")
        print(f"[TRACKING] q_SMR = {q_smr}")
        print(f"[TRACKING] p_HEIDI = {p_heidi}")

        # original pQTL and GWAS cis-region
        target_cis_dir = pqtl_dir / protein
        gwas_file = target_cis_dir / "gwas.parquet"
        pqtl_file = target_cis_dir / "pqtl.parquet"

        # SingleBrain files stored as dat/sc-eQTL/SingleBrain/MG.parquet
        eqtl_file = eqtl_dir / f"{cell_type}.parquet"
        gwas = pl.read_parquet(gwas_file)
        pqtl = pl.read_parquet(pqtl_file)
        eqtl = pl.read_parquet(eqtl_file)
        original_gwas_n = gwas.height
        original_pqtl_n = pqtl.height

        # isolate target gene from the cell-specific eQTL dataset
        # strip Ensembl version so ENSG00000095585.20 == ENSG00000095585
        if "GENE" in eqtl.columns:
            eqtl_gene_col = "GENE"
        elif "Gene" in eqtl.columns:
            eqtl_gene_col = "Gene"
        else:
            print(f"[CONCERN] Gene column not found in {eqtl_file.name}")
            print(f"[CONCERN] Columns found: {eqtl.columns}")
            continue

        if eqtl.get_column(eqtl_gene_col).dtype != pl.String:
            print(f"[CONCERN] {eqtl_gene_col} is not a string column in {eqtl_file.name}")
            continue

        eqtl_target = (
            eqtl
            .with_columns(
                pl.col(eqtl_gene_col).str.split(".").list.first().alias("_gene_no_version")
            )
            .filter(
                pl.col("_gene_no_version").is_in([gene, probe_id])
            )
            .drop("_gene_no_version")
        )

        if eqtl_target.height == 0:
            print(f"[CONCERN] {gene} was not found in the original {cell_type} eQTL parquet")
            continue

        # check SNP columns
        if "SNP" not in gwas.columns:
            print(f"[CONCERN] SNP column not found in {gwas_file}")
            continue

        if "SNP" not in pqtl.columns:
            print(f"[CONCERN] SNP column not found in {pqtl_file}")
            continue

        if "SNP" not in eqtl_target.columns:
            print(f"[CONCERN] SNP column not found in {eqtl_file}")
            continue

        # remove duplicated SNPs before overlap
        gwas = gwas.unique(subset=["SNP"], keep="first")
        pqtl = pqtl.unique(subset=["SNP"], keep="first")
        eqtl_target = eqtl_target.unique(subset=["SNP"], keep="first")

        # get SNPs shared by all three datasets
        shared_snps = (
            gwas.select("SNP")
            .join(pqtl.select("SNP"), on="SNP", how="inner")
            .join(eqtl_target.select("SNP"), on="SNP", how="inner")
            .unique()
        )

        n_shared = shared_snps.height

        if n_shared == 0:
            print(f"[CONCERN] No SNP overlap found between GWAS, pQTL and {cell_type} eQTL for {protein}")
            continue

        # now subset each dataset to exactly the same SNP set
        gwas_matched = gwas.join(shared_snps, on="SNP", how="inner").sort("SNP")
        pqtl_matched = pqtl.join(shared_snps, on="SNP", how="inner").sort("SNP")
        eqtl_matched = eqtl_target.join(shared_snps, on="SNP", how="inner").sort("SNP")
        
        # SNP overlap percentages vs original cis-locus
        gwas_overlap = n_shared / original_gwas_n * 100 if original_gwas_n > 0 else 0
        pqtl_overlap = n_shared / original_pqtl_n * 100 if original_pqtl_n > 0 else 0
        eqtl_overlap = n_shared / eqtl_target.height * 100 if eqtl_target.height > 0 else 0

        print(f"[TRACKING] Number of SNPs previously used for pQTL -> GWAS coloc: {original_gwas_n}")
        print(f"[TRACKING] Number of {cell_type} eQTL SNPs for {gene}: {eqtl_target.height}")
        print(f"[TRACKING] Number of SNPs shared across GWAS + pQTL + eQTL: {n_shared}")
        print(f"[TRACKING] GWAS overlap retained: {gwas_overlap:.2f}%")
        print(f"[TRACKING] pQTL overlap retained: {pqtl_overlap:.2f}%")
        print(f"[TRACKING] eQTL overlap retained: {eqtl_overlap:.2f}%")

        # save:
        # dat/multi_omics_targets/pqtl_dataset/gene_protein/pheno_id/cell_type/
        # sc_eqtl.parquet
        # pqtl.parquet
        # gwas.parquet
        target_out_dir = multi_omics_dir / protein / pheno_id / cell_type
        os.makedirs(target_out_dir, exist_ok=True)
        gwas_out = target_out_dir / "gwas.parquet"
        pqtl_out = target_out_dir / "pqtl.parquet"
        eqtl_out = target_out_dir / "sc_eqtl.parquet"
        gwas_matched.write_parquet(gwas_out)
        pqtl_matched.write_parquet(pqtl_out)
        eqtl_matched.write_parquet(eqtl_out)

        # save basic metadata for R coloc / moloc scripts
        metadata = pl.DataFrame({
            "protein": [protein],
            "gene": [gene],
            "phenotype": [pheno_id],
            "cell_type": [cell_type],
            "pqtl_dataset": [pqtl_dataset],
            "eqtl_dataset": [eqtl_dataset],
            "q_SMR": [q_smr],
            "p_HEIDI": [p_heidi],
            "n_gwas_original": [original_gwas_n],
            "n_pqtl_original": [original_pqtl_n],
            "n_eqtl_gene": [eqtl_target.height],
            "n_shared": [n_shared],
            "gwas_overlap_percent": [gwas_overlap],
            "pqtl_overlap_percent": [pqtl_overlap],
            "eqtl_overlap_percent": [eqtl_overlap]
        })

        metadata_file = target_out_dir / "metadata.tsv"
        metadata.write_csv(metadata_file, separator="\t")
        print(f"[TRACKING] Matched multi-omics files saved for {protein} in {cell_type}")
        print(f"[TRACKING] {gwas_out}")
        print(f"[TRACKING] {pqtl_out}")
        print(f"[TRACKING] {eqtl_out}")

        all_prepared_targets.append({
            "protein": protein,
            "gene": gene,
            "phenotype": pheno_id,
            "cell_type": cell_type,
            "pqtl_dataset": pqtl_dataset,
            "eqtl_dataset": eqtl_dataset,
            "n_shared": n_shared,
            "target_dir": str(target_out_dir)
        })

    if len(all_prepared_targets) == 0:
        print("[CONCERN] No target x cell-type datasets were prepared for coloc / moloc")
        return []

    prepared_targets_df = pl.DataFrame(all_prepared_targets)
    prepared_targets_file = res_dir / f"{pqtl_dataset}_{pheno_id}_prepared_multi_omics_targets.tsv"
    prepared_targets_df.write_csv(prepared_targets_file, separator="\t")
    print("------------------------------------------------------------")
    print(f"[TRACKING] {prepared_targets_df.height} target x cell-type datasets prepared for coloc / moloc")
    print(f"[TRACKING] Prepared target manifest saved to {prepared_targets_file}")
    return prepared_targets_df


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--pheno_id", required=True)
    p.add_argument("--pqtl_dataset", required=True)
    p.add_argument("--eqtl_dataset", required=True)
    p.add_argument("--n_cases", type=int, required=True)
    p.add_argument("--n_controls", type=int, required=True)
    args = p.parse_args()
    
    prepared_targets = target_multi_omics_moloc(
        args.pheno_id,
        args.eqtl_dataset,
        args.pqtl_dataset
    )

    if isinstance(prepared_targets, pl.DataFrame) and prepared_targets.height > 0:
        gwas_eqtl_coloc(
            pheno_id=args.pheno_id,
            pqtl_dataset=args.pqtl_dataset,
            eqtl_dataset=args.eqtl_dataset,
            n_cases=args.n_cases,
            n_controls=args.n_controls
        )


if __name__ == "__main__":
    main()