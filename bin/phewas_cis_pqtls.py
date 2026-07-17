#!/usr/bin/env python3
import argparse
import requests
import os 
import polars as pl 
import pandas as pd 
from drugmr import PheWAS
from statsmodels.stats.multitest import fdrcorrection

# -----------------------------------
# THIS SCRIPT SHALL NOT BE RAN IN HPC
# -----------------------------------

# For any targets on teh final file - focus strictly on cis-pQTL -> store that as a string within variable X
# Query phewas for FinnGen r13 using the lead cis-pQTL SNP
# FinnGen BETA == effect of ALT allele
# Make sure pQTL beta == aligned to the FinnGen ALT allele before running Wald ratio
# For each pheno -> ensure ICD-10 coded
# For each pheno (from p_nominal and beta) -> compute SE


# TO DO'S
# ONLY KEEP THE ONES WHICH FOLLOW ICD-10 CODING
# ICD-10 disease endpoints
# ONLY KEEP ICD-10 disease endpoint pertaining to the defined 2019 chapters (see Supplementary)

# clean pheWAS file to only include the phenos which == present in FinnGen
def clean_phewas_hit(snp: str, rsid: str):
    # snp = "10-96304051-A-G" # chromosome-position-reference-alternative format
    # query SNP using the FinnGen API
    response = requests.get(f"https://r13.finngen.fi/api/variant/{snp}", timeout=60)
    response.raise_for_status()
    data = response.json()
    results = data["results"]
    df = pd.DataFrame(results)
    df = df.dropna(subset=["beta", "sebeta", "pval", "phenocode", "phenostring", "category"])
    print(f"[TRACKING] Number of phenotypes in PheWAS for {rsid}: {len(df)} across {df['category'].nunique()} categories...")
    # retain only ICD disease chapters
    icd_chapters = (
        "I ", "II ", "III ", "IV ", "V ", "VI ", "VII ", "VIII ", "IX ",
        "X ", "XI ", "XII ", "XIII ", "XIV ", "XV ", "XVI ", "XVII "
    )
    df = df[df["category"].str.startswith(icd_chapters)]
    print(f"[TRACKING] Number of ICD disease endpoints retained for {rsid}: {len(df)}")
    df = df[["beta", "sebeta", "pval", "phenocode", "phenostring", "category"]]
    df = df.rename(columns={
        "beta": "BETA",
        "sebeta": "SE",
        "pval": "P",
        "phenocode": "PHENOCODE",
        "phenostring": "PHENOSTRING",
        "category": "CATEGORY",
    })
    return df


# this script runs AFTER THE final TARGET file 
def phewas_for_compelling_targets(pheno_id: str, pqtl_dataset: str, eqtl_dataset: str):
    # ukb_ppp_AD_multi_omics_snp_evidence.tsv
    # for any target.unique()
    # Lead pQTL SNP
    top_snp_file = f"./results/SMR/{eqtl_dataset}/{pheno_id}/{pqtl_dataset}_{pheno_id}_multi_omics_snp_evidence.tsv"
    df = pl.read_csv(top_snp_file, separator="\t")
    # temp_dir
    temp_dir = f"./work/PheWAS/{pqtl_dataset}_{pheno_id}"
    os.makedirs(temp_dir, exist_ok=True)
    out_dir = f"./results/PheWAS/{pqtl_dataset}/{pheno_id}"
    os.makedirs(out_dir, exist_ok=True)
    # only retain the Lead pQTL SNP for each protein
    df = (df.filter(pl.col("snp_type") == "Lead pQTL SNP").unique(subset=["protein", "SNP"], keep="first"))
    print(f"[TRACKING] Number of compelling targets for PheWAS: {df['protein'].n_unique()}...")
    results = []
    for row in df.iter_rows(named=True):
        protein = row["protein"]
        rsid = row["SNP"]
        # A1 == effect allele for the harmonised pQTL beta
        pqtl_effect_allele = str(row["pQTL_A1"]).upper()
        pqtl_other_allele = str(row["pQTL_A2"]).upper()
        chromosome = str(row["CHR"])
        position = int(row["POS"])
        # FinnGen requires chromosome-position-reference-alternative
        # A1/A2 in our file == effect/other rather than necessarily REF/ALT
        # therefore try both allele orders and retain whichever FinnGen accepts
        snp_first_order = f"{chromosome}-{position}-{pqtl_effect_allele}-{pqtl_other_allele}"
        snp_second_order = f"{chromosome}-{position}-{pqtl_other_allele}-{pqtl_effect_allele}"
        response = requests.get(f"https://r13.finngen.fi/api/variant/{snp_first_order}", timeout=60)

        if response.status_code == 200 and len(response.json().get("results", [])) > 0:
            snp = snp_first_order
            finngen_ref = pqtl_effect_allele
            finngen_alt = pqtl_other_allele
        else:
            response = requests.get(f"https://r13.finngen.fi/api/variant/{snp_second_order}", timeout=60)

            if response.status_code == 200 and len(response.json().get("results", [])) > 0:
                snp = snp_second_order
                finngen_ref = pqtl_other_allele
                finngen_alt = pqtl_effect_allele
            else:
                raise ValueError(
                    f"FinnGen variant could not be resolved for {rsid}. "
                    f"Tried {snp_first_order} and {snp_second_order}..."
                )

        print(
            f"[TRACKING] FinnGen REF/ALT resolved for {rsid}: "
            f"REF={finngen_ref}, ALT={finngen_alt}..."
        )

        beta_original = row["pQTL_beta"]
        beta = beta_original
        se = row["pQTL_SE"]
        P = row["pQTL_P"]
        # FinnGen BETA == effect of ALT allele
        # therefore pQTL beta must also correspond to ALT before Wald ratio
        if pqtl_effect_allele == finngen_alt and pqtl_other_allele == finngen_ref:
            pqtl_finngen_flipped = False
        elif pqtl_effect_allele == finngen_ref and pqtl_other_allele == finngen_alt:
            beta = -beta
            pqtl_finngen_flipped = True
        else:
            raise ValueError(
                f"pQTL alleles {pqtl_effect_allele}/{pqtl_other_allele} "
                f"do not match FinnGen REF/ALT {finngen_ref}/{finngen_alt} "
                f"for {rsid}..."
            )
        print(
            f"[TRACKING] pQTL alignment for {rsid}: "
            f"original effect allele={pqtl_effect_allele}, "
            f"FinnGen ALT={finngen_alt}, "
            f"original beta={beta_original}, "
            f"ALT-aligned beta={beta}, "
            f"flipped={pqtl_finngen_flipped}..."
        )
        print(f"[TRACKING] Running FinnGen PheWAS for {protein}: {rsid} ({snp})...")

        # query to phewas db and clean
        df_phewas = clean_phewas_hit(snp=snp, rsid=rsid)
        df_phewas.to_csv(os.path.join(temp_dir, f"{rsid}_raw_hits.csv"), index=False)
        df_phewas = pl.from_pandas(df_phewas)
        # run in-house phewas approach 
        # pQTL = denominator -> B_X
        # GWAS = numerator -> B_Y
        # both betas == aligned to FinnGen ALT allele
        protein_results = []
        for pheno in df_phewas.iter_rows(named=True):
            b = pheno["BETA"]
            s = pheno["SE"]
            res = PheWAS(
                B_Y=b,
                SE_Y=s,
                B_X=beta,
                SE_X=se
            )
            protein_results.append({
                "protein": protein,
                "rsid": rsid,
                "snp": snp,
                "finngen_ref": finngen_ref,
                "finngen_alt": finngen_alt,
                "pqtl_effect_allele_original": pqtl_effect_allele,
                "pqtl_other_allele_original": pqtl_other_allele,
                "beta_pqtl_original": beta_original,
                "beta_pqtl": beta,
                "pqtl_finngen_flipped": pqtl_finngen_flipped,
                "se_pqtl": se,
                "p_pqtl": P,
                "wald_ratio": res["wald_ratio"],
                "se_wald_ratio": res["se_wald_ratio"],
                "p_wald_ratio": res["P_nominal"],
                "beta_phewas": b,
                "se_phewas": s,
                "p_phewas": pheno["P"],
                "PHENOCODE": pheno["PHENOCODE"],
                "PHENOSTRING": pheno["PHENOSTRING"],
                "CATEGORY": pheno["CATEGORY"],
            })

        df_protein_results = pl.DataFrame(protein_results)
        # FDR correct across all phenotypes tested for this protein / SNP
        # rejected, q = fdrcorrection(df_protein_results["p_wald_ratio"].to_numpy(), alpha=0.05, method="indep")
        rejected, q = fdrcorrection(df_protein_results["p_wald_ratio"].to_numpy(), alpha=0.05)
        df_protein_results = df_protein_results.with_columns([pl.Series("q_fdr_wald_ratio", q), pl.Series("fdr_significant", rejected)])
        results.extend(df_protein_results.to_dicts())
    df_results = pl.DataFrame(results)
    df_results = df_results.sort(["protein", "p_wald_ratio"])
    df_results.write_csv(os.path.join(out_dir, f"{pqtl_dataset}_{pheno_id}_PheWAS.tsv"), separator="\t")
    print(f"[TRACKING] PheWAS completed: {df_results.height} associations saved...")


# pheno_id: str, pqtl_dataset: str, eqtl_dataset: str
def main():
    p = argparse.ArgumentParser()
    p.add_argument("--pheno_id", required=True)
    p.add_argument("--pqtl_dataset", required=True)
    p.add_argument("--eqtl_dataset", required=True)
    args = p.parse_args()
    phewas_for_compelling_targets(
        pheno_id=args.pheno_id,
        pqtl_dataset=args.pqtl_dataset,
        eqtl_dataset=args.eqtl_dataset
    )


if __name__ == "__main__":
    main()