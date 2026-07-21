#!/usr/bin/env python3
import argparse
import requests
import os 
from pathlib import Path
import polars as pl 
import pandas as pd 
from drugmr import PheWAS
from statsmodels.stats.multitest import fdrcorrection
from drugmr import PyTwoSampleMR

# -----------------------------------
# THIS SCRIPT SHALL NOT BE RAN IN HPC
# -----------------------------------

# For any targets on teh final file - focus strictly on cis-pQTL -> store that as a string within variable X
# Query phewas for FinnGen r13 using the lead cis-pQTL SNP
# FinnGen BETA == effect of ALT allele
# Make sure pQTL beta + FinnGen beta == aligned to the original outcome GWAS A1 before running Wald ratio / IVW
# For each pheno -> ensure ICD-10 coded
# For each pheno (from p_nominal and beta) -> compute SE
# TO DO'S
# ONLY KEEP THE ONES WHICH FOLLOW ICD-10 CODING
# ICD-10 disease endpoints
# ONLY KEEP ICD-10 disease endpoint pertaining to the defined 2019 chapters (see Supplementary)
# clean pheWAS file to only include the phenos which == present in FinnGen

# ------ ADD COJO STUFF
# ------ IF SIGNAL == 1 SNP THEN WALD AND CARRY ON - 
# ------ ELSE - RUN IVW 
# ------ THEN FDR CORRECT ACROSS ALL ESTIMATES REGARDLESS OF METHOD

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
        # cojo results for any given target 
        mr = PyTwoSampleMR()
        cojo_file = Path(f"./results/COJO/{pqtl_dataset}/{pheno_id}/{protein}/{protein}.jma.cojo")

        if not cojo_file.exists():
            print(f"[TRACKING] No COJO results found for {protein}: {cojo_file}...")
            continue

        df_cojo = pl.read_csv(cojo_file, separator="\t")
        cojo_snps = df_cojo["SNP"].cast(pl.Utf8).to_list()
        n_cojo_snps = len(cojo_snps)

        print(
            f"[TRACKING] Number of independent COJO signals for "
            f"{protein}: {n_cojo_snps}..."
        )
        print(f"[TRACKING] COJO SNPs for {protein}: {cojo_snps}...")


        # read the original pQTL cis-region
        # COJO is used to identify which signals are independent
        # but for MR use the original marginal pQTL beta + SE
        pqtl_file = Path(f"./dat/cis_regions/{pqtl_dataset}/{protein}/pqtl.parquet")
        df_pqtl = pl.read_parquet(pqtl_file)
        df_pqtl = (
            df_pqtl
            .with_columns(
                pl.col("SNP").cast(pl.Utf8),
                pl.col("A1").cast(pl.Utf8).str.to_uppercase(),
                pl.col("A2").cast(pl.Utf8).str.to_uppercase(),
            )
            .filter(pl.col("SNP").is_in(cojo_snps))
            .unique(subset=["SNP"], keep="first")
        )

        if df_pqtl.height != n_cojo_snps:
            missing_snps = list(set(cojo_snps) - set(df_pqtl["SNP"].to_list()))
            print(
                f"[TRACKING] Could not recover all COJO SNPs from the "
                f"original pQTL file for {protein}: {missing_snps}..."
            )
            continue


        # read the original outcome GWAS cis-region
        # A1 == outcome GWAS effect allele used across all downstream analyses
        gwas_file = Path(f"./dat/cis_regions/{pqtl_dataset}/{protein}/gwas.parquet")

        if not gwas_file.exists():
            print(f"[TRACKING] No outcome GWAS cis-region found for {protein}: {gwas_file}...")
            continue

        df_gwas = pl.read_parquet(gwas_file)
        df_gwas = (
            df_gwas
            .select(
                pl.col("SNP").cast(pl.Utf8),
                pl.col("A1").cast(pl.Utf8).str.to_uppercase().alias("GWAS_A1"),
                pl.col("A2").cast(pl.Utf8).str.to_uppercase().alias("GWAS_A2"),
                pl.col("BETA").alias("GWAS_BETA"),
            )
            .filter(pl.col("SNP").is_in(cojo_snps))
            .unique(subset=["SNP"], keep="first")
        )

        if df_gwas.height != n_cojo_snps:
            missing_snps = list(set(cojo_snps) - set(df_gwas["SNP"].to_list()))
            print(
                f"[TRACKING] Could not recover all COJO SNPs from the "
                f"original outcome GWAS file for {protein}: {missing_snps}..."
            )
            continue


        # match original pQTL effects against outcome GWAS A1/A2
        df_pqtl = df_pqtl.join(
            df_gwas,
            on="SNP",
            how="inner"
        )

        if df_pqtl.height != n_cojo_snps:
            missing_snps = list(set(cojo_snps) - set(df_pqtl["SNP"].to_list()))
            print(
                f"[TRACKING] Could not match all COJO SNPs between the "
                f"pQTL and outcome GWAS files for {protein}: {missing_snps}..."
            )
            continue


        # make sure cis-pQTL SNP in protein == Lead pQTL SNP == the same one as the COJO one
        if n_cojo_snps == 1 and cojo_snps[0] != rsid:
            print(
                f"[TRACKING] Lead pQTL SNP does not match the single "
                f"COJO signal for {protein}: lead={rsid}, COJO={cojo_snps[0]}..."
            )
            continue

        if rsid not in cojo_snps:
            print(
                f"[TRACKING] Lead pQTL SNP {rsid} is not present in the "
                f"COJO-selected signals for {protein}..."
            )

        protein_results = []
        signal_phewas = {}
        # query every COJO-selected independent SNP in FinnGen
        for cojo_row in df_pqtl.iter_rows(named=True):
            cojo_rsid = str(cojo_row["SNP"])

            # A1/A2 == effect/other allele for the original pQTL beta
            pqtl_effect_allele = str(cojo_row["A1"]).upper()
            pqtl_other_allele = str(cojo_row["A2"]).upper()

            # GWAS_A1/GWAS_A2 == outcome GWAS effect/other allele
            A1 = str(cojo_row["GWAS_A1"]).upper()
            A2 = str(cojo_row["GWAS_A2"]).upper()
            gwas_beta = float(cojo_row["GWAS_BETA"])

            chromosome = str(cojo_row["CHR"])
            position = int(cojo_row["BP"])

            # FinnGen requires chromosome-position-reference-alternative
            # try both outcome GWAS allele orders and retain whichever FinnGen accepts
            snp_first_order = f"{chromosome}-{position}-{A1}-{A2}"
            snp_second_order = f"{chromosome}-{position}-{A2}-{A1}"
            response = requests.get(f"https://r13.finngen.fi/api/variant/{snp_first_order}", timeout=60)

            if response.status_code == 200 and len(response.json().get("results", [])) > 0:
                snp = snp_first_order
                finngen_ref = A1
                finngen_alt = A2
            else:
                response = requests.get(f"https://r13.finngen.fi/api/variant/{snp_second_order}", timeout=60)

                if response.status_code == 200 and len(response.json().get("results", [])) > 0:
                    snp = snp_second_order
                    finngen_ref = A2
                    finngen_alt = A1
                else:
                    print(
                        f"[TRACKING] FinnGen variant could not be resolved for {cojo_rsid}. "
                        f"Tried {snp_first_order} and {snp_second_order}..."
                    )
                    continue

            print(
                f"[TRACKING] FinnGen REF/ALT resolved for {cojo_rsid}: "
                f"REF={finngen_ref}, ALT={finngen_alt}..."
            )

            beta_original = float(cojo_row["BETA"])
            beta = beta_original
            se = float(cojo_row["SE"])
            P = float(cojo_row["P"])

            # pQTL BETA originally corresponds to pQTL A1
            # align pQTL beta to outcome GWAS A1
            if pqtl_effect_allele == A1 and pqtl_other_allele == A2:
                pqtl_A1_flipped = False
            elif pqtl_effect_allele == A2 and pqtl_other_allele == A1:
                beta = -beta
                pqtl_A1_flipped = True
            else:
                print(
                    f"[TRACKING] pQTL alleles {pqtl_effect_allele}/{pqtl_other_allele} "
                    f"do not match outcome GWAS A1/A2 {A1}/{A2} "
                    f"for {cojo_rsid}..."
                )
                continue

            print(
                f"[TRACKING] pQTL alignment for {cojo_rsid}: "
                f"original effect allele={pqtl_effect_allele}, "
                f"outcome GWAS A1={A1}, "
                f"original beta={beta_original}, "
                f"A1-aligned beta={beta}, "
                f"flipped={pqtl_A1_flipped}..."
            )
            print(f"[TRACKING] Running FinnGen PheWAS for {protein}: {cojo_rsid} ({snp})...")

            # query to phewas db and clean
            df_phewas = clean_phewas_hit(snp=snp, rsid=cojo_rsid)
            df_phewas.to_csv(os.path.join(temp_dir, f"{protein}_{cojo_rsid}_raw_hits.csv"), index=False)
            df_phewas = pl.from_pandas(df_phewas)

            # FinnGen BETA == effect of ALT allele
            # align every FinnGen phenotype beta to outcome GWAS A1
            if finngen_alt == A1 and finngen_ref == A2:
                phewas_A1_flipped = False
            elif finngen_alt == A2 and finngen_ref == A1:
                df_phewas = df_phewas.with_columns(
                    (-pl.col("BETA")).alias("BETA")
                )
                phewas_A1_flipped = True
            else:
                print(
                    f"[TRACKING] FinnGen REF/ALT {finngen_ref}/{finngen_alt} "
                    f"do not match outcome GWAS A1/A2 {A1}/{A2} "
                    f"for {cojo_rsid}..."
                )
                continue

            print(
                f"[TRACKING] FinnGen alignment for {cojo_rsid}: "
                f"FinnGen ALT={finngen_alt}, "
                f"outcome GWAS A1={A1}, "
                f"flipped={phewas_A1_flipped}..."
            )

            signal_phewas[cojo_rsid] = {
                "rsid": cojo_rsid,
                "snp": snp,
                "A1": A1,
                "A2": A2,
                "gwas_beta": gwas_beta,
                "finngen_ref": finngen_ref,
                "finngen_alt": finngen_alt,
                "pqtl_effect_allele_original": pqtl_effect_allele,
                "pqtl_other_allele_original": pqtl_other_allele,
                "beta_pqtl_original": beta_original,
                "beta_pqtl": beta,
                "pqtl_A1_flipped": pqtl_A1_flipped,
                "phewas_A1_flipped": phewas_A1_flipped,
                "se_pqtl": se,
                "p_pqtl": P,
                "df_phewas": df_phewas,
            }


        if len(signal_phewas) == 0:
            print(f"[TRACKING] No COJO-selected SNPs could be queried for {protein}...")
            continue


        # THIS IS IF SNP == 1 FOR COJO RESULT AT A GIVEN TARGET
        if len(signal_phewas) == 1:
            signal = list(signal_phewas.values())[0]
            rsid = signal["rsid"]
            snp = signal["snp"]
            A1 = signal["A1"]
            A2 = signal["A2"]
            gwas_beta = signal["gwas_beta"]
            finngen_ref = signal["finngen_ref"]
            finngen_alt = signal["finngen_alt"]
            pqtl_effect_allele = signal["pqtl_effect_allele_original"]
            pqtl_other_allele = signal["pqtl_other_allele_original"]
            beta_original = signal["beta_pqtl_original"]
            beta = signal["beta_pqtl"]
            pqtl_A1_flipped = signal["pqtl_A1_flipped"]
            phewas_A1_flipped = signal["phewas_A1_flipped"]
            se = signal["se_pqtl"]
            P = signal["p_pqtl"]
            df_phewas = signal["df_phewas"]
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
                    "method": "Wald ratio",
                    "n_instruments": 1,
                    "rsid": rsid,
                    "snp": snp,
                    "A1": A1,
                    "A2": A2,
                    "gwas_beta": gwas_beta,
                    "finngen_ref": finngen_ref,
                    "finngen_alt": finngen_alt,
                    "pqtl_effect_allele_original": pqtl_effect_allele,
                    "pqtl_other_allele_original": pqtl_other_allele,
                    "beta_pqtl_original": beta_original,
                    "beta_pqtl": beta,
                    "pqtl_A1_flipped": pqtl_A1_flipped,
                    "phewas_A1_flipped": phewas_A1_flipped,
                    "se_pqtl": se,
                    "p_pqtl": P,
                    "beta_mr": res["wald_ratio"],
                    "se_mr": res["se_wald_ratio"],
                    "p_mr": res["P_nominal"],
                    "beta_phewas": b,
                    "se_phewas": s,
                    "p_phewas": pheno["P"],
                    "PHENOCODE": pheno["PHENOCODE"],
                    "PHENOSTRING": pheno["PHENOSTRING"],
                    "CATEGORY": pheno["CATEGORY"],
                })


        # THIS IS IF SNP > 1 FOR COJO RESULT AT A GIVEN TARGET
        else:
            common_phenocodes = None
            # only run IVW on phenotypes which appear for every SNP
            for signal in signal_phewas.values():
                phenocodes = set(signal["df_phewas"]["PHENOCODE"].to_list())
                if common_phenocodes is None:
                    common_phenocodes = phenocodes
                else:
                    common_phenocodes = common_phenocodes.intersection(phenocodes)
            common_phenocodes = sorted(common_phenocodes)

            print(
                f"[TRACKING] Number of common PheWAS outcomes across "
                f"{len(signal_phewas)} independent SNPs for "
                f"{protein}: {len(common_phenocodes)}..."
            )

            for phenocode in common_phenocodes:
                exposure_rows = []
                outcome_rows = []
                instrument_rsids = []
                instrument_snps = []
                A1_values = []
                A2_values = []
                gwas_beta_values = []
                finngen_refs = []
                finngen_alts = []
                pqtl_effect_alleles = []
                pqtl_other_alleles = []
                beta_pqtl_original_values = []
                beta_pqtl_values = []
                pqtl_A1_flipped_values = []
                phewas_A1_flipped_values = []
                se_pqtl_values = []
                p_pqtl_values = []
                beta_phewas_values = []
                se_phewas_values = []
                p_phewas_values = []
                phenotype_info = None
                for signal in signal_phewas.values():
                    pheno = (
                        signal["df_phewas"]
                        .filter(pl.col("PHENOCODE") == phenocode)
                        .unique(subset=["PHENOCODE"], keep="first")
                    )
                    if pheno.height == 0:
                        continue

                    pheno = pheno.row(0, named=True)

                    exposure_rows.append({
                        "SNP": signal["rsid"],
                        "BETA": signal["beta_pqtl"],
                        "SE": signal["se_pqtl"],
                    })

                    outcome_rows.append({
                        "SNP": signal["rsid"],
                        "BETA": pheno["BETA"],
                        "SE": pheno["SE"],
                    })

                    instrument_rsids.append(signal["rsid"])
                    instrument_snps.append(signal["snp"])
                    A1_values.append(signal["A1"])
                    A2_values.append(signal["A2"])
                    gwas_beta_values.append(signal["gwas_beta"])
                    finngen_refs.append(signal["finngen_ref"])
                    finngen_alts.append(signal["finngen_alt"])
                    pqtl_effect_alleles.append(signal["pqtl_effect_allele_original"])
                    pqtl_other_alleles.append(signal["pqtl_other_allele_original"])
                    beta_pqtl_original_values.append(signal["beta_pqtl_original"])
                    beta_pqtl_values.append(signal["beta_pqtl"])
                    pqtl_A1_flipped_values.append(signal["pqtl_A1_flipped"])
                    phewas_A1_flipped_values.append(signal["phewas_A1_flipped"])
                    se_pqtl_values.append(signal["se_pqtl"])
                    p_pqtl_values.append(signal["p_pqtl"])
                    beta_phewas_values.append(pheno["BETA"])
                    se_phewas_values.append(pheno["SE"])
                    p_phewas_values.append(pheno["P"])
                    if phenotype_info is None:
                        phenotype_info = {
                            "PHENOCODE": pheno["PHENOCODE"],
                            "PHENOSTRING": pheno["PHENOSTRING"],
                            "CATEGORY": pheno["CATEGORY"],
                        }

                exposure_df = pl.DataFrame(exposure_rows)
                outcome_df = pl.DataFrame(outcome_rows)

                if exposure_df.height < 2:
                    continue

                ivw_result, snp_results = mr.IVW(
                    exposure_df=exposure_df,
                    outcome_df=outcome_df,
                    exposure_snp_col="SNP",
                    exposure_beta_col="BETA",
                    exposure_se_col="SE",
                    outcome_snp_col="SNP",
                    outcome_beta_col="BETA",
                    outcome_se_col="SE"
                )
                ivw_result = ivw_result.row(0, named=True)
                protein_results.append({
                    "protein": protein,
                    "method": "IVW delta",
                    "n_instruments": ivw_result["N_SNPS"],
                    "rsid": ",".join(instrument_rsids),
                    "snp": ",".join(instrument_snps),
                    "A1": ",".join(A1_values),
                    "A2": ",".join(A2_values),
                    "gwas_beta": ",".join(
                        str(value) for value in gwas_beta_values
                    ),
                    "finngen_ref": ",".join(finngen_refs),
                    "finngen_alt": ",".join(finngen_alts),
                    "pqtl_effect_allele_original": ",".join(pqtl_effect_alleles),
                    "pqtl_other_allele_original": ",".join(pqtl_other_alleles),
                    "beta_pqtl_original": ",".join(
                        str(value) for value in beta_pqtl_original_values
                    ),
                    "beta_pqtl": ",".join(
                        str(value) for value in beta_pqtl_values
                    ),
                    "pqtl_A1_flipped": ",".join(
                        str(value) for value in pqtl_A1_flipped_values
                    ),
                    "phewas_A1_flipped": ",".join(
                        str(value) for value in phewas_A1_flipped_values
                    ),
                    "se_pqtl": ",".join(
                        str(value) for value in se_pqtl_values
                    ),
                    "p_pqtl": ",".join(
                        str(value) for value in p_pqtl_values
                    ),
                    "beta_mr": ivw_result["BETA"],
                    "se_mr": ivw_result["SE"],
                    "p_mr": ivw_result["P"],
                    "beta_phewas": ",".join(
                        str(value) for value in beta_phewas_values
                    ),
                    "se_phewas": ",".join(
                        str(value) for value in se_phewas_values
                    ),
                    "p_phewas": ",".join(
                        str(value) for value in p_phewas_values
                    ),
                    "PHENOCODE": phenotype_info["PHENOCODE"],
                    "PHENOSTRING": phenotype_info["PHENOSTRING"],
                    "CATEGORY": phenotype_info["CATEGORY"],
                })


        if len(protein_results) == 0:
            print(f"[TRACKING] No Wald ratio / IVW estimates generated for {protein}...")
            continue

        df_protein_results = pl.DataFrame(protein_results)
        # FDR correct across all phenotypes tested for this protein regardless of method
        # rejected, q = fdrcorrection(df_protein_results["p_mr"].to_numpy(), alpha=0.05, method="indep")
        rejected, q = fdrcorrection(df_protein_results["p_mr"].to_numpy(), alpha=0.05)
        df_protein_results = df_protein_results.with_columns([
            pl.Series("q_fdr_mr", q),
            pl.Series("fdr_significant", rejected)
        ])

        results.extend(df_protein_results.to_dicts())

    if len(results) == 0:
        print("[TRACKING] No PheWAS associations were generated...")
        return

    df_results = pl.DataFrame(results)
    df_results = df_results.sort(["protein", "p_mr"])
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