#!/usr/bin/env python3
import argparse
import requests
import os 
from pathlib import Path
import polars as pl 
import pandas as pd 
from drugmr import PheWAS
from drugmr import paths
# from statsmodels.stats.multitest import fdrcorrection
from drugmr.twosamplemr import PyTwoSampleMR

# -----------------------------------
# THIS SCRIPT SHALL NOT BE RAN IN HPC
# -----------------------------------
# 2,511 ICD coded endpoints vs total 2,755 
finngen_icd_endpoints = 2511
# The COLOC output is used strictly to define which protein targets go into PheWAS
# For every selected target - grab the exact SNPs which were used in cis-MR
# Make A1 the AD risk allele for every cis-MR instrument
# Align pQTL beta + FinnGen beta to the AD risk allele before Wald ratio / IVW
# FinnGen BETA == effect of ALT allele
# For each pheno -> ensure ICD-10 coded
# TO DO'S
# ONLY KEEP THE ONES WHICH FOLLOW ICD-10 CODING
# ICD-10 disease endpoints
# ONLY KEEP ICD-10 disease endpoint pertaining to the defined 2019 chapters (see Supplementary)
# clean pheWAS file to only include the phenos which == present in FinnGen

# ------ NO COJO HERE
# ------ IF 1 cis-MR INSTRUMENT IS AVAILABLE IN FINNGEN THEN WALD AND CARRY ON - 
# ------ IF >1 cis-MR INSTRUMENT IS AVAILABLE IN FINNGEN THEN RUN IVW
# ------ THEN BONFERRONI CORRECT ACROSS ALL ESTIMATES REGARDLESS OF METHOD


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


# this script runs AFTER cis-MR + pairwise pQTL-GWAS COLOC
def phewas_for_compelling_targets(pheno_id: str, pqtl_dataset: str, local_results_dir: str = "results"):
    # COLOC defines the targets only
    # each protein_id == its own protein assay / aptamer
    coloc_file = paths.coloc_out(pqtl_dataset, pheno_id, local_results_dir)
    df_coloc = pl.read_csv(coloc_file, separator="\t")

    if "protein_id" in df_coloc.columns:
        df_coloc = df_coloc.rename({"protein_id": "protein"})

    if "protein" not in df_coloc.columns:
        raise ValueError(
            f"Could not find protein_id or protein in COLOC file: {coloc_file}"
        )
    compelling_targets = (df_coloc.select(pl.col("protein").cast(pl.Utf8)).drop_nulls().unique(maintain_order=True))
    # exact harmonised instruments which were used in the original cis-MR
    instruments_file = paths.mr_instruments_out(pqtl_dataset, pheno_id, local_results_dir)
    df_instruments = pl.read_csv(instruments_file, separator="\t")
    required_instrument_cols = [
        "protein",
        "pqtl_dataset",
        "outcome_trait",
        "SNP",
        "effect_allele.exposure",
        "other_allele.exposure",
        "beta.exposure",
        "se.exposure",
        "pval.exposure",
        "effect_allele.outcome",
        "other_allele.outcome",
        "beta.outcome",
        "used_in_mr",
    ]

    df_instruments = (
        df_instruments
        .with_columns(
            pl.col("protein").cast(pl.Utf8),
            pl.col("pqtl_dataset").cast(pl.Utf8),
            pl.col("outcome_trait").cast(pl.Utf8),
            pl.col("SNP").cast(pl.Utf8),
            pl.col("effect_allele.exposure").cast(pl.Utf8).str.to_uppercase(),
            pl.col("other_allele.exposure").cast(pl.Utf8).str.to_uppercase(),
            pl.col("beta.exposure").cast(pl.Float64),
            pl.col("se.exposure").cast(pl.Float64),
            pl.col("pval.exposure").cast(pl.Float64),
            pl.col("effect_allele.outcome").cast(pl.Utf8).str.to_uppercase(),
            pl.col("other_allele.outcome").cast(pl.Utf8).str.to_uppercase(),
            pl.col("beta.outcome").cast(pl.Float64),
            (
                pl.col("used_in_mr")
                .cast(pl.Utf8)
                .str.to_uppercase()
                .eq("TRUE")
            ).alias("used_in_mr"),
        )
        .filter(
            (pl.col("pqtl_dataset") == pqtl_dataset) &
            (pl.col("outcome_trait") == pheno_id) &
            (pl.col("used_in_mr") == True)
        )
        .join(compelling_targets, on="protein", how="inner")
        .unique(subset=["protein", "SNP"], keep="first")
    )

    # temp_dir
    temp_dir = f"./work/PheWAS-FinnGen/{pqtl_dataset}_{pheno_id}"
    os.makedirs(temp_dir, exist_ok=True)
    phewas_out_file = paths.phewas_out(pqtl_dataset, pheno_id, local_results_dir)
    os.makedirs(phewas_out_file.parent, exist_ok=True)
    print(f"[TRACKING] Number of compelling COLOC targets for PheWAS: {compelling_targets.height}...")
    print(f"[TRACKING] Number of retained cis-MR instrument rows: {df_instruments.height}...")
    results = []
    mr = PyTwoSampleMR()
    # run every protein assay / aptamer within the COLOC file
    for target_row in compelling_targets.iter_rows(named=True):
        protein = target_row["protein"]
        df_target_instruments = (
            df_instruments
            .filter(pl.col("protein") == protein)
            .sort("instrument_rank")
            if "instrument_rank" in df_instruments.columns
            else df_instruments.filter(pl.col("protein") == protein)
        )

        if df_target_instruments.height == 0:
            print(
                f"[TRACKING] No used_in_mr cis-MR instruments found for "
                f"{protein} in {instruments_file}..."
            )
            continue

        instrument_snps = df_target_instruments["SNP"].to_list()
        n_cis_mr_instruments = len(instrument_snps)
        print(
            f"[TRACKING] Number of cis-MR instruments for "
            f"{protein}: {n_cis_mr_instruments}..."
        )
        print(f"[TRACKING] cis-MR SNPs for {protein}: {instrument_snps}...")

        # CHR + BP are only recovered here to construct the FinnGen API variant
        # exposure beta + SE + alleles remain strictly from the cis-MR instruments file
        pqtl_file = Path(f"./dat/cis_regions/{pqtl_dataset}/{protein}/pqtl.parquet")
        if not pqtl_file.exists():
            print(f"[TRACKING] Original pQTL cis-region not found for {protein}: {pqtl_file}...")
            continue

        df_pqtl_positions = (
            pl.read_parquet(pqtl_file)
            .select(
                pl.col("SNP").cast(pl.Utf8),
                pl.col("CHR").cast(pl.Utf8),
                pl.col("BP").cast(pl.Int64),
            )
            .filter(pl.col("SNP").is_in(instrument_snps))
            .unique(subset=["SNP"], keep="first")
        )

        df_target_instruments = df_target_instruments.join(
            df_pqtl_positions,
            on="SNP",
            how="inner"
        )

        if df_target_instruments.height != len(instrument_snps):
            recovered_snps = set(df_target_instruments["SNP"].to_list())
            missing_snps = sorted(set(instrument_snps) - recovered_snps)
            print(
                f"[TRACKING] Could not recover CHR/BP for all cis-MR "
                f"instruments for {protein}: {missing_snps}..."
            )
            continue

        signal_phewas = {}
        # query every SNP which was actually used in cis-MR
        for instrument_row in df_target_instruments.iter_rows(named=True):
            rsid = str(instrument_row["SNP"])
            exposure_effect_allele = str(instrument_row["effect_allele.exposure"]).upper()
            exposure_other_allele = str(instrument_row["other_allele.exposure"]).upper()
            beta_exposure = float(instrument_row["beta.exposure"])
            se_exposure = float(instrument_row["se.exposure"])
            p_exposure = float(instrument_row["pval.exposure"])
            ad_effect_allele = str(instrument_row["effect_allele.outcome"]).upper()
            ad_other_allele = str(instrument_row["other_allele.outcome"]).upper()
            beta_ad_original = float(instrument_row["beta.outcome"])
            chromosome = str(instrument_row["CHR"]).replace("chr", "")
            position = int(instrument_row["BP"])
            # make A1 the AD risk-increasing allele
            # this guarantees beta_ad > 0 relative to A1
            if beta_ad_original > 0:
                A1 = ad_effect_allele
                A2 = ad_other_allele
                beta_ad = beta_ad_original
                ad_A1_flipped = False
            elif beta_ad_original < 0:
                A1 = ad_other_allele
                A2 = ad_effect_allele
                beta_ad = -beta_ad_original
                ad_A1_flipped = True
            else:
                print(f"[TRACKING] AD beta is zero for {rsid}...")
                continue

            # pQTL BETA originally corresponds to effect_allele.exposure
            # align pQTL beta to the AD risk allele A1
            beta_exposure_original = beta_exposure
            if (
                exposure_effect_allele == A1 and
                exposure_other_allele == A2
            ):
                exposure_A1_flipped = False

            elif (
                exposure_effect_allele == A2 and
                exposure_other_allele == A1
            ):
                beta_exposure = -beta_exposure
                exposure_A1_flipped = True

            else:
                print(
                    f"[TRACKING] Exposure alleles "
                    f"{exposure_effect_allele}/{exposure_other_allele} "
                    f"do not match AD risk alleles {A1}/{A2} "
                    f"for {rsid}..."
                )
                continue

            print(
                f"[TRACKING] AD risk alignment for {rsid}: "
                f"A1={A1}, A2={A2}, "
                f"AD beta={beta_ad}, "
                f"original pQTL beta={beta_exposure_original}, "
                f"A1-aligned pQTL beta={beta_exposure}, "
                f"pQTL flipped={exposure_A1_flipped}..."
            )

            # FinnGen requires chromosome-position-reference-alternative
            # try both allele orders and retain whichever FinnGen accepts
            snp_first_order = (
                f"{chromosome}-{position}-"
                f"{A1}-{A2}"
            )
            snp_second_order = (
                f"{chromosome}-{position}-"
                f"{A2}-{A1}"
            )
            response = requests.get(
                f"https://r13.finngen.fi/api/variant/{snp_first_order}",
                timeout=60
            )

            if response.status_code == 200 and len(response.json().get("results", [])) > 0:
                snp = snp_first_order
                finngen_ref = A1
                finngen_alt = A2
            else:
                response = requests.get(
                    f"https://r13.finngen.fi/api/variant/{snp_second_order}",
                    timeout=60
                )

                if response.status_code == 200 and len(response.json().get("results", [])) > 0:
                    snp = snp_second_order
                    finngen_ref = A2
                    finngen_alt = A1
                else:
                    print(
                        f"[TRACKING] FinnGen variant could not be resolved for {rsid}. "
                        f"Tried {snp_first_order} and {snp_second_order}..."
                    )
                    continue

            print(
                f"[TRACKING] FinnGen REF/ALT resolved for {rsid}: "
                f"REF={finngen_ref}, ALT={finngen_alt}..."
            )
            print(f"[TRACKING] Running FinnGen PheWAS for {protein}: {rsid} ({snp})...")

            # query to phewas db and clean
            df_phewas = clean_phewas_hit(snp=snp, rsid=rsid)
            df_phewas.to_csv(os.path.join(temp_dir, f"{protein}_{rsid}_raw_hits.csv"), index=False)
            df_phewas = pl.from_pandas(df_phewas)

            # FinnGen BETA == effect of ALT allele
            # align FinnGen beta directly to the AD risk allele A1
            if (
                finngen_alt == A1 and
                finngen_ref == A2
            ):
                phewas_A1_flipped = False
            elif (
                finngen_alt == A2 and
                finngen_ref == A1
            ):
                df_phewas = df_phewas.with_columns(
                    (-pl.col("BETA")).alias("BETA")
                )
                phewas_A1_flipped = True
            else:
                print(
                    f"[TRACKING] FinnGen REF/ALT {finngen_ref}/{finngen_alt} "
                    f"do not match AD risk alleles {A1}/{A2} "
                    f"for {rsid}..."
                )
                continue

            print(
                f"[TRACKING] FinnGen alignment for {rsid}: "
                f"FinnGen ALT={finngen_alt}, "
                f"AD risk A1={A1}, "
                f"flipped={phewas_A1_flipped}..."
            )

            signal_phewas[rsid] = {
                "rsid": rsid,
                "snp": snp,
                "A1": A1,
                "A2": A2,
                "ad_effect_allele_original": ad_effect_allele,
                "ad_other_allele_original": ad_other_allele,
                "beta_ad_original": beta_ad_original,
                "beta_ad": beta_ad,
                "ad_A1_flipped": ad_A1_flipped,
                "exposure_effect_allele": exposure_effect_allele,
                "exposure_other_allele": exposure_other_allele,
                "finngen_ref": finngen_ref,
                "finngen_alt": finngen_alt,
                "beta_exposure_original": beta_exposure_original,
                "beta_exposure": beta_exposure,
                "exposure_A1_flipped": exposure_A1_flipped,
                "se_exposure": se_exposure,
                "p_exposure": p_exposure,
                "phewas_A1_flipped": phewas_A1_flipped,
                "df_phewas": df_phewas,
            }

        # define how many of the original cis-MR instruments were found in FinnGen
        n_instruments_original = n_cis_mr_instruments
        n_instruments_available = len(signal_phewas)
        missing_snps = sorted(set(instrument_snps) - set(signal_phewas.keys()))
        instrument_completeness = (
            n_instruments_available / n_instruments_original
        )

        if n_instruments_available == 0:
            print(
                f"[TRACKING] No cis-MR instruments found in FinnGen "
                f"for {protein}..."
            )
            continue

        if n_instruments_available < n_instruments_original:
            print(
                f"[TRACKING] Using available FinnGen instruments for "
                f"{protein}: {n_instruments_available}/"
                f"{n_instruments_original}; "
                f"missing={missing_snps}..."
            )

        protein_results = []

        # THIS IS IF SNP == 1 FOR AVAILABLE cis-MR INSTRUMENTS AT A GIVEN TARGET
        if n_instruments_available == 1:
            signal = list(signal_phewas.values())[0]
            for pheno in signal["df_phewas"].iter_rows(named=True):
                res = PheWAS(
                    B_Y=pheno["BETA"],
                    SE_Y=pheno["SE"],
                    B_X=signal["beta_exposure"],
                    SE_X=signal["se_exposure"]
                )
                protein_results.append({
                    "protein": protein,
                    "pqtl_dataset": pqtl_dataset,
                    "method": "Wald ratio",
                    "n_instruments_original": n_instruments_original,
                    "n_instruments_available": n_instruments_available,
                    "n_instruments": 1,
                    "missing_instruments": ",".join(missing_snps),
                    "instrument_completeness": instrument_completeness,
                    "rsid": signal["rsid"],
                    "snp": signal["snp"],
                    "A1": signal["A1"],
                    "A2": signal["A2"],
                    "ad_effect_allele_original": signal["ad_effect_allele_original"],
                    "ad_other_allele_original": signal["ad_other_allele_original"],
                    "beta_ad_original": str(signal["beta_ad_original"]),
                    "beta_ad": str(signal["beta_ad"]),
                    "ad_A1_flipped": str(signal["ad_A1_flipped"]),
                    "exposure_effect_allele": signal["exposure_effect_allele"],
                    "exposure_other_allele": signal["exposure_other_allele"],
                    "finngen_ref": signal["finngen_ref"],
                    "finngen_alt": signal["finngen_alt"],
                    "beta_exposure_original": str(signal["beta_exposure_original"]),
                    "beta_exposure": str(signal["beta_exposure"]),
                    "exposure_A1_flipped": str(signal["exposure_A1_flipped"]),
                    "se_exposure": str(signal["se_exposure"]),
                    "p_exposure": str(signal["p_exposure"]),
                    "phewas_A1_flipped": str(
                        signal["phewas_A1_flipped"]
                    ),
                    "beta_mr": res["wald_ratio"],
                    "se_mr": res["se_wald_ratio"],
                    "p_mr": res["P_nominal"],
                    "beta_phewas": str(pheno["BETA"]),
                    "se_phewas": str(pheno["SE"]),
                    "p_phewas": str(pheno["P"]),
                    "PHENOCODE": pheno["PHENOCODE"],
                    "PHENOSTRING": pheno["PHENOSTRING"],
                    "CATEGORY": pheno["CATEGORY"],
                })

        # THIS IS IF SNP > 1 FOR AVAILABLE cis-MR INSTRUMENTS AT A GIVEN TARGET
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
                f"{len(signal_phewas)} cis-MR instruments for "
                f"{protein}: {len(common_phenocodes)}..."
            )

            for phenocode in common_phenocodes:
                exposure_rows = []
                outcome_rows = []
                instrument_rsids = []
                instrument_snps = []
                A1_values = []
                A2_values = []
                ad_effect_alleles = []
                ad_other_alleles = []
                beta_ad_original_values = []
                beta_ad_values = []
                ad_A1_flipped_values = []
                exposure_effect_alleles = []
                exposure_other_alleles = []
                finngen_refs = []
                finngen_alts = []
                beta_exposure_original_values = []
                beta_exposure_values = []
                exposure_A1_flipped_values = []
                se_exposure_values = []
                p_exposure_values = []
                phewas_flipped_values = []
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
                        "BETA": signal["beta_exposure"],
                        "SE": signal["se_exposure"],
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
                    ad_effect_alleles.append(signal["ad_effect_allele_original"])
                    ad_other_alleles.append(signal["ad_other_allele_original"])
                    beta_ad_original_values.append(signal["beta_ad_original"])
                    beta_ad_values.append(signal["beta_ad"])
                    ad_A1_flipped_values.append(signal["ad_A1_flipped"])
                    exposure_effect_alleles.append(signal["exposure_effect_allele"])
                    exposure_other_alleles.append(signal["exposure_other_allele"])
                    finngen_refs.append(signal["finngen_ref"])
                    finngen_alts.append(signal["finngen_alt"])
                    beta_exposure_original_values.append(
                        signal["beta_exposure_original"]
                    )
                    beta_exposure_values.append(signal["beta_exposure"])
                    exposure_A1_flipped_values.append(
                        signal["exposure_A1_flipped"]
                    )
                    se_exposure_values.append(signal["se_exposure"])
                    p_exposure_values.append(signal["p_exposure"])
                    phewas_flipped_values.append(
                        signal["phewas_A1_flipped"]
                    )
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
                    print(
                        f"[TRACKING] Incomplete IVW instrument set for "
                        f"{protein} and {phenocode}..."
                    )
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
                    "pqtl_dataset": pqtl_dataset,
                    "method": "IVW delta",
                    "n_instruments_original": n_instruments_original,
                    "n_instruments_available": n_instruments_available,
                    "n_instruments": ivw_result["N_SNPS"],
                    "missing_instruments": ",".join(missing_snps),
                    "instrument_completeness": instrument_completeness,
                    "rsid": ",".join(instrument_rsids),
                    "snp": ",".join(instrument_snps),
                    "A1": ",".join(A1_values),
                    "A2": ",".join(A2_values),
                    "ad_effect_allele_original": ",".join(ad_effect_alleles),
                    "ad_other_allele_original": ",".join(ad_other_alleles),
                    "beta_ad_original": ",".join(
                        str(value) for value in beta_ad_original_values
                    ),
                    "beta_ad": ",".join(
                        str(value) for value in beta_ad_values
                    ),
                    "ad_A1_flipped": ",".join(
                        str(value) for value in ad_A1_flipped_values
                    ),
                    "exposure_effect_allele": ",".join(exposure_effect_alleles),
                    "exposure_other_allele": ",".join(exposure_other_alleles),
                    "finngen_ref": ",".join(finngen_refs),
                    "finngen_alt": ",".join(finngen_alts),
                    "beta_exposure_original": ",".join(
                        str(value) for value in beta_exposure_original_values
                    ),
                    "beta_exposure": ",".join(
                        str(value) for value in beta_exposure_values
                    ),
                    "exposure_A1_flipped": ",".join(
                        str(value) for value in exposure_A1_flipped_values
                    ),
                    "se_exposure": ",".join(
                        str(value) for value in se_exposure_values
                    ),
                    "p_exposure": ",".join(
                        str(value) for value in p_exposure_values
                    ),
                    "phewas_A1_flipped": ",".join(
                        str(value) for value in phewas_flipped_values
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
        # Bonferroni correct across all predefined FinnGen ICD endpoints for this protein
        df_protein_results = df_protein_results.with_columns([
            pl.min_horizontal(
                pl.col("p_mr") * finngen_icd_endpoints,
                pl.lit(1.0)
            ).alias("p_bonferroni"),
            (
                pl.col("p_mr") < (0.05 / finngen_icd_endpoints)
            ).alias("bonferroni_significant")
        ])

        results.extend(df_protein_results.to_dicts())

    if len(results) == 0:
        print("[TRACKING] No PheWAS associations were generated...")
        return

    df_results = pl.DataFrame(results, infer_schema_length=None)
    df_results = df_results.sort(["protein", "p_mr"])
    df_results.write_csv(phewas_out_file, separator="\t")
    print(f"[TRACKING] PheWAS completed: {df_results.height} associations saved...")


# pheno_id: str, pqtl_dataset: str
def main():
    p = argparse.ArgumentParser()
    p.add_argument("--pheno_id", required=True)
    p.add_argument("--pqtl_dataset", required=True)
    p.add_argument("--local_results_dir", default="results")
    args = p.parse_args()
    phewas_for_compelling_targets(
        pheno_id=args.pheno_id,
        pqtl_dataset=args.pqtl_dataset,
        local_results_dir=args.local_results_dir
    )


if __name__ == "__main__":
    main()