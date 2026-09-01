#!/usr/bin/env python3
import argparse
import os
from pathlib import Path
import polars as pl
import requests
from drugmr.phewas import PheWAS
from drugmr import paths
from drugmr.twosamplemr import PyTwoSampleMR

# Bonferroni correction is applied per-protein, across however many endpoints
# were actually tested for that protein - NOT a fixed constant (see
# n_endpoints_tested below). Kept only as a documentation reference.
UKB_TOPMED_TOTAL_PHENOS = 1419
coloc_threshold = 0

def grab_phewas_info(snp: str, rsid: str):
    snp = snp.replace(":", "-")
    url = f"https://pheweb.org/UKB-TOPMed/api/variant/{snp}"
    response = requests.get(url, timeout=60)
    response.raise_for_status()
    data = response.json()
    if not isinstance(data, dict):
        raise ValueError(f"Invalid UKBB PheWAS response for {rsid}: {snp}")
    phenos = data.get("phenos", [])
    if len(phenos) == 0:
        raise ValueError(f"No UKBB PheWAS phenotypes returned for {rsid}: {snp}")
    df = pl.DataFrame(phenos)
    print(f"[TRACKING] Number of phenotypes in TOPMed-imputed UKBB PheWAS for {rsid}: {df.height}")
    required_cols = ["beta", "af", "sebeta", "pval", "phenocode", "phenostring", "category"]
    missing_cols = [col for col in required_cols if col not in df.columns]
    if missing_cols:
        raise ValueError(f"Missing UKBB PheWAS columns for {rsid}: {missing_cols}")
    df = (
        df
        .select(required_cols)
        .rename({
            "beta": "BETA",
            "sebeta": "SE",
            "af": "FRQ",
            "pval": "P",
            "phenocode": "PHENOCODE",
            "phenostring": "PHENOSTRING",
            "category": "CATEGORY"
        })
        .with_columns(
            pl.col("BETA").cast(pl.Float64),
            pl.col("SE").cast(pl.Float64),
            pl.col("FRQ").cast(pl.Float64),
            pl.col("P").cast(pl.Float64),
            pl.col("PHENOCODE").cast(pl.Utf8),
            pl.col("PHENOSTRING").cast(pl.Utf8),
            pl.col("CATEGORY").cast(pl.Utf8)
        )
        .drop_nulls(["BETA", "SE", "P", "PHENOCODE"])
    )
    print(f"[TRACKING] Dataframe shape after keeping only relevant stats: {df.shape}")
    return df


def resolve_ukbb_variant(chromosome: str, position: int, A1: str, A2: str, rsid: str):
    snp_first_order = f"{chromosome}-{position}-{A1}-{A2}"
    snp_second_order = f"{chromosome}-{position}-{A2}-{A1}"
    for snp, ukb_ref, ukb_alt in [(snp_first_order, A1, A2), (snp_second_order, A2, A1)]:
        try:
            response = requests.get(f"https://pheweb.org/UKB-TOPMed/api/variant/{snp}", timeout=60)
            if response.status_code != 200:
                continue
            data = response.json()
            if not isinstance(data, dict):
                continue
            if len(data.get("phenos", [])) > 0:
                return snp, ukb_ref, ukb_alt
        except (requests.RequestException, ValueError):
            continue
    print(f"[TRACKING] UKBB variant could not be resolved for {rsid}. Tried {snp_first_order} and {snp_second_order}...")
    return None, None, None


def phewas_mr_on_ukbb(pqtl_dataset: str, pheno_id: str, local_results_dir: str = "results"):
    coloc_file = paths.coloc_out(pqtl_dataset, pheno_id, local_results_dir)
    df_coloc = pl.read_csv(coloc_file, separator="\t")
    if "protein_id" in df_coloc.columns:
        df_coloc = df_coloc.rename({"protein_id": "protein"})
    if "protein" not in df_coloc.columns:
        raise ValueError(f"Could not find protein_id or protein in COLOC file: {coloc_file}")
    pp_h4_col = next((col for col in ["PP.H4.abf", "PP.H4", "pp_h4", "PPH4"] if col in df_coloc.columns), None)
    if pp_h4_col is not None:
        df_coloc = df_coloc.filter(pl.col(pp_h4_col).cast(pl.Float64, strict=False) >= coloc_threshold)
    else:
        print(f"[TRACKING] No PP.H4 column found in {coloc_file}; using every protein in the COLOC file...")
    compelling_targets = df_coloc.select(pl.col("protein").cast(pl.Utf8)).drop_nulls().unique(maintain_order=True)

    # PWCoCo (conditional coloc) - complementary to standard COLOC above, not a
    # replacement (see project_pwcoco_wiring memory): a target that colocalises
    # under EITHER method should reach UKB PheWAS, so PWCoCo-passing proteins are
    # unioned in below. pwcoco_out() may not exist - PWCoCo runs non-fatally in
    # local.py/hpc.py, so a failed or not-yet-run PWCoCo step must not break this.
    pwcoco_file = paths.pwcoco_out(pqtl_dataset, pheno_id, local_results_dir)
    if Path(pwcoco_file).exists():
        df_pwcoco = pl.read_csv(pwcoco_file, separator="\t")
        if "H4" in df_pwcoco.columns and "protein" in df_pwcoco.columns:
            pwcoco_targets = (
                df_pwcoco
                .filter(pl.col("H4").cast(pl.Float64, strict=False) >= coloc_threshold)
                .select(pl.col("protein").cast(pl.Utf8))
                .drop_nulls()
                .unique(maintain_order=True)
            )
            compelling_targets = compelling_targets.vstack(pwcoco_targets).unique(maintain_order=True)
    else:
        print(f"[TRACKING] No PWCoCo output found at {pwcoco_file}; using standard COLOC targets only...")

    # UKB PheWAS is a FALLBACK, only run for a target when NONE of its retained
    # cis-MR instruments were available in FinnGen - it is not a second parallel
    # PheWAS source for every compelling target. Requires bin/phewas_cis_pqtls.py
    # to have already run and written its coverage manifest.
    finngen_coverage_file = paths.phewas_finngen_coverage_out(pqtl_dataset, pheno_id, local_results_dir)
    if not Path(finngen_coverage_file).exists():
        raise FileNotFoundError(
            f"FinnGen PheWAS coverage manifest not found: {finngen_coverage_file}. "
            "UKB PheWAS fallback requires bin/phewas_cis_pqtls.py to run first."
        )
    df_finngen_coverage = pl.read_csv(finngen_coverage_file, separator="\t")
    finngen_uncovered_targets = (
        df_finngen_coverage
        .filter(~pl.col("finngen_covered"))
        .select(pl.col("protein").cast(pl.Utf8))
        .unique(maintain_order=True)
    )
    n_before_fallback_filter = compelling_targets.height
    compelling_targets = compelling_targets.join(finngen_uncovered_targets, on="protein", how="inner")
    print(
        f"[TRACKING] UKB PheWAS fallback targets (zero FinnGen instrument coverage): "
        f"{compelling_targets.height}/{n_before_fallback_filter}..."
    )
    if compelling_targets.height == 0:
        print("[TRACKING] Every compelling target already has FinnGen PheWAS coverage - nothing to run in UKB...")
        return

    instruments_file = paths.mr_instruments_out(pqtl_dataset, pheno_id, local_results_dir)
    df_instruments = pl.read_csv(instruments_file, separator="\t")
    required_instrument_cols = [
        "protein", "pqtl_dataset", "outcome_trait", "SNP",
        "effect_allele.exposure", "other_allele.exposure",
        "beta.exposure", "se.exposure", "pval.exposure",
        "effect_allele.outcome", "other_allele.outcome",
        "beta.outcome", "used_in_mr"
    ]
    missing_cols = [col for col in required_instrument_cols if col not in df_instruments.columns]
    if missing_cols:
        raise ValueError(f"Missing cis-MR instrument columns in {instruments_file}: {missing_cols}")
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
            pl.col("used_in_mr").cast(pl.Utf8).str.to_uppercase().eq("TRUE").alias("used_in_mr")
        )
        .filter(
            (pl.col("pqtl_dataset") == pqtl_dataset) &
            (pl.col("outcome_trait") == pheno_id) &
            pl.col("used_in_mr")
        )
        .join(compelling_targets, on="protein", how="inner")
        .unique(subset=["protein", "SNP"], keep="first")
    )
    temp_dir = f"./work/PheWAS_UKBB/{pqtl_dataset}_{pheno_id}"
    os.makedirs(temp_dir, exist_ok=True)
    phewas_ukbb_out_file = paths.phewas_ukbb_out(pqtl_dataset, pheno_id, local_results_dir)
    os.makedirs(phewas_ukbb_out_file.parent, exist_ok=True)
    print(f"[TRACKING] Number of compelling COLOC targets for PheWAS: {compelling_targets.height}...")
    print(f"[TRACKING] Number of retained cis-MR instrument rows: {df_instruments.height}...")
    results = []
    mr = PyTwoSampleMR()
    # run every protein assay / aptamer within the COLOC file
    for target_row in compelling_targets.iter_rows(named=True):
        protein = target_row["protein"]
        df_target_instruments = df_instruments.filter(pl.col("protein") == protein)
        if "instrument_rank" in df_target_instruments.columns:
            df_target_instruments = df_target_instruments.sort("instrument_rank")
        instrument_rsids_original = df_target_instruments["SNP"].to_list()
        n_cis_mr_instruments = len(instrument_rsids_original)
        print(f"[TRACKING] Number of cis-MR instruments for {protein}: {n_cis_mr_instruments}...")
        print(f"[TRACKING] cis-MR SNPs for {protein}: {instrument_rsids_original}...")
        if n_cis_mr_instruments == 0:
            print(f"[TRACKING] No retained cis-MR instruments for {protein}...")
            continue
        # recover pos and chr
        pqtl_file = Path(f"./dat/cis_regions/{pqtl_dataset}/{protein}/pqtl.parquet")
        if not pqtl_file.exists():
            print(f"[TRACKING] Missing pQTL file for {protein}: {pqtl_file}...")
            continue
        df_pqtl = pl.read_parquet(pqtl_file)
        chr_col = next((col for col in ["CHR", "chr", "chrom", "chromosome"] if col in df_pqtl.columns), None)
        bp_col = next((col for col in ["BP", "bp", "POS", "pos", "position"] if col in df_pqtl.columns), None)
        if chr_col is None or bp_col is None or "SNP" not in df_pqtl.columns:
            print(f"[TRACKING] Could not find SNP/CHR/BP columns in {pqtl_file}...")
            continue
        df_pqtl_positions = (
            df_pqtl
            .select(
                pl.col("SNP").cast(pl.Utf8),
                pl.col(chr_col).cast(pl.Utf8).alias("CHR"),
                pl.col(bp_col).cast(pl.Int64).alias("BP")
            )
            .filter(pl.col("SNP").is_in(instrument_rsids_original))
            .unique(subset=["SNP"], keep="first")
        )
        df_target_instruments = df_target_instruments.join(df_pqtl_positions, on="SNP", how="inner")
        # check if snp == missing
        if df_target_instruments.height != len(instrument_rsids_original):
            recovered_snps = set(df_target_instruments["SNP"].to_list())
            missing_position_snps = sorted(set(instrument_rsids_original) - recovered_snps)
            print(f"[TRACKING] Could not recover CHR/BP for all cis-MR instruments for {protein}: {missing_position_snps}...")
        if df_target_instruments.height == 0:
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
            chromosome = str(instrument_row["CHR"]).lower().replace("chr", "")
            position = int(instrument_row["BP"])
            # make A1 the pQTL exposure-increasing (protein abundance-increasing)
            # allele - this guarantees beta_exposure > 0 relative to A1, matching
            # bin/phewas_cis_pqtls.py's FinnGen convention exactly, so beta_mr is
            # on the same axis (effect of higher protein abundance) whichever
            # source produced it. This used to align to the AD risk allele
            # instead, which put FinnGen and UKB beta_mr on different axes and
            # made the additional-indication / adverse-effect classification
            # (which compares beta_mr's sign to the primary protein->AD beta)
            # silently wrong for UKB-sourced hits.
            beta_exposure_original = beta_exposure
            if beta_exposure_original > 0:
                A1 = exposure_effect_allele
                A2 = exposure_other_allele
                beta_exposure = beta_exposure_original
                exposure_A1_flipped = False
            elif beta_exposure_original < 0:
                A1 = exposure_other_allele
                A2 = exposure_effect_allele
                beta_exposure = -beta_exposure_original
                exposure_A1_flipped = True
            else:
                print(f"[TRACKING] pQTL beta is zero for {rsid}...")
                continue
            # align AD beta to the pQTL risk allele A1
            if ad_effect_allele == A1 and ad_other_allele == A2:
                ad_A1_flipped = False
                beta_ad = beta_ad_original
            elif ad_effect_allele == A2 and ad_other_allele == A1:
                beta_ad = -beta_ad_original
                ad_A1_flipped = True
            else:
                print(f"[TRACKING] AD alleles {ad_effect_allele}/{ad_other_allele} do not match pQTL risk alleles {A1}/{A2} for {rsid}...")
                continue
            print(f"[TRACKING] pQTL risk alignment for {rsid}: A1={A1}, A2={A2}, pQTL beta={beta_exposure}, original AD beta={beta_ad_original}, A1-aligned AD beta={beta_ad}, AD flipped={ad_A1_flipped}...")
            snp, ukb_ref, ukb_alt = resolve_ukbb_variant(chromosome, position, A1, A2, rsid)
            if snp is None:
                continue
            print(f"[TRACKING] UKBB REF/ALT resolved for {rsid}: REF={ukb_ref}, ALT={ukb_alt}...")
            print(f"[TRACKING] Running UKBB PheWAS for {protein}: {rsid} ({snp})...")
            # query to phewas db and clean
            try:
                df_phewas = grab_phewas_info(snp=snp, rsid=rsid)
            except (requests.RequestException, ValueError) as error:
                print(f"[TRACKING] UKBB PheWAS query failed for {protein}: {rsid}: {error}...")
                continue
            df_phewas.write_csv(os.path.join(temp_dir, f"{protein}_{rsid}_raw_hits.csv"))
            # align UKBB beta directly to the pQTL risk (protein-increasing) allele A1
            # PheWeb beta is relative to ALT
            if ukb_alt == A1 and ukb_ref == A2:
                phewas_A1_flipped = False
            elif ukb_alt == A2 and ukb_ref == A1:
                df_phewas = df_phewas.with_columns((-pl.col("BETA")).alias("BETA"))
                phewas_A1_flipped = True
            else:
                print(f"[TRACKING] UKBB REF/ALT {ukb_ref}/{ukb_alt} do not match AD risk alleles {A1}/{A2} for {rsid}...")
                continue
            print(f"[TRACKING] UKBB alignment for {rsid}: UKBB ALT={ukb_alt}, AD risk A1={A1}, flipped={phewas_A1_flipped}...")
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
                "ukb_ref": ukb_ref,
                "ukb_alt": ukb_alt,
                "beta_exposure_original": beta_exposure_original,
                "beta_exposure": beta_exposure,
                "exposure_A1_flipped": exposure_A1_flipped,
                "se_exposure": se_exposure,
                "p_exposure": p_exposure,
                "phewas_A1_flipped": phewas_A1_flipped,
                "df_phewas": df_phewas
            }
        # define how many of the original cis-MR instruments were found in UKBB
        n_instruments_original = n_cis_mr_instruments
        n_instruments_available = len(signal_phewas)
        missing_snps = sorted(set(instrument_rsids_original) - set(signal_phewas.keys()))
        instrument_completeness = n_instruments_available / n_instruments_original
        if n_instruments_available == 0:
            print(f"[TRACKING] No cis-MR instruments found in UKBB for {protein}...")
            continue
        if n_instruments_available < n_instruments_original:
            print(f"[TRACKING] Using available UKBB instruments for {protein}: {n_instruments_available}/{n_instruments_original}; missing={missing_snps}...")
        protein_results = []
        # THIS IS IF SNP == 1 FOR AVAILABLE cis-MR INSTRUMENTS AT A GIVEN TARGET
        if n_instruments_available == 1:
            signal = list(signal_phewas.values())[0]
            for pheno in signal["df_phewas"].iter_rows(named=True):
                res = PheWAS(B_Y=pheno["BETA"], SE_Y=pheno["SE"], B_X=signal["beta_exposure"], SE_X=signal["se_exposure"])
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
                    "ukb_ref": signal["ukb_ref"],
                    "ukb_alt": signal["ukb_alt"],
                    "beta_exposure_original": str(signal["beta_exposure_original"]),
                    "beta_exposure": str(signal["beta_exposure"]),
                    "exposure_A1_flipped": str(signal["exposure_A1_flipped"]),
                    "se_exposure": str(signal["se_exposure"]),
                    "p_exposure": str(signal["p_exposure"]),
                    "phewas_A1_flipped": str(signal["phewas_A1_flipped"]),
                    "beta_mr": res["wald_ratio"],
                    "se_mr": res["se_wald_ratio"],
                    "p_mr": res["P_nominal"],
                    "beta_phewas": str(pheno["BETA"]),
                    "se_phewas": str(pheno["SE"]),
                    "p_phewas": str(pheno["P"]),
                    "PHENOCODE": pheno["PHENOCODE"],
                    "PHENOSTRING": pheno["PHENOSTRING"],
                    "CATEGORY": pheno["CATEGORY"]
                })
        # THIS IS IF SNP > 1 FOR AVAILABLE cis-MR INSTRUMENTS AT A GIVEN TARGET
        else:
            common_phenocodes = None
            # only run IVW on phenotypes which appear for every SNP
            for signal in signal_phewas.values():
                phenocodes = set(signal["df_phewas"]["PHENOCODE"].to_list())
                common_phenocodes = phenocodes if common_phenocodes is None else common_phenocodes.intersection(phenocodes)
            common_phenocodes = sorted(common_phenocodes or [])
            print(f"[TRACKING] Number of common PheWAS outcomes across {len(signal_phewas)} cis-MR instruments for {protein}: {len(common_phenocodes)}...")
            for phenocode in common_phenocodes:
                exposure_rows = []
                outcome_rows = []
                instrument_rsids = []
                instrument_variants = []
                A1_values = []
                A2_values = []
                ad_effect_alleles = []
                ad_other_alleles = []
                beta_ad_original_values = []
                beta_ad_values = []
                ad_A1_flipped_values = []
                exposure_effect_alleles = []
                exposure_other_alleles = []
                ukb_refs = []
                ukb_alts = []
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
                    pheno = signal["df_phewas"].filter(pl.col("PHENOCODE") == phenocode).unique(subset=["PHENOCODE"], keep="first")
                    if pheno.height == 0:
                        continue
                    pheno = pheno.row(0, named=True)
                    exposure_rows.append({"SNP": signal["rsid"], "BETA": signal["beta_exposure"], "SE": signal["se_exposure"]})
                    outcome_rows.append({"SNP": signal["rsid"], "BETA": pheno["BETA"], "SE": pheno["SE"]})
                    instrument_rsids.append(signal["rsid"])
                    instrument_variants.append(signal["snp"])
                    A1_values.append(signal["A1"])
                    A2_values.append(signal["A2"])
                    ad_effect_alleles.append(signal["ad_effect_allele_original"])
                    ad_other_alleles.append(signal["ad_other_allele_original"])
                    beta_ad_original_values.append(signal["beta_ad_original"])
                    beta_ad_values.append(signal["beta_ad"])
                    ad_A1_flipped_values.append(signal["ad_A1_flipped"])
                    exposure_effect_alleles.append(signal["exposure_effect_allele"])
                    exposure_other_alleles.append(signal["exposure_other_allele"])
                    ukb_refs.append(signal["ukb_ref"])
                    ukb_alts.append(signal["ukb_alt"])
                    beta_exposure_original_values.append(signal["beta_exposure_original"])
                    beta_exposure_values.append(signal["beta_exposure"])
                    exposure_A1_flipped_values.append(signal["exposure_A1_flipped"])
                    se_exposure_values.append(signal["se_exposure"])
                    p_exposure_values.append(signal["p_exposure"])
                    phewas_flipped_values.append(signal["phewas_A1_flipped"])
                    beta_phewas_values.append(pheno["BETA"])
                    se_phewas_values.append(pheno["SE"])
                    p_phewas_values.append(pheno["P"])
                    if phenotype_info is None:
                        phenotype_info = {"PHENOCODE": pheno["PHENOCODE"], "PHENOSTRING": pheno["PHENOSTRING"], "CATEGORY": pheno["CATEGORY"]}
                exposure_df = pl.DataFrame(exposure_rows)
                outcome_df = pl.DataFrame(outcome_rows)
                if exposure_df.height < 2:
                    print(f"[TRACKING] Incomplete IVW instrument set for {protein} and {phenocode}...")
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
                    "snp": ",".join(instrument_variants),
                    "A1": ",".join(A1_values),
                    "A2": ",".join(A2_values),
                    "ad_effect_allele_original": ",".join(ad_effect_alleles),
                    "ad_other_allele_original": ",".join(ad_other_alleles),
                    "beta_ad_original": ",".join(str(value) for value in beta_ad_original_values),
                    "beta_ad": ",".join(str(value) for value in beta_ad_values),
                    "ad_A1_flipped": ",".join(str(value) for value in ad_A1_flipped_values),
                    "exposure_effect_allele": ",".join(exposure_effect_alleles),
                    "exposure_other_allele": ",".join(exposure_other_alleles),
                    "ukb_ref": ",".join(ukb_refs),
                    "ukb_alt": ",".join(ukb_alts),
                    "beta_exposure_original": ",".join(str(value) for value in beta_exposure_original_values),
                    "beta_exposure": ",".join(str(value) for value in beta_exposure_values),
                    "exposure_A1_flipped": ",".join(str(value) for value in exposure_A1_flipped_values),
                    "se_exposure": ",".join(str(value) for value in se_exposure_values),
                    "p_exposure": ",".join(str(value) for value in p_exposure_values),
                    "phewas_A1_flipped": ",".join(str(value) for value in phewas_flipped_values),
                    "beta_mr": ivw_result["BETA"],
                    "se_mr": ivw_result["SE"],
                    "p_mr": ivw_result["P"],
                    "beta_phewas": ",".join(str(value) for value in beta_phewas_values),
                    "se_phewas": ",".join(str(value) for value in se_phewas_values),
                    "p_phewas": ",".join(str(value) for value in p_phewas_values),
                    "PHENOCODE": phenotype_info["PHENOCODE"],
                    "PHENOSTRING": phenotype_info["PHENOSTRING"],
                    "CATEGORY": phenotype_info["CATEGORY"]
                })
        if len(protein_results) == 0:
            print(f"[TRACKING] No Wald ratio / IVW estimates generated for {protein}...")
            continue
        df_protein_results = pl.DataFrame(protein_results, infer_schema_length=None).with_columns(pl.col("p_mr").cast(pl.Float64))
        # Bonferroni correct across the endpoints actually tested for THIS protein
        # (this protein's own row count here), not a fixed global constant.
        n_endpoints_tested = df_protein_results.height
        df_protein_results = df_protein_results.with_columns(
            pl.lit(n_endpoints_tested).alias("n_endpoints_tested"),
            pl.min_horizontal(pl.col("p_mr") * n_endpoints_tested, pl.lit(1.0)).alias("p_bonferroni"),
            (pl.col("p_mr") < (0.05 / n_endpoints_tested)).alias("bonferroni_significant")
        )
        results.extend(df_protein_results.to_dicts())
    if len(results) == 0:
        print("[TRACKING] No PheWAS associations were generated...")
        return
    df_results = pl.DataFrame(results, infer_schema_length=None).sort(["protein", "p_mr"])
    output_file = phewas_ukbb_out_file
    df_results.write_csv(output_file, separator="\t")
    print(f"[TRACKING] PheWAS completed: {df_results.height} associations saved to {output_file}...")


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--pheno_id", required=True)
    p.add_argument("--pqtl_dataset", required=True)
    p.add_argument("--local_results_dir", default="results")
    args = p.parse_args()
    phewas_mr_on_ukbb(pheno_id=args.pheno_id, pqtl_dataset=args.pqtl_dataset, local_results_dir=args.local_results_dir)


if __name__ == "__main__":
    main()