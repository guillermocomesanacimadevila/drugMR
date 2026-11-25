import argparse
import os
import shutil
import subprocess
from pathlib import Path

import pandas as pd
import polars as pl
from drugmr.methods.cojo import COJO
from drugmr.methods.twosamplemr import PyTwoSampleMR
from drugmr.paths import (
    coloc_out,
    locus_conditioning_out,
    locus_conditioning_work_dir,
    mr_instruments_out,
    mr_out,
    pwcoco_out,
)
from drugmr.utils import (
    compute_ld_to_lead,
    detect_relevant_snps,
    extract_coloc_or_pwcoco_targets,
    grab_cis_mr_hits,
    ld_clumping,
)


def new_record(protein):

    """ Every target gets the same columns, whatever happens to it """

    return {
        "protein": protein,
        "status": None,
        "n_gws_snps": 0,
        "n_gws_snps_not_in_panel": 0,
        "top_gws_snp_in_panel": None,
        "own_signal_leads": "",
        "rivals": "",
        "flagged_rivals": "",
        "n_instruments_used": 0,
        "mr_beta": None,
        "mr_se": None,
        "mr_p": None,
        "mr_beta_cond": None,
        "mr_se_cond": None,
        "mr_p_cond": None,
        "mr_beta_ratio": None,
        "coloc_pp4": None,
        "coloc_pp4_cond": None,
        "coloc_pp3": None,
        "coloc_pp3_cond": None,
    }

def condition_on_nearby_targets(
        pqtl_dataset,
        pheno_id,
        ref_bfile,
        n_cases,
        n_controls,
        local_results_dir="results",
        cis_regions_dir=None,
        wald_fdr_q=0.05,
        ivw_fdr_q=0.05,
        cochran_q_pval=0.05,
        egger_intercept_pval_min=0,
        min_instruments_for_ivw=3,
        pp4_threshold=0.7,
        p1=1e-4,
        p2=1e-4,
        p12=1e-5,
        gws_pval=5e-8,
        clump_r2=0.1,
        clump_kb=1000,
        own_signal_r2=0.8,
        flag_r2=0.3,
        cojo_collinear=0.99,
        target_build="GRCh38"
):


    """
    Inputs:
    -- hits (proteins that passes cis-MR & coloc and/or PWCoCo)
    -- outcome_gwas (.parquet file for that protein's cis-region from outcome GWAS)
    -- exposure_pqtl (the same as above but the pQTL)
    -- ref_bfile (self-explanatory)
    -- window -> fixed at 1Mb

    Logic:
    1. -> for protein in hits:
        -> extract cis-instrument SNPs
        -> grab outcome_gwas cis-region

    2. Grab SNPs where P<5E-08 (from the GWAS cis-region for protein X)
    -> if len(snps) == 0:
        continue
        
    3. Clump ALL GWS SNPs (instruments included) -> (r2 0.1, 1000 kb) -> leads
    for each lead:
       if lead in instruments            -> drop
       r2 = max r2 with any instrument
       if r2 >= 0.8                      -> drop (X's own signal)
       elif 0.3 <= r2 < 0.8              -> rivals + flagged
       else                              -> rivals
    if rivals empty -> next protein

    4. Convert outcome GWAS cis-region to .ma and store in tmp/
    -> Add leads in rivals: list -> txt cond-file
    -> Check also if rivals fall due to != coverage with reference panel
    -> output -> conditioned outcome GWAS

    5. Re-run (with unchanged pQTL cis-region)
    -> Re-run cis-MR for protein X
    -> Re-run coloc for protein X
    -> Re-run PWCoCo for protein X
    """

    cis_mr_res = mr_out(pheno_id=pheno_id, pqtl_dataset=pqtl_dataset, out_dir=local_results_dir)
    coloc_res = coloc_out(pheno_id=pheno_id, pqtl_dataset=pqtl_dataset, out_dir=local_results_dir)
    pwcoco_res = pwcoco_out(pheno_id=pheno_id, pqtl_dataset=pqtl_dataset, out_dir=local_results_dir)
    instruments_df = pl.read_csv(mr_instruments_out(pheno_id=pheno_id, pqtl_dataset=pqtl_dataset, out_dir=local_results_dir), separator="\t")

    # grab cis-MR hits first!
    cis_mr_hits = grab_cis_mr_hits(
        csv_file=cis_mr_res,
        wald_fdr_q=wald_fdr_q,
        ivw_fdr_q=ivw_fdr_q,
        cochran_q_pval=cochran_q_pval,
        egger_intercept_pval_min=egger_intercept_pval_min,
        min_instruments_for_ivw=min_instruments_for_ivw
    )

    coloc_hits = extract_coloc_or_pwcoco_targets(
        coloc_csv_file=coloc_res,
        pwcoco_csv_file=pwcoco_res,
        pp4_thresh=pp4_threshold,
        method=("coloc",)
    )

    pwcoco_hits = extract_coloc_or_pwcoco_targets(
        coloc_csv_file=coloc_res,
        pwcoco_csv_file=pwcoco_res,
        pp4_thresh=pp4_threshold,
        method=("pwcoco",)
    )

    pgen_targets = [p for p in cis_mr_hits if p in coloc_hits or p in pwcoco_hits]

    # SNPs present in the LD reference panel
    # a GWS AD SNP missing here can never become a lead, so never a rival
    # .bim columns: CHR SNP CM BP A1 A2, only SNP needed (read as text)
    bim = pl.read_csv(
        f"{ref_bfile}.bim",
        separator="\t",
        has_header=False,
        columns=[1],
        new_columns=["SNP"],
        infer_schema=False,
    )
    panel_snps = set(bim["SNP"].to_list())
    del bim

    # APOE window, same coordinates as bin/qc_gwas.py --remove_apoe
    if target_build == "GRCh37":
        apoe_start = 44000000
        apoe_end = 46500000
    elif target_build == "GRCh38":
        apoe_start = 43495848
        apoe_end = 45996742
    else:
        raise ValueError(f"Unsupported target_build: {target_build}")

    # original cis-MR results (instrument count per protein)
    original_mr = pl.read_csv(cis_mr_res, separator="\t", infer_schema_length=None)

    # original coloc results (baseline PP4/PP3 per protein)
    original_coloc = pl.read_csv(coloc_res, separator="\t", infer_schema_length=None)

    work_root = locus_conditioning_work_dir(local_results_dir)

    records = []

    for protein in pgen_targets:

        # save in record
        record = new_record(protein)

        work = work_root / protein
        os.makedirs(work, exist_ok=True)

        if cis_regions_dir:
            cis_regions = Path(cis_regions_dir) / str(protein)
        else:
            cis_regions = Path(f"./dat/cis_regions/{pqtl_dataset}/{protein}")
        gwas_parquet = pl.read_parquet(cis_regions / "gwas.parquet")

        # whole cis window sits on one chromosome
        chromosome = gwas_parquet["CHR"][0]

        # APOE: too strong to condition on with 1000G LD, so not attempted
        window_start = gwas_parquet["BP"].min()
        window_end = gwas_parquet["BP"].max()
        on_chr19 = str(chromosome).removeprefix("chr") == "19"
        if on_chr19 and window_start <= apoe_end and window_end >= apoe_start:
            print(f"{protein}: cis window overlaps APOE, not conditioned")
            record["status"] = "apoe_locus"
            records.append(record)
            continue

        # extract cis_instruments
        cis_instruments = []
        for row in instruments_df.iter_rows(named=True):
            p = row["protein"]
            i = row["SNP"]
            if str(p) == str(protein) and row["mr_keep"]:
                cis_instruments.append(str(i))

        # continue here!
        candidate_snps = detect_relevant_snps(
            dataframe=gwas_parquet,
            p_col="P",
            p_thresh=gws_pval
        )

        if len(candidate_snps) == 0:
            print(f"No GWS SNPs in protein {protein} cis-region for the outcome {pheno_id} GWAS...")
            record["status"] = "no_gws_signal"
            records.append(record)
            continue

        record["n_gws_snps"] = len(candidate_snps)

        # how much of the GWS signal the panel can see
        n_missing = 0
        for snp in candidate_snps["SNP"].to_list():
            if str(snp) not in panel_snps:
                n_missing = n_missing + 1
        record["n_gws_snps_not_in_panel"] = n_missing

        # the strongest AD SNP in the window must be visible, or a rival can hide
        top_gws_snp = str(candidate_snps.sort("P")["SNP"][0])
        record["top_gws_snp_in_panel"] = top_gws_snp in panel_snps
        if top_gws_snp not in panel_snps:
            print(f"[CONCERN] {protein}: top GWS SNP {top_gws_snp} not in LD panel")

        # LD of each cis-MR instrument against the window 
        instrument_ld = {} # dictionary -> mapping lead snp: df with LD matrix
        cannot_asess = False
        for instrument in cis_instruments:
            try:
                ld = compute_ld_to_lead(
                    ref_bfile=ref_bfile,
                    lead_snp=instrument,
                    chromosome=chromosome,
                    out_file=work / f"ld_{instrument}.parquet"
                )
                instrument_ld[instrument] = ld
            except RuntimeError as error:
                print(f"[CONCERN] {protein}: no LD for instrument {instrument}: {error}")
                cannot_asess = True
        if cannot_asess:
            print(f"[CONCERN] {protein}: cannot assess, skipping")
            record["status"] = "ld_failed"
            records.append(record)
            continue

        # save candidate SNPs to temp work directory
        gws_file = work / "gwas_snps.txt"
        candidate_snps.select(["SNP", "P"]).write_csv(gws_file, separator="\t")

        # plink writes no .clumped when nothing clumps: never read last run's file
        clumped_file = work / "clump.clumped"
        if clumped_file.exists():
            clumped_file.unlink()

        # LD clumping here!!!!!!
        try:
            ld_clumping(
                ref_bfile=ref_bfile,
                sumstats=gws_file,
                clump_pval=gws_pval,
                clump_r2=clump_r2,
                clump_kb=clump_kb,
                out=work / "clump",
                chromosome=chromosome
            )
        except subprocess.CalledProcessError as error:
            print(f"[CONCERN] {protein}: Cooked! Clumping failed: {error}")
            record["status"] = "clump_failed"
            records.append(record)
            continue

        # read clumped file
        if not clumped_file.exists():
            print(f"{protein}: no leads after clumping")
            record["status"] = "no_leads"
            records.append(record)
            continue

        clumped = pd.read_csv(clumped_file, sep=r"\s+")
        leads = clumped["SNP"].tolist()
        clumped = pl.from_pandas(clumped) # back to polars!

        # leads without the cis_instruments
        final_leads = [str(i) for i in leads if i not in cis_instruments]

        if len(final_leads) == 0:
            print(f"{protein}: only X's own signal, no rivals")
            record["status"] = "own_signal_only"
            records.append(record)
            continue

        # split leads by max r2 with any cis-instrument
        own_signal_leads = [] # leads in high LD with an instrument (X's own signal)
        rivals = [] # leads to condition on (cond file)
        flagged = [] # rivals partly correlated with X's signal
        for lead in final_leads: # for each lead!!!
            # check SNP in highest LD with each i
            r2_values = []
            for instrument in instrument_ld:
                ld = instrument_ld[instrument]
                match = ld.filter(pl.col("SNP") == lead)
                if len(match) > 0:
                    r2_values.append(match["R2"][0])

            # highest r2 = closest to X's signal; not in any LD table -> 0
            if len(r2_values) > 0:
                max_r2 = max(r2_values)
            else:
                max_r2 = 0

            # X's own signal: do not condition on it
            if max_r2 >= own_signal_r2:
                own_signal_leads.append(lead)
                print(f"{protein}: {lead} r2 {max_r2:.2f} with an instrument, X's own signal")

            # partly correlated: condition on it, but flag it
            elif max_r2 >= flag_r2:
                rivals.append(lead)
                flagged.append(lead)
                print(f"{protein}: {lead} r2 {max_r2:.2f}, rival (flagged)")

            # clean rival
            else:
                rivals.append(lead)
                print(f"{protein}: {lead} r2 {max_r2:.2f}, rival")

        # store the split whatever happens next
        record["own_signal_leads"] = ",".join(own_signal_leads)
        record["rivals"] = ",".join(rivals)
        record["flagged_rivals"] = ",".join(flagged)

        if len(rivals) == 0:
            print(f"{protein}: no rival AD signal")
            record["status"] = "own_signal_only"
            records.append(record)
            continue

        # now apply drugmr/methods/cojo.py functionality!!!
        cojo = COJO(
            ref_bfile=ref_bfile,
            work_dir=work,
            collinear=cojo_collinear
        )

        ma = cojo.to_ma(sumstats=gwas_parquet)

        # condition the AD GWAS on rivals
        # Chr  SNP  bp  refA  freq  b  se  p  n  freq_geno  bC  bC_se  pC
        conditioned = cojo.conditional(
            ma=ma,
            condition_on=rivals,
            chromosome=chromosome,
            tag=protein
        )

        if conditioned is None:
            print(f"[CONCERN] {protein}: COJO gave no result, skipping")
            record["status"] = "cojo_failed"
            records.append(record)
            continue

        # re-run cis-MR (IVW + Cochran's Q and wald ratio)
        mr = PyTwoSampleMR()

        # instruments_df 
        # from this df -> filter SNP/s + pQTL beta + se + n_instruments, etc...
        # add one new col -> f"{pheno_id}_cond..."
        # check n_instruments -> wald - IVW - IVW + Cochran's Q?

        exposure = instruments_df.filter((pl.col("protein") == protein) & pl.col("mr_keep"))

        # filter from conditioned file the required cols whereby instruments in exposure and rename cols
        aligned_rows = []
        for row in exposure.iter_rows(named=True):
            snp = str(row["SNP"])

            # match on the conditioned df
            matched = conditioned.filter(pl.col("SNP") == snp)
            if len(matched) == 0:
                print(f"[CONCERN] {protein}: instrument {snp} missing from COJO output")
                continue

            # capture alleles
            ref_allele = str(matched["refA"][0]).upper()
            effect_allele = str(row["effect_allele.exposure"]).upper()
            other_allele = str(row["other_allele.exposure"]).upper()

            if ref_allele == effect_allele:
                sign = 1
            elif ref_allele == other_allele:
                sign = -1
            else:
                print(f"[CONCERN] {protein}: {snp} alleles do not match, dropped")
                continue

            b = sign * matched["b"][0]
            bC = sign * matched["bC"][0]

            # aligned unconditioned b must equal the harmonised outcome beta from cis-MR
            if abs(b - row["beta.outcome"]) > 1e-4:
                print(f"[CONCERN] {protein}: {snp} b {b:.4f} != beta.outcome {row['beta.outcome']:.4f}")

            aligned_rows.append({
                "SNP": snp,
                "b": b,
                "se": matched["se"][0],
                "bC": bC,
                "bC_se": matched["bC_se"][0],
            })

        if len(aligned_rows) == 0:
            print(f"[CONCERN] {protein}: no instruments left after COJO")
            record["status"] = "no_instruments_after_cojo"
            records.append(record)
            continue

        outcome = pl.DataFrame(aligned_rows)
        record["n_instruments_used"] = len(aligned_rows)

        # now cis-MR!!!!
        original_row = original_mr.filter(pl.col("protein") == protein)

        # instruments lost to COJO or alignment: fewer SNPs than the original cis-MR
        if len(aligned_rows) < original_row["n_instruments"][0]:
            print(f"[CONCERN] {protein}: {len(aligned_rows)} of {original_row['n_instruments'][0]} instruments left after COJO")

        # unconditioned baseline with the same Python MR, so SEs are comparable
        raw_outcome = outcome.select(
            pl.col("SNP"),
            pl.col("b").alias("beta.outcome"),
            pl.col("se").alias("se.outcome"),
        )

        if len(raw_outcome) == 1:
            raw_result = mr.wald_ratio(
                exposure_df=exposure,
                outcome_df=raw_outcome,
                exposure_snp_col="SNP",
                exposure_beta_col="beta.exposure",
                exposure_se_col="se.exposure",
                outcome_snp_col="SNP",
                outcome_beta_col="beta.outcome",
                outcome_se_col="se.outcome",
            )
        else:
            raw_result, raw_snp_results = mr.IVW(
                exposure_df=exposure,
                outcome_df=raw_outcome,
                exposure_snp_col="SNP",
                exposure_beta_col="beta.exposure",
                exposure_se_col="se.exposure",
                outcome_snp_col="SNP",
                outcome_beta_col="beta.outcome",
                outcome_se_col="se.outcome",
            )

        record["mr_beta"] = raw_result["BETA"][0]
        record["mr_se"] = raw_result["SE"][0]
        record["mr_p"] = raw_result["P"][0]

        # conditioned sumstats in the same column names as the exposure file
        cond_outcome = outcome.select(
            pl.col("SNP"),
            pl.col("bC").alias("beta.outcome"),
            pl.col("bC_se").alias("se.outcome"),
        )

        if len(cond_outcome) == 1:
            result = mr.wald_ratio(
                exposure_df=exposure,
                outcome_df=cond_outcome,
                exposure_snp_col="SNP",
                exposure_beta_col="beta.exposure",
                exposure_se_col="se.exposure",
                outcome_snp_col="SNP",
                outcome_beta_col="beta.outcome",
                outcome_se_col="se.outcome",
            )
        else:
            result, snp_results = mr.IVW(
                exposure_df=exposure,
                outcome_df=cond_outcome,
                exposure_snp_col="SNP",
                exposure_beta_col="beta.exposure",
                exposure_se_col="se.exposure",
                outcome_snp_col="SNP",
                outcome_beta_col="beta.outcome",
                outcome_se_col="se.outcome",
            ) 

        record["mr_beta_cond"] = result["BETA"][0]
        record["mr_se_cond"] = result["SE"][0]
        record["mr_p_cond"] = result["P"][0]

        # near 1: conditioning on the rivals did not move X's estimate
        record["mr_beta_ratio"] = record["mr_beta_cond"] / record["mr_beta"]


        # COLOC
        original_coloc_row = original_coloc.filter(pl.col("protein") == protein)
        if len(original_coloc_row) > 0: # which it should NOT
            record["coloc_pp4"] = original_coloc_row["PP.H4.abf"][0]
            record["coloc_pp3"] = original_coloc_row["PP.H3.abf"][0]

        # grab relevant cols from conditioned df
        cond_cols = conditioned.select([
            pl.col("SNP").cast(pl.String),
            pl.col("bC"),
            pl.col("bC_se"),
            pl.col("pC"),
        ])

        gwas_cond = gwas_parquet.with_columns(pl.col("SNP").cast(pl.String))
        gwas_cond = gwas_cond.join(cond_cols, on="SNP", how="inner")

        # rivals themselves and SNPs COJO left out come back as NA
        gwas_cond = gwas_cond.drop_nulls(subset=["bC", "bC_se"])

        # swap in conditioned sumstats (bC is on GWAS A1, no flip)
        gwas_cond = gwas_cond.drop(["BETA", "SE", "P"])
        gwas_cond = gwas_cond.rename({"bC": "BETA", "bC_se": "SE", "pC": "P"})

        gwas_cond_file = work / "gwas_conditioned.parquet"
        gwas_cond.write_parquet(gwas_cond_file)

        # coloc.R on the conditioned GWAS, pQTL unchanged, written into the work folder
        coloc_script = str(Path(os.environ.get("PYTHONPATH", ".")) / "bin" / "coloc.R")
        coloc_root = work / "conditioned"
        cmd = [
            "Rscript", coloc_script,
            pqtl_dataset, str(protein), pheno_id,
            str(gwas_cond_file), str(cis_regions / "pqtl.parquet"),
            str(n_cases), str(n_controls),
            str(coloc_root), str(pp4_threshold),
            str(p1), str(p2), str(p12),
        ]
        try:
            subprocess.run(cmd, check=True)

            # coloc.R writes <coloc_root>/coloc/<pheno_id>_<protein>_coloc.tsv
            coloc_file = coloc_root / "coloc" / f"{pheno_id}_{protein}_coloc.tsv"
            coloc_row = pl.read_csv(coloc_file, separator="\t")
            record["coloc_pp4_cond"] = coloc_row["PP.H4.abf"][0]
            record["coloc_pp3_cond"] = coloc_row["PP.H3.abf"][0]
        except subprocess.CalledProcessError as error:
            print(f"[CONCERN] {protein}: conditioned coloc failed: {error}")

        record["status"] = "conditioned"
        records.append(record)

    # one row per target, whatever happened to it
    if len(records) > 0:
        table = pl.DataFrame(records)

        # hit lists come from sets, so fix the row order for reproducible output
        table = table.sort("protein")
    else:
        # no targets: still write the header so readers do not break
        column_names = list(new_record("").keys())
        table = pl.DataFrame(schema=column_names)
    out_file = locus_conditioning_out(pqtl_dataset=pqtl_dataset, pheno_id=pheno_id, out_dir=local_results_dir)
    os.makedirs(out_file.parent, exist_ok=True)
    table.write_csv(out_file, separator="\t")
    print(f"[DONE] Saved locus conditioning table: {out_file}")

    # scratch (ld, clump, .ma, .cma.cojo, conditioned parquet, coloc) not needed past this step
    if work_root.exists():
        shutil.rmtree(work_root)

    return table


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--pqtl_dataset", required=True)
    p.add_argument("--pheno_id", required=True)
    p.add_argument("--ref_bfile", required=True)
    p.add_argument("--n_cases", required=True, type=int)
    p.add_argument("--n_controls", required=True, type=int)
    p.add_argument("--local_results_dir", default="results")
    p.add_argument("--cis_regions_dir", default=None)
    p.add_argument("--wald_fdr_q", type=float, default=0.05)
    p.add_argument("--ivw_fdr_q", type=float, default=0.05)
    p.add_argument("--cochran_q_pval", type=float, default=0.05)
    p.add_argument("--egger_intercept_pval_min", type=float, default=0)
    p.add_argument("--min_instruments_for_ivw", type=int, default=3)
    p.add_argument("--pp4_threshold", type=float, default=0.7)
    p.add_argument("--p1", type=float, default=1e-4)
    p.add_argument("--p2", type=float, default=1e-4)
    p.add_argument("--p12", type=float, default=1e-5)
    p.add_argument("--gws_pval", type=float, default=5e-8)
    p.add_argument("--clump_r2", type=float, default=0.1)
    p.add_argument("--clump_kb", type=int, default=1000)
    p.add_argument("--own_signal_r2", type=float, default=0.8)
    p.add_argument("--flag_r2", type=float, default=0.3)
    p.add_argument("--cojo_collinear", type=float, default=0.99)
    p.add_argument("--target_build", default="GRCh38", choices=["GRCh37", "GRCh38"])
    args = p.parse_args()

    condition_on_nearby_targets(
        pqtl_dataset=args.pqtl_dataset,
        pheno_id=args.pheno_id,
        ref_bfile=args.ref_bfile,
        n_cases=args.n_cases,
        n_controls=args.n_controls,
        local_results_dir=args.local_results_dir,
        cis_regions_dir=args.cis_regions_dir,
        wald_fdr_q=args.wald_fdr_q,
        ivw_fdr_q=args.ivw_fdr_q,
        cochran_q_pval=args.cochran_q_pval,
        egger_intercept_pval_min=args.egger_intercept_pval_min,
        min_instruments_for_ivw=args.min_instruments_for_ivw,
        pp4_threshold=args.pp4_threshold,
        p1=args.p1,
        p2=args.p2,
        p12=args.p12,
        gws_pval=args.gws_pval,
        clump_r2=args.clump_r2,
        clump_kb=args.clump_kb,
        own_signal_r2=args.own_signal_r2,
        flag_r2=args.flag_r2,
        cojo_collinear=args.cojo_collinear,
        target_build=args.target_build,
    )

if __name__ == "__main__":
    main()