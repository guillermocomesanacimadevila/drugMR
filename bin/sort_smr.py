#!/usr/bin/env python3
import argparse
import polars as pl
from pathlib import Path
import subprocess
import shutil
from drugmr import SMR
from drugmr import paths
import pandas as pd
import os
from statsmodels.stats.multitest import fdrcorrection

# ------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------
# MAIN TO DO'S
# -> SLAP FUNCTION 3 (MAYBE 1 ONTO A DIFFERNT SCRIPT -> MAY CRASH SMR IF != RESULTS)
# -> CONSEQUENTLY UPDATE drugmr/local.py and drugmr/hpc.py
# ------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------

# need probs another function with mediators
# -> load cis-MR results for pQTL dataset X  -> same for coloc -> check whether gene which passes coloc thresh and MR estimate FDR
# -> check on all cells - save the SMR output for that in 1 or > 1 cells onto results/...

# ADD TO DOCKER IMAGE
# - STATSMODELS
# - SMR PACKAGE (as part of ref/)

# COLOC and cis-MR filter for promising targets
def extract_promising_targets(pqtl_dataset: str, pheno_id: str, local_results_dir: str = "results"):
    # extract stuff from here
    cis_mr_res = paths.mr_out(pqtl_dataset, pheno_id, local_results_dir)
    cis_mr_df = pl.read_csv(cis_mr_res, separator="\t")
    coloc_res = paths.coloc_out(pqtl_dataset, pheno_id, local_results_dir)
    coloc_df = pl.read_csv(coloc_res, separator="\t")

    # base parameters
    # wald_fdr = 0.05
    # ivw_fdr = 0.05
    # cochran_q = 0.05
    # coloc_thresh = 0.75 #### subject to change

    wald_hits = []
    ivw_hits = []
    coloc_hits = []

    for row in cis_mr_df.iter_rows(named=True):
        # separate where n_instruments == 1 or > 1
        n_instruments = row["n_instruments"]
        wald_fdr = row["Wald_FDR_q"]
        ivw_fdr = row["IVW_FDR_q"]
        cochran_q = row["Q_pval"]
        if n_instruments == 1 and wald_fdr is not None and wald_fdr < 0.05:
            wald_hits.append(row["protein"])
        elif (n_instruments > 1 and ivw_fdr is not None and ivw_fdr < 0.05 and cochran_q is not None and cochran_q > 0.05):
            ivw_hits.append(row["protein"])

    for row in coloc_df.iter_rows(named=True):
        pp4 = row["PP.H4.abf"]
        if pp4 is not None and pp4 > 0.75:
            coloc_hits.append(row["protein_id"])

    # compile final hits
    mr_hits = wald_hits + ivw_hits
    final_hits = [i for i in mr_hits if i in coloc_hits]
    print(f"[TRACKING] {len(wald_hits)} Wald ratio hits found")
    print(f"[TRACKING] {len(ivw_hits)} IVW hits found")
    print(f"[TRACKING] {len(coloc_hits)} coloc hits found")
    print(f"[TRACKING] {len(final_hits)} final promising targets found")
    return final_hits


# shared GWAS -> .ma prep (SMR GWAS format == SNP A1 A2 freq b se p N)
# used by both the bulk and single-cell SMR runs - column logic is dataset-agnostic
def prepare_smr_gwas(sumstats: str, pheno_id: str):
    # temp dir to store .ma file per pheno
    temp_dir = "./work/SMR/"
    temp_dir = Path(temp_dir)
    os.makedirs(temp_dir, exist_ok=True)

    # Store sumstats (temp) within work/ as a .ma file
    # which then delete after all cell types are ran -> just a temp file
    # Store sumstats (temp) within work/ as a .ma file
    # which then delete after all cell types are ran -> just a temp file
    df = pl.read_csv(sumstats, separator="\t")
    n_before = df.height
    # rename and save as .ma
    # SMR GWAS format == SNP A1 A2 freq b se p N
    df = (
        df
        .select([
            pl.col("SNP"),
            pl.col("A1"),
            pl.col("A2"),
            pl.col("FRQ").alias("freq"),
            pl.col("BETA").alias("b"),
            pl.col("SE").alias("se"),
            pl.col("P").alias("p"),
            pl.col("N")
        ])
        .filter(
            pl.col("SNP").is_not_null(),
            pl.col("A1").is_not_null(),
            pl.col("A2").is_not_null(),
            pl.col("freq").is_not_null(),
            pl.col("b").is_not_null(),
            pl.col("se").is_not_null(),
            pl.col("p").is_not_null(),
            pl.col("N").is_not_null(),
            ~pl.col("SNP").str.contains(","),
            ~pl.col("SNP").str.contains(";"),
            ~pl.col("SNP").str.contains(" ")
        )
        .with_columns(
            pl.col("N").round(0).cast(pl.Int64)
        )
    )

    print(f"[TRACKING] Removed {n_before - df.height} invalid / incomplete GWAS rows for SMR")
    print(f"[TRACKING] {df.height} GWAS variants retained for SMR")
    temp_sumstats = temp_dir / f"{pheno_id}.ma"
    df.write_csv(temp_sumstats, separator="\t")
    return temp_sumstats


# some SMR output rows have literal "NA" text (e.g. HEIDI not computable for that probe)
# instead of an empty field - treat it as null so numeric columns like p_HEIDI stay numeric
# rather than getting inferred as a string column
def read_smr_tsv(path: Path):
    return pl.read_csv(path, separator="\t", null_values=["NA"])


# shared FDR correction on a single .smr output file (kept as-is across every probe
# in the file, not just the promising targets) - used by both bulk and single-cell runs
def fdr_correct_smr_file(f: Path, pheno_id: str, label: str):
    smr_df = read_smr_tsv(f)
    # now FDR correct p_SMR and add q_SMR col
    if "p_SMR" not in smr_df.columns:
        print(f"[CONCERN] p_SMR not found in {f.name}")
        return

    # this runs every time (even when the smr binary itself was skipped), so drop any
    # q_SMR (or a stray q_SMR_right left over from an earlier un-idempotent rerun)
    # before recomputing, otherwise the join below collides with itself on rerun
    stale_cols = [col for col in smr_df.columns if col == "q_SMR" or col.startswith("q_SMR_right")]
    if stale_cols:
        smr_df = smr_df.drop(stale_cols)

    valid_p = (
        smr_df
        .filter(pl.col("p_SMR").is_not_null())
        .get_column("p_SMR")
        .to_numpy()
    )

    if len(valid_p) == 0:
        print(f"[CONCERN] No valid p_SMR values found in {f.name}")
        return

    _, q_values = fdrcorrection(
        valid_p,
        alpha=0.05,
        method="indep"
    )
    q_df = (
        smr_df
        .filter(pl.col("p_SMR").is_not_null())
        .select(
            pl.int_range(pl.len()).alias("_row_id")
        )
        .with_columns(
            pl.Series("q_SMR", q_values)
        )
    )
    smr_df = (
        smr_df
        .with_row_index("_row_id")
        .join(q_df, on="_row_id", how="left")
        .drop("_row_id")
    )

    # overwrite same SMR result with q_SMR added
    smr_df.write_csv(f, separator="\t")
    print(f"[TRACKING] FDR corrected SMR results saved for {pheno_id} in {label}")


# SMR reports A1 arbitrarily (whichever allele the bfile / eQTL reference happened to
# assign) - flip so A1 is always the AD risk allele (b_GWAS > 0), same convention as
# the "make A1 the GWAS risk allele" block in bin/compile_cis_hit_info.py.
# b_SMR itself is untouched (ratio of the two negated betas is unchanged).
def align_to_risk_allele(df: pl.DataFrame):
    if "b_GWAS" not in df.columns or "A1" not in df.columns or "A2" not in df.columns:
        print("[CONCERN] Cannot align alleles to the AD risk allele - missing A1 / A2 / b_GWAS")
        return df

    flip = pl.col("b_GWAS") < 0

    df = df.with_columns(
        pl.when(flip).then(pl.col("A2")).otherwise(pl.col("A1")).alias("A1"),
        pl.when(flip).then(pl.col("A1")).otherwise(pl.col("A2")).alias("A2"),
        pl.when(flip).then(-pl.col("b_GWAS")).otherwise(pl.col("b_GWAS")).alias("b_GWAS")
    )

    if "b_eQTL" in df.columns:
        df = df.with_columns(
            pl.when(flip).then(-pl.col("b_eQTL")).otherwise(pl.col("b_eQTL")).alias("b_eQTL")
        )

    if "Freq" in df.columns:
        df = df.with_columns(
            pl.when(flip).then(1 - pl.col("Freq")).otherwise(pl.col("Freq")).alias("Freq")
        )

    return df


# for single-cell targets, replace the SMR-reported eQTL beta with the value from the
# original per-cell-type eQTL file (source of truth) - re-signed to match the SMR
# file's own A1 (never touching A1/A2 themselves, so b_SMR / align_to_risk_allele stay
# internally consistent). Not used for bulk - those files are ingested as pre-computed.
def pull_original_sc_eqtl_beta(target_smr: pl.DataFrame, eqtl_dataset: str, cell: str):
    if target_smr.height == 0:
        return target_smr

    parquet_path = Path(f"./dat/sc-eQTL/{eqtl_dataset}/{cell}.parquet")
    if not parquet_path.exists():
        print(f"[CONCERN] Original eQTL file not found for {cell}: {parquet_path}")
        return target_smr

    needed = target_smr.select(["probeID", "topSNP"]).unique()
    needed_genes = needed.get_column("probeID").to_list()
    needed_snps = needed.get_column("topSNP").to_list()

    original = (
        pl.scan_parquet(parquet_path)
        .filter(pl.col("GENE").is_in(needed_genes) & pl.col("SNP").is_in(needed_snps))
        .select(["GENE", "SNP", "EA", "BETA"])
        .collect()
        .rename({"GENE": "probeID", "SNP": "topSNP", "EA": "orig_EA", "BETA": "orig_BETA"})
        .join(needed, on=["probeID", "topSNP"], how="inner")
        .unique(subset=["probeID", "topSNP"])
    )

    if original.height == 0:
        print(f"[CONCERN] No matching original eQTL rows found for {cell}")
        return target_smr

    target_smr = target_smr.join(original, on=["probeID", "topSNP"], how="left")

    target_smr = target_smr.with_columns(
        pl.when(pl.col("orig_EA") == pl.col("A1")).then(pl.col("orig_BETA"))
        .when(pl.col("orig_EA") == pl.col("A2")).then(-pl.col("orig_BETA"))
        .otherwise(None)
        .alias("b_eqtl_from_source")
    )

    n_mismatched = target_smr.filter(
        pl.col("orig_EA").is_not_null() & pl.col("b_eqtl_from_source").is_null()
    ).height
    if n_mismatched > 0:
        print(
            f"[CONCERN] {n_mismatched} row(s) in {cell} had an allele mismatch between "
            f"the SMR A1/A2 and the original eQTL file - b_eQTL left as SMR-reported"
        )

    target_smr = target_smr.with_columns(
        pl.when(pl.col("b_eqtl_from_source").is_not_null())
        .then(pl.col("b_eqtl_from_source"))
        .otherwise(pl.col("b_eQTL"))
        .alias("b_eQTL")
    ).drop(["orig_EA", "orig_BETA", "b_eqtl_from_source"])

    # b_SMR is a ratio of b_GWAS/b_eQTL - recompute so it stays consistent with the
    # (possibly now different-signed) b_eQTL rather than leaving a stale ratio
    target_smr = target_smr.with_columns(
        pl.when(pl.col("b_eQTL") != 0)
        .then(pl.col("b_GWAS") / pl.col("b_eQTL"))
        .otherwise(pl.col("b_SMR"))
        .alias("b_SMR")
    )

    return target_smr


def run_single_cell_smr(pqtl_dataset: str, eqtl_dataset: str, pheno_id: str, sumstats: str, ref_bfile: str, maf: float, local_results_dir: str = "results"):
    ref_bfile = Path(ref_bfile)
    eqtl_temp = eqtl_dataset.lower()

    if eqtl_temp == "singlebrain":
        temp_sumstats = prepare_smr_gwas(sumstats, pheno_id)
        cell_types = ["Ast", "Ext", "MG", "OD", "OPC", "End", "IN"]
        # cell_types = ["MG"]
        eqtls = f"./dat/sc-eQTL/{eqtl_dataset}/SMR_ready"
        eqtls = Path(eqtls)
        for cell in cell_types:
            cell_dir = eqtls / cell
            besd_file = cell_dir / f"{cell}.besd"
            esi_file = cell_dir / f"{cell}.esi"
            epi_file = cell_dir / f"{cell}.epi"

            if not cell_dir.exists():
                print(f"[CONCERN] Cell type directory {cell_dir} not found")
                continue

            if not besd_file.exists():
                print(f"[CONCERN] {besd_file} not found")
                continue

            if not esi_file.exists():
                print(f"[CONCERN] {esi_file} not found")
                continue

            if not epi_file.exists():
                print(f"[CONCERN] {epi_file} not found")
                continue

            print(f"[TRACKING] Cell type {cell} found!")

            # use prefix without .besd / .esi / .epi for SMR

            #####
            #####
            #####
            beqtl_summary = cell_dir / cell
            #####
            #####
            #####

            # check whether SMR has already been ran for trait X in cell type Y
            # single-cell results live under results/SMR/sc/... (sibling to results/SMR/bulk/...)
            smr_res = paths.smr_raw_dir(f"sc/{eqtl_dataset}/{cell}", pheno_id, local_results_dir)
            existing_smr = [f for f in smr_res.glob(f"*{pheno_id}*.smr") if f.stat().st_size > 0]

            if len(existing_smr) > 0:
                print(f"[TRACKING] SMR already completed for {pheno_id} in {cell} - skipping SMR")
            else:
                SMR(
                    pheno_id=pheno_id,
                    sumstats=temp_sumstats,
                    ref_bfile=ref_bfile,
                    beqtl_summary=beqtl_summary,
                    eqtl_dataset=f"sc/{eqtl_dataset}/{cell}",
                    peqtl_smr=5.0e-8, #### change to default one
                    peqtl_heidi=1.57e-3, ###### change to real default
                    thread_num=8,
                    maf=maf
                )

            # load SMR results
            # saving into out_dir 1 results file per cell type for trait X
            # results/SMR/sc/SingleBrain/{cell}/{pheno_id}/...
            smr_res = paths.smr_raw_dir(f"sc/{eqtl_dataset}/{cell}", pheno_id, local_results_dir)
            for f in smr_res.glob("*.smr"):
                if pheno_id in f.name:
                    fdr_correct_smr_file(f, pheno_id, cell)

        # delete temp GWAS .ma only after all cell types are ran
        if temp_sumstats.exists():
            temp_sumstats.unlink()

        hits = extract_promising_targets(pheno_id=pheno_id, pqtl_dataset=pqtl_dataset, local_results_dir=local_results_dir)
        # now extract all of the SMR data from the results for each cell type pertaining to those targets and store as a dataframe in results/SMR/dataset
        # rows == 1 SMR result for target X on cell-type Y
        # so 7 cell types x X targets in terms of rows
        all_target_smr = []
        for cell in cell_types:
            smr_res = paths.smr_raw_dir(f"sc/{eqtl_dataset}/{cell}", pheno_id, local_results_dir)
            for f in smr_res.glob("*.smr"):
                if pheno_id not in f.name:
                    continue
                smr_df = read_smr_tsv(f)

                # SMR usually calls the gene / probe column Probe
                # match the gene part of GENE_UNIPROT targets to the SMR Probe column
                if "Gene" not in smr_df.columns:
                    print(f"[CONCERN] Gene column not found in {f.name}")
                    continue

                target_map = {
                    target.split("_")[0]: target for target in hits
                }

                target_genes = list(target_map.keys())

                target_smr = (
                    smr_df
                    .filter(pl.col("Gene").is_in(target_genes))
                    .with_columns(
                        pl.col("Gene").replace(target_map).alias("protein"),
                        pl.lit(cell).alias("cell_type"),
                        pl.lit("single_cell").alias("data_type"),
                        pl.lit(pheno_id).alias("phenotype"),
                        pl.lit(eqtl_dataset).alias("eqtl_dataset"),
                        pl.lit(pqtl_dataset).alias("pqtl_dataset")
                    )
                )

                target_smr = pull_original_sc_eqtl_beta(target_smr, eqtl_dataset, cell)

                if target_smr.height > 0:
                    all_target_smr.append(target_smr)

        out_file = paths.smr_sc_out(pqtl_dataset, pheno_id, eqtl_dataset, local_results_dir)
        os.makedirs(out_file.parent, exist_ok=True)

        if len(all_target_smr) > 0:
            final_smr_df = pl.concat(all_target_smr, how="diagonal_relaxed")
            final_smr_df = align_to_risk_allele(final_smr_df)
            final_smr_df.write_csv(out_file, separator="\t")
            print(f"[TRACKING] Compiled promising target SMR results saved to {out_file}")
        else:
            print(f"[CONCERN] No SMR results found for the promising {pqtl_dataset} targets")


# raw bulk eQTL besd/esi/epi under dat/bulk-eQTL/{eqtl_dataset}/ are split per chromosome
# (unlike single-cell's genome-wide SMR_ready sets), so 1 tissue == 1 smr call per
# chromosome rather than 1 call overall. Returns {tissue_label: {chr_num: besd_prefix}},
# or None if eqtl_dataset has no raw dat/bulk-eQTL directory (e.g. it's only available as
# a pre-computed .smr, like eQTLGen) - in which case there's nothing to freshly run.
def bulk_tissue_prefixes(eqtl_dataset: str):
    base_dir = Path(f"./dat/bulk-eQTL/{eqtl_dataset}")
    if not base_dir.exists():
        return None

    tissues = {}

    if eqtl_dataset == "GTEx_v10":
        for tissue_dir in sorted(base_dir.iterdir()):
            if not tissue_dir.is_dir():
                continue
            tissue = tissue_dir.name
            label = f"GTEx_{tissue}_v10"
            chr_prefixes = {}
            for chr_num in range(1, 23):
                prefix = tissue_dir / f"{tissue}.v10.eQTL.cis_qtl_pairs.{chr_num}"
                if Path(f"{prefix}.besd").exists():
                    chr_prefixes[chr_num] = prefix
            if chr_prefixes:
                tissues[label] = chr_prefixes
    elif eqtl_dataset == "MetaBrain":
        label = "BrainMeta"
        chr_prefixes = {}
        for chr_num in range(1, 23):
            prefix = base_dir / f"BrainMeta_cis_eQTL_chr{chr_num}"
            if Path(f"{prefix}.besd").exists():
                chr_prefixes[chr_num] = prefix
        if chr_prefixes:
            tissues[label] = chr_prefixes
    else:
        print(f"[CONCERN] No raw bulk eQTL layout known for {eqtl_dataset} - cannot run SMR from scratch")
        return None

    return tissues


# run the smr binary against the raw dat/bulk-eQTL besd/esi/epi files, 1 chromosome at a
# time per tissue (each chromosome's besd/esi/epi only covers that chromosome's SNPs/probes,
# so per-chromosome .smr outputs concatenate cleanly into 1 genome-wide file - no besd
# merging needed). Skips a tissue entirely if its final .smr already exists (same
# idempotency convention as run_single_cell_smr). A literal qtl_name column is stamped onto
# the concatenated output and it's FDR-corrected via the same helper single-cell SMR uses,
# so the result is indistinguishable from a "pre-computed" bulk file to ingest_bulk_smr.
def run_bulk_smr(pqtl_dataset: str, eqtl_dataset: str, pheno_id: str, sumstats: str, ref_bfile: str, maf: float, local_results_dir: str = "results"):
    tissues = bulk_tissue_prefixes(eqtl_dataset)

    if not tissues:
        print(f"[TRACKING] No raw bulk eQTL besd/esi/epi files found under ./dat/bulk-eQTL/{eqtl_dataset} - nothing to run SMR on")
        return

    ref_bfile = Path(ref_bfile)
    temp_sumstats = prepare_smr_gwas(sumstats, pheno_id)

    for label, chr_prefixes in tissues.items():
        bulk_dir = paths.smr_bulk_dir(eqtl_dataset, local_results_dir)
        final_dir = bulk_dir / f"eQTL_{label}"
        final_file = final_dir / f"{pheno_id}_{label}.smr"

        # pre-existing bulk results don't all follow the same directory convention (GTEx's
        # are nested under an eQTL_<tissue>/ subdirectory, MetaBrain's legacy file sits flat
        # at the top level) - search for either rather than assuming one, so a legacy file
        # is recognised as "already done" instead of silently duplicated alongside a fresh one
        existing = [
            f for f in bulk_dir.rglob(f"*{pheno_id}*{label}*.smr")
            if f.stat().st_size > 0
        ] if bulk_dir.exists() else []

        if existing:
            print(f"[TRACKING] Bulk SMR already completed for {pheno_id} in {label} - skipping SMR ({existing[0]})")
            continue

        print(f"[TRACKING] Running SMR for {label} ({len(chr_prefixes)} chromosome(s))")
        chr_smr_files = []
        for chr_num, prefix in sorted(chr_prefixes.items()):
            SMR(
                pheno_id=pheno_id,
                sumstats=temp_sumstats,
                ref_bfile=ref_bfile,
                beqtl_summary=prefix,
                eqtl_dataset=f"bulk_raw/{eqtl_dataset}/{label}/chr{chr_num}",
                peqtl_smr=5.0e-8, #### change to default one
                peqtl_heidi=1.57e-3, ###### change to real default
                thread_num=8,
                maf=maf
            )
            chr_out = Path(f"{paths.smr_raw_prefix(f'bulk_raw/{eqtl_dataset}/{label}/chr{chr_num}', pheno_id, local_results_dir)}.smr")
            if chr_out.exists() and chr_out.stat().st_size > 0:
                chr_smr_files.append(chr_out)
            else:
                print(f"[CONCERN] SMR produced no output for {label} chr{chr_num}")

        if not chr_smr_files:
            print(f"[CONCERN] No per-chromosome SMR output generated for {label} - skipping")
            continue

        combined = pl.concat(
            [read_smr_tsv(f) for f in chr_smr_files],
            how="diagonal_relaxed"
        ).with_columns(pl.lit(f"eQTL_{label}").alias("qtl_name"))

        final_dir.mkdir(parents=True, exist_ok=True)
        combined.write_csv(final_file, separator="\t")
        fdr_correct_smr_file(final_file, pheno_id, label)
        print(f"[DONE] Saved genome-wide bulk SMR results for {label}: {final_file}")

        scratch_dir = Path(local_results_dir) / "SMR" / "bulk_raw" / eqtl_dataset / label
        if scratch_dir.exists():
            shutil.rmtree(scratch_dir)

    if temp_sumstats.exists():
        temp_sumstats.unlink()


# bulk eQTL SMR (eQTLGen / MetaBrain / GTEx_v10) - loads whatever .smr-shaped files sit
# under results/SMR/bulk/{eqtl_dataset}/, harmonises the column names onto the same shape
# single-cell uses, and extracts the promising targets. Those files are either produced
# just above by run_bulk_smr() (from the raw dat/bulk-eQTL besd/esi/epi, when available) or
# pre-computed elsewhere (e.g. eQTLGen, which has no raw dat/bulk-eQTL directory). GTEx_v10
# is tissue-resolved (1 file per tissue via rglob, same idea as single-cell's per-cell
# loop); eQTLGen / MetaBrain are flat (1 file for the dataset).
def ingest_bulk_smr(pqtl_dataset: str, eqtl_dataset: str, pheno_id: str, local_results_dir: str = "results"):
    bulk_dir = paths.smr_bulk_dir(eqtl_dataset, local_results_dir)
    smr_files = sorted(bulk_dir.rglob(f"*{pheno_id}*.smr"))

    if len(smr_files) == 0:
        print(f"[CONCERN] No pre-computed bulk SMR files found under {bulk_dir} for {pheno_id}")
        return

    print(f"[TRACKING] Found {len(smr_files)} pre-computed bulk SMR file(s) for {eqtl_dataset}")

    hits = extract_promising_targets(pheno_id=pheno_id, pqtl_dataset=pqtl_dataset, local_results_dir=local_results_dir)
    target_map = {
        target.split("_")[0]: target for target in hits
    }
    target_genes = list(target_map.keys())

    all_target_smr = []
    for f in smr_files:
        smr_df = read_smr_tsv(f)

        # pre-computed bulk files use "index" (gene symbol) + "p_SMR_FDR" where the raw
        # smr binary output (single-cell) uses "Gene" + a freshly-computed "q_SMR"
        rename_map = {}
        if "index" in smr_df.columns and "Gene" not in smr_df.columns:
            rename_map["index"] = "Gene"
        if "p_SMR_FDR" in smr_df.columns and "q_SMR" not in smr_df.columns:
            rename_map["p_SMR_FDR"] = "q_SMR"
        if rename_map:
            smr_df = smr_df.rename(rename_map)

        if "Gene" not in smr_df.columns:
            print(f"[CONCERN] Gene column not found in {f.name}")
            continue

        # sub-dataset / tissue label straight from qtl_name
        # (e.g. eQTL_GTEx_Brain_Cortex_v10 -> GTEx_Brain_Cortex_v10, eQTL_eQTLGen -> eQTLGen)
        label = eqtl_dataset
        if "qtl_name" in smr_df.columns:
            qtl_names = smr_df.get_column("qtl_name").drop_nulls().unique().to_list()
            if len(qtl_names) == 1:
                label = qtl_names[0].removeprefix("eQTL_")

        target_smr = (
            smr_df
            .filter(pl.col("Gene").is_in(target_genes))
            .with_columns(
                pl.col("Gene").replace(target_map).alias("protein"),
                pl.lit(label).alias("cell_type"),
                pl.lit("bulk").alias("data_type"),
                pl.lit(pheno_id).alias("phenotype"),
                pl.lit(eqtl_dataset).alias("eqtl_dataset"),
                pl.lit(pqtl_dataset).alias("pqtl_dataset")
            )
        )

        if target_smr.height > 0:
            all_target_smr.append(target_smr)

    out_file = paths.smr_bulk_out(pqtl_dataset, pheno_id, eqtl_dataset, local_results_dir)
    os.makedirs(out_file.parent, exist_ok=True)

    if len(all_target_smr) > 0:
        final_smr_df = pl.concat(all_target_smr, how="diagonal_relaxed")
        final_smr_df = align_to_risk_allele(final_smr_df)
        final_smr_df.write_csv(out_file, separator="\t")
        print(f"[TRACKING] Compiled promising target SMR results saved to {out_file}")
    else:
        print(f"[CONCERN] No pre-computed bulk SMR results found for the promising {pqtl_dataset} targets")


def compile_multi_omics_targets(pheno_id: str, pqtl_dataset: str, eqtl_dataset: str, eqtl_mode: str, local_results_dir: str = "results"):
    # single-cell results live under results/SMR/sc/..., bulk under results/SMR/bulk/...
    targets_path = (
        paths.smr_sc_out(pqtl_dataset, pheno_id, eqtl_dataset, local_results_dir)
        if eqtl_mode == "single_cell"
        else paths.smr_bulk_out(pqtl_dataset, pheno_id, eqtl_dataset, local_results_dir)
    )
    df = read_smr_tsv(targets_path)

    # HEIDI CUT OFF = 0.01
    # SMR CUT OFF = 0.05
    # q_SMR - P_HEIDI

    heidi_col = "P_HEIDI"
    if "p_HEIDI" in df.columns:
        heidi_col = "p_HEIDI"

    final_targets_df = (
        df
        .filter(
            pl.col("q_SMR").is_not_null(),
            pl.col(heidi_col).is_not_null(),
            pl.col("q_SMR") < 0.05,
            pl.col(heidi_col) > 0.01
        )
        .sort(["protein", "cell_type", "q_SMR"])
    )

    if final_targets_df.height == 0:
        print(f"[CONCERN] No drug targets passed cis-MR (pQTLs) + COLOC + {eqtl_dataset} eQTL SMR")
        return []

    # canonical combined output (bulk + single-cell hits together) that the dashboard reads
    # upsert: drop any stale rows for this eqtl_dataset, then append the fresh ones,
    # so bulk and single-cell runs (in either order) compose instead of overwriting each other
    combined_file = paths.smr_final_targets_out(pqtl_dataset, pheno_id, local_results_dir)
    os.makedirs(combined_file.parent, exist_ok=True)

    if combined_file.exists() and combined_file.stat().st_size > 0:
        existing_df = read_smr_tsv(combined_file)
        existing_df = existing_df.filter(pl.col("eqtl_dataset") != eqtl_dataset)
        combined_df = pl.concat([existing_df, final_targets_df], how="diagonal_relaxed")
    else:
        combined_df = final_targets_df

    combined_df.write_csv(combined_file, separator="\t")

    # unique targets only for the next steps
    targets = (
        final_targets_df
        .select("protein")
        .unique()
        .sort("protein")
        .get_column("protein")
        .to_list()
    )

    print(f"[TRACKING] {final_targets_df.height} target x cell-type SMR hits found for {eqtl_dataset}")
    print(f"[TRACKING] {len(targets)} unique drug targets passed cis-MR (pQTLs) + COLOC + {eqtl_dataset} eQTL SMR")
    print(f"[TRACKING] Combined multi-omics target results saved to {combined_file}")
    print(f"[TRACKING] Drug targets: {targets}")
    return targets


# THEN
# -> For each prioritised target
# -> Check original cis-region (matched with GWAS)
# -> Match cis-region with sc-eQTL for cell type X
# -> RUN MOLOC / pairwise coloc

# sumstats: str, ref_bfile: str, maf: float
def main():
    p = argparse.ArgumentParser()
    p.add_argument("--pheno_id", required=True)
    p.add_argument("--sumstats", required=True)
    p.add_argument("--pqtl_dataset", required=True)
    p.add_argument("--eqtl_dataset", required=True)
    p.add_argument("--eqtl_mode", required=True, choices=["bulk", "single_cell"])
    p.add_argument("--ref_bfile", required=True)
    p.add_argument("--maf", type=float, default=0.01)
    p.add_argument("--local_results_dir", default="results")
    args = p.parse_args()

    # running SMR (bulk or single-cell, depending on --eqtl_mode)
    # bulk: run fresh SMR from raw dat/bulk-eQTL besd/esi/epi where available (no-op if
    # already done, or if this eqtl_dataset only exists as a pre-computed .smr, e.g.
    # eQTLGen), then ingest whatever ends up under results/SMR/bulk/{eqtl_dataset}/
    if args.eqtl_mode == "bulk":
        run_bulk_smr(
            pqtl_dataset=args.pqtl_dataset,
            eqtl_dataset=args.eqtl_dataset,
            pheno_id=args.pheno_id,
            sumstats=args.sumstats,
            ref_bfile=args.ref_bfile,
            maf=args.maf,
            local_results_dir=args.local_results_dir
        )
        ingest_bulk_smr(
            pqtl_dataset=args.pqtl_dataset,
            eqtl_dataset=args.eqtl_dataset,
            pheno_id=args.pheno_id,
            local_results_dir=args.local_results_dir
        )
    else:
        run_single_cell_smr(
            pqtl_dataset=args.pqtl_dataset,
            eqtl_dataset=args.eqtl_dataset,
            pheno_id=args.pheno_id,
            sumstats=args.sumstats,
            maf=args.maf,
            ref_bfile=args.ref_bfile,
            local_results_dir=args.local_results_dir
        )

    # final hits
    compile_multi_omics_targets(
        pheno_id=args.pheno_id,
        pqtl_dataset=args.pqtl_dataset,
        eqtl_dataset=args.eqtl_dataset,
        eqtl_mode=args.eqtl_mode,
        local_results_dir=args.local_results_dir
    )

if __name__ == "__main__":
    main()
