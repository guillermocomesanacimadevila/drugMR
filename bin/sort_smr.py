#!/usr/bin/env python3
import argparse
import os
import re
import shutil
from pathlib import Path

import polars as pl
from statsmodels.stats.multitest import fdrcorrection

from drugmr import paths
from drugmr.smr import SMRUtils

_smr_utils = SMRUtils(manifest_path=paths.DEFAULT_QTL_MANIFEST_PATH)  # ncbi_ref_path not needed - BESD is always pre-built for registered datasets, ETL-from-scratch never triggers


def resolve_qtl_type(qtl_dataset: str) -> str:
    try:
        row = _smr_utils.qtl_manifest.get_row(qtl_dataset.lower())
    except ValueError:
        return "eqtl"
    return row.get("qtl_type") or "eqtl"

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
def extract_promising_targets(
        pqtl_dataset: str,
        pheno_id: str,
        local_results_dir: str = "results",
        coloc_file: str | None = None,
        wald_fdr_q: float = 0.05,
        ivw_fdr_q: float = 0.05,
        cochran_q_pval: float = 0.05,
        pp4_threshold: float = 0.75,
):
    # extract stuff from here
    cis_mr_res = paths.mr_out(pqtl_dataset, pheno_id, local_results_dir)
    cis_mr_df = pl.read_csv(cis_mr_res, separator="\t")
    coloc_res = coloc_file or paths.coloc_out(pqtl_dataset, pheno_id, local_results_dir)
    coloc_df = pl.read_csv(coloc_res, separator="\t")

    wald_hits = []
    ivw_hits = []
    coloc_hits = set()

    for row in cis_mr_df.iter_rows(named=True):
        # separate where n_instruments == 1 or > 1
        n_instruments = row["n_instruments"]
        wald_fdr = row["Wald_FDR_q"]
        ivw_fdr = row["IVW_FDR_q"]
        cochran_q = row["Q_pval"]
        if n_instruments == 1 and wald_fdr is not None and wald_fdr < wald_fdr_q:
            wald_hits.append(row["protein"])
        elif (n_instruments > 1 and ivw_fdr is not None and ivw_fdr < ivw_fdr_q and cochran_q is not None and cochran_q > cochran_q_pval):
            ivw_hits.append(row["protein"])

    for row in coloc_df.iter_rows(named=True):
        pp4 = row["PP.H4.abf"]
        if pp4 is not None and pp4 > pp4_threshold:
            coloc_hits.add(row["protein"])

    # PWCoCo (conditional coloc) - complementary to standard COLOC above, not a
    # replacement (see project_pwcoco_wiring memory): a target that colocalises
    # under EITHER method should reach SMR, so PWCoCo-passing proteins are unioned
    # into coloc_hits here rather than requiring standard COLOC specifically.
    # pwcoco_out() may not exist - PWCoCo runs non-fatally in local.py/hpc.py, so a
    # failed or not-yet-run PWCoCo step must not break SMR eligibility.
    pwcoco_res = paths.pwcoco_out(pqtl_dataset, pheno_id, local_results_dir)
    if Path(pwcoco_res).exists():
        pwcoco_df = pl.read_csv(pwcoco_res, separator="\t")
        for row in pwcoco_df.iter_rows(named=True):
            h4 = row["H4"]
            if h4 is not None and h4 > pp4_threshold:
                coloc_hits.add(row["protein"])
    else:
        print(f"[TRACKING] No PWCoCo output found at {pwcoco_res}; using standard COLOC hits only...")

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
def prepare_smr_gwas(sumstats: str, pheno_id: str, local_results_dir: str = "results"):
    # temp dir to store .ma file per pheno
    temp_dir = paths.work_dir_for_results_dir(local_results_dir) / "smr"
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

    if "b_QTL" in df.columns:
        df = df.with_columns(
            pl.when(flip).then(-pl.col("b_QTL")).otherwise(pl.col("b_QTL")).alias("b_QTL")
        )

    if "Freq" in df.columns:
        df = df.with_columns(
            pl.when(flip).then(1 - pl.col("Freq")).otherwise(pl.col("Freq")).alias("Freq")
        )

    return df


# for single-cell targets, replace the SMR-reported QTL beta with the value from the
# original per-cell-type QTL file (source of truth) - re-signed to match the SMR
# file's own A1 (never touching A1/A2 themselves, so b_SMR / align_to_risk_allele stay
# internally consistent). Not used for bulk - those files are ingested as pre-computed.
def pull_original_sc_qtl_beta(target_smr: pl.DataFrame, qtl_dataset: str, cell: str):
    if target_smr.height == 0:
        return target_smr

    parquet_path = _smr_utils.resolve_sc_qtl_file(qtl_dataset, cell)
    if parquet_path is None or not parquet_path.exists():
        print(f"[CONCERN] Original QTL file not found for {qtl_dataset}/{cell}")
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
        print(f"[CONCERN] No matching original QTL rows found for {cell}")
        return target_smr

    target_smr = target_smr.join(original, on=["probeID", "topSNP"], how="left")

    target_smr = target_smr.with_columns(
        pl.when(pl.col("orig_EA") == pl.col("A1")).then(pl.col("orig_BETA"))
        .when(pl.col("orig_EA") == pl.col("A2")).then(-pl.col("orig_BETA"))
        .otherwise(None)
        .alias("b_qtl_from_source")
    )

    n_mismatched = target_smr.filter(
        pl.col("orig_EA").is_not_null() & pl.col("b_qtl_from_source").is_null()
    ).height
    if n_mismatched > 0:
        print(
            f"[CONCERN] {n_mismatched} row(s) in {cell} had an allele mismatch between "
            f"the SMR A1/A2 and the original QTL file - b_QTL left as SMR-reported"
        )

    target_smr = target_smr.with_columns(
        pl.when(pl.col("b_qtl_from_source").is_not_null())
        .then(pl.col("b_qtl_from_source"))
        .otherwise(pl.col("b_QTL"))
        .alias("b_QTL")
    ).drop(["orig_EA", "orig_BETA", "b_qtl_from_source"])

    # b_SMR is a ratio of b_GWAS/b_QTL - recompute so it stays consistent with the
    # (possibly now different-signed) b_QTL rather than leaving a stale ratio
    target_smr = target_smr.with_columns(
        pl.when(pl.col("b_QTL") != 0)
        .then(pl.col("b_GWAS") / pl.col("b_QTL"))
        .otherwise(pl.col("b_SMR"))
        .alias("b_SMR")
    )

    return target_smr


def run_single_cell_smr(pqtl_dataset: str, qtl_dataset: str, pheno_id: str, sumstats: str, ref_bfile: str, maf: float, local_results_dir: str = "results", synthesis_dir: str = "synthesis", coloc_file: str | None = None, wald_fdr_q: float = 0.05, ivw_fdr_q: float = 0.05, cochran_q_pval: float = 0.05, pp4_threshold: float = 0.75, p_qtl_smr: float = 5.0e-8, p_qtl_heidi: float = 1.57e-3):
    ref_bfile = Path(ref_bfile)
    qtl_temp = qtl_dataset.lower()

    if qtl_temp == "singlebrain":
        temp_sumstats = prepare_smr_gwas(sumstats, pheno_id, local_results_dir)

        # manifest-driven BESD discovery (detect_qtl_split finds the real
        # SMR_ready/<cell>/<cell>.besd triples) instead of a hardcoded cell
        # list + manual path/existence checks - verified byte-identical
        # against the old hardcoded discovery for all 7 real cell types
        besd_prefixes = _smr_utils.ensure_besd(qtl_temp, esd_dir=Path(synthesis_dir) / "qtl_esd")
        cell_types = sorted(prefix.name for prefix in besd_prefixes.values())

        for cell in cell_types:
            beqtl_summary = next(p for p in besd_prefixes.values() if p.name == cell)
            print(f"[TRACKING] Cell type {cell} found!")

            # check whether SMR has already been ran for trait X in cell type Y - this
            # depends only on (pheno_id, qtl_dataset, cell), never on pqtl_dataset, so
            # it's looked up under the shared synthesis/ tree (not local_results_dir)
            # and reused by every pqtl_dataset run instead of being recomputed per run
            smr_res = paths.smr_raw_dir(f"sc/{qtl_dataset}/{cell}", pheno_id, synthesis_dir)
            existing_smr = [f for f in smr_res.glob(f"*{pheno_id}*.smr") if f.stat().st_size > 0]

            if len(existing_smr) > 0:
                print(f"[TRACKING] SMR already completed for {pheno_id} in {cell} - skipping SMR")
            else:
                _smr_utils.run_smr(
                    pheno_id=pheno_id,
                    sumstats=temp_sumstats,
                    ref_bfile=ref_bfile,
                    beqtl_summary=beqtl_summary,
                    qtl_dataset=f"sc/{qtl_dataset}/{cell}",
                    p_qtl_smr=p_qtl_smr,
                    p_qtl_heidi=p_qtl_heidi,
                    thread_num=8,
                    maf=maf,
                    out_dir=synthesis_dir
                )

            # load SMR results
            # saving into out_dir 1 results file per cell type for trait X
            # synthesis/SMR/sc/SingleBrain/{cell}/{pheno_id}/...
            smr_res = paths.smr_raw_dir(f"sc/{qtl_dataset}/{cell}", pheno_id, synthesis_dir)
            for f in smr_res.glob("*.smr"):
                if pheno_id in f.name:
                    fdr_correct_smr_file(f, pheno_id, cell)

        # delete temp GWAS .ma only after all cell types are ran
        if temp_sumstats.exists():
            temp_sumstats.unlink()

        hits = extract_promising_targets(pheno_id=pheno_id, pqtl_dataset=pqtl_dataset, local_results_dir=local_results_dir, coloc_file=coloc_file, wald_fdr_q=wald_fdr_q, ivw_fdr_q=ivw_fdr_q, cochran_q_pval=cochran_q_pval, pp4_threshold=pp4_threshold)
        # now extract all of the SMR data from the results for each cell type pertaining to those targets and store as a dataframe in results/SMR/dataset
        # rows == 1 SMR result for target X on cell-type Y
        # so 7 cell types x X targets in terms of rows
        all_target_smr = []
        for cell in cell_types:
            smr_res = paths.smr_raw_dir(f"sc/{qtl_dataset}/{cell}", pheno_id, synthesis_dir)
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
                        pl.lit(qtl_dataset).alias("qtl_dataset"),
                        pl.lit(pqtl_dataset).alias("pqtl_dataset"),
                        pl.lit(resolve_qtl_type(qtl_dataset)).alias("qtl_type")
                    )
                )

                target_smr = pull_original_sc_qtl_beta(target_smr, qtl_dataset, cell)

                if target_smr.height > 0:
                    all_target_smr.append(target_smr)

        out_file = paths.smr_sc_out(pqtl_dataset, pheno_id, qtl_dataset, local_results_dir)
        os.makedirs(out_file.parent, exist_ok=True)

        if len(all_target_smr) > 0:
            final_smr_df = pl.concat(all_target_smr, how="diagonal_relaxed")
            final_smr_df = align_to_risk_allele(final_smr_df)
            final_smr_df.write_csv(out_file, separator="\t")
            print(f"[TRACKING] Compiled promising target SMR results saved to {out_file}")
        else:
            print(f"[CONCERN] No SMR results found for the promising {pqtl_dataset} targets")


# manifest-driven replacement for the old directory-scanning version - verified
# byte-identical against it for both real GTEx_v10 (14 tissues) and MetaBrain
# (1 label). A dataset registered as several per-region manifest rows (e.g.
# GTEx_v10's 14 real tissues, sharing parent_dataset=GTEx_v10) is looked up via
# that grouping; a single-row dataset (e.g. MetaBrain) is looked up directly.
# Chromosome number is the trailing digits of ensure_besd()'s own per-prefix
# label, regardless of whether the underlying naming uses "_chr4" (MetaBrain)
# or ".4" (GTEx_v10) - both real conventions. Returns {tissue_label: {chr_num:
# besd_prefix}}, or None if nothing is registered/found for qtl_dataset.
def bulk_tissue_prefixes(qtl_dataset: str, synthesis_dir: str = "synthesis"):
    manifest = _smr_utils.qtl_manifest
    rows = manifest.get_rows_by_parent(qtl_dataset)
    if not rows:
        try:
            rows = [manifest.get_row(qtl_dataset.lower())]
        except ValueError:
            print(f"[CONCERN] No raw bulk QTL layout known for {qtl_dataset} - cannot run SMR from scratch")
            return None

    tissues = {}
    for row in rows:
        if row.get("parent_dataset"):
            tissue = Path(row["path"]).parent.name
            label = f"{qtl_dataset.removesuffix('_v10')}_{tissue}_v10"
        else:
            label = Path(row["path"]).stem.split("_")[0]

        besd_prefixes = _smr_utils.ensure_besd(row["dataset"], esd_dir=Path(synthesis_dir) / "qtl_esd")
        if not besd_prefixes:
            continue

        chr_prefixes = {}
        for prefix_label, prefix in besd_prefixes.items():
            match = re.search(r"(\d+)$", prefix_label)
            if match:
                chr_prefixes[int(match.group(1))] = prefix
        if chr_prefixes:
            tissues[label] = chr_prefixes

    return tissues if tissues else None


# run the smr binary against the raw dat/bulk-eQTL besd/esi/epi files, 1 chromosome at a
# time per tissue (each chromosome's besd/esi/epi only covers that chromosome's SNPs/probes,
# so per-chromosome .smr outputs concatenate cleanly into 1 genome-wide file - no besd
# merging needed). Skips a tissue entirely if its final .smr already exists (same
# idempotency convention as run_single_cell_smr). A literal qtl_name column is stamped onto
# the concatenated output and it's FDR-corrected via the same helper single-cell SMR uses,
# so the result is indistinguishable from a "pre-computed" bulk file to ingest_bulk_smr.
def run_bulk_smr(pqtl_dataset: str, qtl_dataset: str, pheno_id: str, sumstats: str, ref_bfile: str, maf: float, local_results_dir: str = "results", synthesis_dir: str = "synthesis", p_qtl_smr: float = 5.0e-8, p_qtl_heidi: float = 1.57e-3):
    tissues = bulk_tissue_prefixes(qtl_dataset, synthesis_dir=synthesis_dir)

    if not tissues:
        print(f"[TRACKING] No raw bulk QTL besd/esi/epi files found under ./dat/bulk-eQTL/{qtl_dataset} - nothing to run SMR on")
        return

    ref_bfile = Path(ref_bfile)
    temp_sumstats = prepare_smr_gwas(sumstats, pheno_id, local_results_dir)

    for label, chr_prefixes in tissues.items():
        # bulk SMR(GWAS x eQTL) for a given (pheno_id, qtl_dataset, label) never depends
        # on pqtl_dataset - always looked up/written under the shared synthesis/ tree
        # (not local_results_dir) so every pqtl_dataset run reuses the same completed
        # genome-wide SMR instead of re-running the smr binary per chromosome each time
        bulk_dir = paths.smr_bulk_dir(qtl_dataset, synthesis_dir)
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
            _smr_utils.run_smr(
                pheno_id=pheno_id,
                sumstats=temp_sumstats,
                ref_bfile=ref_bfile,
                beqtl_summary=prefix,
                qtl_dataset=f"bulk_raw/{qtl_dataset}/{label}/chr{chr_num}",
                p_qtl_smr=p_qtl_smr,
                p_qtl_heidi=p_qtl_heidi,
                thread_num=8,
                maf=maf,
                out_dir=synthesis_dir
            )
            chr_out = Path(f"{paths.smr_raw_prefix(f'bulk_raw/{qtl_dataset}/{label}/chr{chr_num}', pheno_id, synthesis_dir)}.smr")
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
        ).with_columns(pl.lit(f"QTL_{label}").alias("qtl_name"))

        final_dir.mkdir(parents=True, exist_ok=True)
        combined.write_csv(final_file, separator="\t")
        fdr_correct_smr_file(final_file, pheno_id, label)
        print(f"[DONE] Saved genome-wide bulk SMR results for {label}: {final_file}")

        scratch_dir = Path(synthesis_dir) / "SMR" / "bulk_raw" / qtl_dataset / label
        if scratch_dir.exists():
            shutil.rmtree(scratch_dir)

    if temp_sumstats.exists():
        temp_sumstats.unlink()


# bulk eQTL SMR (MetaBrain / GTEx_v10) - loads whatever .smr-shaped files sit
# under results/SMR/bulk/{qtl_dataset}/, harmonises the column names onto the same shape
# single-cell uses, and extracts the promising targets. Those files are either produced
# just above by run_bulk_smr() (from the raw dat/bulk-eQTL besd/esi/epi, when available) or
# pre-computed elsewhere (a dataset with no raw dat/bulk-eQTL directory). GTEx_v10
# is tissue-resolved (1 file per tissue via rglob, same idea as single-cell's per-cell
# loop); MetaBrain is flat (1 file for the dataset).
def ingest_bulk_smr(pqtl_dataset: str, qtl_dataset: str, pheno_id: str, local_results_dir: str = "results", synthesis_dir: str = "synthesis", coloc_file: str | None = None, wald_fdr_q: float = 0.05, ivw_fdr_q: float = 0.05, cochran_q_pval: float = 0.05, pp4_threshold: float = 0.75):
    # same shared synthesis/ tree run_bulk_smr() writes to/checks - independent of
    # local_results_dir (this dataset's run-scoped out_dir) since the underlying SMR
    # computation is keyed only by (pheno_id, qtl_dataset), not pqtl_dataset
    bulk_dir = paths.smr_bulk_dir(qtl_dataset, synthesis_dir)
    smr_files = sorted(bulk_dir.rglob(f"*{pheno_id}*.smr"))

    if len(smr_files) == 0:
        print(f"[CONCERN] No pre-computed bulk SMR files found under {bulk_dir} for {pheno_id}")
        return

    print(f"[TRACKING] Found {len(smr_files)} pre-computed bulk SMR file(s) for {qtl_dataset}")

    hits = extract_promising_targets(pheno_id=pheno_id, pqtl_dataset=pqtl_dataset, local_results_dir=local_results_dir, coloc_file=coloc_file, wald_fdr_q=wald_fdr_q, ivw_fdr_q=ivw_fdr_q, cochran_q_pval=cochran_q_pval, pp4_threshold=pp4_threshold)
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
        # (e.g. QTL_GTEx_Brain_Cortex_v10 -> GTEx_Brain_Cortex_v10)
        label = qtl_dataset
        if "qtl_name" in smr_df.columns:
            qtl_names = smr_df.get_column("qtl_name").drop_nulls().unique().to_list()
            if len(qtl_names) == 1:
                label = qtl_names[0].removeprefix("QTL_")

        target_smr = (
            smr_df
            .filter(pl.col("Gene").is_in(target_genes))
            .with_columns(
                pl.col("Gene").replace(target_map).alias("protein"),
                pl.lit(label).alias("cell_type"),
                pl.lit("bulk").alias("data_type"),
                pl.lit(pheno_id).alias("phenotype"),
                pl.lit(qtl_dataset).alias("qtl_dataset"),
                pl.lit(pqtl_dataset).alias("pqtl_dataset"),
                pl.lit(resolve_qtl_type(qtl_dataset)).alias("qtl_type")
            )
        )

        if target_smr.height > 0:
            all_target_smr.append(target_smr)

    out_file = paths.smr_bulk_out(pqtl_dataset, pheno_id, qtl_dataset, local_results_dir)
    os.makedirs(out_file.parent, exist_ok=True)

    if len(all_target_smr) > 0:
        final_smr_df = pl.concat(all_target_smr, how="diagonal_relaxed")
        final_smr_df = align_to_risk_allele(final_smr_df)
        final_smr_df.write_csv(out_file, separator="\t")
        print(f"[TRACKING] Compiled promising target SMR results saved to {out_file}")
    else:
        print(f"[CONCERN] No pre-computed bulk SMR results found for the promising {pqtl_dataset} targets")


# match sql/schema.sql's smr_results column names - the SMR tool's own output
# uses these exact spellings (probeID, topSNP, ...) and plain lowercasing
# doesn't produce the schema's snake_case names for them. Applied to each side
# before the upsert concat in compile_multi_omics_targets() (rather than once
# after) so a freshly-computed frame (raw names) and an existing_df read back
# from an already-migrated file (already-renamed names) don't collide into
# duplicate half-populated columns.
def rename_smr_to_schema(frame: pl.DataFrame) -> pl.DataFrame:
    rename_map = {
        "probeID": "probe_id",
        "ProbeChr": "probe_chr",
        "topSNP": "top_snp",
        "topSNP_chr": "top_snp_chr",
        "topSNP_bp": "top_snp_bp",
        "start": "start_bp",
        "end": "end_bp",
    }
    rename_map = {k: v for k, v in rename_map.items() if k in frame.columns}
    return frame.rename(rename_map) if rename_map else frame


# shared by compile_multi_omics_targets() (per-call upsert, local.py/hpc.py) and
# merge_multi_omics_targets_batch() (fan-in, Nextflow) - HEIDI CUT OFF = 0.01,
# SMR CUT OFF = 0.05. Extracted so both call sites apply the identical gate
# rather than risking the two drifting apart.
def filter_smr_targets(df: pl.DataFrame, qtl_dataset: str, p_smr_threshold: float = 0.05, p_heidi_threshold: float = 0.01) -> pl.DataFrame:
    heidi_col = "P_HEIDI"
    if "p_HEIDI" in df.columns:
        heidi_col = "p_HEIDI"

    final_targets_df = (
        df
        .filter(
            pl.col("q_SMR").is_not_null(),
            pl.col(heidi_col).is_not_null(),
            pl.col("q_SMR") < p_smr_threshold,
            pl.col(heidi_col) > p_heidi_threshold
        )
        .sort(["protein", "cell_type", "q_SMR"])
    )

    if final_targets_df.height == 0:
        print(f"[CONCERN] No drug targets passed cis-MR (pQTLs) + COLOC + {qtl_dataset} QTL SMR")
        return final_targets_df

    return rename_smr_to_schema(final_targets_df)


def compile_multi_omics_targets(pheno_id: str, pqtl_dataset: str, qtl_dataset: str, qtl_mode: str, local_results_dir: str = "results", p_smr_threshold: float = 0.05, p_heidi_threshold: float = 0.01):
    # single-cell results live under results/SMR/sc/..., bulk under results/SMR/bulk/...
    targets_path = (
        paths.smr_sc_out(pqtl_dataset, pheno_id, qtl_dataset, local_results_dir)
        if qtl_mode == "single_cell"
        else paths.smr_bulk_out(pqtl_dataset, pheno_id, qtl_dataset, local_results_dir)
    )
    df = read_smr_tsv(targets_path)
    final_targets_df = filter_smr_targets(df, qtl_dataset, p_smr_threshold, p_heidi_threshold)

    if final_targets_df.height == 0:
        return []

    # canonical combined output (bulk + single-cell hits together) that the dashboard reads
    # upsert: drop any stale rows for this qtl_dataset, then append the fresh ones,
    # so bulk and single-cell runs (in either order) compose instead of overwriting each other
    combined_file = paths.smr_final_targets_out(pqtl_dataset, pheno_id, local_results_dir)
    os.makedirs(combined_file.parent, exist_ok=True)

    if combined_file.exists() and combined_file.stat().st_size > 0:
        existing_df = rename_smr_to_schema(read_smr_tsv(combined_file))
        existing_df = existing_df.filter(pl.col("qtl_dataset").str.to_lowercase() != qtl_dataset.lower())
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

    print(f"[TRACKING] {final_targets_df.height} target x cell-type SMR hits found for {qtl_dataset}")
    print(f"[TRACKING] {len(targets)} unique drug targets passed cis-MR (pQTLs) + COLOC + {qtl_dataset} QTL SMR")
    print(f"[TRACKING] Combined multi-omics target results saved to {combined_file}")
    print(f"[TRACKING] Drug targets: {targets}")
    return targets


# Nextflow fan-in counterpart to compile_multi_omics_targets(): instead of
# upserting one qtl_dataset's rows into a shared file across N sequential
# calls (only safe when those calls share one real persistent file, as in
# local.py/hpc.py), this takes every qtl_dataset's already-computed
# promising_targets.tsv at once (`inputs`: (qtl_dataset, qtl_mode,
# promising_targets_path) tuples) and writes the combined file in a single
# shot - correct regardless of whether the producing tasks ran in parallel or
# in isolated sandboxes, unlike the upsert pattern above (see
# project_nextflow_migration memory for why the upsert can't run per-task).
def merge_multi_omics_targets_batch(pheno_id: str, pqtl_dataset: str, inputs: list[tuple[str, str, str]], local_results_dir: str = "results", p_smr_threshold: float = 0.05, p_heidi_threshold: float = 0.01):
    filtered_frames = []
    for qtl_dataset, qtl_mode, targets_path in inputs:
        df = read_smr_tsv(Path(targets_path))
        filtered = filter_smr_targets(df, qtl_dataset, p_smr_threshold, p_heidi_threshold)
        if filtered.height > 0:
            filtered_frames.append(filtered)

    combined_file = paths.smr_final_targets_out(pqtl_dataset, pheno_id, local_results_dir)
    os.makedirs(combined_file.parent, exist_ok=True)

    if not filtered_frames:
        print(f"[CONCERN] No drug targets passed cis-MR (pQTLs) + COLOC + QTL SMR across any of {len(inputs)} QTL dataset(s)")
        pl.DataFrame().write_csv(combined_file, separator="\t")
        return []

    combined_df = pl.concat(filtered_frames, how="diagonal_relaxed")
    combined_df.write_csv(combined_file, separator="\t")

    targets = (
        combined_df
        .select("protein")
        .unique()
        .sort("protein")
        .get_column("protein")
        .to_list()
    )

    print(f"[TRACKING] {combined_df.height} target x cell-type SMR hits found across {len(inputs)} QTL dataset(s)")
    print(f"[TRACKING] {len(targets)} unique drug targets passed cis-MR (pQTLs) + COLOC + QTL SMR")
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
    p.add_argument("--sumstats", default=None)
    p.add_argument("--pqtl_dataset", required=True)
    p.add_argument("--qtl_dataset", default=None)
    p.add_argument("--qtl_mode", required=True, choices=["bulk", "single_cell", "merge"])
    p.add_argument("--ref_bfile", default=None)
    p.add_argument("--maf", type=float, default=0.01)
    p.add_argument("--local_results_dir", default="results")
    p.add_argument("--repo_root", default=None)
    p.add_argument("--synthesis_dir", default="synthesis")
    p.add_argument("--manifest_path", default=paths.DEFAULT_QTL_MANIFEST_PATH)
    p.add_argument("--gene_annotation", default=None)
    p.add_argument("--coloc_file", default=None)
    p.add_argument("--wald_fdr_q", type=float, default=0.05)
    p.add_argument("--ivw_fdr_q", type=float, default=0.05)
    p.add_argument("--cochran_q_pval", type=float, default=0.05)
    p.add_argument("--pp4_threshold", type=float, default=0.75)
    p.add_argument("--p_qtl_smr", type=float, default=5.0e-8)
    p.add_argument("--p_qtl_heidi", type=float, default=1.57e-3)
    p.add_argument("--p_smr_threshold", type=float, default=0.05)
    p.add_argument("--p_heidi_threshold", type=float, default=0.01)
    p.add_argument("--skip_merge", action="store_true")
    # merge mode only - repeatable "qtl_dataset:qtl_mode:promising_targets_path",
    # one per upstream SMR_BULK/SMR_SC task being fanned in
    p.add_argument("--input", action="append", default=[])
    args = p.parse_args()

    _smr_utils.manifest_path = args.manifest_path
    gene_annotation = Path(args.gene_annotation) if args.gene_annotation else None
    if gene_annotation is not None and args.repo_root and not gene_annotation.is_absolute():
        gene_annotation = Path(args.repo_root) / gene_annotation
    _smr_utils.ncbi_ref_path = str(gene_annotation) if gene_annotation is not None else None
    if args.repo_root:
        _smr_utils.base_dir = args.repo_root

    if args.qtl_mode == "merge":
        inputs = []
        for raw in args.input:
            qtl_dataset, qtl_mode, targets_path = raw.split(":", 2)
            inputs.append((qtl_dataset, qtl_mode, targets_path))
        merge_multi_omics_targets_batch(
            pheno_id=args.pheno_id,
            pqtl_dataset=args.pqtl_dataset,
            inputs=inputs,
            local_results_dir=args.local_results_dir,
            p_smr_threshold=args.p_smr_threshold,
            p_heidi_threshold=args.p_heidi_threshold,
        )
        return

    missing = [
        name for name, val in [("--sumstats", args.sumstats), ("--qtl_dataset", args.qtl_dataset), ("--ref_bfile", args.ref_bfile)]
        if val is None
    ]
    if missing:
        p.error(f"--qtl_mode {args.qtl_mode} requires: {', '.join(missing)}")

    if args.qtl_mode == "bulk":
        run_bulk_smr(
            pqtl_dataset=args.pqtl_dataset,
            qtl_dataset=args.qtl_dataset,
            pheno_id=args.pheno_id,
            sumstats=args.sumstats,
            ref_bfile=args.ref_bfile,
            maf=args.maf,
            local_results_dir=args.local_results_dir,
            synthesis_dir=args.synthesis_dir,
            p_qtl_smr=args.p_qtl_smr,
            p_qtl_heidi=args.p_qtl_heidi
        )
        ingest_bulk_smr(
            pqtl_dataset=args.pqtl_dataset,
            qtl_dataset=args.qtl_dataset,
            pheno_id=args.pheno_id,
            local_results_dir=args.local_results_dir,
            synthesis_dir=args.synthesis_dir,
            coloc_file=args.coloc_file,
            wald_fdr_q=args.wald_fdr_q,
            ivw_fdr_q=args.ivw_fdr_q,
            cochran_q_pval=args.cochran_q_pval,
            pp4_threshold=args.pp4_threshold
        )
    else:
        run_single_cell_smr(
            pqtl_dataset=args.pqtl_dataset,
            qtl_dataset=args.qtl_dataset,
            pheno_id=args.pheno_id,
            sumstats=args.sumstats,
            maf=args.maf,
            ref_bfile=args.ref_bfile,
            local_results_dir=args.local_results_dir,
            synthesis_dir=args.synthesis_dir,
            coloc_file=args.coloc_file,
            wald_fdr_q=args.wald_fdr_q,
            ivw_fdr_q=args.ivw_fdr_q,
            cochran_q_pval=args.cochran_q_pval,
            pp4_threshold=args.pp4_threshold,
            p_qtl_smr=args.p_qtl_smr,
            p_qtl_heidi=args.p_qtl_heidi
        )

    if not args.skip_merge:
        compile_multi_omics_targets(
            pheno_id=args.pheno_id,
            pqtl_dataset=args.pqtl_dataset,
            qtl_dataset=args.qtl_dataset,
            qtl_mode=args.qtl_mode,
            local_results_dir=args.local_results_dir,
            p_smr_threshold=args.p_smr_threshold,
            p_heidi_threshold=args.p_heidi_threshold
        )

if __name__ == "__main__":
    main()
