#!/usr/bin/env python3
import argparse
import subprocess
import tempfile
import time
from pathlib import Path

import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
import polars as pl
import requests
import streamlit as st
from liftover import ChainFile
from plotly.subplots import make_subplots

from drugmr import paths, registry
from bin.load_db_into_postgres import PostgresLoader, PostgresReader

# shared plotting conventions so charts look consistent across tabs rather than each
# px.* call picking its own default palette
SIGNIFICANCE_COLOR_MAP = {True: "#d62728", False: "#7f7f7f"}  # red = significant, grey = not
SEQUENTIAL_SCALE = "Viridis"  # continuous significance / intensity scale

# status colors reused across the prioritisation Sankey - fixed, never themed
STATUS_GOOD = "#0ca30c"
STATUS_CRITICAL = "#d03b3b"
STATUS_MUTED = "#898781"
SANKEY_BULK_COLOR = "#2a78d6"
SANKEY_SC_COLOR = "#eb6834"
SANKEY_BOTH_COLOR = "#1baf7a"

# These 3 hex values get reused for 2 UNRELATED tri-state groupings (deliberately,
# not by accident): the Final Targets Sankey's bulk/sc/both eQTL-support column,
# and its coloc_support (COLOC-only/PWCoCo-only/both) column - both happen to
# want a blue/orange/teal "1 of 2 methods / both" pattern, and this also matches
# the green/blue/orange badge convention used for coloc_support on the Overview
# cards and Target Profile (st.badge colors "green"/"blue"/"orange"). Aliased
# below purely so a reader of coloc_support code isn't confused by "BULK"/"SC"
# naming that has nothing to do with COLOC/PWCoCo - same colors, on purpose.
COLOC_SUPPORT_BOTH_COLOR = SANKEY_BOTH_COLOR
COLOC_SUPPORT_COLOC_ONLY_COLOR = SANKEY_BULK_COLOR
COLOC_SUPPORT_PWCOCO_ONLY_COLOR = SANKEY_SC_COLOR

# single source of truth for "where am I in the pipeline" - the Overview tab's step
# map and every downstream tab's stage caption are both built from this list, so the
# two can never drift out of sync with each other or with the st.tabs() labels below
PIPELINE_STAGES = [
    dict(title="cis-MR", blurb="Mendelian randomisation of cis-instrumented protein abundance on the outcome."),
    dict(title="pQTL–GWAS COLOC", blurb="Pairwise colocalisation confirming the pQTL and GWAS signals share one causal variant."),
    dict(title="FinnGen PheWAS", blurb="Phenome-wide MR classifying Bonferroni-significant hits as potential additional indications or adverse effects."),
    dict(title="UKB PheWAS", blurb="Fallback phenome-wide MR in UK Biobank EHR-derived phenotypes, for targets uncovered by FinnGen."),
    dict(title="SMR (bulk/sc eQTL)", blurb="SMR + HEIDI test that the pQTL signal also acts through transcription."),
    dict(title="HyPrColoc (bulk/sc eQTL)", blurb="pQTL + GWAS + eQTL signals sharing one causal variant - via HyPrColoc's clustering, or via PWCoCo-QTL's SNP-level triangulation."),
    dict(title="Final Targets", blurb="Targets surviving every stage above, each reported at its correct SNP."),
]


def stage_caption(stage_number: int):
    """Small 'you are here' overline shown above a tab's main header.

    Kept as 1 short call site per tab (rather than baking the text into every
    st.subheader) so the wording can't drift between tabs and PIPELINE_STAGES.
    """
    total = len(PIPELINE_STAGES)
    title = PIPELINE_STAGES[stage_number - 1]["title"]
    st.caption(f"STAGE {stage_number} OF {total} · {title.upper()}")


def hex_to_rgba(hex_color: str, alpha: float):
    hex_color = hex_color.lstrip("#")
    r, g, b = (int(hex_color[i:i + 2], 16) for i in (0, 2, 4))
    return f"rgba({r},{g},{b},{alpha})"


def format_protein_list_html(proteins, per_line: int = 6, limit: int = 48):
    proteins = sorted(proteins)

    if not proteins:
        return "(no targets)"

    truncated = len(proteins) > limit
    shown = proteins[:limit]
    lines = [", ".join(shown[i:i + per_line]) for i in range(0, len(shown), per_line)]
    text = "<br>".join(lines)

    if truncated:
        text += f"<br>+{len(proteins) - limit} more"

    return text


def layout_sankey_columns(column_indices: list, node_values: list, n_nodes: int, margin: float = 0.06, min_gap: float = 0.08):
    """Vertical layout for a multi-column Sankey, solved one column at a time.

    Laying the flow out as a strict tree (each node nested inside its parent's
    vertical band) collapses once the pass lane narrows: with ~130 proteins in
    at cis-MR and ~12 surviving to SMR, that lane's four SMR children have a
    twelfth of the height to share and their labels overlap. Positioning each
    column independently instead keeps nodes in pass-lane-first order (so no
    ribbon crosses another) while guaranteeing every node a `min_gap` slot of
    its own. Node *heights* stay value-proportional - that's Plotly's shared
    scale and it's what makes the drop-off readable - only the centres move.
    """
    node_y = [0.5] * n_nodes
    low, high = margin, 1 - margin

    # Plotly sizes node heights on one scale shared by the whole diagram (the
    # fullest column fills the plot), so spans have to be measured against that
    # same scale rather than each column's own total. Measuring per column would
    # stretch a stage's nodes to fill the height even after most targets have
    # dropped out, pushing small nodes to the bottom edge on long swooping
    # ribbons that no longer line up with the blocks they connect.
    scale = max((sum(node_values[index] for index in indices) for indices in column_indices if indices), default=1) or 1

    for indices in column_indices:
        if not indices:
            continue

        # stack from the top so the surviving lane stays a near-horizontal band
        # and drop-outs peel off underneath it
        centers = []
        cursor = low

        for index in indices:
            span = (high - low) * node_values[index] / scale
            centers.append(cursor + span / 2)
            cursor += span

        for position in range(1, len(centers)):
            centers[position] = max(centers[position], centers[position - 1] + min_gap)

        overflow = centers[-1] - high

        if overflow > 0:
            centers = [center - overflow for center in centers]

            for position in range(len(centers) - 1, 0, -1):
                centers[position - 1] = min(centers[position - 1], centers[position] - min_gap)

        for index, center in zip(indices, centers):
            node_y[index] = min(max(center, 0.02), 0.98)

    return node_y


# KEY CHANGES DOWN THE LINE WITH MORE PQTL DATASETS
# -> CHANGE THE DASHBOARD FUNCT TO ADD MORE PQTL DATASETS
# biomarker meta analysis: https://pmc.ncbi.nlm.nih.gov/articles/instance/12136742/pdf/nihpp-rs6597595v1.pdf

def retention(current: int, previous: int):
    return 0.0 if previous == 0 else 100 * current / previous


def available_cols(df: pd.DataFrame, cols: list[str]):
    return [col for col in cols if col in df.columns]


def load_required_tsv(file: Path, label: str):
    if not file.exists():
        st.error(f"{label} result file not found: {file}")
        st.stop()

    try:
        df = pd.read_csv(file, sep="\t", low_memory=False)
    except (OSError, pd.errors.ParserError, UnicodeDecodeError) as error:
        st.error(f"{label} result file could not be read: {file}")
        st.exception(error)
        st.stop()

    if df.empty:
        st.error(f"{label} result file is empty: {file}")
        st.stop()

    return df


def load_optional_tsv(file: Path, label: str, warn_if_missing: bool = True):
    # warn_if_missing=False is for results that are legitimately, routinely absent
    # (e.g. 0 targets triangulated) rather than "this pipeline step never ran" -
    # showing a raw filesystem path in a warning banner on every tab for an
    # expected empty result reads as a crash to a standard user, not a status
    if not file.exists():
        if warn_if_missing:
            st.warning(f"{label} result file not found: {file}")
        return pd.DataFrame()

    try:
        df = pd.read_csv(file, sep="\t", low_memory=False)
    except (OSError, pd.errors.ParserError, UnicodeDecodeError) as error:
        if warn_if_missing:
            st.warning(f"{label} result file could not be read: {file}")
            st.exception(error)
        return pd.DataFrame()

    if df.empty:
        if warn_if_missing:
            st.warning(f"{label} result file is empty: {file}")
        return pd.DataFrame()

    return df


def find_result_file(project_dir: Path, candidate_files: list[Path], candidate_names: list[str]):
    for file in candidate_files:
        if file.exists():
            return file

    matches = []

    for candidate_name in candidate_names:
        matches.extend(list((project_dir / "results").rglob(candidate_name)))

    matches = sorted(set(matches))

    if len(matches) == 1:
        return matches[0]

    return None


def legacy_resolve_dataset_files(project_dir: Path, phenotype: str, dataset_id: str):
    """Pre-Phase-3 candidate-path guessing, kept as a fallback for any
    (phenotype, dataset_id) with no runs/registry.json entry yet - e.g. a
    dataset that predates the runs/ migration and was never backfilled."""
    mr_file = find_result_file(
        project_dir,
        [
            project_dir / f"results/cis-MR/{dataset_id}_{phenotype}_all_MR.tsv",
            project_dir / f"results/cis_MR/{dataset_id}_{phenotype}_all_MR.tsv",
            project_dir / f"results/cis_MR/{dataset_id}/{phenotype}/{dataset_id}_{phenotype}_all_MR.tsv",
            project_dir / f"results/MR/{dataset_id}/{phenotype}/{dataset_id}_{phenotype}_all_MR.tsv",
            project_dir / f"results/MR/{dataset_id}_{phenotype}_all_MR.tsv"
        ],
        [
            f"{dataset_id}_{phenotype}_all_MR.tsv",
            f"{dataset_id}_{phenotype}_MR.tsv"
        ]
    )

    coloc_file = find_result_file(
        project_dir,
        [
            project_dir / f"results/coloc/{dataset_id}/{dataset_id}_{phenotype}_all_coloc.tsv",
            project_dir / f"results/coloc/{dataset_id}/{phenotype}/{dataset_id}_{phenotype}_all_coloc.tsv",
            project_dir / f"results/COLOC/{dataset_id}/{phenotype}/{dataset_id}_{phenotype}_all_coloc.tsv",
            project_dir / f"results/COLOC/{dataset_id}_{phenotype}_all_coloc.tsv",
            project_dir / f"results/coloc/{dataset_id}_{phenotype}_all_coloc.tsv"
        ],
        [
            f"{dataset_id}_{phenotype}_all_coloc.tsv",
            f"{dataset_id}_{phenotype}_coloc.tsv",
            f"{dataset_id}_{phenotype}_COLOC.tsv"
        ]
    )

    return {
        "mr": mr_file,
        "coloc": coloc_file,
        "finngen_phewas": project_dir / "results" / "PheWAS-FinnGen" / dataset_id / phenotype / f"{dataset_id}_{phenotype}_PheWAS-FinnGen.tsv",
        "ukb_phewas": project_dir / "results" / "PheWAS_UKBB" / dataset_id / phenotype / f"{dataset_id}_{phenotype}_PheWAS.tsv",
        "target_info": project_dir / "results" / "target_stats" / dataset_id / phenotype / f"{dataset_id}_{phenotype}_top_cis_hits.tsv",
        "smr": project_dir / "results" / "SMR" / f"{dataset_id}_{phenotype}_final_multi_omics_targets.tsv",
        "hyprcoloc": project_dir / "results" / "hyprcoloc" / f"{dataset_id}_{phenotype}_all_hyprcoloc.tsv",
        "pwcoco": project_dir / "results" / "pwcoco" / dataset_id / f"{dataset_id}_{phenotype}_all_pwcoco.tsv",
        "pwcoco_eqtl_pqtl": project_dir / "results" / "pwcoco_qtl" / dataset_id / f"{dataset_id}_{phenotype}_all_pwcoco_eqtl_pqtl.tsv",
        "pwcoco_eqtl_gwas": project_dir / "results" / "pwcoco_qtl" / dataset_id / f"{dataset_id}_{phenotype}_all_pwcoco_eqtl_gwas.tsv",
    }


def resolve_dataset_files(project_dir: Path, phenotype: str, dataset_id: str, run_id: str = "latest"):
    """Resolve this dataset's 6 dashboard files via runs/registry.json first -
    falling back to legacy_resolve_dataset_files() when this (phenotype,
    dataset_id) has no registry entry (never run through the migrated
    pipeline). Returns (run_id_used_or_None, {file_key: Path}).
    """
    runs_root = str(project_dir / "runs")
    resolved_run_id = run_id
    if run_id == "latest":
        resolved_run_id = registry.get_latest_run_id(phenotype, dataset_id, root=runs_root)

    if resolved_run_id is not None:
        out_dir = str(project_dir / paths.run_results_dir(resolved_run_id))
        return resolved_run_id, {
            "mr": project_dir / paths.mr_out(dataset_id, phenotype, out_dir),
            "coloc": project_dir / paths.coloc_out(dataset_id, phenotype, out_dir),
            "finngen_phewas": project_dir / paths.phewas_out(dataset_id, phenotype, out_dir),
            "ukb_phewas": project_dir / paths.phewas_ukbb_out(dataset_id, phenotype, out_dir),
            "target_info": project_dir / paths.target_stats_out(dataset_id, phenotype, out_dir),
            "smr": project_dir / paths.smr_final_targets_out(dataset_id, phenotype, out_dir),
            "hyprcoloc": project_dir / paths.hyprcoloc_out(dataset_id, phenotype, out_dir),
            "pwcoco": project_dir / paths.pwcoco_out(dataset_id, phenotype, out_dir),
            "pwcoco_eqtl_pqtl": project_dir / paths.pwcoco_eqtl_pqtl_out(dataset_id, phenotype, out_dir),
            "pwcoco_eqtl_gwas": project_dir / paths.pwcoco_eqtl_gwas_out(dataset_id, phenotype, out_dir),
        }

    return None, legacy_resolve_dataset_files(project_dir, phenotype, dataset_id)


def filter_protein(df: pd.DataFrame, protein: str):
    if df.empty or not protein or "protein" not in df.columns:
        return df

    return df[
        df["protein"]
        .astype(str)
        .str.contains(protein, case=False, na=False, regex=False)
    ].copy()


def safe_nunique(df: pd.DataFrame, col: str):
    if df.empty or col not in df.columns:
        return 0

    return df[col].nunique()


def safe_median(df: pd.DataFrame, col: str, scientific: bool = False):
    if df.empty or col not in df.columns:
        return "NA"

    values = pd.to_numeric(df[col], errors="coerce").dropna()

    if values.empty:
        return "NA"

    if scientific:
        return f"{values.median():.3e}"

    return f"{values.median():.3f}"


def standardise_columns(df: pd.DataFrame):
    df = df.copy()

    df.columns = (
        pd.Index(df.columns)
        .astype(str)
        .str.strip()
        .str.lower()
        .str.replace(r"[^a-z0-9]+", "_", regex=True)
        .str.strip("_")
    )

    if df.columns.duplicated().any():
        duplicated_cols = df.columns[df.columns.duplicated()].unique().tolist()
        st.error(f"Duplicated columns after standardisation: {duplicated_cols}")
        st.stop()

    return df


# 1 PWCoCo output table -> {protein: {snp: h4}}, keeping only rows that clear
# pp4_thresh - mirrors bin/pwcoco_qtl_wrapper.py's snp_h4_map(), but computed live
# here against the sidebar's own PP.H4 slider (pp4) rather than trusting a
# pre-computed file frozen at whatever threshold the pipeline last used. SNP1/SNP2
# is literally "unconditioned" for PWCoCo's own unconditioned row, and a
# conditioned row's SNP carries a trailing "*" (PWCoCo's conditioning-SNP marker) -
# stripped here so the same variant matches across combos regardless of source.
def compute_snp_h4_map(df: pd.DataFrame, pp4_thresh: float):
    m = {}
    if df.empty or "protein" not in df.columns or "h4" not in df.columns:
        return m

    h4 = pd.to_numeric(df["h4"], errors="coerce")
    passing = df[h4.fillna(0) >= pp4_thresh]

    for _, row in passing.iterrows():
        protein = str(row["protein"])
        row_h4 = float(row["h4"])
        for col in ("snp1", "snp2"):
            if col not in row.index or pd.isna(row[col]):
                continue
            snp = str(row[col])
            if snp and snp != "unconditioned":
                snp = snp.rstrip("*")
                protein_map = m.setdefault(protein, {})
                protein_map[snp] = max(protein_map.get(snp, 0), row_h4)

    return m


# 1 row per (target, SNP) pair where the SAME SNP clears pp4_thresh in ALL 3
# combos (pQTL-GWAS, eQTL-pQTL, eQTL-GWAS) at once - the live equivalent of
# bin/pwcoco_qtl_wrapper.py's shared_rows. Every row here is a triangulated
# target; the "protein" column's unique values are triangulated_proteins.
def compute_shared_snp_table(pqtl_gwas_df: pd.DataFrame, eqtl_pqtl_df: pd.DataFrame, eqtl_gwas_df: pd.DataFrame, pp4_thresh: float):
    pg_map = compute_snp_h4_map(pqtl_gwas_df, pp4_thresh)
    ep_map = compute_snp_h4_map(eqtl_pqtl_df, pp4_thresh)
    eg_map = compute_snp_h4_map(eqtl_gwas_df, pp4_thresh)

    rows = []
    for protein in set(pg_map) & set(ep_map) & set(eg_map):
        shared_snps = set(pg_map[protein]) & set(ep_map[protein]) & set(eg_map[protein])
        for snp in shared_snps:
            rows.append({
                "protein": protein,
                "snp": snp,
                "pqtl_gwas_h4": pg_map[protein][snp],
                "eqtl_pqtl_h4": ep_map[protein][snp],
                "eqtl_gwas_h4": eg_map[protein][snp],
            })

    return pd.DataFrame(rows)



def prepare_phewas(df: pd.DataFrame):
    if df.empty:
        return df

    df = standardise_columns(df)

    df = df.rename(columns={
        "protein_id": "protein",
    })

    # A1/A2 already come from the original outcome GWAS
    # do not overwrite them with PheWAS REF/ALT
    for col in ["a1", "a2", "ukb_ref", "ukb_alt", "finngen_ref", "finngen_alt"]:
        if col in df.columns:
            df[col] = df[col].astype(str).str.upper()

    for col in [
        "n_instruments",
        "n_instruments_original",
        "n_instruments_available",
        "instrument_completeness",
        "beta_mr",
        "se_mr",
        "p_mr",
        "p_fdr",
        "fdr_q",
        "p_bonferroni"
    ]:
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors="coerce")

    if "fdr_significant" in df.columns:
        df["fdr_significant"] = (
            df["fdr_significant"]
            .astype(str)
            .str.lower()
            .isin(["true", "1", "yes"])
        )
    elif "fdr_q" in df.columns:
        df["fdr_significant"] = df["fdr_q"].fillna(np.inf) <= 0.05
    elif "p_fdr" in df.columns:
        df["fdr_significant"] = df["p_fdr"].fillna(np.inf) <= 0.05

    # deliberately NOT computing bonferroni_significant/p_bonferroni here (this
    # used to trust whatever the raw file had, or a naive <=0.05 cutoff) -
    # classify_phewas_associations() always recomputes both fresh, per protein,
    # from p_mr, before anything reaches the UI on every real call path
    # (render_phewas_section, compute_phewas_classification_status) - keeping a
    # 2nd, different (and never-actually-used) computation here was dead code
    # that misleadingly implied this function's own correction was load-bearing.

    return df


def classify_phewas_associations(phewas_df: pd.DataFrame, mr_outcome_df: pd.DataFrame) -> pd.DataFrame:
    """Row-level PheWAS MR classification against the primary protein->AD effect.

    Recomputes Bonferroni correction itself, per protein, across the number of
    endpoints actually loaded for that protein in this source (rather than
    trusting whatever fixed-constant correction may be baked into an older
    output file) - this is what makes the correction self-correcting even
    against results/ files written before the per-protein Bonferroni fix in
    bin/phewas_cis_pqtls.py / bin/ukb_phewas.py.

    Classification rule (Bonferroni-significant associations only):
      - same sign as the primary cis-MR beta -> "additional_indication"
        (i increased AD risk -> inhibition indicated -> a positive PheWAS
        estimate is a same-direction repurposing signal; ii decreased AD risk
        -> augmentation indicated -> a negative PheWAS estimate is the
        same-direction signal - both collapse to "same sign as primary beta")
      - opposite sign -> "adverse_effect"
    Non-significant associations are "not_classified". Requires beta_mr to be
    on the protein-abundance-increasing-allele axis for both PheWAS sources
    (true as of the ukb_phewas.py alignment fix - historic UKB rows written
    before that fix are on the wrong axis and will misclassify).
    """
    df = phewas_df.copy()

    if df.empty or not {"protein", "p_mr", "beta_mr"}.issubset(df.columns):
        df["phewas_classification"] = pd.Series(dtype=object)
        return df

    df["protein"] = df["protein"].astype(str)
    df["p_mr"] = pd.to_numeric(df["p_mr"], errors="coerce")
    df["beta_mr"] = pd.to_numeric(df["beta_mr"], errors="coerce")
    if "se_mr" in df.columns:
        df["se_mr"] = pd.to_numeric(df["se_mr"], errors="coerce")

    valid = df["p_mr"].notna() & (df["p_mr"] > 0) & df["beta_mr"].notna()
    df["n_endpoints_tested"] = df.loc[valid].groupby("protein")["p_mr"].transform("count")
    n_endpoints_tested = df["n_endpoints_tested"]
    df["p_bonferroni"] = np.minimum(df["p_mr"] * n_endpoints_tested, 1.0)
    df["bonferroni_significant"] = valid & (df["p_mr"] < (0.05 / n_endpoints_tested))

    # odds ratio view - both FinnGen (R13 endpoints) and UKB-TOPMed (PheCodes) here
    # are binary disease/case-control phenotypes, so beta_mr is a log-odds ratio
    if "se_mr" in df.columns:
        df["or_mr"] = np.exp(df["beta_mr"])
        df["or_ci_low"] = np.exp(df["beta_mr"] - 1.96 * df["se_mr"])
        df["or_ci_high"] = np.exp(df["beta_mr"] + 1.96 * df["se_mr"])
        df["beta_ci_low"] = df["beta_mr"] - 1.96 * df["se_mr"]
        df["beta_ci_high"] = df["beta_mr"] + 1.96 * df["se_mr"]

    primary_beta = (
        mr_outcome_df.drop_duplicates("protein").set_index("protein")["mr_beta"]
        if not mr_outcome_df.empty and {"protein", "mr_beta"}.issubset(mr_outcome_df.columns)
        else pd.Series(dtype=float)
    )
    df["primary_mr_beta"] = df["protein"].map(primary_beta)
    df["primary_mr_direction"] = np.select(
        [df["primary_mr_beta"] > 0, df["primary_mr_beta"] < 0],
        ["increases_ad_risk", "decreases_ad_risk"],
        default=None,
    )

    has_primary = df["primary_mr_beta"].notna() & (df["primary_mr_beta"] != 0)
    same_direction = np.sign(df["beta_mr"]) == np.sign(df["primary_mr_beta"])

    df["phewas_classification"] = "not_classified"
    classifiable = df["bonferroni_significant"].fillna(False) & has_primary
    df.loc[classifiable & same_direction, "phewas_classification"] = "additional_indication"
    df.loc[classifiable & ~same_direction, "phewas_classification"] = "adverse_effect"

    return df


def compute_phewas_classification_status(phewas_df: pd.DataFrame, mr_outcome_df: pd.DataFrame, proteins):
    """Per-protein worst-case PheWAS classification for prioritisation gates.

    'adverse_effect' if the protein has >=1 opposite-direction Bonferroni-
    significant hit; else 'additional_indication' if it has >=1 same-direction
    one (informational - NOT a gate failure); else 'none'. Only 'adverse_effect'
    should exclude a target from Prioritised Targets / the Final Targets Sankey.
    """
    proteins = list(proteins)
    status = {protein: "none" for protein in proteins}

    classified = classify_phewas_associations(phewas_df, mr_outcome_df)
    if classified.empty or "protein" not in classified.columns:
        return status

    adverse_proteins = set(
        classified.loc[classified["phewas_classification"] == "adverse_effect", "protein"].unique()
    )
    indication_proteins = set(
        classified.loc[classified["phewas_classification"] == "additional_indication", "protein"].unique()
    )

    for protein in proteins:
        if protein in adverse_proteins:
            status[protein] = "adverse_effect"
        elif protein in indication_proteins:
            status[protein] = "additional_indication"

    return status


def compute_hyprcoloc_pass_status(hyprcoloc_df: pd.DataFrame, proteins, threshold: float):
    """Per-protein HyPrColoc pass flag for the prioritisation Sankey.

    A protein PASSES only if it has a HyPrColoc cluster whose traits include the
    pQTL, GWAS and eQTL trait together (i.e. all 3 traits share a single causal
    variant, not just 2 of the 3) with posterior_prob >= threshold. Proteins with
    no HyPrColoc row, or only 2-trait / low-probability clusters, count as FAIL.
    """
    proteins = list(proteins)
    status = {protein: False for protein in proteins}

    if hyprcoloc_df.empty or not {"protein", "traits", "posterior_prob"}.issubset(hyprcoloc_df.columns):
        return status

    df = hyprcoloc_df.copy()
    df["protein"] = df["protein"].astype(str)
    df["posterior_prob"] = pd.to_numeric(df["posterior_prob"], errors="coerce")
    traits_lower = df["traits"].astype(str).str.lower()

    all_three_traits = (
        traits_lower.str.contains("pqtl_") &
        traits_lower.str.contains("gwas_") &
        traits_lower.str.contains("eqtl_")
    )

    passing = all_three_traits & (df["posterior_prob"].fillna(0) >= threshold)
    passing_proteins = set(df.loc[passing, "protein"].unique())

    for protein in proteins:
        if protein in passing_proteins:
            status[protein] = True

    return status


# HyPrColoc rows (1 per protein x cell_type/tissue, occasionally more when a cluster
# couldn't be resolved in 1 shot) that actually PASS - same "all 3 traits together,
# posterior_prob >= threshold" rule as compute_hyprcoloc_pass_status, but returns the
# rows themselves (not just a per-protein bool) so the Final Targets table can read off
# each row's candidate_snp and its aligned alleles/betas (a1/a2/gwas_beta/gwas_p/
# pqtl_beta/pqtl_p/eqtl_beta/eqtl_p - written by bin/hyprcoloc_targets.py's
# attach_candidate_snp_stats, p == 0 already replaced with 1e-300 at the source).
# Older HyPrColoc runs predating that change won't have those columns yet - callers
# should handle their absence rather than assume they're always there.
def select_hyprcoloc_candidate_rows(hyprcoloc_df: pd.DataFrame, threshold: float):
    required = {"protein", "cell_type", "traits", "posterior_prob", "candidate_snp"}

    if hyprcoloc_df.empty or not required.issubset(hyprcoloc_df.columns):
        return pd.DataFrame()

    df = hyprcoloc_df.copy()
    df["protein"] = df["protein"].astype(str)
    df["posterior_prob"] = pd.to_numeric(df["posterior_prob"], errors="coerce")
    traits_lower = df["traits"].astype(str).str.lower()

    all_three_traits = (
        traits_lower.str.contains("pqtl_") &
        traits_lower.str.contains("gwas_") &
        traits_lower.str.contains("eqtl_")
    )

    passing = df[all_three_traits & (df["posterior_prob"].fillna(0) >= threshold)].copy()

    if passing.empty:
        return passing

    # 1 row per protein x cell_type (x data_type, when present) - keep the
    # highest-probability cluster if more than 1 qualifying row exists
    group_cols = [col for col in ["protein", "cell_type", "data_type"] if col in passing.columns]
    passing = passing.sort_values("posterior_prob", ascending=False).drop_duplicates(subset=group_cols, keep="first")
    return passing


# ---------------------------------------------------------------------------
# Regional locus plot (stacked GWAS / pQTL / optional eQTL + gene track)
#
# GWAS + pQTL regional summary stats are read straight from the per-target
# cis-region parquets bin/pwcoco_wrapper.py and bin/hyprcoloc_targets.py
# already produce (dat/cis_regions/{pqtl_dataset}/{protein}/{pqtl,gwas}.parquet)
# - the FULL cis window, not just the handful of MR instrument SNPs. eQTL
# regional data mirrors bin/hyprcoloc_targets.py's load_eqtl_table() exactly,
# but keeps every SNP (not collapsed to the lead one) and keeps CHR/BP for
# plotting. Gene body coordinates (chr/start/end/strand) come for free from
# smr_final_targets_out's own gene annotation - no external gene reference
# or live API call needed.
# ---------------------------------------------------------------------------

# fallback flat colours for the regional locus plot's GWAS/pQTL/eQTL panels,
# used ONLY when LD (r²) colouring is unavailable (see ld_available below) -
# reuses the same 3 hex values as SANKEY_BULK/SC/BOTH_COLOR, which is a 3rd,
# also-unrelated use of this palette (that pair means bulk/sc/both eQTL there);
# not a naming/semantics link, just a shared "3 blue/orange/teal" palette
REGIONAL_TRACK_COLORS = {
    "GWAS": SANKEY_BULK_COLOR,
    "pQTL": SANKEY_BOTH_COLOR,
    "eQTL": SANKEY_SC_COLOR,
}

# standard LocusZoom-style LD binning - red (tight LD with the candidate) through
# blue (independent), grey for SNPs absent from the reference panel, purple
# diamond for the candidate itself
LD_BIN_EDGES = [-0.01, 0.2, 0.4, 0.6, 0.8, 1.0]
LD_BIN_LABELS = ["0.0–0.2", "0.2–0.4", "0.4–0.6", "0.6–0.8", "0.8–1.0"]
LD_BIN_COLORS = ["#4575b4", "#91bfdb", "#fee090", "#fc8d59", "#d73027"]
LD_NO_DATA_COLOR = "#bdbdbd"
CANDIDATE_COLOR = "#762a83"
GWAS_SIGNIFICANCE_P = 5e-8


@st.cache_data(show_spinner=False)
def load_regional_cis_data(pqtl_dataset: str, protein: str):
    """Full regional pQTL + GWAS summary stats for 1 target's cis window."""
    project_dir = Path(__file__).resolve().parent.parent
    cis_dir = project_dir / "dat" / "cis_regions" / pqtl_dataset / protein
    pqtl_file = cis_dir / "pqtl.parquet"
    gwas_file = cis_dir / "gwas.parquet"

    if not pqtl_file.exists() or not gwas_file.exists():
        return pd.DataFrame(), pd.DataFrame()

    cols = ["SNP", "CHR", "BP", "A1", "A2", "BETA", "SE", "P"]
    pqtl = pl.read_parquet(pqtl_file).select(cols).to_pandas()
    gwas = pl.read_parquet(gwas_file).select(cols).to_pandas()
    pqtl.columns = pqtl.columns.str.lower()
    gwas.columns = gwas.columns.str.lower()
    return pqtl, gwas


def _liftover_hg19_to_hg38(df: pd.DataFrame, chr_col: str = "chr", bp_col: str = "bp") -> pd.DataFrame:
    """In-place-style liftover of chr_col/bp_col from GRCh37/hg19 to GRCh38, using
    the same chain file bin/qc_gwas.py already uses elsewhere in this pipeline.
    Rows that fail to lift, or whose bp_col is null, are dropped. Updates BOTH
    columns from the lift hit (not just bp_col) - a hg19->hg38 liftover can move
    a coordinate onto a different chromosome, same as bin/qc_gwas.py's own
    liftover_df_to_hg38 accounts for."""
    project_dir = Path(__file__).resolve().parent.parent
    chain_file = project_dir / "dat" / "ref" / "liftover" / "hg19ToHg38.over.chain"
    if not chain_file.exists() or df.empty:
        return df

    df = df.dropna(subset=[bp_col]).copy()
    if df.empty:
        return df

    converter = ChainFile(str(chain_file), one_based=True)
    lifted_chr, lifted_bp = [], []
    for chrom, pos in zip(df[chr_col].astype(str), df[bp_col].astype(int)):
        query_chrom = chrom if chrom.startswith("chr") else f"chr{chrom}"
        hits = converter[query_chrom][int(pos)]
        if hits:
            lifted_chr.append(hits[0][0].replace("chr", ""))
            lifted_bp.append(hits[0][1])
        else:
            lifted_chr.append(None)
            lifted_bp.append(None)

    df[chr_col] = lifted_chr
    df[bp_col] = lifted_bp
    return df.dropna(subset=[chr_col, bp_col])


@st.cache_data(show_spinner=False)
def load_regional_eqtl_data(data_type: str, eqtl_dataset: str, cell_type: str, base_gene_id: str):
    """Full regional eQTL summary stats for 1 gene - same file layout and allele
    convention as bin/hyprcoloc_targets.py's load_eqtl_table(), but every SNP in
    the window is kept (not collapsed to the lead one) and CHR/BP are retained.

    IMPORTANT: MetaBrain's own BP column is GRCh37/hg19, not GRCh38 like every
    other file this pipeline reads (confirmed directly - e.g. rs5848 sits at
    chr17:42,430,244 in MetaBrain vs chr17:44,352,876 everywhere else, a ~1.92 Mb
    offset). GTEx_v10 and SingleBrain are already GRCh38 (verified against the
    same SNP) and are left untouched. Only MetaBrain gets lifted over here."""
    project_dir = Path(__file__).resolve().parent.parent

    if data_type == "single_cell":
        eqtl_file = project_dir / "dat" / "sc-eQTL" / eqtl_dataset / f"{cell_type}.parquet"
        if not eqtl_file.exists():
            return pd.DataFrame()
        df = (
            pl.scan_parquet(eqtl_file)
            .filter(pl.col("GENE").str.split(".").list.first() == base_gene_id)
            .select(["SNP", "CHR", "BP", "A1", "A2", "EA", "BETA", "SE", "P"])
            .with_columns(
                pl.when(pl.col("EA") == pl.col("A2")).then(pl.col("A1")).otherwise(pl.col("A2")).alias("eqtl_a2")
            )
            .select(["SNP", "CHR", "BP", pl.col("EA").alias("A1"), pl.col("eqtl_a2").alias("A2"), "BETA", "SE", "P"])
            .collect()
            .to_pandas()
        )
    elif data_type == "bulk":
        if eqtl_dataset == "GTEx_v10":
            tissue = cell_type.removeprefix("GTEx_").removesuffix("_v10")
            eqtl_file = project_dir / "dat" / "bulk-eQTL" / "GTEx_v10" / tissue / f"{tissue}.parquet"
        elif eqtl_dataset == "MetaBrain":
            eqtl_file = project_dir / "dat" / "bulk-eQTL" / "MetaBrain" / "BrainMeta_cis_eQTL.parquet"
        else:
            return pd.DataFrame()

        if not eqtl_file.exists():
            return pd.DataFrame()

        df = (
            pl.scan_parquet(eqtl_file)
            .filter(pl.col("Probe").str.split(".").list.first() == base_gene_id)
            .select(["SNP", pl.col("Chr").alias("CHR"), "BP", "A1", "A2", pl.col("b").alias("BETA"), "SE", pl.col("p").alias("P")])
            .collect()
            .to_pandas()
        )
        df.columns = df.columns.str.lower()
        if eqtl_dataset == "MetaBrain":
            df = _liftover_hg19_to_hg38(df)
        return df
    else:
        return pd.DataFrame()

    df.columns = df.columns.str.lower()
    return df


@st.cache_data(show_spinner=False)
def load_regional_ld(candidate_snp: str, chrom, window_kb: int = 5000):
    """r² between candidate_snp and every other SNP on the same chromosome within
    window_kb, computed live from the pipeline's own 1000G EUR Phase 3 reference
    panel (dat/ref/1000G_EUR_Phase3_plink/1000G.EUR.QC.ALL) - the SAME reference
    used for LD clumping and PWCoCo elsewhere in this pipeline, via plink --r2.

    Matched onto plotted SNPs by rsID ONLY, never by position: this reference
    panel is GRCh37/hg19 (confirmed - e.g. rs4292 sits at chr17:61,554,341 here
    vs chr17:63,476,980 in this pipeline's GRCh38 pQTL/GWAS/eQTL data, ~1.9 Mb
    off, consistent with a build offset) while every other file this dashboard
    reads is GRCh38. BP values are therefore NOT comparable across the two -
    only rsIDs are build-stable enough to join on.
    """
    project_dir = Path(__file__).resolve().parent.parent
    ref_bfile = project_dir / "dat" / "ref" / "1000G_EUR_Phase3_plink" / "1000G.EUR.QC.ALL"

    # NOTE: Path.with_suffix() only replaces text after the LAST dot, which would
    # mangle this filename (it has dots in the basename itself) - plain string
    # concatenation is the correct way to append a fixed suffix here.
    if not Path(f"{ref_bfile}.bed").exists():
        return pd.DataFrame()

    try:
        chrom_int = int(chrom)
    except (TypeError, ValueError):
        return pd.DataFrame()

    with tempfile.TemporaryDirectory() as tmp_dir:
        out_prefix = Path(tmp_dir) / "ld"
        cmd = [
            "plink", "--bfile", str(ref_bfile),
            "--chr", str(chrom_int),
            "--r2", "--ld-snp", str(candidate_snp),
            "--ld-window-kb", str(window_kb),
            "--ld-window", "999999",
            "--ld-window-r2", "0",
            "--out", str(out_prefix),
        ]
        try:
            with st.spinner("Computing LD (r²) against the 1000G reference panel..."):
                subprocess.run(cmd, check=True, capture_output=True, text=True, timeout=120)
        except (subprocess.CalledProcessError, subprocess.TimeoutExpired, FileNotFoundError):
            return pd.DataFrame()

        ld_file = out_prefix.with_suffix(".ld")
        if not ld_file.exists():
            return pd.DataFrame()

        ld = pd.read_csv(ld_file, sep=r"\s+")
        ld.columns = ld.columns.str.lower()
        if "snp_b" not in ld.columns or "r2" not in ld.columns:
            return pd.DataFrame()
        return ld[["snp_b", "r2"]].rename(columns={"snp_b": "snp"})


@st.cache_data(show_spinner=False)
def load_genes_in_region(chrom, start_bp: int, end_bp: int):
    """Every GRCh38 protein-coding gene overlapping [start_bp, end_bp] on chrom,
    from Ensembl's REST API (rest.ensembl.org - GRCh38 by default), so the gene
    track shows the whole local gene neighbourhood, not just the 1 target gene.
    Cached per-region, so this is a live call only the first time a given window
    is viewed - not on every rerender.

    Ensembl's overlap/region endpoint is a shared, rate-limited public service -
    confirmed empirically it can 503 transiently even on modest (~2 Mb) queries
    with no consistent size threshold (500kb/1Mb/3Mb all succeeded, 2Mb failed
    once then a similar-sized query succeeded), so this retries a couple of
    times with backoff before giving up and falling back to the single-gene
    SMR annotation (handled by the caller)."""
    start_bp, end_bp = int(start_bp), int(end_bp)
    max_span = 4_900_000
    if end_bp - start_bp > max_span:
        mid = (start_bp + end_bp) // 2
        start_bp, end_bp = mid - max_span // 2, mid + max_span // 2

    chrom_str = str(chrom).replace("chr", "")
    genes = None
    with st.spinner("Loading gene annotations for this region..."):
        for attempt, backoff in enumerate((0, 3, 8)):
            if backoff:
                time.sleep(backoff)
            try:
                response = requests.get(
                    f"https://rest.ensembl.org/overlap/region/human/{chrom_str}:{start_bp}-{end_bp}",
                    params={"feature": "gene", "content-type": "application/json"},
                    timeout=45,
                )
                response.raise_for_status()
                genes = response.json()
                break
            except (requests.RequestException, ValueError):
                continue

    if genes is None:
        return pd.DataFrame()

    if not genes:
        return pd.DataFrame()

    df = pd.DataFrame(genes)
    required = {"start", "end", "strand", "gene_id"}
    if not required.issubset(df.columns):
        return pd.DataFrame()

    if "biotype" in df.columns:
        df = df[df["biotype"] == "protein_coding"]

    if df.empty:
        return pd.DataFrame()

    df["symbol"] = df["external_name"] if "external_name" in df.columns else df["gene_id"]
    df["symbol"] = df["symbol"].fillna(df["gene_id"])
    df["strand"] = df["strand"].map({1: "+", -1: "-"}).fillna("")
    return df[["gene_id", "symbol", "start", "end", "strand"]].sort_values("start").reset_index(drop=True)


def _pack_gene_lanes(genes_df: pd.DataFrame, window_bp: float):
    """Greedy interval packing so overlapping/nearby genes stack into separate
    lanes instead of colliding - same idea IGV/LocusZoom gene tracks use. Adds
    a 0-indexed 'lane' column; returns (genes_df, n_lanes)."""
    if genes_df.empty:
        return genes_df, 1

    padding = max(window_bp * 0.03, 1000)
    genes_df = genes_df.sort_values("start").copy()
    lane_ends = []
    lanes = []

    for _, gene in genes_df.iterrows():
        placed = False
        for lane_idx, lane_end in enumerate(lane_ends):
            if gene["start"] > lane_end + padding:
                lane_ends[lane_idx] = gene["end"]
                lanes.append(lane_idx)
                placed = True
                break
        if not placed:
            lane_ends.append(gene["end"])
            lanes.append(len(lane_ends) - 1)

    genes_df["lane"] = lanes
    return genes_df, len(lane_ends)


def _select_labelled_genes(genes_df: pd.DataFrame, window_bp: float, target_symbol: str):
    """Collision-free static-label placement, per lane: the target gene always
    gets a label; every other gene gets one ONLY if it fits without overlapping
    an already-placed label in its lane (checked against every reserved
    interval, not just the previous one, since the target is placed first and
    is not necessarily left-to-right in position). Genes that lose out on a
    static label still get a hover-visible marker at render time - nothing
    becomes undiscoverable, it just isn't drawn as text when the window is too
    gene-dense to fit every name without overlapping. Adds a 'show_label' bool
    column."""
    if genes_df.empty:
        genes_df["show_label"] = pd.Series(dtype=bool)
        return genes_df

    genes_df = genes_df.copy()
    genes_df["show_label"] = False

    # empirical per-character label footprint, in bp - deliberately generous
    # (favours skipping a label over risking a collision)
    char_width_bp = window_bp * 0.011
    min_gap_bp = window_bp * 0.006
    min_label_width_bp = window_bp * 0.02

    def label_span(gene):
        mid = (gene["start"] + gene["end"]) / 2
        half_width = max(len(str(gene["symbol"])) * char_width_bp, min_label_width_bp) / 2
        return mid - half_width, mid + half_width

    is_target = genes_df["symbol"].astype(str).str.upper() == target_symbol
    # target(s) placed first so they always win a label; the rest follow in
    # genomic order for a stable, deterministic result
    ordered_index = list(genes_df[is_target].index) + list(genes_df[~is_target].sort_values("start").index)

    reserved_by_lane = {}
    for idx in ordered_index:
        gene = genes_df.loc[idx]
        lane = gene["lane"]
        lo, hi = label_span(gene)
        occupied = reserved_by_lane.setdefault(lane, [])
        collides = any(not (hi + min_gap_bp < r_lo or lo - min_gap_bp > r_hi) for r_lo, r_hi in occupied)
        if not collides:
            occupied.append((lo, hi))
            genes_df.loc[idx, "show_label"] = True

    return genes_df


def _regional_eqtl_options(smr_rows: pd.DataFrame):
    """1 selectable option per unique (eqtl_dataset, data_type, cell_type,
    probeid) combination this target has SMR support in, with a friendly label.
    Returns (options, default_label) - default is whichever combo has the
    strongest SMR evidence (lowest q_SMR, falling back to p_SMR), so the eQTL
    panel is populated automatically rather than defaulting to hidden."""
    required = {"eqtl_dataset", "data_type", "cell_type", "probeid"}
    if smr_rows.empty or not required.issubset(smr_rows.columns):
        return [], None

    rank_col = "q_smr" if "q_smr" in smr_rows.columns else ("p_smr" if "p_smr" in smr_rows.columns else None)
    group_cols = list(required)
    df = smr_rows.dropna(subset=group_cols).copy()

    if rank_col is not None:
        df[rank_col] = pd.to_numeric(df[rank_col], errors="coerce")
        df = df.sort_values(rank_col, na_position="last")

    combos = df.drop_duplicates(subset=group_cols, keep="first")

    options = []
    for _, row in combos.iterrows():
        data_type = str(row["data_type"])
        eqtl_dataset = str(row["eqtl_dataset"])
        cell_type = str(row["cell_type"])
        probeid = str(row["probeid"])

        if data_type == "bulk" and eqtl_dataset == "GTEx_v10":
            label = f"GTEx v10 · {cell_type.removeprefix('GTEx_').removesuffix('_v10')}"
        elif data_type == "bulk":
            label = f"{eqtl_dataset} (bulk)"
        else:
            label = f"{eqtl_dataset} · {cell_type} (single-cell)"

        options.append({
            "label": label,
            "data_type": data_type,
            "eqtl_dataset": eqtl_dataset,
            "cell_type": cell_type,
            "base_gene_id": probeid.split(".")[0],
            "rank": row.get(rank_col) if rank_col is not None else None,
        })

    default_label = options[0]["label"] if options else None  # already best-ranked, before alpha sort
    options.sort(key=lambda o: o["label"])
    return options, default_label


def render_regional_locus_plot(protein: str, pqtl_dataset: str, smr_rows: pd.DataFrame, hypr_rows: pd.DataFrame, key_prefix: str):
    """Stacked regional association plot (GWAS + pQTL, optional eQTL, optional
    gene track) for 1 target - the visual counterpart to the PP.H4/H4 badges
    shown above it: do these signals actually overlap at this locus?"""
    pqtl_df, gwas_df = load_regional_cis_data(pqtl_dataset, protein)

    if pqtl_df.empty or gwas_df.empty:
        st.info(
            "No regional cis-window summary statistics found on disk for this target - "
            "the plot needs `dat/cis_regions/.../{pqtl,gwas}.parquet`, produced alongside PWCoCo."
        )
        return

    # candidate SNP: prefer HyPrColoc's colocalisation-informed candidate (Stage 6),
    # fall back to this target's own top cis-pQTL hit (always available once COLOC
    # has run) - same "best evidence currently available, not a hard requirement"
    # pattern used throughout render_target_profile
    candidate_snp = None
    candidate_source = None
    if not hypr_rows.empty and "candidate_snp" in hypr_rows.columns:
        candidate_values = hypr_rows["candidate_snp"].dropna()
        if not candidate_values.empty:
            candidate_snp = str(candidate_values.iloc[0])
            candidate_source = "HyPrColoc cluster's candidate SNP"

    if candidate_snp is None:
        top_pqtl = pqtl_df.sort_values("p").iloc[0]
        candidate_snp = str(top_pqtl["snp"])
        candidate_source = "top cis-pQTL hit"

    candidate_bp = None
    candidate_chr = None
    for df in (pqtl_df, gwas_df):
        match = df[df["snp"] == candidate_snp]
        if not match.empty:
            candidate_bp = float(match.iloc[0]["bp"])
            candidate_chr = match.iloc[0]["chr"]
            break

    eqtl_options, default_eqtl_label = _regional_eqtl_options(smr_rows)
    with st.container(border=True):
        option_labels = ["None"] + [o["label"] for o in eqtl_options]
        default_index = option_labels.index(default_eqtl_label) if default_eqtl_label in option_labels else 0
        selected_label = st.selectbox(
            "eQTL track",
            option_labels,
            index=default_index,
            key=f"{key_prefix}_regional_eqtl_select",
            help=(
                "Adds a 3rd panel showing this gene's eQTL association in the chosen "
                "tissue/cell type. Defaults to whichever SMR dataset gave this target "
                "its strongest transcriptional support."
            ),
        )

    # GWAS/pQTL's own cis-window defines the region shown throughout - eQTL cis-
    # windows are typically wider (confirmed - MetaBrain/GTEx routinely span
    # several Mb vs pQTL/GWAS's ~1-2 Mb), and letting that wider window drive
    # the shared x-axis stretched the whole plot and squashed the real peak
    # into a sliver in the middle. eQTL is cropped to this boundary instead.
    window_start_bp = float(min(gwas_df["bp"].min(), pqtl_df["bp"].min()))
    window_end_bp = float(max(gwas_df["bp"].max(), pqtl_df["bp"].max()))

    eqtl_df = pd.DataFrame()
    eqtl_label = None
    if selected_label != "None":
        selected = next(o for o in eqtl_options if o["label"] == selected_label)
        eqtl_df_full = load_regional_eqtl_data(
            data_type=selected["data_type"],
            eqtl_dataset=selected["eqtl_dataset"],
            cell_type=selected["cell_type"],
            base_gene_id=selected["base_gene_id"],
        )
        eqtl_label = selected_label
        if eqtl_df_full.empty:
            st.caption(f"No eQTL rows found for this gene in {eqtl_label} within the cached parquet - showing GWAS/pQTL only.")
        else:
            eqtl_df = eqtl_df_full[
                (eqtl_df_full["bp"] >= window_start_bp) & (eqtl_df_full["bp"] <= window_end_bp)
            ]
            n_cropped = len(eqtl_df_full) - len(eqtl_df)
            if n_cropped > 0:
                st.caption(
                    f"{eqtl_label}'s eQTL cis-window extends beyond the GWAS/pQTL region shown here - "
                    f"{n_cropped:,} SNP(s) outside {window_start_bp/1e6:.2f}-{window_end_bp/1e6:.2f} Mb are not plotted."
                )
            if eqtl_df.empty:
                st.caption(f"No {eqtl_label} eQTL SNPs fall within the GWAS/pQTL window - showing GWAS/pQTL only.")

    # every protein-coding gene in the GWAS/pQTL window (GRCh38) - matches the
    # region actually shown, since eQTL is now cropped to the same boundary
    window_chr = candidate_chr if candidate_chr is not None else gwas_df["chr"].iloc[0]
    genes_df = load_genes_in_region(window_chr, window_start_bp, window_end_bp)

    if genes_df.empty:
        # Ensembl unreachable, or nothing protein-coding overlapped the window -
        # fall back to just this target's own gene, from its SMR annotation
        # (bin/sort_smr.py's own gene reference), rather than showing no track at all
        smr_gene_cols = {"chr", "start", "end", "strand"}
        if not smr_rows.empty and smr_gene_cols.issubset(smr_rows.columns):
            fallback_row = smr_rows.dropna(subset=list(smr_gene_cols)).head(1)
            if not fallback_row.empty:
                fallback = fallback_row.iloc[0]
                genes_df = pd.DataFrame([{
                    "gene_id": None,
                    "symbol": protein.split("_")[0],
                    "start": fallback["start"],
                    "end": fallback["end"],
                    "strand": fallback["strand"],
                }])

    n_gene_lanes = 0
    target_symbol = protein.split("_")[0].upper()
    if not genes_df.empty:
        genes_df, n_gene_lanes = _pack_gene_lanes(genes_df, window_end_bp - window_start_bp)
        genes_df = _select_labelled_genes(genes_df, window_end_bp - window_start_bp, target_symbol)

    # r² against the candidate SNP, from the pipeline's own 1000G EUR reference -
    # graceful no-op (flat single colour per track) if plink/the reference panel
    # is unavailable, or the candidate itself isn't in the panel
    ld_df = pd.DataFrame()
    if candidate_chr is not None:
        ld_df = load_regional_ld(candidate_snp, candidate_chr)
    ld_available = not ld_df.empty

    if candidate_bp is not None:
        ld_note = " Points are coloured by LD (r²) with this variant." if ld_available else " LD (r²) colouring unavailable for this variant - showing uncoloured points."
        st.caption(f"Candidate variant **{candidate_snp}** ({candidate_source}) marked in every panel.{ld_note}")
    else:
        st.caption(f"Candidate variant **{candidate_snp}** ({candidate_source}) - position not found in the plotted window, marker omitted.")

    show_gene_track = not genes_df.empty
    if not show_gene_track:
        st.caption("Gene track unavailable - couldn't reach Ensembl and no SMR gene annotation exists yet for this target.")
    row_specs = [("GWAS", gwas_df), ("pQTL", pqtl_df)]
    if not eqtl_df.empty:
        row_specs.append((f"eQTL · {eqtl_label}", eqtl_df))
    n_data_rows = len(row_specs)
    if show_gene_track:
        row_specs.append(("Gene", None))

    n_rows = len(row_specs)
    # kept proportional to the actual pixel height formula below (230px per
    # data row, 60+35px/lane for the gene row) so the 2 never drift apart as
    # n_gene_lanes grows - they used to use unrelated formulas, which visually
    # over- or under-allocated the gene track for gene-dense windows
    gene_row_px = 60 + 35 * n_gene_lanes
    gene_row_height = gene_row_px / 230
    row_heights = [1.0] * n_data_rows + ([gene_row_height] if show_gene_track else [])

    fig = make_subplots(
        rows=n_rows, cols=1, shared_xaxes=True, vertical_spacing=0.05,
        row_heights=row_heights,
        subplot_titles=[label for label, _ in row_specs],
    )

    # a legend entry is shown once, on whichever panel first actually has a
    # point for it - hardcoding "only row 1 (GWAS)" would silently drop legend
    # entries for an LD bin/the candidate SNP that's absent from GWAS but
    # present in pQTL/eQTL
    legend_shown = set()

    def _show_legend_once(legend_key: str) -> bool:
        if legend_key in legend_shown:
            return False
        legend_shown.add(legend_key)
        return True

    for i, (label, df) in enumerate(row_specs, start=1):
        if label == "Gene":
            gene_mid_mb, gene_y, gene_hover = [], [], []
            for _, gene in genes_df.iterrows():
                start_mb = float(gene["start"]) / 1_000_000
                end_mb = float(gene["end"]) / 1_000_000
                lane = int(gene["lane"]) if "lane" in gene else 0
                y = -lane
                is_target = str(gene["symbol"]).upper() == target_symbol
                color = CANDIDATE_COLOR if is_target else "#90a4ae"
                strand = str(gene.get("strand", ""))
                arrow = "→" if strand == "+" else ("←" if strand == "-" else "")
                fig.add_shape(
                    type="line", x0=start_mb, x1=end_mb, y0=y, y1=y,
                    line=dict(color=color, width=10 if is_target else 6), row=i, col=1,
                )
                # only genes that survived _select_labelled_genes' collision check
                # get a static text label - every gene still gets a hover point
                # below, so an unlabelled gene's name is still 1 hover away
                if gene.get("show_label", False):
                    label_text = f"<b>{gene['symbol']}</b> {arrow}" if is_target else f"{gene['symbol']} {arrow}"
                    fig.add_annotation(
                        x=(start_mb + end_mb) / 2, y=y + 0.38, showarrow=False,
                        text=label_text, row=i, col=1,
                        font=dict(size=11 if is_target else 9, color=color),
                    )
                gene_mid_mb.append((start_mb + end_mb) / 2)
                gene_y.append(y)
                gene_hover.append(f"{gene['symbol']} ({strand})" if strand else str(gene["symbol"]))

            fig.add_trace(
                go.Scatter(
                    x=gene_mid_mb, y=gene_y, mode="markers",
                    marker=dict(size=6, color="rgba(0,0,0,0)"), showlegend=False,
                    hovertext=gene_hover, hoverinfo="text",
                ),
                row=i, col=1,
            )
            fig.update_yaxes(visible=False, range=[-n_gene_lanes + 0.4, 1], row=i, col=1)
            fig.update_xaxes(showgrid=False, row=i, col=1)
            continue

        plot_df = df.copy()
        plot_df["bp_mb"] = plot_df["bp"] / 1_000_000
        plot_df["minus_log10_p"] = -np.log10(plot_df["p"].clip(lower=1e-300))
        is_candidate = plot_df["snp"] == candidate_snp
        candidate_points = plot_df[is_candidate]
        other_points = plot_df[~is_candidate]

        if ld_available:
            merged = other_points.merge(ld_df, on="snp", how="left")
            merged["ld_bin"] = pd.cut(merged["r2"], bins=LD_BIN_EDGES, labels=LD_BIN_LABELS)

            for color, bin_label in zip(LD_BIN_COLORS, LD_BIN_LABELS):
                bucket = merged[merged["ld_bin"] == bin_label]
                if bucket.empty:
                    continue
                fig.add_trace(
                    go.Scatter(
                        x=bucket["bp_mb"], y=bucket["minus_log10_p"], mode="markers",
                        marker=dict(size=7, color=color, line=dict(width=0.4, color="white")),
                        name=f"r² {bin_label}", legendgroup=f"ld_{bin_label}",
                        showlegend=_show_legend_once(f"ld_{bin_label}"),
                        customdata=np.stack([bucket["snp"], bucket["r2"]], axis=-1),
                        hovertemplate="%{customdata[0]}<br>pos=%{x:.3f} Mb<br>-log10(p)=%{y:.2f}<br>r²=%{customdata[1]:.2f}<extra>" + label + "</extra>",
                    ),
                    row=i, col=1,
                )

            no_ld = merged[merged["r2"].isna()]
            if not no_ld.empty:
                fig.add_trace(
                    go.Scatter(
                        x=no_ld["bp_mb"], y=no_ld["minus_log10_p"], mode="markers",
                        marker=dict(size=6, color=LD_NO_DATA_COLOR, opacity=0.6),
                        name="No LD data", legendgroup="ld_none",
                        showlegend=_show_legend_once("ld_none"),
                        customdata=no_ld[["snp"]],
                        hovertemplate="%{customdata[0]}<br>pos=%{x:.3f} Mb<br>-log10(p)=%{y:.2f}<extra>" + label + "</extra>",
                    ),
                    row=i, col=1,
                )
        else:
            track_color = REGIONAL_TRACK_COLORS.get(label.split(" ")[0], "#4c78a8")
            fig.add_trace(
                go.Scatter(
                    x=other_points["bp_mb"], y=other_points["minus_log10_p"], mode="markers",
                    marker=dict(size=6, color=track_color, opacity=0.75),
                    name=label, showlegend=False,
                    customdata=other_points[["snp"]],
                    hovertemplate="%{customdata[0]}<br>pos=%{x:.3f} Mb<br>-log10(p)=%{y:.2f}<extra>" + label + "</extra>",
                ),
                row=i, col=1,
            )

        if not candidate_points.empty:
            fig.add_trace(
                go.Scatter(
                    x=candidate_points["bp_mb"], y=candidate_points["minus_log10_p"], mode="markers",
                    marker=dict(symbol="diamond", size=14, color=CANDIDATE_COLOR, line=dict(width=1.5, color="white")),
                    name="Candidate SNP", legendgroup="candidate",
                    showlegend=_show_legend_once("candidate"),
                    customdata=candidate_points[["snp"]],
                    hovertemplate="%{customdata[0]} (candidate)<br>pos=%{x:.3f} Mb<br>-log10(p)=%{y:.2f}<extra>" + label + "</extra>",
                ),
                row=i, col=1,
            )

        fig.update_yaxes(title_text="−log₁₀(p)", title_font=dict(size=11), row=i, col=1)
        fig.update_xaxes(showgrid=True, gridcolor="#eef1f4", row=i, col=1)

        if candidate_bp is not None:
            fig.add_vline(x=candidate_bp / 1_000_000, line_dash="dot", line_color="#9e9e9e", line_width=1, row=i, col=1)

        if label == "GWAS":
            fig.add_hline(
                y=-np.log10(GWAS_SIGNIFICANCE_P), line_dash="dash", line_color="#d62728", line_width=1,
                opacity=0.7, row=i, col=1,
                annotation_text="genome-wide significance (5×10⁻⁸)", annotation_position="top left",
                annotation_font=dict(size=9, color="#d62728"),
            )

    fig.update_xaxes(title_text="Position (Mb, GRCh38)", title_font=dict(size=12), row=n_rows, col=1)
    # restyle ONLY the auto-generated subplot titles (their text is exactly a
    # row_specs label) - a selector-less for_each_annotation would also clobber
    # the gene-track labels and the GWAS significance-line annotation added
    # above, overriding their deliberate size/colour and, for gene labels,
    # invalidating _select_labelled_genes' collision-free placement (which
    # assumes the smaller font it was sized for)
    subplot_title_texts = {label for label, _ in row_specs}
    fig.for_each_annotation(
        lambda a: a.update(font=dict(size=13, color="#37474f")) if a.text in subplot_title_texts else None
    )
    fig.update_layout(
        height=230 * n_data_rows + (gene_row_px if show_gene_track else 0),
        template="plotly_white",
        font=dict(family="-apple-system, Helvetica Neue, Arial, sans-serif", size=12),
        margin=dict(t=70, b=50, l=60, r=30),
        legend=dict(orientation="h", yanchor="bottom", y=1.0, xanchor="right", x=1.0, font=dict(size=10)),
        plot_bgcolor="white",
        hovermode="closest",
    )

    st.plotly_chart(fig, width="stretch", key=f"{key_prefix}_regional_plot")

    with st.expander("View regional summary statistics"):
        data_specs = [(label, df) for label, df in row_specs if label != "Gene"]
        tabs = st.tabs([label for label, _ in data_specs])
        for idx, (tab, (label, df)) in enumerate(zip(tabs, data_specs)):
            with tab:
                display_df = df[["snp", "chr", "bp", "a1", "a2", "beta", "se", "p"]].sort_values("p")
                st.dataframe(display_df, width="stretch", hide_index=True)
                st.download_button(
                    label=f"Download {label} regional summary stats",
                    data=display_df.to_csv(index=False, sep="\t"),
                    file_name=f"{protein}_{pqtl_dataset}_{label.replace(' ', '_').replace('·', '')}_regional.tsv",
                    mime="text/tab-separated-values",
                    key=f"{key_prefix}_download_regional_{idx}",
                )


def subset_phewas_outcome(df: pd.DataFrame, outcome: str):
    if df.empty:
        return df

    if "outcome_trait" in df.columns:
        return df[df["outcome_trait"] == outcome].copy()

    if "pheno_id" in df.columns:
        return df[df["pheno_id"] == outcome].copy()

    return df.copy()


def probability_strength_color(value, threshold: float = 0.70):
    """Tiering for an already-passing probability (COLOC PP.H4, HyPrColoc posterior).

    Every row shown has already cleared `threshold` (the sidebar's configured
    gate, e.g. PP.H4 >= pp4), so this isn't a pass/fail color - it's how
    comfortably a row cleared it. Tiers are computed relative to `threshold`
    (not hardcoded 0.80/0.90) so they never contradict a value that a plain
    pass/fail badge elsewhere (e.g. Target Profile's `color="green" if
    coloc_pp4 >= pp4 else "red"`) would already call green.
    """
    if pd.isna(value):
        return "gray"

    span = max(1.0 - threshold, 1e-6)
    if value >= threshold + 0.8 * span:
        return "green"

    if value >= threshold + 0.5 * span:
        return "blue"

    return "orange"


def significance_strength_color(value):
    """Same tiering idea as probability_strength_color, for a p-value/FDR column."""
    if pd.isna(value):
        return "gray"

    if value <= 1e-4:
        return "green"

    if value <= 1e-2:
        return "blue"

    return "orange"


def render_prioritised_target_cards(df: pd.DataFrame, outcome: str, pp4: float = 0.70, n_columns: int = 2):
    """1 scientific 'card' per prioritised protein instead of a dense numeric table.

    The Overview tab's raw dataframe of ~20 columns was the single most confusing part
    of the dashboard for a first-time reader - too dense to scan for "which targets look
    good" at a glance. The full table (every column) + a TSV download are still available
    in the expander directly below this grid for anyone who wants to inspect or export it.
    """
    columns = st.columns(n_columns)

    for position, (_, row) in enumerate(df.iterrows()):
        protein_value = row.get("protein")
        protein = protein_value if isinstance(protein_value, str) and protein_value else "Unknown protein"
        mr_beta = row.get("mr_beta")
        pp_h4 = row.get("pp_h4_abf")
        mr_fdr = row.get("mr_fdr_q")
        snp = row.get("snp")
        a1 = row.get("a1")
        a2 = row.get("a2")
        mr_method = row.get("mr_method")
        n_instruments = row.get("n_instruments")

        with columns[position % n_columns]:
            with st.container(border=True):
                st.markdown(f"**{protein}**")

                if pd.notna(mr_beta) and mr_beta > 0:
                    # st.badge(f"Raises {outcome} risk", icon="⬆️", color="red")
                    st.badge(f"Raises {outcome} risk", color="red")
                elif pd.notna(mr_beta) and mr_beta < 0:
                    # st.badge(f"Lowers {outcome} risk", icon="⬇️", color="green")
                    st.badge(f"Lowers {outcome} risk", color="green")
                else:
                    st.badge("Direction unavailable", color="gray")

                # stacked, not side-by-side st.columns(2) - 2 badges sharing a
                # card's half-width truncates ("PP.H4 0.9..." / "MR FD...") at
                # anything narrower than a maximised wide monitor; full card
                # width per badge is the only layout that's robust to that
                pp_label = f"PP.H4 = {pp_h4:.2f}" if pd.notna(pp_h4) else "PP.H4 = NA"
                st.badge(pp_label, color=probability_strength_color(pp_h4, threshold=pp4))

                fdr_label = f"MR FDR = {mr_fdr:.1e}" if pd.notna(mr_fdr) else "MR FDR = NA"
                st.badge(fdr_label, color=significance_strength_color(mr_fdr))

                coloc_support = row.get("coloc_support")

                if isinstance(coloc_support, str) and coloc_support:
                    support_label, support_color = {
                        "both": ("COLOC + PWCoCo", "green"),
                        "coloc_only": ("COLOC only", "blue"),
                        "pwcoco_only": ("PWCoCo only", "orange"),
                    }.get(coloc_support, (coloc_support, "gray"))
                    st.badge(support_label, color=support_color)

                method_bits = []

                if isinstance(mr_method, str) and mr_method:
                    method_bits.append(mr_method)

                if pd.notna(n_instruments):
                    method_bits.append(f"{int(n_instruments)} instrument(s)")

                if pd.notna(mr_beta):
                    method_bits.append(f"beta = {mr_beta:.3f}")

                if method_bits:
                    st.caption(" · ".join(method_bits))

                snp_label = snp if isinstance(snp, str) and snp else "NA"
                allele_label = f"{a1}/{a2}" if isinstance(a1, str) and isinstance(a2, str) else "?/?"
                st.caption(f"Top SNP: **{snp_label}** ({allele_label})")


def render_phewas_section(
    phewas_outcome: pd.DataFrame,
    mr_coloc_pass: pd.DataFrame,
    mr_outcome: pd.DataFrame,
    source_name: str,
    source_description: str,
    pqtl_dataset: str,
    outcome: str,
    key_prefix: str,
    stage_number: int,
    is_fallback: bool = False,
):
    stage_caption(stage_number)
    st.subheader(f"{source_name} phenome-wide MR · indications & adverse effects")

    if is_fallback:
        st.caption(
            f"{source_name} phenome-wide MR is a FALLBACK: it only runs for targets where "
            "none of the retained cis-MR instruments were available in FinnGen."
        )

    with st.expander("Classification methodology", expanded=False):
        st.markdown(
            "- Cis-pQTL instruments are aligned to the **protein abundance-increasing "
            "allele**; the Wald ratio is used for 1 retained instrument, IVW for ≥ 2 "
            "(delta-method SEs).\n"
            "- Bonferroni correction is applied **per protein**, across the number of "
            "endpoints actually tested for that protein in this source "
            "(Pₐₑᵢ < 0.05 defines significance).\n"
            "- Only Bonferroni-significant binary-disease associations are classified:\n"
            "  - Higher protein abundance **increases** AD risk → inhibition indicated: "
            "a **positive** PheWAS estimate is a potential **additional indication**; a "
            "**negative** estimate is a potential **adverse effect**.\n"
            "  - Higher protein abundance **decreases** AD risk → augmentation "
            "indicated: a **negative** PheWAS estimate is a potential **additional "
            "indication**; a **positive** estimate is a potential **adverse effect**.\n"
            "- Non-significant associations are **not classified**."
        )

    if phewas_outcome.empty:
        fallback_note = " This is expected if every target had >=1 retained instrument in FinnGen already." if is_fallback else ""
        st.info(f"No local {source_name} PheWAS results are available for this outcome.{fallback_note}")
        return

    if "protein" not in phewas_outcome.columns:
        st.error(f"The {source_name} PheWAS result file does not contain a protein column.")
        return

    classified_outcome = classify_phewas_associations(phewas_outcome, mr_outcome)

    phewas_targets = sorted(classified_outcome["protein"].dropna().astype(str).unique())

    if len(phewas_targets) == 0:
        fallback_note = " This is expected if every target had >=1 retained instrument in FinnGen already." if is_fallback else ""
        st.info(f"No proteins were found in the {source_name} PheWAS table.{fallback_note}")
        return

    default_phewas_target = 0
    prioritised_target_names = mr_coloc_pass["protein"].dropna().astype(str).unique().tolist()

    for target in prioritised_target_names:
        if target in phewas_targets:
            default_phewas_target = phewas_targets.index(target)
            break

    with st.container(border=True):
        target_col, target_info_col = st.columns([2, 1])

        with target_col:
            selected_phewas_target = st.selectbox(
                "Target",
                phewas_targets,
                index=default_phewas_target,
                key=f"{key_prefix}_selected_phewas_target"
            )

        with target_info_col:
            st.metric("Targets with PheWAS results", len(phewas_targets))

    target_phewas = classified_outcome[
        classified_outcome["protein"].astype(str) == selected_phewas_target
    ].copy()

    p_col = "p_mr" if "p_mr" in target_phewas.columns else None
    beta_col = "beta_mr" if "beta_mr" in target_phewas.columns else None

    if p_col is None or beta_col is None:
        st.error(
            f"The {source_name} PheWAS result file needs the MR effect column "
            "(beta_mr) and the MR p-value column (p_mr)."
        )
        return

    target_phewas = target_phewas[
        target_phewas[p_col].notna() &
        target_phewas[beta_col].notna() &
        (target_phewas[p_col] > 0)
    ].copy()

    if target_phewas.empty:
        st.info(f"No valid {source_name} PheWAS associations were found for {selected_phewas_target}.")
        return

    target_phewas["minus_log10_p"] = -np.log10(target_phewas[p_col])

    phenotype_col = "phenostring" if "phenostring" in target_phewas.columns else "phenocode"

    if phenotype_col not in target_phewas.columns:
        st.error(f"The {source_name} PheWAS result file does not contain a phenotype column.")
        return

    n_endpoints_target = (
        int(target_phewas["n_endpoints_tested"].dropna().iloc[0])
        if "n_endpoints_tested" in target_phewas.columns and target_phewas["n_endpoints_tested"].notna().any()
        else len(target_phewas)
    )
    n_nominal = int((target_phewas[p_col] < 0.05).sum())
    n_bonferroni = int(target_phewas["bonferroni_significant"].fillna(False).sum())
    n_indications = int((target_phewas["phewas_classification"] == "additional_indication").sum())
    n_adverse = int((target_phewas["phewas_classification"] == "adverse_effect").sum())

    # primary protein->AD direction, for framing every hit below
    primary_beta = (
        target_phewas["primary_mr_beta"].dropna().iloc[0]
        if "primary_mr_beta" in target_phewas.columns and target_phewas["primary_mr_beta"].notna().any()
        else None
    )

    if primary_beta is not None and primary_beta > 0:
        st.info(
            f"**{selected_phewas_target}**: higher protein abundance **increases** AD risk "
            "→ **inhibition** indicated. A same-direction (positive) PheWAS hit below is "
            "a potential **additional indication**; an opposite-direction (negative) hit is a "
            "potential **adverse effect**."
        )
    elif primary_beta is not None and primary_beta < 0:
        st.info(
            f"**{selected_phewas_target}**: higher protein abundance **decreases** AD risk "
            "→ **augmentation** indicated. A same-direction (negative) PheWAS hit below "
            "is a potential **additional indication**; an opposite-direction (positive) hit "
            "is a potential **adverse effect**."
        )
    else:
        st.warning(
            f"No primary cis-MR beta was found for {selected_phewas_target} in this outcome, "
            "so associations below cannot be classified as indication vs. adverse effect."
        )

    with st.container(border=True):
        metric1, metric2, metric3, metric4, metric5 = st.columns(5)
        metric1.metric(f"{source_name} endpoints tested", n_endpoints_target)
        metric2.metric("Nominal, p<0.05 (incl. Bonferroni-sig.)", n_nominal)
        metric3.metric("↳ Bonferroni-significant", n_bonferroni)
        metric4.metric("Additional indications", n_indications)
        metric5.metric("Adverse effects", n_adverse)

    st.caption(
        f"PheWAS MR estimates show the effect of genetically predicted protein abundance "
        f"on each {source_description}, aligned to the protein abundance-increasing allele "
        "throughout. Wald ratio is used for targets with one available cis-MR instrument, "
        "IVW for > 1. Both PheWAS sources here are binary disease/case-control endpoints, so "
        "effects are also reported as odds ratios with 95% CIs."
    )

    classification_labels = {
        "additional_indication": "Additional indication",
        "adverse_effect": "Adverse effect",
        "not_classified": "Not classified",
    }
    classification_colors = {
        "Additional indication": "#2e7d32",
        "Adverse effect": "#c62828",
        "Not classified": "#9e9e9e",
    }
    target_phewas["classification_label"] = target_phewas["phewas_classification"].map(classification_labels).fillna("Not classified")

    if {"or_mr", "or_ci_low", "or_ci_high"}.issubset(target_phewas.columns):
        target_phewas["or_display"] = target_phewas.apply(
            lambda row: (
                f"{row['or_mr']:.2f} ({row['or_ci_low']:.2f}–{row['or_ci_high']:.2f})"
                if pd.notna(row["or_mr"]) and pd.notna(row["or_ci_low"]) and pd.notna(row["or_ci_high"])
                else "NA"
            ),
            axis=1,
        )

    st.divider()
    st.subheader("Phenome-wide association landscape")

    plot_kwargs = {
        "data_frame": target_phewas,
        "x": beta_col,
        "y": "minus_log10_p",
        "hover_name": phenotype_col,
        "color": "classification_label",
        "color_discrete_map": classification_colors,
        "category_orders": {"classification_label": list(classification_colors.keys())},
        "hover_data": {
            beta_col: ":.4f",
            p_col: ":.3e",
            "minus_log10_p": False,
            "classification_label": True,
        },
        "labels": {
            beta_col: "PheWAS MR beta",
            "minus_log10_p": "-log10(PheWAS p-value)",
            "classification_label": "Classification",
        },
        "title": f"{source_name} PheWAS profile: {selected_phewas_target}",
        "height": 600,
        "template": "plotly_white"
    }

    if "phenocode" in target_phewas.columns:
        plot_kwargs["hover_data"]["phenocode"] = True

    phewas_fig = px.scatter(**plot_kwargs)
    phewas_fig.add_hline(y=-np.log10(0.05 / n_endpoints_target), line_dash="dash", line_color="grey")
    phewas_fig.add_vline(x=0, line_dash="dash", line_color="grey")
    st.plotly_chart(phewas_fig, width="stretch")

    classification_display_cols = {
        "phenostring": "Phenotype",
        "phenocode": "Code",
        "category": "Category",
        "beta_mr": "MR beta",
        "or_display": "OR (95% CI)",
        "p_mr": "MR p-value",
        "p_bonferroni": "Bonferroni p-value",
    }

    def _render_classification_table(df_subset, empty_message):
        if df_subset.empty:
            st.caption(empty_message)
            return
        cols = [c for c in classification_display_cols if c in df_subset.columns]
        st.dataframe(
            df_subset.sort_values("p_mr")[cols].rename(columns=classification_display_cols),
            width="stretch",
            hide_index=True,
        )

    st.divider()
    ind_col, adv_col = st.columns(2)

    with ind_col:
        with st.container(border=True):
            st.badge(f"{n_indications} potential additional indication(s)", color="green")
            _render_classification_table(
                target_phewas[target_phewas["phewas_classification"] == "additional_indication"],
                f"No same-direction Bonferroni-significant {source_name} hits for this target.",
            )

    with adv_col:
        with st.container(border=True):
            st.badge(f"{n_adverse} potential adverse effect(s)", color="red")
            _render_classification_table(
                target_phewas[target_phewas["phewas_classification"] == "adverse_effect"],
                f"No opposite-direction Bonferroni-significant {source_name} hits for this target.",
            )

    full_cols = [
        "protein",
        "method",
        "n_instruments_original",
        "n_instruments_available",
        "n_instruments",
        "instrument_completeness",
        "missing_instruments",
        "rsid",
        "snp",
        "a1",
        "a2",
        "ukb_ref",
        "ukb_alt",
        "phenocode",
        "phenostring",
        "category",
        "beta_mr",
        "se_mr",
        "p_mr",
        "n_endpoints_tested",
        "p_bonferroni",
        "bonferroni_significant",
        "or_display",
        "primary_mr_beta",
        "phewas_classification",
    ]
    full_cols = available_cols(target_phewas, full_cols)

    full_column_names = {
        "protein": "Protein",
        "method": "Method",
        "n_instruments_original": "N instruments (original)",
        "n_instruments_available": "N instruments (available)",
        "n_instruments": "N instruments used",
        "instrument_completeness": "Instrument completeness",
        "missing_instruments": "Missing instruments",
        "rsid": "RSID",
        "snp": "SNP",
        "a1": "Effect allele",
        "a2": "Other allele",
        "ukb_ref": "UKB reference allele",
        "ukb_alt": "UKB alternate allele",
        "phenocode": "Phenotype code",
        "phenostring": "Phenotype",
        "category": "Category",
        "beta_mr": "MR beta",
        "se_mr": "MR SE",
        "p_mr": "MR p-value",
        "n_endpoints_tested": "Endpoints tested (this protein)",
        "p_bonferroni": "Bonferroni p-value",
        "bonferroni_significant": "Bonferroni significant",
        "or_display": "OR (95% CI)",
        "primary_mr_beta": "Primary cis-MR beta (protein→AD)",
        "phewas_classification": "Classification",
    }

    with st.expander(f"View all {source_name} PheWAS associations ({len(target_phewas)} endpoints)"):
        st.dataframe(
            target_phewas[full_cols].sort_values(p_col, ascending=True).rename(columns=full_column_names),
            width="stretch",
            hide_index=True
        )

    # same curated, renamed columns as the "View all" table above - not the raw
    # frame, which still carries internal helper columns (minus_log10_p,
    # classification_label, beta_ci_low/high) the on-screen table never shows
    st.download_button(
        label=f"Download {selected_phewas_target} {source_name} PheWAS results",
        data=target_phewas[full_cols].sort_values(p_col, ascending=True).rename(columns=full_column_names).to_csv(index=False, sep="\t"),
        file_name=f"{selected_phewas_target}_{outcome}_{source_name.replace(' ', '_')}_PheWAS.tsv",
        mime="text/tab-separated-values",
        key=f"{key_prefix}_download_phewas_{pqtl_dataset}_{outcome}_{selected_phewas_target}",
        width="stretch"
    )


def render_target_profile(
    protein: str,
    pqtl_dataset: str,
    mr_outcome: pd.DataFrame,
    mr_pass_proteins: set,
    coloc_outcome: pd.DataFrame,
    pwcoco_outcome: pd.DataFrame,
    coloc_support_status: dict,
    pp4: float,
    finngen_phewas_outcome: pd.DataFrame,
    ukb_phewas_outcome: pd.DataFrame,
    smr_display: pd.DataFrame,
    hyprcoloc_display: pd.DataFrame,
    smr_fdr_threshold: float,
    heidi_p_threshold: float,
    hyprcoloc_pp_threshold: float,
):
    """1 target's complete evidence trail, stage by stage, in a single vertical
    read instead of checking each stage's tab separately - the dashboard's other
    tabs are organised by pipeline stage (1 table per stage, every protein at
    once), which is how the pipeline is built but not how a geneticist actually
    wants to browse it ("show me everything about this protein")."""

    support_labels = {
        "both": ("Supported by both COLOC + PWCoCo", "green"),
        "coloc_only": ("Supported by standard COLOC only", "blue"),
        "pwcoco_only": ("Supported by PWCoCo only", "orange"),
    }

    # --- resolve every stage's pass/fail up front, so the verdict banner at the
    # top can be computed before the stage-by-stage detail below it ---
    mr_rows = mr_outcome[mr_outcome["protein"].astype(str) == protein] if "protein" in mr_outcome.columns else mr_outcome.iloc[0:0]
    passed_mr = protein in mr_pass_proteins

    coloc_row = coloc_outcome[coloc_outcome["protein"].astype(str) == protein] if "protein" in coloc_outcome.columns else coloc_outcome.iloc[0:0]
    pwcoco_rows = pwcoco_outcome[pwcoco_outcome["protein"].astype(str) == protein] if "protein" in pwcoco_outcome.columns else pwcoco_outcome.iloc[0:0]
    support = coloc_support_status.get(protein)
    passed_coloc_stage = support in support_labels

    finngen_status = compute_phewas_classification_status(finngen_phewas_outcome, mr_outcome, [protein]).get(protein, "none")
    ukb_status = compute_phewas_classification_status(ukb_phewas_outcome, mr_outcome, [protein]).get(protein, "none")
    passed_safety = passed_coloc_stage and finngen_status != "adverse_effect" and ukb_status != "adverse_effect"

    smr_rows = smr_display[smr_display["protein"].astype(str) == protein] if "protein" in smr_display.columns else smr_display.iloc[0:0]
    smr_pass_rows = (
        smr_rows[
            smr_rows["q_smr"].notna() & (smr_rows["q_smr"] < smr_fdr_threshold) &
            smr_rows["p_heidi"].notna() & (smr_rows["p_heidi"] > heidi_p_threshold)
        ]
        if {"q_smr", "p_heidi"}.issubset(smr_rows.columns) else smr_rows.iloc[0:0]
    )
    has_smr_support = passed_safety and not smr_pass_rows.empty

    hypr_rows = hyprcoloc_display[hyprcoloc_display["protein"].astype(str) == protein] if "protein" in hyprcoloc_display.columns else hyprcoloc_display.iloc[0:0]
    passed_hyprcoloc = has_smr_support and compute_hyprcoloc_pass_status(hypr_rows, [protein], hyprcoloc_pp_threshold).get(protein, False)

    # --- verdict banner ---
    if passed_hyprcoloc:
        st.success("FINAL TARGET · passed every stage, including HyPrColoc's 3-trait colocalisation check.")
    elif has_smr_support:
        st.info("Reached SMR support, but did not clear HyPrColoc's 3-trait colocalisation threshold.")
    elif passed_safety:
        st.success("PRIORITISED · passed cis-MR, colocalisation and PheWAS safety. No SMR/eQTL support found or tested yet.")
    elif passed_coloc_stage:
        st.warning("ADVERSE EFFECT FLAG · a Bonferroni-significant FinnGen or UKB PheWAS hit runs opposite to the primary protein→AD effect direction - excluded from Prioritised Targets.")
    elif passed_mr:
        st.error("STOPPED AT COLOCALISATION · passed cis-MR, but neither standard COLOC nor PWCoCo cleared the PP.H4 threshold.")
    else:
        st.error("STOPPED AT cis-MR · did not clear the MR FDR / Cochran Q thresholds.")

    st.divider()

    # --- Stage 1: cis-MR ---
    st.markdown("#### Stage 1 · cis-MR")
    with st.container(border=True):
        st.badge("PASSED" if passed_mr else "FAILED", color="green" if passed_mr else "red")
        if mr_rows.empty:
            st.caption("No cis-MR row found for this target.")
        else:
            row = mr_rows.iloc[0]
            c1, c2, c3, c4 = st.columns(4)
            c1.metric("Method", row.get("mr_method", "NA") if pd.notna(row.get("mr_method")) else "NA")
            c2.metric("Instruments", int(row["n_instruments"]) if pd.notna(row.get("n_instruments")) else "NA")
            c3.metric("MR beta", f"{row['mr_beta']:.3f}" if pd.notna(row.get("mr_beta")) else "NA")
            c4.metric("MR FDR q", f"{row['mr_fdr_q']:.2e}" if pd.notna(row.get("mr_fdr_q")) else "NA")

    # --- Stage 2: standard COLOC + PWCoCo ---
    st.markdown("#### Stage 2 · pQTL–GWAS colocalisation")
    with st.container(border=True):
        if not passed_mr:
            st.caption("Not reached - target did not pass cis-MR.")
        else:
            col1, col2 = st.columns(2)

            with col1:
                st.caption("Standard COLOC")
                if coloc_row.empty or "pp_h4_abf" not in coloc_row.columns or pd.isna(coloc_row.iloc[0].get("pp_h4_abf")):
                    st.badge("No COLOC result", color="gray")
                else:
                    coloc_pp4 = coloc_row.iloc[0]["pp_h4_abf"]
                    st.badge(f"PP.H4 = {coloc_pp4:.3f}", color="green" if coloc_pp4 >= pp4 else "red")

            with col2:
                st.caption("PWCoCo (conditional coloc)")
                if pwcoco_rows.empty or "h4" not in pwcoco_rows.columns:
                    st.badge("No PWCoCo result", color="gray")
                else:
                    max_h4 = pd.to_numeric(pwcoco_rows["h4"], errors="coerce").max()
                    if pd.isna(max_h4):
                        st.badge("No PWCoCo result", color="gray")
                    else:
                        st.badge(
                            f"Best PP.H4 = {max_h4:.3f} ({len(pwcoco_rows)} signal(s) tested)",
                            color="green" if max_h4 >= pp4 else "red"
                        )

            st.divider()
            label, color = support_labels.get(support, ("Neither method colocalised", "red"))
            st.badge(label, color=color)

    if passed_mr:
        st.markdown("##### Regional association plot")
        render_regional_locus_plot(
            protein=protein,
            pqtl_dataset=pqtl_dataset,
            smr_rows=smr_rows,
            hypr_rows=hypr_rows,
            key_prefix="target_profile",
        )

    # --- Stage 3 / 4: phenome-wide MR (FinnGen primary, UKB fallback) ---
    st.markdown("#### Stage 3 · FinnGen phenome-wide MR")
    with st.container(border=True):
        if not passed_coloc_stage:
            st.caption("Not reached.")
        elif finngen_status == "adverse_effect":
            st.badge("ADVERSE EFFECT", color="red")
            st.caption("Bonferroni-significant FinnGen hit running opposite to the primary protein→AD direction.")
        elif finngen_status == "additional_indication":
            st.badge("ADDITIONAL INDICATION", color="green")
            st.caption("Bonferroni-significant FinnGen hit running the same direction as the primary protein→AD effect - a potential repurposing signal, not a safety concern.")
        else:
            st.badge("NO SIGNIFICANT SIGNAL", color="green")

    st.markdown("#### Stage 4 · UKB phenome-wide MR (fallback)")
    with st.container(border=True):
        if not passed_coloc_stage:
            st.caption("Not reached.")
        elif ukb_phewas_outcome.empty or ("protein" in ukb_phewas_outcome.columns and protein not in ukb_phewas_outcome["protein"].astype(str).unique()):
            st.badge("NOT RUN", color="gray")
            st.caption("UKB is only run when none of this target's retained cis-MR instruments were available in FinnGen.")
        elif ukb_status == "adverse_effect":
            st.badge("ADVERSE EFFECT", color="red")
            st.caption("Bonferroni-significant UKB hit running opposite to the primary protein→AD direction.")
        elif ukb_status == "additional_indication":
            st.badge("ADDITIONAL INDICATION", color="green")
            st.caption("Bonferroni-significant UKB hit running the same direction as the primary protein→AD effect - a potential repurposing signal, not a safety concern.")
        else:
            st.badge("NO SIGNIFICANT SIGNAL", color="green")

    # --- Stage 5: SMR ---
    st.markdown("#### Stage 5 · SMR (bulk/sc eQTL)")
    with st.container(border=True):
        if not passed_safety:
            st.caption("Not reached.")
        elif smr_rows.empty:
            st.badge("No SMR support", color="gray")
        else:
            display_cols = available_cols(
                smr_rows, ["eqtl_dataset", "data_type", "cell_type", "b_smr", "p_smr", "q_smr", "p_heidi"]
            )
            st.dataframe(smr_rows[display_cols], hide_index=True, width="stretch")
            if smr_pass_rows.empty:
                st.badge("No dataset clears the SMR FDR / HEIDI thresholds", color="red")
            else:
                n_datasets = smr_pass_rows["eqtl_dataset"].nunique() if "eqtl_dataset" in smr_pass_rows.columns else len(smr_pass_rows)
                st.badge(f"SMR support in {n_datasets} dataset(s)", color="green")

    # --- Stage 6: HyPrColoc ---
    st.markdown("#### Stage 6 · HyPrColoc")
    with st.container(border=True):
        if not has_smr_support:
            st.caption("Not reached - no SMR/eQTL support to test.")
        elif hypr_rows.empty:
            st.badge("No HyPrColoc result", color="gray")
        else:
            st.badge("PASSED" if passed_hyprcoloc else "FAILED", color="green" if passed_hyprcoloc else "red")
            display_cols = available_cols(hypr_rows, ["cell_type", "traits", "posterior_prob", "candidate_snp"])
            st.dataframe(hypr_rows[display_cols], hide_index=True, width="stretch")


def dashboard(db_name: str, port_number: str, phenotype: str, pqtl_dataset: str):
    mr_table = "cis_mr_results"
    coloc_table = "coloc_results"
    finngen_phewas_table = "finngen_phewas_safety"
    ukb_phewas_table = "ukb_phewas_safety"

    # main aesthetics
    st.set_page_config(
        page_title=f"{db_name}",
        page_icon="",
        layout="wide",
        initial_sidebar_state="expanded"
    )

    # custom stylesheet on top of .streamlit/config.toml's base theme - only
    # targets data-testid/data-baseweb attributes and other stable ARIA/role
    # attributes Streamlit keeps stable across versions for theming/testing
    # (not the auto-generated st-emotion-cache-* build hashes, which change per
    # build and would silently stop matching on any upgrade). Note: st.dataframe
    # renders its grid on an HTML canvas (glide-data-grid), not real <table>/<th>
    # DOM - its internal header/cell styling can't be reached from CSS at all;
    # it follows .streamlit/config.toml's theme directly instead, which is the
    # real reason getting that file actually picked up (see drugmr/local.py's
    # results() / drugmr/hpc.py's run_dashboard_local(), both now pin cwd to
    # project_root for exactly this) matters more here than any CSS rule could.
    st.markdown(
        """
        <style>
        :root { --accent: #0E7C86; }

        /* tighter, calmer vertical rhythm - the default stacks widgets with a
           lot of dead air, which is what makes a data-dense page feel busy */
        div[data-testid="stVerticalBlock"] { gap: 0.6rem; }

        /* headers: a little more weight for a cleaner hierarchy */
        h1, h2, h3 { font-weight: 650; letter-spacing: -0.01em; }

        /* bordered containers (st.container(border=True)) - used throughout
           for stage cards - get real depth instead of a flat grey outline */
        div[data-testid="stVerticalBlockBorderWrapper"] {
            border-radius: 12px !important;
            box-shadow: 0 1px 4px rgba(16, 24, 32, 0.06);
        }

        /* metrics: bigger numbers read faster at a glance than the default size.
           Labels wrap onto a 2nd line instead of Streamlit's default silent
           ellipsis-truncation - a label that's too long for a narrow metric
           column (e.g. in a 4-up row) stays fully readable rather than cutting
           off mid-word. */
        div[data-testid="stMetricValue"] { font-size: 1.65rem; font-weight: 650; }
        div[data-testid="stMetricLabel"], div[data-testid="stMetricDelta"] {
            font-size: 0.82rem;
            opacity: 0.75;
            white-space: normal !important;
            overflow: visible !important;
            text-overflow: unset !important;
        }

        /* multiselect chips (e.g. the HyPrColoc cell-type/tissue filter) - the
           default size reads as a wall of pills when 10+ are selected at once;
           tighter padding and a smaller font let more fit per line */
        div[data-baseweb="tag"] {
            font-size: 0.78rem;
            padding-top: 0.05rem;
            padding-bottom: 0.05rem;
        }

        /* tabs: bolder labels, a clear accent-coloured underline on whichever
           tab is active, and more breathing room so a 5-wide tab bar doesn't
           feel cramped */
        button[data-baseweb="tab"] {
            font-weight: 600;
            padding-top: 0.55rem;
            padding-bottom: 0.55rem;
        }
        button[data-baseweb="tab"][aria-selected="true"] { color: var(--accent); }
        div[data-baseweb="tab-highlight"] { background-color: var(--accent); height: 3px; }

        /* sidebar: a visible seam from the main content, and tighter expander
           spacing so the 6 stage-grouped threshold sections read as 1 coherent
           panel rather than 6 disconnected boxes */
        section[data-testid="stSidebar"] {
            border-right: 1px solid rgba(16, 24, 32, 0.08);
        }
        section[data-testid="stSidebar"] div[data-testid="stExpander"] {
            margin-bottom: 0.35rem;
        }

        /* alert boxes (st.info/success/warning/error) - slightly rounder to
           match the card styling above instead of Streamlit's sharper default */
        div[data-testid="stAlert"] { border-radius: 10px; }

        /* dataframes/tables: only the outer wrapper is real DOM (the grid
           itself is a canvas - glide-data-grid - so its internal header/cell
           styling can't be reached from CSS; it follows the app theme
           directly once .streamlit/config.toml actually loads). Rounding +
           a thin border here just makes the wrapper match the card language
           used everywhere else instead of Streamlit's flat default edge. */
        div[data-testid="stDataFrame"], div[data-testid="stTable"] {
            border-radius: 10px;
            overflow: hidden;
            border: 1px solid rgba(16, 24, 32, 0.08);
        }

        /* slider thumb - previously confirmed live (getComputedStyle) that this
           rendered Streamlit's hardcoded default red (#FF4B4B) instead of
           .streamlit/config.toml's primaryColor, because the theme file wasn't
           being picked up at all (a CWD issue - Streamlit only finds
           .streamlit/config.toml relative to the directory `streamlit run` is
           invoked FROM, and drugmr/local.py's results() / drugmr/hpc.py's
           run_dashboard_local() launched it without pinning that directory).
           Now fixed at the source (both launchers pin cwd=project_root), so
           this CSS block is redundant defense-in-depth, not the real fix -
           safe to keep since it just reinforces the same color the theme
           itself now sets. role="slider" is a stable ARIA attribute (unlike
           the st-emotion-cache-* build-hash classes elsewhere on this
           element), so it survives Streamlit upgrades either way. */
        div[data-testid="stSlider"] [role="slider"] {
            background-color: var(--accent) !important;
            border-color: var(--accent) !important;
        }
        div[data-testid="stSliderThumbValue"] { color: var(--accent) !important; }
        </style>
        """,
        unsafe_allow_html=True,
    )

    # pQTL dataset selection schema
    # CLI pQTL dataset is used as the default dashboard selection
    dataset_names = {
        "ukb_ppp": "UKB-PPP",
        "decode": "deCODE",
        "wu_csf": "WU-CSF",
        "wingo_brain": "Wingo_Brain"
    }

    dataset_ns = {
        "ukb_ppp": 54219,
        "decode": 35559,
        "wu_csf": 3506,
        "wingo_brain": 1013
    }

    project_dir = Path(__file__).resolve().parent.parent

    # check which datasets have the required dashboard files - resolved via
    # runs/registry.json first (see resolve_dataset_files()), falling back to
    # legacy candidate-path guessing for any dataset never run through runs/
    dataset_result_files = {}
    dataset_run_ids = {}
    available_datasets = []

    for dataset_id in dataset_names:
        run_id_used, files = resolve_dataset_files(project_dir, phenotype, dataset_id, run_id="latest")

        required_files = [files["mr"], files["coloc"]]

        if all(file is not None and file.exists() for file in required_files):
            available_datasets.append(dataset_id)
            dataset_result_files[dataset_id] = files
            dataset_run_ids[dataset_id] = run_id_used

    if len(available_datasets) == 0:
        st.error(f"No dataset has a complete set of cis-MR and COLOC dashboard files for {phenotype}.")
        st.stop()

    # use the CLI dataset as default
    # otherwise use the first complete dataset which was found
    if pqtl_dataset not in available_datasets:
        pqtl_dataset = available_datasets[0]

    available_dataset_names = [dataset_names[dataset_id] for dataset_id in available_datasets]
    default_dataset_name = dataset_names[pqtl_dataset]

    st.title(f"{db_name}")
    st.caption("Genetically supported drug target discovery and clinical safety dashboard")

    with st.container(border=True):
        dataset_col, run_col, dataset_info_col = st.columns([2, 1, 1])

        with dataset_col:
            selected_dataset_name = st.segmented_control(
                "pQTL dataset",
                available_dataset_names,
                default=default_dataset_name,
                selection_mode="single",
                key="pqtl_dataset_selector"
            )

        dataset_ids = {dataset_name: dataset_id for dataset_id, dataset_name in dataset_names.items()}
        if selected_dataset_name is None:
            selected_dataset_name = default_dataset_name

        pqtl_dataset = dataset_ids[selected_dataset_name]
        dataset_name = dataset_names[pqtl_dataset]
        dataset_n = dataset_ns[pqtl_dataset]

        # run selector - history comes from runs/registry.json; a dataset resolved
        # via legacy_resolve_dataset_files() (no registry entry) has no history at all
        run_history = registry.load_registry(root=str(project_dir / "runs")).get(
            f"{phenotype}__{pqtl_dataset}", {}
        ).get("history", [])

        with run_col:
            if run_history:
                run_options = ["latest"] + list(reversed(run_history))
                selected_run = st.selectbox("Run", run_options, index=0, key="run_selector")
            else:
                selected_run = "latest"
                st.caption("No run history (legacy path)")

        if selected_run != "latest":
            run_id_used, dataset_result_files[pqtl_dataset] = resolve_dataset_files(
                project_dir, phenotype, pqtl_dataset, run_id=selected_run
            )
            dataset_run_ids[pqtl_dataset] = run_id_used

        with dataset_info_col:
            st.metric("pQTL sample size", f"{dataset_n:,}")

    st.divider()

    # corresponding selected dataset result files
    mr_file = dataset_result_files[pqtl_dataset]["mr"]
    coloc_file = dataset_result_files[pqtl_dataset]["coloc"]
    finngen_phewas_file = dataset_result_files[pqtl_dataset]["finngen_phewas"]
    ukb_phewas_file = dataset_result_files[pqtl_dataset]["ukb_phewas"]
    target_info_file = dataset_result_files[pqtl_dataset]["target_info"]
    smr_file = dataset_result_files[pqtl_dataset]["smr"]
    hyprcoloc_file = dataset_result_files[pqtl_dataset]["hyprcoloc"]
    pwcoco_file = dataset_result_files[pqtl_dataset]["pwcoco"]
    pwcoco_eqtl_pqtl_file = dataset_result_files[pqtl_dataset]["pwcoco_eqtl_pqtl"]
    pwcoco_eqtl_gwas_file = dataset_result_files[pqtl_dataset]["pwcoco_eqtl_gwas"]

    # push this run's result files straight into PostgreSQL (schema-matching
    # since bin/coloc_targets.py, bin/sort_smr.py, bin/pwcoco_wrapper.py,
    # bin/pwcoco_qtl_wrapper.py and bin/compile_cis_hit_info.py all emit
    # sql/schema.sql-shaped columns now) - scoped to this dataset's actual
    # run_id, not the dashboard's own in-memory transforms further below
    run_id = dataset_run_ids[pqtl_dataset]

    postgres_tables = [
        ("cis_mr_results", mr_file, True),
        ("coloc_results", coloc_file, True),
        ("finngen_phewas_safety", finngen_phewas_file, False),
        ("ukb_phewas_safety", ukb_phewas_file, False),
        ("target_stats", target_info_file, False),
        ("smr_results", smr_file, False),
        ("hyprcoloc_results", hyprcoloc_file, False),
        ("pwcoco_results", pwcoco_file, False),
        ("pwcoco_eqtl_pqtl_results", pwcoco_eqtl_pqtl_file, False),
        ("pwcoco_eqtl_gwas_results", pwcoco_eqtl_gwas_file, False),
    ]

    loader = PostgresLoader(run_id=run_id, db_id=db_name)
    postgres_table_available = {}

    for table, file_path, required in postgres_tables:
        if file_path is None or not file_path.exists():
            postgres_table_available[table] = False
            continue
        try:
            loader.load_table(results_file=file_path, pqtl_dataset=pqtl_dataset, table=table)
            postgres_table_available[table] = True
        except Exception as error:
            postgres_table_available[table] = False
            if required:
                st.error(f"The {table} dashboard table could not be refreshed.")
                st.exception(error)
                st.stop()
            else:
                st.warning(f"The {table} dashboard table could not be refreshed.")
                st.exception(error)

    # load local result files into PostgreSQL for the dashboard
    mr = load_required_tsv(mr_file, "cis-MR")
    coloc = load_required_tsv(coloc_file, "pQTL–GWAS COLOC")
    finngen_phewas = load_optional_tsv(finngen_phewas_file, "FinnGen PheWAS safety")
    ukb_phewas = load_optional_tsv(ukb_phewas_file, "UKB PheWAS safety")
    target_info = load_optional_tsv(target_info_file, "Harmonised target information")
    smr = load_optional_tsv(smr_file, "SMR (bulk/sc eQTL)")
    hyprcoloc = load_optional_tsv(hyprcoloc_file, "HyPrColoc (bulk/sc eQTL)")
    pwcoco = load_optional_tsv(pwcoco_file, "PWCoCo (conditional coloc)")
    pwcoco_eqtl_pqtl = load_optional_tsv(pwcoco_eqtl_pqtl_file, "PWCoCo (eQTL-pQTL)")
    pwcoco_eqtl_gwas = load_optional_tsv(pwcoco_eqtl_gwas_file, "PWCoCo (eQTL-GWAS)")

    # standardise MR + pQTL COLOC columns before loading into PostgreSQL
    # avoids dataset-specific differences such as Wald_beta vs wald_beta
    mr = standardise_columns(mr)
    coloc = standardise_columns(coloc)

    if not target_info.empty:
        target_info = standardise_columns(target_info)

    if not smr.empty:
        smr = standardise_columns(smr)

    if not hyprcoloc.empty:
        hyprcoloc = standardise_columns(hyprcoloc)

    if not pwcoco.empty:
        pwcoco = standardise_columns(pwcoco)

    if not pwcoco_eqtl_pqtl.empty:
        pwcoco_eqtl_pqtl = standardise_columns(pwcoco_eqtl_pqtl)

    if not pwcoco_eqtl_gwas.empty:
        pwcoco_eqtl_gwas = standardise_columns(pwcoco_eqtl_gwas)

    # make protein column consistent before loading into PostgreSQL
    if "protein_id" in mr.columns:
        mr = mr.rename(columns={"protein_id": "protein"})

    if "protein_id" in coloc.columns:
        coloc = coloc.rename(columns={"protein_id": "protein"})

    if not target_info.empty and "protein_id" in target_info.columns:
        target_info = target_info.rename(columns={"protein_id": "protein"})

    if not smr.empty and "protein_id" in smr.columns:
        smr = smr.rename(columns={"protein_id": "protein"})

    if not hyprcoloc.empty and "protein_id" in hyprcoloc.columns:
        hyprcoloc = hyprcoloc.rename(columns={"protein_id": "protein"})

    # make sure the selected dataset is always recorded
    if "pqtl_dataset" not in mr.columns:
        mr["pqtl_dataset"] = pqtl_dataset

    if "pqtl_dataset" not in coloc.columns:
        coloc["pqtl_dataset"] = pqtl_dataset


    finngen_phewas_available = postgres_table_available["finngen_phewas_safety"]
    ukb_phewas_available = postgres_table_available["ukb_phewas_safety"]

    # values captured here (this is where mr/coloc/etc. hold the freshly-loaded
    # TSV row counts, before later stages filter/transform them) - actual
    # rendering is deferred to the bottom of the sidebar (near the Legend), since
    # this is pure engineering telemetry with zero value to the geneticist this
    # dashboard is for; it doesn't need to be the first thing anyone sees
    tracking_info = {
        "mr_rows": len(mr),
        "coloc_rows": len(coloc),
        "finngen_phewas_rows": len(finngen_phewas) if finngen_phewas_available else None,
        "ukb_phewas_rows": len(ukb_phewas) if ukb_phewas_available else None,
        "target_info_rows": len(target_info) if not target_info.empty else None,
        "smr_rows": len(smr) if not smr.empty else None,
        "hyprcoloc_rows": len(hyprcoloc) if not hyprcoloc.empty else None,
        "pwcoco_rows": len(pwcoco) if not pwcoco.empty else None,
        "mr_table": mr_table,
        "coloc_table": coloc_table,
        "finngen_phewas_table": finngen_phewas_table,
        "ukb_phewas_table": ukb_phewas_table,
    }

    # load MR + COLOC results
    reader = PostgresReader(run_id=run_id, db_id=db_name)
    mr = reader.get_table(mr_table)
    coloc = reader.get_table(coloc_table)

    if finngen_phewas_available:
        finngen_phewas = prepare_phewas(reader.get_table(finngen_phewas_table))
    else:
        finngen_phewas = pd.DataFrame()

    if ukb_phewas_available:
        ukb_phewas = prepare_phewas(reader.get_table(ukb_phewas_table))
    else:
        ukb_phewas = pd.DataFrame()

    # MR ammenities
    # standardise numeric MR columns
    mr_numeric_cols = [
        "n_instruments",
        "ivw_beta",
        "ivw_se",
        "ivw_pval",
        "ivw_fdr_q",
        "wald_beta",
        "wald_se",
        "wald_pval",
        "wald_fdr_q",
        "q_pval",
        "egger_intercept_pval"
    ]

    for col in mr_numeric_cols:
        if col not in mr.columns:
            mr[col] = np.nan

        mr[col] = pd.to_numeric(mr[col], errors="coerce")

    coloc_numeric_cols = [
        "pp_h0_abf",
        "pp_h1_abf",
        "pp_h2_abf",
        "pp_h3_abf",
        "pp_h4_abf"
    ]

    for col in coloc_numeric_cols:
        if col in coloc.columns:
            coloc[col] = pd.to_numeric(coloc[col], errors="coerce")

    if "protein_id" in mr.columns:
        mr = mr.rename(columns={"protein_id": "protein"})

    if "protein_id" in coloc.columns:
        coloc = coloc.rename(columns={"protein_id": "protein"})

    if not target_info.empty:
        target_numeric_cols = [
            "frq",
            "gwas_beta",
            "gwas_se",
            "gwas_p",
            "pqtl_beta",
            "pqtl_se",
            "pqtl_p"
        ]

        for col in target_numeric_cols:
            if col in target_info.columns:
                target_info[col] = pd.to_numeric(target_info[col], errors="coerce")

        for col in ["a1", "a2"]:
            if col in target_info.columns:
                target_info[col] = target_info[col].astype(str).str.upper()

        if "protein" in target_info.columns:
            target_info = target_info.drop_duplicates(subset=["protein"])

    required_mr_cols = ["protein", "outcome_trait", "n_instruments"]
    missing_mr_cols = [col for col in required_mr_cols if col not in mr.columns]

    if len(missing_mr_cols) > 0:
        st.error(
            f"cis-MR result file is missing required columns: {missing_mr_cols}. "
            f"File: {mr_file}"
        )
        st.stop()

    required_coloc_cols = ["protein", "outcome_trait"]
    missing_coloc_cols = [col for col in required_coloc_cols if col not in coloc.columns]

    if len(missing_coloc_cols) > 0:
        st.error(
            f"COLOC result file is missing required columns: {missing_coloc_cols}. "
            f"File: {coloc_file}"
        )
        st.stop()

    # if 1 instrument -> use Wald
    # otherwise -> use IVW
    mr["mr_method"] = np.where(mr["n_instruments"] == 1, "Wald ratio", "IVW")
    mr["mr_beta"] = np.where(mr["n_instruments"] == 1, mr["wald_beta"], mr["ivw_beta"])
    mr["mr_se"] = np.where(mr["n_instruments"] == 1, mr["wald_se"], mr["ivw_se"])
    mr["mr_pval"] = np.where(mr["n_instruments"] == 1, mr["wald_pval"], mr["ivw_pval"])
    mr["mr_fdr_q"] = np.where(mr["n_instruments"] == 1, mr["wald_fdr_q"], mr["ivw_fdr_q"])

    # standardise selected pQTL dataset
    selected_pqtl_dataset = (
        str(pqtl_dataset)
        .strip()
        .lower()
        .replace("-", "_")
        .replace(" ", "_")
    )

    # subset database tables to the selected pQTL dataset where possible
    for dataset_col in ["pqtl_dataset", "dataset"]:
        if dataset_col in mr.columns:
            mr_dataset_values = (
                mr[dataset_col]
                .astype(str)
                .str.strip()
                .str.lower()
                .str.replace("-", "_", regex=False)
                .str.replace(" ", "_", regex=False)
            )

            mr = mr[mr_dataset_values == selected_pqtl_dataset].copy()
            break

    for dataset_col in ["pqtl_dataset", "dataset"]:
        if dataset_col in coloc.columns:
            coloc_dataset_values = (
                coloc[dataset_col]
                .astype(str)
                .str.strip()
                .str.lower()
                .str.replace("-", "_", regex=False)
                .str.replace(" ", "_", regex=False)
            )

            coloc = coloc[coloc_dataset_values == selected_pqtl_dataset].copy()
            break

    # available outcomes and default CLI phenotype
    outcomes = sorted(mr["outcome_trait"].dropna().unique())

    if len(outcomes) == 0:
        st.error(f"No cis-MR results were found in {mr_file} for {dataset_name}.")
        st.stop()

    default_outcome = outcomes.index(phenotype) if phenotype in outcomes else 0

    # sidebar filters
    st.sidebar.title("Dashboard controls")
    st.sidebar.caption(f"{dataset_name} → {phenotype}")
    st.sidebar.divider()

    with st.sidebar.expander("Analysis selection", expanded=True):
        outcome = st.selectbox("Outcome", outcomes, index=default_outcome)

    with st.sidebar.expander("Stage 1 · cis-MR", expanded=True):
        fdr = st.slider(
            "MR FDR threshold", 0.0, 1.0, 0.05, 0.01,
            help="A target passes cis-MR when its FDR-adjusted MR p-value is at or below this."
        )
        q_pval = st.slider(
            "Minimum Cochran Q p-value", 0.0, 1.0, 0.05, 0.01,
            help="IVW targets only (2+ instruments) - guards against instrument heterogeneity. Wald-ratio (1-instrument) targets are unaffected."
        )

    with st.sidebar.expander("Stage 2 · Colocalisation", expanded=True):
        pp4 = st.slider(
            "PP.H4 threshold", 0.0, 1.0, 0.70, 0.01,
            help="Shared by BOTH standard COLOC and PWCoCo - a target clears this stage if EITHER method's PP.H4 is at or above this bar."
        )

    st.sidebar.caption("Stage 3-4 · PheWAS (FinnGen/UKB) has no threshold here - it uses a fixed per-protein Bonferroni cutoff, not a configurable slider.")

    with st.sidebar.expander("Stage 5 · SMR / HEIDI", expanded=False):
        smr_fdr_threshold = st.slider("SMR FDR (q_SMR) threshold", 0.0, 1.0, 0.05, 0.01)
        heidi_p_threshold = st.slider(
            "Minimum HEIDI p-value", 0.0, 1.0, 0.01, 0.01,
            help="Guards against SMR support being driven by linkage rather than a truly shared causal variant."
        )

    with st.sidebar.expander("Stage 6 · HyPrColoc", expanded=False):
        hyprcoloc_pp_threshold = st.slider("HyPrColoc posterior probability threshold", 0.0, 1.0, 0.70, 0.01)

    with st.sidebar.expander("Protein filter", expanded=True):
        protein = st.text_input(
            "Protein search",
            placeholder="e.g. PILRA, APOE or Q9UKJ1",
            help="Filters every tab at once. For a single protein's full story without affecting other tabs, use the Target Profile tab instead."
        )

    with st.sidebar.expander("Legend", expanded=False):
        st.badge("Passed / safe / supported / both methods agree", color="green")
        st.badge("Failed / safety flag", color="red")
        st.badge("Not reached / no data", color="gray")
        st.caption(
            "Blue and orange repeat at 3 different \"which of 2 methods supported "
            "this\" splits in the pipeline - same 2 colours each time, but a "
            "different pair of methods depending on where you see them: standard "
            "COLOC vs. PWCoCo (Stage 2), bulk vs. single-cell eQTL (SMR, Stage 5), "
            "or HyPrColoc vs. PWCoCo-QTL (Stage 6). Hover a Sankey node's tooltip "
            "to see exactly which pair applies."
        )
        st.badge("1st method only (blue)", color="blue")
        st.badge("2nd method only (orange)", color="orange")

    with st.sidebar.expander("Tracking", expanded=False):
        st.write(f"Loaded {tracking_info['mr_rows']} rows into {tracking_info['mr_table']}")
        st.write(f"Loaded {tracking_info['coloc_rows']} rows into {tracking_info['coloc_table']}")

        if tracking_info["finngen_phewas_rows"] is not None:
            st.write(f"Loaded {tracking_info['finngen_phewas_rows']} rows into {tracking_info['finngen_phewas_table']}")

        if tracking_info["ukb_phewas_rows"] is not None:
            st.write(f"Loaded {tracking_info['ukb_phewas_rows']} rows into {tracking_info['ukb_phewas_table']}")

        if tracking_info["target_info_rows"] is not None:
            st.write(f"Loaded {tracking_info['target_info_rows']} harmonised top cis-hit rows")

        if tracking_info["smr_rows"] is not None:
            st.write(f"Loaded {tracking_info['smr_rows']} SMR (bulk/sc eQTL) rows")

        if tracking_info["hyprcoloc_rows"] is not None:
            st.write(f"Loaded {tracking_info['hyprcoloc_rows']} HyPrColoc (bulk/sc eQTL) rows")

        if tracking_info["pwcoco_rows"] is not None:
            st.write(f"Loaded {tracking_info['pwcoco_rows']} PWCoCo (conditional coloc) rows")

    st.header(f"{dataset_name} → {outcome}")
    st.caption(
        f"N = {dataset_n:,} | MR FDR ≤ {fdr:.2f} | Q p ≥ {q_pval:.2f} | "
        f"pQTL–GWAS PP.H4 ≥ {pp4:.2f} | SMR FDR ≤ {smr_fdr_threshold:.2f} | "
        f"HEIDI p ≥ {heidi_p_threshold:.2f} | HyPrColoc PP ≥ {hyprcoloc_pp_threshold:.2f}"
    )

    # subset everything to selected outcome
    mr_outcome = mr[mr["outcome_trait"] == outcome].copy()
    coloc_outcome = coloc[coloc["outcome_trait"] == outcome].copy()

    finngen_phewas_outcome = subset_phewas_outcome(finngen_phewas, outcome)
    ukb_phewas_outcome = subset_phewas_outcome(ukb_phewas, outcome)

    # STAGE 1
    # cis-MR supported proteins
    mr_pass = mr_outcome.copy()

    if "mr_fdr_q" in mr_pass.columns:
        mr_pass = mr_pass[mr_pass["mr_fdr_q"].fillna(np.inf) <= fdr]

    # apply Cochran Q only to IVW proteins
    # Wald proteins have no Cochran Q so keep them
    if "q_pval" in mr_pass.columns:
        mr_pass = mr_pass[
            ((mr_pass["mr_method"] == "IVW") & (mr_pass["q_pval"].fillna(-np.inf) >= q_pval))
            |
            (mr_pass["mr_method"] == "Wald ratio")
        ]


    # STAGE 2
    # pQTL - GWAS COLOC
    coloc_pass = coloc_outcome.copy()

    if "pp_h4_abf" not in coloc_pass.columns:
        st.error(
            "The COLOC results do not contain the required PP.H4 column. "
            f"Available columns: {list(coloc_pass.columns)}"
        )
        st.stop()

    coloc_pass["pp_h4_abf"] = pd.to_numeric(
        coloc_pass["pp_h4_abf"],
        errors="coerce"
    )

    coloc_pass = coloc_pass[
        coloc_pass["pp_h4_abf"].fillna(0) >= pp4
    ].copy()

    # STAGE 2b
    # PWCoCo (conditional coloc via GCTA-COJO-style stepwise selection) - runs
    # alongside standard pairwise COLOC above, not instead of it (see
    # project_pwcoco_wiring memory). A protein that colocalises under EITHER
    # method is carried forward as a prioritised target; coloc_support below
    # records which method(s) actually supported it rather than dropping a
    # discordant hit. PWCoCo's own .coloc output uses the same H0-H4 coloc.abf
    # posterior-probability scheme as standard COLOC, so it reuses the same
    # PP.H4 threshold (pp4) instead of a separate PWCoCo-specific control.
    # PWCoCo isn't split by outcome_trait (no such column), same as SMR/HyPrColoc.
    pwcoco_outcome = pwcoco.copy()
    pwcoco_pass = pwcoco_outcome.copy()

    if not pwcoco_pass.empty and "h4" in pwcoco_pass.columns:
        pwcoco_pass["h4"] = pd.to_numeric(pwcoco_pass["h4"], errors="coerce")
        pwcoco_pass = pwcoco_pass[pwcoco_pass["h4"].fillna(0) >= pp4].copy()
    else:
        pwcoco_pass = pwcoco_pass.iloc[0:0]

    coloc_pass_proteins = (
        set(coloc_pass["protein"].dropna().astype(str)) if "protein" in coloc_pass.columns else set()
    )
    pwcoco_pass_proteins = (
        set(pwcoco_pass["protein"].dropna().astype(str)) if "protein" in pwcoco_pass.columns else set()
    )
    coloc_support_pass_proteins = coloc_pass_proteins | pwcoco_pass_proteins

    def _coloc_support_label(protein_id):
        in_coloc = protein_id in coloc_pass_proteins
        in_pwcoco = protein_id in pwcoco_pass_proteins
        if in_coloc and in_pwcoco:
            return "both"
        return "coloc_only" if in_coloc else "pwcoco_only"

    # proteins which pass cis-MR AND (standard COLOC OR PWCoCo) - a left merge
    # (not inner) so a pwcoco_only protein keeps its mr_pass columns even though
    # it has no coloc_pass row to join against (its pp_h4_abf etc. come through
    # as NaN, which is the correct signal that standard COLOC didn't support it)
    mr_coloc_pass = mr_pass[
        mr_pass["protein"].astype(str).isin(coloc_support_pass_proteins)
    ].merge(
        coloc_pass,
        on="protein",
        how="left",
        suffixes=("_mr", "_pqtl_coloc")
    )

    mr_coloc_pass["coloc_support"] = mr_coloc_pass["protein"].astype(str).map(_coloc_support_label)

    # preserve assay-specific protein IDs
    # only remove fully duplicated merged rows
    mr_coloc_pass = mr_coloc_pass.drop_duplicates()

    # single source of truth for "which method(s) supported this protein" - computed
    # once here (pre-protein-search-filter, so it's stable regardless of the sidebar
    # search box) and reused by both the Final Targets Sankey and the Target Profile tab
    coloc_support_status = (
        mr_coloc_pass.drop_duplicates("protein").set_index("protein")["coloc_support"].to_dict()
        if "protein" in mr_coloc_pass.columns and "coloc_support" in mr_coloc_pass.columns
        else {}
    )

    # add harmonised top cis-hit information
    if not target_info.empty and "protein" in target_info.columns:
        target_cols = [
            "protein",
            "snp",
            "a1",
            "a2",
            "frq",
            "gwas_beta",
            "gwas_se",
            "gwas_p",
            "pqtl_beta",
            "pqtl_se",
            "pqtl_p"
        ]

        target_cols = available_cols(target_info, target_cols)

        mr_coloc_pass = mr_coloc_pass.merge(
            target_info[target_cols],
            on="protein",
            how="left"
        )


    # SMR (bulk/sc eQTL) and HyPrColoc (bulk/sc eQTL) hits are not split by
    # outcome_trait like mr/coloc are, so they're carried through as their own
    # tables rather than an "_outcome" subset
    smr_display = smr.copy()
    hyprcoloc_display = hyprcoloc.copy()
    pwcoco_display = pwcoco_outcome.copy()
    pwcoco_eqtl_pqtl_display = pwcoco_eqtl_pqtl.copy()
    pwcoco_eqtl_gwas_display = pwcoco_eqtl_gwas.copy()

    # computed live against the sidebar's pp4 slider (see compute_shared_snp_table) -
    # NOT read from the pipeline's pre-computed pwcoco_shared_snps.tsv, so this
    # always agrees with coloc_pass/pwcoco_pass above regardless of where pp4 is
    # set, rather than silently disagreeing with them at any threshold other than
    # whatever the pipeline used when it last wrote that file. Reused by the Final
    # Targets Sankey (three_trait_* below) and the PWCoCo-QTL tab's own table.
    shared_snp_table = compute_shared_snp_table(pwcoco_outcome, pwcoco_eqtl_pqtl, pwcoco_eqtl_gwas, pp4)
    triangulated_proteins = (
        set(shared_snp_table["protein"].dropna().astype(str))
        if "protein" in shared_snp_table.columns
        else set()
    )

    if "posterior_prob" in hyprcoloc_display.columns:
        hyprcoloc_display["posterior_prob"] = pd.to_numeric(hyprcoloc_display["posterior_prob"], errors="coerce")

    # protein search
    if protein:
        mr_outcome = filter_protein(mr_outcome, protein)
        mr_pass = filter_protein(mr_pass, protein)
        coloc_pass = filter_protein(coloc_pass, protein)
        mr_coloc_pass = filter_protein(mr_coloc_pass, protein)
        finngen_phewas_outcome = filter_protein(finngen_phewas_outcome, protein)
        ukb_phewas_outcome = filter_protein(ukb_phewas_outcome, protein)
        smr_display = filter_protein(smr_display, protein)
        hyprcoloc_display = filter_protein(hyprcoloc_display, protein)
        pwcoco_display = filter_protein(pwcoco_display, protein)

    # round coloc posterior probs
    for col in coloc_numeric_cols:
        if col in coloc_pass.columns:
            coloc_pass[col] = coloc_pass[col].round(3)

        if col in mr_coloc_pass.columns:
            mr_coloc_pass[col] = mr_coloc_pass[col].round(3)

    # SINGLE SOURCE OF TRUTH for every per-protein stage set used anywhere in this
    # dashboard (Overview's headline metrics + funnel chart, the Final Targets
    # Sankey, and the Final target list) - computed once, here, before any tab
    # renders, so a protein's stage-N pass/fail status can never drift between
    # 2 different displays that happened to compute it separately. If you need
    # a new "how many targets reached stage X" number anywhere in the app, read
    # it from one of these sets/counts rather than re-deriving it locally.
    mr_pass_proteins = set(mr_pass["protein"].dropna().astype(str)) if "protein" in mr_pass.columns else set()

    # coloc_support_status (computed in STAGE 2b above) already records, per
    # protein, which of standard COLOC / PWCoCo (or both) supported it
    coloc_support_pass_set = mr_pass_proteins & set(coloc_support_status.keys())
    coloc_support_both_set = {p for p in coloc_support_pass_set if coloc_support_status.get(p) == "both"}
    coloc_support_coloc_only_set = {p for p in coloc_support_pass_set if coloc_support_status.get(p) == "coloc_only"}
    coloc_support_pwcoco_only_set = {p for p in coloc_support_pass_set if coloc_support_status.get(p) == "pwcoco_only"}

    # everything downstream flows from the combined (COLOC OR PWCoCo) pass set
    coloc_pass_set = coloc_support_pass_set

    finngen_status = compute_phewas_classification_status(finngen_phewas_outcome, mr_outcome, coloc_pass_set)
    finngen_pass_set = {protein for protein in coloc_pass_set if finngen_status.get(protein) != "adverse_effect"}

    ukb_status = compute_phewas_classification_status(ukb_phewas_outcome, mr_outcome, finngen_pass_set)
    ukb_pass_set = {protein for protein in finngen_pass_set if ukb_status.get(protein) != "adverse_effect"}

    smr_pass_rows = smr_display.copy()

    if not smr_pass_rows.empty and {"q_smr", "p_heidi"}.issubset(smr_pass_rows.columns):
        smr_pass_rows = smr_pass_rows[
            smr_pass_rows["q_smr"].notna() & (smr_pass_rows["q_smr"] < smr_fdr_threshold) &
            smr_pass_rows["p_heidi"].notna() & (smr_pass_rows["p_heidi"] > heidi_p_threshold)
        ]
    else:
        smr_pass_rows = smr_pass_rows.iloc[0:0]

    if "data_type" in smr_pass_rows.columns and "protein" in smr_pass_rows.columns:
        bulk_pass_set = set(smr_pass_rows.loc[smr_pass_rows["data_type"] == "bulk", "protein"].dropna().astype(str))
        sc_pass_set = set(smr_pass_rows.loc[smr_pass_rows["data_type"] != "bulk", "protein"].dropna().astype(str))
    else:
        bulk_pass_set = set()
        sc_pass_set = set()

    coloc_fail_set = mr_pass_proteins - coloc_pass_set
    finngen_fail_set = coloc_pass_set - finngen_pass_set
    ukb_fail_set = finngen_pass_set - ukb_pass_set
    both_set = ukb_pass_set & bulk_pass_set & sc_pass_set
    bulk_only_set = (ukb_pass_set & bulk_pass_set) - sc_pass_set
    sc_only_set = (ukb_pass_set & sc_pass_set) - bulk_pass_set
    neither_set = ukb_pass_set - bulk_pass_set - sc_pass_set

    # HyPrColoc runs for any target with SMR support, whether that came from bulk
    # eQTL, single-cell eQTL, or both - no-SMR-support targets are terminal at the
    # SMR stage. HyPrColoc (bin/hyprcoloc_targets.py) runs on the FULL coloc_support
    # union (standard-COLOC-only, PWCoCo-only, or both - see project_pwcoco_wiring
    # memory), so there's no separate "not run" carve-out here - a missing
    # HyPrColoc row (e.g. <2 shared SNPs) falls into the ordinary "failed" bucket.
    smr_eligible_set = both_set | sc_only_set | bulk_only_set
    hyprcoloc_testable_set = smr_eligible_set
    hyprcoloc_status = compute_hyprcoloc_pass_status(hyprcoloc_display, hyprcoloc_testable_set, hyprcoloc_pp_threshold)
    hyprcoloc_pass_set = {protein for protein in hyprcoloc_testable_set if hyprcoloc_status.get(protein)}
    hyprcoloc_fail_set = hyprcoloc_testable_set - hyprcoloc_pass_set

    # PWCoCo-QTL triangulation - co-equal alternative to HyPrColoc above, not a
    # downstream refinement of it (same relationship PWCoCo(pQTL-GWAS) has to
    # standard pairwise COLOC(pQTL-GWAS) - see coloc_support above). Both ask
    # whether pQTL, eQTL and GWAS share 1 causal variant: HyPrColoc via a single
    # shared-cluster assumption, PWCoCo-QTL via conditioning (pQTL-GWAS, eQTL-pQTL
    # and eQTL-GWAS PWCoCo all converging on the SAME SNP - see triangulated_proteins
    # above, computed live via compute_shared_snp_table()). A target reaches Final
    # Targets if EITHER method supports it.
    pwcoco_qtl_pass_set = hyprcoloc_testable_set & triangulated_proteins
    three_trait_both_set = hyprcoloc_pass_set & pwcoco_qtl_pass_set
    three_trait_hyprcoloc_only_set = hyprcoloc_pass_set - pwcoco_qtl_pass_set
    three_trait_pwcoco_qtl_only_set = pwcoco_qtl_pass_set - hyprcoloc_pass_set
    three_trait_fail_set = hyprcoloc_testable_set - (hyprcoloc_pass_set | pwcoco_qtl_pass_set)

    # "Prioritised" (Overview headline metric) = passed cis-MR + COLOC-or-PWCoCo +
    # safety-cleared on FinnGen/UKB PheWAS - i.e. ukb_pass_set above, same set the
    # Sankey's "UKB passed" node uses. Kept as its own dataframe (not just a count)
    # since the Overview tab's prioritised-target cards below need the full rows.
    mr_coloc_safe_pass = (
        mr_coloc_pass[mr_coloc_pass["protein"].astype(str).isin(ukb_pass_set)].copy()
        if "protein" in mr_coloc_pass.columns else mr_coloc_pass.copy()
    )

    # main staged target counts - every number here reads directly off the shared
    # sets computed above (mr_pass_proteins, coloc_pass_set, ukb_pass_set,
    # smr_eligible_set, hyprcoloc_pass_set), so this bar, the Overview funnel chart,
    # the Final Targets Sankey and the Final target list can never disagree - all 4
    # displays are reading the same underlying sets, not 4 separate computations of
    # "how many targets are at stage N".
    n_tested = safe_nunique(mr_outcome, "protein")
    n_mr = safe_nunique(mr_pass, "protein")
    n_mr_coloc = safe_nunique(mr_coloc_pass, "protein")
    n_mr_coloc_safe = safe_nunique(mr_coloc_safe_pass, "protein")
    n_finngen_phewas = safe_nunique(finngen_phewas_outcome, "protein")
    n_ukb_phewas = safe_nunique(ukb_phewas_outcome, "protein")
    n_smr_eligible = len(smr_eligible_set)
    # "Multi-omics" headline = EITHER HyPrColoc OR PWCoCo-QTL triangulation, not
    # HyPrColoc alone - see three_trait_* above (same union philosophy as
    # coloc_support). Matches the Final Targets Sankey's own three_trait_* nodes.
    n_multi_omics_pass = len(hyprcoloc_pass_set | pwcoco_qtl_pass_set)

    # persistent context bar - stays visible above every tab, so "how many targets
    # survived so far, and what does this dashboard even show" never requires
    # hunting through tabs to answer. Runs the complete pipeline through to
    # Multi-omics (HyPrColoc) - the same final number as the Final Targets tab's
    # "Multi-omics" unique-target count - so nowhere in the app stops the funnel
    # early and leaves a reader wondering why a later tab shows a smaller number.
    with st.container(border=True):
        funnel_col1, funnel_col2, funnel_col3, funnel_col4, funnel_col5, funnel_col6 = st.columns(6)
        # delta text kept short (no "of <previous stage name>") - the stage name is
        # already the label right above it, and a longer string just gets clipped
        # by Streamlit's fixed-width delta pill at 6-across
        funnel_col1.metric("Tested", n_tested)
        funnel_col2.metric(
            "Passed cis-MR", n_mr,
            f"{retention(n_mr, n_tested):.0f}% retained"
        )
        funnel_col3.metric(
            "Colocalised", n_mr_coloc,
            f"{retention(n_mr_coloc, n_mr):.0f}% retained"
        )
        funnel_col4.metric(
            "Safety-cleared", n_mr_coloc_safe,
            f"{retention(n_mr_coloc_safe, n_mr_coloc):.0f}% retained"
        )
        funnel_col5.metric(
            "SMR-supported", n_smr_eligible,
            f"{retention(n_smr_eligible, n_mr_coloc_safe):.0f}% retained"
        )
        funnel_col6.metric(
            "Multi-omics", n_multi_omics_pass,
            f"{retention(n_multi_omics_pass, n_smr_eligible):.0f}% retained"
        )
        st.caption(
            f"{dataset_name} pQTLs → **{outcome}**. Colocalised = passed standard COLOC or "
            "PWCoCo. Safety-cleared = also no adverse FinnGen/UKB PheWAS hit ('Prioritised' "
            "on the Overview tab's target cards below). SMR-supported = also cleared SMR "
            "FDR/HEIDI. Multi-omics = pQTL+GWAS+eQTL share 1 causal variant, via **either** "
            "HyPrColoc's clustering **or** PWCoCo-QTL's SNP-level triangulation - this is the "
            "exact same count as the **Multi-omics** view on the **7. Final Targets** tab. "
            "Full branching detail (COLOC-vs-PWCoCo, bulk-vs-single-cell, HyPrColoc-vs-"
            "PWCoCo-QTL) is on that tab's Sankey; colour legend is in the sidebar."
        )

    tab1, tab_profile, tab_evidence, tab8, tab9, tab10 = st.tabs([
        "Overview",
        "Target Profile",
        "Evidence by Stage",
        "7. Final Targets",
        "PWCoCo (conditional coloc)",
        "PWCoCo-QTL (eQTL triangulation)"
    ])

    with tab1:
        st.caption("OVERVIEW · START HERE")
        st.subheader("How this dashboard is organised")
        st.caption(
            "Every protein below moves through the same 7-stage pipeline, in order. A "
            "target only reaches a later stage once it has already passed everything "
            "before it. Each stage's full results table is under **Evidence by Stage**; "
            "for 1 target's complete story instead, use **Target Profile**."
        )

        # 2 rows of 4 rather than 1 row of 7 - at a typical laptop viewport (not
        # maximised, or a split-screen window), 7 equal st.columns() go narrow
        # enough that Streamlit wraps every word onto its own line and the cards
        # become unreadable; a max-4-wide grid stays legible at realistic widths
        pipeline_rows = [PIPELINE_STAGES[i:i + 4] for i in range(0, len(PIPELINE_STAGES), 4)]
        position = 1

        for row_stages in pipeline_rows:
            # size the row to how many cards it actually has (not always 4) -
            # a trailing row of 3 previously still reserved a 4th empty slot,
            # leaving a visible gap that read as "something's missing"
            pipeline_columns = st.columns(len(row_stages))
            for pipeline_column, stage in zip(pipeline_columns, row_stages):
                with pipeline_column:
                    with st.container(border=True):
                        st.markdown(f"**{position}. {stage['title']}**")
                        st.caption(stage["blurb"])
                position += 1

        st.divider()
        st.subheader("Target prioritisation")
        st.caption(
            "The same 6 numbers as the bar above the tabs, shown here at scale - the "
            "complete flow through every stage, all the way to the same **Multi-omics** "
            "count on the **7. Final Targets** tab. That tab's Sankey diagram shows the "
            "same funnel with the extra branching detail (COLOC-vs-PWCoCo, bulk-vs-"
            "single-cell) this chart collapses for readability."
        )

        funnel_df = pd.DataFrame({
            "stage": [
                "Proteins tested by cis-MR",
                "cis-MR supported",
                "+ COLOC or PWCoCo",
                "+ no adverse PheWAS hit (FinnGen/UKB)",
                "+ SMR FDR/HEIDI support",
                "+ multi-omics confirmed (HyPrColoc or PWCoCo-QTL)"
            ],
            "n_targets": [
                n_tested,
                n_mr,
                n_mr_coloc,
                n_mr_coloc_safe,
                n_smr_eligible,
                n_multi_omics_pass
            ]
        })

        # log-scale x-axis - the drop-off from stage 1 to stage 2 is routinely
        # >10x, which on a linear axis squashes every later bar into a sliver
        # a few pixels wide (all visually identical regardless of whether it's
        # 39, 22 or 10) - log scale keeps each stage's own drop-off legible
        # while still showing the overall taper. Each bar's exact count is
        # printed on it regardless, so this doesn't cost any precision.
        funnel_fig = px.bar(
            funnel_df,
            x="n_targets",
            y="stage",
            orientation="h",
            text="n_targets",
            title="Progressive target prioritisation",
            labels={"n_targets": "Number of unique proteins (log scale)", "stage": ""},
            height=520,
            template="plotly_white",
            log_x=True
        )

        funnel_fig.update_yaxes(
            categoryorder="array",
            categoryarray=funnel_df["stage"][::-1]
        )
        funnel_fig.update_traces(textposition="outside", textfont=dict(size=13))
        # dtick=1 restricts the log axis to whole powers of 10 (1, 10, 100, 1000) -
        # Plotly's default log-axis tick spacing otherwise interleaves minor ticks
        # (2, 5) between them, which reads as a broken sequence ("10 2 5 100...")
        # to anyone not used to log axes
        funnel_fig.update_xaxes(dtick=1)
        funnel_fig.update_layout(showlegend=False, margin=dict(l=20, r=60, t=60, b=20))
        st.plotly_chart(funnel_fig, width="stretch")

        st.divider()
        st.subheader("Prioritised targets")
        st.caption(
            "Every protein that has passed cis-MR and colocalisation (standard COLOC or PWCoCo) at the "
            "thresholds set in the sidebar, and has no Bonferroni-significant FinnGen or "
            "UKB phenome-wide MR hit classified as a potential adverse effect (i.e. running "
            "opposite to the primary protein→AD effect direction), gets 1 card below. A same-"
            "direction hit is a potential additional indication and does not exclude a "
            "target. Betas and alleles are harmonised to the outcome GWAS risk allele."
        )

        if not mr_coloc_safe_pass.empty:
            st.success(
                f"{n_mr_coloc_safe} unique target(s) passed the selected cis-MR and pairwise COLOC "
                "thresholds with no Bonferroni-significant adverse-effect signal in FinnGen/UKB phenome-wide MR."
            )

            if n_mr_coloc_safe < n_mr_coloc:
                st.caption(
                    f"{n_mr_coloc - n_mr_coloc_safe} additional target(s) passed cis-MR + COLOC but were "
                    "excluded here for a Bonferroni-significant, opposite-direction (adverse-effect) FinnGen/UKB PheWAS hit - "
                    "see the Sankey diagram on the **7. Final Targets** tab for the full breakdown."
                )

            cards_df = mr_coloc_safe_pass.sort_values(
                ["pp_h4_abf", "mr_fdr_q"],
                ascending=[False, True],
                na_position="last"
            ) if "pp_h4_abf" in mr_coloc_safe_pass.columns else mr_coloc_safe_pass

            render_prioritised_target_cards(cards_df, outcome, pp4=pp4)

            st.divider()

            prioritised_cols = [
                "protein",
                "snp",
                "a1",
                "a2",
                "frq",
                "gwas_beta",
                "gwas_se",
                "gwas_p",
                "pqtl_beta",
                "pqtl_se",
                "pqtl_p",
                "mr_method",
                "n_instruments",
                "mr_beta",
                "mr_se",
                "mr_pval",
                "mr_fdr_q",
                "q_pval",
                "egger_intercept_pval",
                "pp_h4_abf",
                "coloc_support"
            ]

            prioritised_cols = available_cols(mr_coloc_safe_pass, prioritised_cols)

            if "pp_h4_abf" in mr_coloc_safe_pass.columns:
                mr_coloc_safe_pass = mr_coloc_safe_pass.sort_values(
                    ["pp_h4_abf", "mr_fdr_q"],
                    ascending=[False, True],
                    na_position="last"
                )

            overview_table = mr_coloc_safe_pass[prioritised_cols].copy()

            overview_column_names = {
                "protein": "Protein",
                "snp": "Top SNP",
                "a1": "Risk allele",
                "a2": "Other allele",
                "frq": "Risk allele frequency",
                "gwas_beta": "GWAS beta",
                "gwas_se": "GWAS SE",
                "gwas_p": "GWAS p-value",
                "pqtl_beta": "pQTL beta",
                "pqtl_se": "pQTL SE",
                "pqtl_p": "pQTL p-value",
                "mr_method": "MR method",
                "n_instruments": "N instruments",
                "mr_beta": "MR beta",
                "mr_se": "MR SE",
                "mr_pval": "MR p-value",
                "mr_fdr_q": "MR FDR",
                "q_pval": "Cochran Q p-value",
                "egger_intercept_pval": "Egger intercept p-value",
                "pp_h4_abf": "COLOC PP.H4",
                "coloc_support": "Coloc support"
            }

            overview_table = overview_table.rename(columns=overview_column_names)

            with st.expander("View full data table & download every column"):
                st.dataframe(
                    overview_table,
                    width="stretch",
                    hide_index=True,
                    column_config={
                        "Risk allele frequency": st.column_config.NumberColumn(format="%.3f"),
                        "GWAS beta": st.column_config.NumberColumn(format="%.4f"),
                        "GWAS SE": st.column_config.NumberColumn(format="%.4f"),
                        "GWAS p-value": st.column_config.NumberColumn(format="%.3e"),
                        "pQTL beta": st.column_config.NumberColumn(format="%.4f"),
                        "pQTL SE": st.column_config.NumberColumn(format="%.4f"),
                        "pQTL p-value": st.column_config.NumberColumn(format="%.3e"),
                        "MR beta": st.column_config.NumberColumn(format="%.4f"),
                        "MR SE": st.column_config.NumberColumn(format="%.4f"),
                        "MR p-value": st.column_config.NumberColumn(format="%.3e"),
                        "MR FDR": st.column_config.NumberColumn(format="%.3e"),
                        "Cochran Q p-value": st.column_config.NumberColumn(format="%.3e"),
                        "Egger intercept p-value": st.column_config.NumberColumn(format="%.3e"),
                        "COLOC PP.H4": st.column_config.NumberColumn(format="%.3f")
                    }
                )

                st.download_button(
                    label="Download prioritised targets",
                    data=overview_table.to_csv(index=False, sep="\t"),
                    file_name=f"{pqtl_dataset}_{outcome}_prioritised_target_overview.tsv",
                    mime="text/tab-separated-values",
                    key="download_prioritised_targets_overview",
                    width="stretch"
                )

        elif n_mr_coloc > 0:
            st.info(
                f"{n_mr_coloc} target(s) passed the selected cis-MR and pQTL COLOC thresholds, but all were "
                "excluded here for a Bonferroni-significant, opposite-direction (adverse-effect) FinnGen/UKB PheWAS hit - see the "
                "Sankey diagram on the **7. Final Targets** tab for the full breakdown."
            )
        else:
            st.info("No proteins currently pass both the selected cis-MR and pQTL COLOC thresholds.")

        if finngen_phewas_available or ukb_phewas_available:
            st.caption(
                f"FinnGen PheWAS results are available for {n_finngen_phewas} unique target(s); "
                f"UKB PheWAS results are available for {n_ukb_phewas} unique target(s)."
            )

    with tab_profile:
        st.caption("SEARCH ANY TARGET · FULL EVIDENCE TRAIL IN ONE PLACE")
        st.subheader("Target Profile")
        st.caption(
            "Every other tab is organised by pipeline stage (1 table per stage, every "
            "protein at once). This tab flips that around: pick 1 target and see its "
            "complete evidence trail - cis-MR, colocalisation, PheWAS safety, SMR and "
            "HyPrColoc - top to bottom, in the order the pipeline actually applies them."
        )

        profile_proteins = (
            sorted(mr_outcome["protein"].dropna().astype(str).unique().tolist())
            if "protein" in mr_outcome.columns else []
        )

        if not profile_proteins:
            st.info("No proteins are available to look up for this outcome / pQTL dataset.")
        else:
            default_index = 0
            if protein:
                matches = [p for p in profile_proteins if protein.lower() in p.lower()]
                if matches:
                    default_index = profile_proteins.index(matches[0])

            selected_target = st.selectbox(
                "Choose a target",
                profile_proteins,
                index=default_index,
                key="target_profile_selector"
            )

            mr_pass_proteins_for_profile = (
                set(mr_pass["protein"].dropna().astype(str)) if "protein" in mr_pass.columns else set()
            )

            render_target_profile(
                protein=selected_target,
                pqtl_dataset=pqtl_dataset,
                mr_outcome=mr_outcome,
                mr_pass_proteins=mr_pass_proteins_for_profile,
                coloc_outcome=coloc_outcome,
                pwcoco_outcome=pwcoco_outcome,
                coloc_support_status=coloc_support_status,
                pp4=pp4,
                finngen_phewas_outcome=finngen_phewas_outcome,
                ukb_phewas_outcome=ukb_phewas_outcome,
                smr_display=smr_display,
                hyprcoloc_display=hyprcoloc_display,
                smr_fdr_threshold=smr_fdr_threshold,
                heidi_p_threshold=heidi_p_threshold,
                hyprcoloc_pp_threshold=hyprcoloc_pp_threshold,
            )

    with tab_evidence:
        st.caption(
            "Detailed per-stage tables, for auditing exact numbers behind a call - most people "
            "want the **Overview** or **Target Profile** tab instead."
        )

        tab2, tab3, tab4, tab5, tab6, tab7 = st.tabs([
            "1. cis-MR",
            "2. pQTL–GWAS COLOC",
            "3. FinnGen PheWAS",
            "4. UKB PheWAS",
            "5. SMR (bulk/sc eQTL)",
            "6. HyPrColoc (bulk/sc eQTL)"
        ])

        with tab2:
            stage_caption(1)
            st.subheader("cis-MR results")
            st.caption(
                "Wald ratio is used for proteins with one instrument. IVW is used for proteins "
                "with more than one instrument."
            )

            show_all_mr = st.toggle("Show all tested cis-MR proteins", value=False)
            mr_display = mr_outcome if show_all_mr else mr_pass
            n_ivw = (mr_display["mr_method"] == "IVW").sum()
            n_wald = (mr_display["mr_method"] == "Wald ratio").sum()

            with st.container(border=True):
                col1, col2, col3 = st.columns(3)
                col1.metric("MR proteins shown", mr_display["protein"].nunique())
                col2.metric("IVW proteins", int(n_ivw))
                col3.metric("Wald proteins", int(n_wald))

            display_cols = [
                "protein",
                "outcome_trait",
                "n_instruments",
                "mr_method",
                "mr_beta",
                "mr_se",
                "mr_pval",
                "mr_fdr_q",
                "q_pval",
                "egger_intercept_pval",
                "ivw_beta",
                "ivw_se",
                "ivw_pval",
                "ivw_fdr_q",
                "wald_beta",
                "wald_se",
                "wald_pval",
                "wald_fdr_q"
            ]

            display_cols = available_cols(mr_display, display_cols)
            remaining_cols = [col for col in mr_display.columns if col not in display_cols]

            st.subheader("MR association table")

            mr_table_column_names = {
                "protein": "Protein",
                "outcome_trait": "Outcome",
                "n_instruments": "N instruments",
                "mr_method": "MR method",
                "mr_beta": "MR beta",
                "mr_se": "MR SE",
                "mr_pval": "MR p-value",
                "mr_fdr_q": "MR FDR",
                "q_pval": "Cochran Q p-value",
                "egger_intercept_pval": "Egger intercept p-value",
                "ivw_beta": "IVW beta",
                "ivw_se": "IVW SE",
                "ivw_pval": "IVW p-value",
                "ivw_fdr_q": "IVW FDR",
                "wald_beta": "Wald beta",
                "wald_se": "Wald SE",
                "wald_pval": "Wald p-value",
                "wald_fdr_q": "Wald FDR",
            }

            st.dataframe(
                mr_display[display_cols + remaining_cols].rename(columns=mr_table_column_names),
                width="stretch",
                hide_index=True
            )

            st.download_button(
                label="Download cis-MR results",
                data=mr_display[display_cols + remaining_cols].rename(columns=mr_table_column_names).to_csv(index=False, sep="\t"),
                file_name=f"{pqtl_dataset}_{outcome}_cis_MR.tsv",
                mime="text/tab-separated-values",
                key="download_cis_mr_results",
                width="stretch",
            )

            st.divider()

            # primary MR volcano plot
            plot_df = mr_display[
                mr_display["mr_pval"].notna() &
                mr_display["mr_beta"].notna() &
                (mr_display["mr_pval"] > 0)
            ].copy()

            if not plot_df.empty:
                st.subheader("MR effect landscape")

                plot_df["minus_log10_mr_pval"] = -np.log10(plot_df["mr_pval"])
                plot_df["significant"] = plot_df["mr_fdr_q"] < 0.05

                fig = px.scatter(
                    plot_df,
                    x="mr_beta",
                    y="minus_log10_mr_pval",
                    hover_name="protein",
                    color="significant",
                    color_discrete_map=SIGNIFICANCE_COLOR_MAP,
                    symbol="mr_method",
                    hover_data={
                        "mr_method": True,
                        "n_instruments": True,
                        "mr_beta": ":.4f",
                        "mr_se": ":.4f",
                        "mr_pval": ":.3e",
                        "mr_fdr_q": ":.3e",
                        "minus_log10_mr_pval": False
                    },
                    labels={
                        "mr_beta": "Primary MR beta",
                        "minus_log10_mr_pval": "-log10(primary MR p-value)",
                        "mr_method": "MR method",
                        "significant": "FDR < 0.05"
                    },
                    title="Primary cis-MR volcano plot",
                    height=600,
                    template="plotly_white"
                )

                fig.add_hline(y=-np.log10(0.05), line_dash="dash", line_color="grey")
                fig.add_vline(x=0, line_dash="dash", line_color="grey")
                st.plotly_chart(fig, width="stretch")

            else:
                st.info("No MR results remain after applying the selected filters.")

        with tab3:
            stage_caption(2)
            st.subheader("cis-MR + pQTL–GWAS COLOC targets")
            st.caption(
                "Targets shown here pass cis-MR and cleared the PP.H4 threshold via standard COLOC, "
                "PWCoCo, or both - check the **Coloc support** column. A `pwcoco_only` row's PP.H0-H4 "
                "columns are blank by design (standard COLOC genuinely didn't support it); see the "
                "**PWCoCo** tab for its conditional-analysis result instead."
            )

            if not mr_coloc_pass.empty:
                with st.container(border=True):
                    col1, col2, col3 = st.columns(3)
                    col1.metric("Colocalised proteins", mr_coloc_pass["protein"].nunique())
                    col2.metric("Median PP.H4", safe_median(mr_coloc_pass, "pp_h4_abf"))
                    col3.metric("Median MR FDR", safe_median(mr_coloc_pass, "mr_fdr_q", scientific=True))

                prioritised_cols = [
                    "protein",
                    "coloc_support",
                    "mr_method",
                    "n_instruments",
                    "mr_beta",
                    "mr_se",
                    "mr_pval",
                    "mr_fdr_q",
                    "q_pval",
                    "egger_intercept_pval",
                    "snp",
                    "pp_h0_abf",
                    "pp_h1_abf",
                    "pp_h2_abf",
                    "pp_h3_abf",
                    "pp_h4_abf"
                ]

                prioritised_cols = available_cols(mr_coloc_pass, prioritised_cols)
                remaining_cols = [col for col in mr_coloc_pass.columns if col not in prioritised_cols]

                if "pp_h4_abf" in mr_coloc_pass.columns:
                    mr_coloc_pass = mr_coloc_pass.sort_values(
                        ["pp_h4_abf", "mr_fdr_q"],
                        ascending=[False, True],
                        na_position="last"
                    )

                st.divider()
                st.subheader("Prioritised target table")

                prioritised_table_column_names = {
                    "protein": "Protein",
                    "coloc_support": "Coloc support",
                    "mr_method": "MR method",
                    "n_instruments": "N instruments",
                    "mr_beta": "MR beta",
                    "mr_se": "MR SE",
                    "mr_pval": "MR p-value",
                    "mr_fdr_q": "MR FDR",
                    "q_pval": "Cochran Q p-value",
                    "egger_intercept_pval": "Egger intercept p-value",
                    "snp": "Top SNP",
                    "pp_h0_abf": "PP.H0",
                    "pp_h1_abf": "PP.H1",
                    "pp_h2_abf": "PP.H2",
                    "pp_h3_abf": "PP.H3",
                    "pp_h4_abf": "PP.H4",
                }

                st.dataframe(
                    mr_coloc_pass[prioritised_cols + remaining_cols].rename(columns=prioritised_table_column_names),
                    width="stretch",
                    hide_index=True
                )

                st.download_button(
                    label="Download prioritised targets",
                    data=mr_coloc_pass.to_csv(index=False, sep="\t"),
                    file_name=f"{outcome}_prioritised_targets.tsv",
                    mime="text/tab-separated-values",
                    key="download_prioritised_targets_coloc",
                    width="stretch"
                )

            else:
                st.info("No proteins currently pass both the selected cis-MR and pQTL COLOC thresholds.")

        with tab4:
            render_phewas_section(
                phewas_outcome=finngen_phewas_outcome,
                mr_coloc_pass=mr_coloc_pass,
                mr_outcome=mr_outcome,
                source_name="FinnGen",
                source_description="FinnGen phenotype (ICD coded)",
                pqtl_dataset=pqtl_dataset,
                outcome=outcome,
                key_prefix="finngen",
                stage_number=3
            )

        with tab5:
            render_phewas_section(
                phewas_outcome=ukb_phewas_outcome,
                mr_coloc_pass=mr_coloc_pass,
                mr_outcome=mr_outcome,
                source_name="UKB TOPMed",
                source_description="UK Biobank EHR-derived PheCode phenotype",
                pqtl_dataset=pqtl_dataset,
                outcome=outcome,
                key_prefix="ukb",
                stage_number=4,
                is_fallback=True,
            )

        with tab6:
            stage_caption(5)
            st.subheader("SMR (bulk / single-cell eQTL)")
            st.caption(
                "Targets shown here already passed cis-MR + pQTL–GWAS COLOC, and additionally "
                "passed SMR (FDR-corrected) + HEIDI in the configured bulk and/or single-cell eQTL "
                "dataset(s). Alleles are aligned to the outcome risk allele, same convention as the "
                "cis-MR/COLOC top-hit table."
            )

            if smr_display.empty:
                st.info("No SMR results are available for this pQTL dataset yet.")
            else:
                # the caption above claims rows here already passed SMR FDR + HEIDI
                # at the sidebar's configured thresholds - that wasn't actually
                # being enforced (the sidebar sliders were a no-op on this tab,
                # only the data-type/cell-type pickers below did any filtering),
                # so apply the same q_smr/p_heidi gate used everywhere else
                # (render_target_profile, the Final Targets Sankey) before anything
                # else narrows the view further
                n_smr_before_threshold = safe_nunique(smr_display, "protein")
                if {"q_smr", "p_heidi"}.issubset(smr_display.columns):
                    smr_filtered = smr_display[
                        smr_display["q_smr"].notna() & (smr_display["q_smr"] < smr_fdr_threshold) &
                        smr_display["p_heidi"].notna() & (smr_display["p_heidi"] > heidi_p_threshold)
                    ].copy()
                else:
                    smr_filtered = smr_display.copy()

                if smr_filtered.empty:
                    st.info(
                        f"No targets currently pass SMR FDR < {smr_fdr_threshold:.2f} and "
                        f"HEIDI p > {heidi_p_threshold:.2f} for this pQTL dataset (out of "
                        f"{n_smr_before_threshold} target(s) with any SMR result). Try relaxing "
                        "the thresholds in the sidebar."
                    )

                if "data_type" in smr_filtered.columns:
                    available_data_types = ["All"] + sorted(smr_filtered["data_type"].dropna().unique().tolist())

                    data_type_choice = st.segmented_control(
                        "eQTL data type",
                        available_data_types,
                        default="All",
                        selection_mode="single",
                        key="smr_data_type_selector"
                    )

                    if data_type_choice and data_type_choice != "All":
                        smr_filtered = smr_filtered[smr_filtered["data_type"] == data_type_choice]

                if "cell_type" in smr_filtered.columns:
                    available_cell_types = sorted(smr_filtered["cell_type"].dropna().unique().tolist())
                    selected_cell_types = st.multiselect(
                        "Cell type ('bulk' rows have no single cell-type)",
                        available_cell_types,
                        default=available_cell_types
                    )

                    if selected_cell_types:
                        smr_filtered = smr_filtered[smr_filtered["cell_type"].isin(selected_cell_types)]

                with st.container(border=True):
                    metric1, metric2 = st.columns(2)
                    metric1.metric("Unique targets shown", safe_nunique(smr_filtered, "protein"))
                    metric2.metric("Target x cell/bulk rows", len(smr_filtered))

                st.divider()
                st.subheader("Target support landscape")
                st.caption(
                    "Each cell is a target that passed cis-MR + COLOC + SMR + HEIDI in that cell "
                    "type or bulk/tissue dataset. Color intensity is -log10(SMR FDR); blank cells "
                    "mean that target did not pass in that context."
                )

                heatmap_df = smr_filtered[
                    smr_filtered["q_smr"].notna() & (smr_filtered["q_smr"] > 0)
                ].copy() if "q_smr" in smr_filtered.columns else pd.DataFrame()

                if not heatmap_df.empty:
                    heatmap_df["neg_log10_q_smr"] = -np.log10(heatmap_df["q_smr"])
                    pivot = heatmap_df.pivot_table(
                        index="protein",
                        columns="cell_type",
                        values="neg_log10_q_smr",
                        aggfunc="max"
                    )
                    # most broadly-supported targets (fewest blank cells) first
                    pivot = pivot.loc[pivot.notna().sum(axis=1).sort_values(ascending=False).index]

                    heatmap_fig = px.imshow(
                        pivot,
                        color_continuous_scale=SEQUENTIAL_SCALE,
                        aspect="auto",
                        labels={"color": "-log10(SMR FDR)", "x": "Cell type / bulk dataset", "y": "Protein"},
                        title="Target x cell-type/tissue SMR support"
                    )
                    heatmap_fig.update_layout(height=max(400, 32 * len(pivot)))
                    st.plotly_chart(heatmap_fig, width="stretch")
                else:
                    st.info("No targets with a valid SMR FDR value to plot.")

                st.divider()

                smr_cols = [
                    "protein",
                    "data_type",
                    "cell_type",
                    "gene",
                    "a1",
                    "a2",
                    "b_gwas",
                    "b_eqtl",
                    "b_smr",
                    "p_smr",
                    "q_smr",
                    "p_heidi",
                    "eqtl_dataset"
                ]

                smr_cols = available_cols(smr_filtered, smr_cols)
                remaining_smr_cols = [col for col in smr_filtered.columns if col not in smr_cols]

                smr_table = smr_filtered[smr_cols + remaining_smr_cols].copy()

                smr_column_names = {
                    "protein": "Protein",
                    "data_type": "Data type",
                    "cell_type": "Cell type",
                    "gene": "Gene",
                    "a1": "Risk allele",
                    "a2": "Other allele",
                    "b_gwas": "GWAS beta (risk allele)",
                    "b_eqtl": "eQTL beta (risk allele)",
                    "b_smr": "SMR beta",
                    "p_smr": "SMR p-value",
                    "q_smr": "SMR FDR",
                    "p_heidi": "HEIDI p-value",
                    "eqtl_dataset": "eQTL dataset"
                }

                smr_table = smr_table.rename(columns=smr_column_names)

                st.dataframe(
                    smr_table,
                    width="stretch",
                    hide_index=True,
                    column_config={
                        "GWAS beta (risk allele)": st.column_config.NumberColumn(format="%.4f"),
                        "eQTL beta (risk allele)": st.column_config.NumberColumn(format="%.4f"),
                        "SMR beta": st.column_config.NumberColumn(format="%.4f"),
                        "SMR p-value": st.column_config.NumberColumn(format="%.3e"),
                        "SMR FDR": st.column_config.NumberColumn(format="%.3e"),
                        "HEIDI p-value": st.column_config.NumberColumn(format="%.3e")
                    }
                )

                st.download_button(
                    label="Download SMR hits",
                    data=smr_table.to_csv(index=False, sep="\t"),
                    file_name=f"{pqtl_dataset}_{phenotype}_SMR_hits.tsv",
                    mime="text/tab-separated-values",
                    key="download_smr_hits",
                    width="stretch"
                )

        with tab7:
            stage_caption(6)
            st.subheader("HyPrColoc (bulk/sc eQTL)")
            st.caption(
                "For every target x cell-type/tissue hit that already passed cis-MR + pQTL–GWAS "
                "COLOC + SMR + HEIDI in a bulk or single-cell eQTL dataset, HyPrColoc jointly tests "
                "the pQTL, GWAS and eQTL association signals in that target's cis-region for a "
                "single shared causal variant, restricted to the SNPs shared across all three and "
                "aligned onto a common effect allele. The table below only shows rows where "
                "HyPrColoc actually put all 3 traits into 1 credible set - that 3-way test is the "
                "entire point of running HyPrColoc, so a cluster missing the eQTL trait (it either "
                "clustered separately or joined no cluster at all) isn't shown as a result here, "
                "regardless of how confident the pQTL+GWAS-only cluster it did find is."
            )

            if hyprcoloc_display.empty:
                st.info("No HyPrColoc results are available for this pQTL dataset yet.")
            else:
                hyprcoloc_filtered = hyprcoloc_display.copy()

                if "data_type" in hyprcoloc_filtered.columns:
                    available_hyprcoloc_data_types = ["All"] + sorted(hyprcoloc_filtered["data_type"].dropna().unique().tolist())

                    hyprcoloc_data_type_choice = st.segmented_control(
                        "eQTL data type",
                        available_hyprcoloc_data_types,
                        default="All",
                        selection_mode="single",
                        key="hyprcoloc_data_type_selector"
                    )

                    if hyprcoloc_data_type_choice and hyprcoloc_data_type_choice != "All":
                        hyprcoloc_filtered = hyprcoloc_filtered[hyprcoloc_filtered["data_type"] == hyprcoloc_data_type_choice]

                if "cell_type" in hyprcoloc_filtered.columns:
                    available_hyprcoloc_cell_types = sorted(hyprcoloc_filtered["cell_type"].dropna().unique().tolist())
                    selected_hyprcoloc_cell_types = st.multiselect(
                        "Cell type / tissue",
                        available_hyprcoloc_cell_types,
                        default=available_hyprcoloc_cell_types,
                        key="hyprcoloc_cell_type_selector"
                    )

                    if selected_hyprcoloc_cell_types:
                        hyprcoloc_filtered = hyprcoloc_filtered[hyprcoloc_filtered["cell_type"].isin(selected_hyprcoloc_cell_types)]

                if "traits" in hyprcoloc_filtered.columns:
                    traits_lower = hyprcoloc_filtered["traits"].astype(str).str.lower()
                    all_3_traits_clustered = (
                        traits_lower.str.contains("pqtl_") &
                        traits_lower.str.contains("gwas_") &
                        traits_lower.str.contains("eqtl_")
                    )
                else:
                    all_3_traits_clustered = pd.Series(False, index=hyprcoloc_filtered.index)

                # HyPrColoc's whole purpose here is testing whether pQTL + GWAS + eQTL
                # share 1 causal variant - a credible set HyPrColoc couldn't put the eQTL
                # trait into (it clustered separately, or didn't join any cluster at all)
                # isn't a 3-way colocalisation result, so it doesn't belong in this table
                # at all, however high its own (pQTL+GWAS-only) posterior_prob is. Drop
                # those rows here rather than just marking them unticked, so every metric
                # and row below is scoped to genuine 3-trait attempts.
                n_before_3trait_filter = len(hyprcoloc_filtered)
                hyprcoloc_filtered = hyprcoloc_filtered[all_3_traits_clustered].copy()
                n_dropped_2trait = n_before_3trait_filter - len(hyprcoloc_filtered)
                if n_dropped_2trait > 0:
                    st.caption(
                        f"{n_dropped_2trait} target x cell-type/tissue row(s) excluded below - "
                        "HyPrColoc could not put the eQTL trait into the same credible set as "
                        "pQTL + GWAS for those (it clustered separately or didn't join any "
                        "cluster), so they aren't a 3-way colocalisation result."
                    )

                # every remaining row is already a 3-trait cluster (2-trait-only rows were
                # dropped above), so the tick/cross column here only needs to check the
                # posterior-probability threshold - matching compute_hyprcoloc_pass_status
                # and the "Passing threshold" metric above.
                if "posterior_prob" in hyprcoloc_filtered.columns:
                    posterior_prob_numeric = pd.to_numeric(hyprcoloc_filtered["posterior_prob"], errors="coerce")
                    hyprcoloc_filtered["hyprcoloc_passed"] = posterior_prob_numeric.fillna(0) >= hyprcoloc_pp_threshold
                else:
                    hyprcoloc_filtered["hyprcoloc_passed"] = False

                hyprcoloc_pass_rows = hyprcoloc_filtered[hyprcoloc_filtered["hyprcoloc_passed"]]

                with st.container(border=True):
                    metric1, metric2, metric3, metric4 = st.columns(4)
                    metric1.metric("Unique targets tested", safe_nunique(hyprcoloc_filtered, "protein"))
                    metric2.metric("Target x cell-type rows", len(hyprcoloc_filtered))
                    metric3.metric("Passing threshold", len(hyprcoloc_pass_rows))
                    metric4.metric("Median posterior probability", safe_median(hyprcoloc_filtered, "posterior_prob"))

                st.divider()
                st.subheader("HyPrColoc results")
                st.caption(
                    f"A row **passes** when its cluster contains the pQTL, GWAS and eQTL trait "
                    f"together with posterior probability ≥ {hyprcoloc_pp_threshold:.2f}. Some "
                    "targets have more than 1 row when HyPrColoc could not put every trait into a "
                    "single cluster (e.g. the eQTL signal clusters separately from pQTL + GWAS)."
                )

                hyprcoloc_cols = [
                    "protein",
                    "data_type",
                    "cell_type",
                    "traits",
                    "hyprcoloc_passed",
                    "posterior_prob",
                    "regional_prob",
                    "candidate_snp",
                    "posterior_explained_by_snp",
                    "dropped_trait",
                    "n_snps"
                ]

                hyprcoloc_cols = available_cols(hyprcoloc_filtered, hyprcoloc_cols)
                remaining_hyprcoloc_cols = [col for col in hyprcoloc_filtered.columns if col not in hyprcoloc_cols]

                hyprcoloc_table = hyprcoloc_filtered[hyprcoloc_cols + remaining_hyprcoloc_cols].copy()

                if "posterior_prob" in hyprcoloc_table.columns:
                    hyprcoloc_table = hyprcoloc_table.sort_values(
                        "posterior_prob",
                        ascending=False,
                        na_position="last"
                    )

                hyprcoloc_column_names = {
                    "protein": "Target",
                    "data_type": "Data type",
                    "cell_type": "Cell type",
                    "traits": "Clustered traits",
                    "hyprcoloc_passed": "HyPrColoc passed",
                    "posterior_prob": "Posterior probability",
                    "regional_prob": "Regional probability",
                    "candidate_snp": "Candidate SNP",
                    "posterior_explained_by_snp": "Posterior explained by SNP",
                    "dropped_trait": "Dropped trait",
                    "n_snps": "N SNPs"
                }

                hyprcoloc_table = hyprcoloc_table.rename(columns=hyprcoloc_column_names)

                st.dataframe(
                    hyprcoloc_table,
                    width="stretch",
                    hide_index=True,
                    column_config={
                        "Posterior probability": st.column_config.NumberColumn(format="%.4f"),
                        "Regional probability": st.column_config.NumberColumn(format="%.4f"),
                        "Posterior explained by SNP": st.column_config.NumberColumn(format="%.4f")
                    }
                )

                st.download_button(
                    label="Download HyPrColoc results",
                    data=hyprcoloc_table.to_csv(index=False, sep="\t"),
                    file_name=f"{pqtl_dataset}_{phenotype}_HyPrColoc.tsv",
                    mime="text/tab-separated-values",
                    key="download_hyprcoloc_hits",
                    width="stretch"
                )

    with tab8:
        stage_caption(7)
        st.subheader("Final Targets")
        st.caption(
            "The complete set of targets which passed cis-MR + pQTL–GWAS COLOC + SMR + HEIDI, "
            "broken down by the cell type or bulk/tissue dataset each was supported in. For "
            "single-cell rows, the eQTL beta is sourced from the original per-cell-type eQTL "
            "file rather than the raw SMR output (SMR's own allele coding doesn't always match "
            "the original file's effect allele)."
        )

        st.divider()
        st.subheader("Target prioritisation flow")
        st.caption(
            "How the proteins that passed cis-MR for this pQTL dataset narrow down to the "
            "final targets below. A target only reaches a downstream check once it has "
            "already passed everything upstream, so each stage is drawn from the survivors "
            "of the one before it. Hover any stage or ribbon to see exactly which genes it holds."
        )

        with st.expander("How each stage is decided"):
            st.markdown(
                "- **cis-MR**: this flow starts from the proteins that already passed; the "
                "much larger screening drop-off across every tested protein is the funnel in "
                "the Overview tab.\n"
                "- **pQTL–GWAS COLOC**: passes on the posterior-probability threshold set in "
                "the sidebar. PWCoCo (a conditional-analysis variant of COLOC, see the "
                "**PWCoCo** tab) runs alongside it on the same targets - passing *either* "
                "method is enough to continue, split into \"Both methods\" / \"COLOC only\" / "
                "\"PWCoCo only\" lanes below so discordant hits stay visible rather than being "
                "silently dropped.\n"
                "- **FinnGen / UKB phenome-wide MR**: UKB only runs as a fallback for targets "
                "with zero retained instruments in FinnGen. A target fails *only* when it has a "
                "Bonferroni-significant hit classified as an **adverse effect** (running opposite "
                "to the primary protein→AD effect direction). A same-direction hit is an "
                "**additional indication** and does not gate; no significant hit, or no PheWAS "
                "coverage at all, also count as passing.\n"
                f"- **SMR support**: requires SMR FDR (`q_SMR`) < {smr_fdr_threshold:.2f} and "
                f"HEIDI p-value > {heidi_p_threshold:.2f}, split by whether that support came "
                "from bulk/tissue eQTL data, single-cell eQTL data, or both.\n"
                "- **HyPrColoc**: runs against whichever eQTL dataset(s) supported the target's "
                "SMR stage (bulk, single-cell, or both); no-SMR-support targets end at the SMR "
                "stage. Runs on targets supported by standard COLOC, PWCoCo, or both (same "
                "**Coloc support** union as the COLOC/PWCoCo stage) - passes when a HyPrColoc "
                "cluster contains the pQTL, GWAS *and* eQTL trait together (not just 2 of the 3) "
                f"with posterior probability ≥ {hyprcoloc_pp_threshold:.2f}."
            )

        # every set used below (mr_pass_proteins, coloc_pass_set, finngen_pass_set,
        # ukb_pass_set, smr_eligible_set, hyprcoloc_pass_set, etc.) is computed once,
        # shared with the Overview funnel and the Final target list further down -
        # see the "SINGLE SOURCE OF TRUTH" block earlier in this function. Nothing
        # is recomputed here, so this diagram can never disagree with those.
        if not mr_pass_proteins:
            st.info("No proteins currently pass the cis-MR threshold, so no flow can be drawn.")
        else:
            # single source of truth for the diagram, its hover text and the
            # drill-down selector below - every count is derived from these sets.
            # `label` is what the chart prints (kept short so five columns fit
            # without wrapping), `name` is the formal stage name used on hover
            # and in the selector.
            sankey_groups = [
                dict(key="mr", column=0, name="cis-MR passed", label="cis-MR passed",
                     proteins=mr_pass_proteins, color=STATUS_MUTED, dropout=False),
                dict(key="coloc_support_both", column=1, name="COLOC + PWCoCo agree", label="Both methods",
                     proteins=coloc_support_both_set, color=COLOC_SUPPORT_BOTH_COLOR, dropout=False),
                dict(key="coloc_support_coloc_only", column=1, name="COLOC passed, PWCoCo did not", label="COLOC only",
                     proteins=coloc_support_coloc_only_set, color=COLOC_SUPPORT_COLOC_ONLY_COLOR, dropout=False),
                dict(key="coloc_support_pwcoco_only", column=1, name="PWCoCo passed, COLOC did not", label="PWCoCo only",
                     proteins=coloc_support_pwcoco_only_set, color=COLOC_SUPPORT_PWCOCO_ONLY_COLOR, dropout=False),
                dict(key="coloc_fail", column=1, name="Neither COLOC nor PWCoCo passed", label="Neither",
                     proteins=coloc_fail_set, color=STATUS_CRITICAL, dropout=True),
                dict(key="finngen_pass", column=2, name="FinnGen safety passed", label="FinnGen passed",
                     proteins=finngen_pass_set, color=STATUS_GOOD, dropout=False),
                dict(key="finngen_fail", column=2, name="FinnGen safety failed", label="FinnGen failed",
                     proteins=finngen_fail_set, color=STATUS_CRITICAL, dropout=True),
                dict(key="ukb_pass", column=3, name="UKB safety passed", label="UKB passed",
                     proteins=ukb_pass_set, color=STATUS_GOOD, dropout=False),
                dict(key="ukb_fail", column=3, name="UKB safety failed", label="UKB failed",
                     proteins=ukb_fail_set, color=STATUS_CRITICAL, dropout=True),
                dict(key="smr_both", column=4, name="SMR bulk + single-cell", label="Bulk + single-cell",
                     proteins=both_set, color=SANKEY_BOTH_COLOR, dropout=False),
                dict(key="smr_bulk", column=4, name="SMR bulk only", label="Bulk only",
                     proteins=bulk_only_set, color=SANKEY_BULK_COLOR, dropout=False),
                dict(key="smr_sc", column=4, name="SMR single-cell only", label="Single-cell only",
                     proteins=sc_only_set, color=SANKEY_SC_COLOR, dropout=False),
                dict(key="smr_none", column=4, name="SMR: no support", label="No SMR support",
                     proteins=neither_set, color=STATUS_MUTED, dropout=True),
                dict(key="three_trait_both", column=5, name="HyPrColoc AND PWCoCo-QTL triangulation agree", label="Both methods",
                     proteins=three_trait_both_set, color=COLOC_SUPPORT_BOTH_COLOR, dropout=False),
                dict(key="three_trait_hyprcoloc_only", column=5, name="HyPrColoc passed, PWCoCo-QTL did not triangulate", label="HyPrColoc only",
                     proteins=three_trait_hyprcoloc_only_set, color=COLOC_SUPPORT_COLOC_ONLY_COLOR, dropout=False),
                dict(key="three_trait_pwcoco_qtl_only", column=5, name="PWCoCo-QTL triangulated, HyPrColoc did not", label="PWCoCo-QTL only",
                     proteins=three_trait_pwcoco_qtl_only_set, color=COLOC_SUPPORT_PWCOCO_ONLY_COLOR, dropout=False),
                dict(key="three_trait_fail", column=5, name="Neither HyPrColoc nor PWCoCo-QTL supported", label="Neither",
                     proteins=three_trait_fail_set, color=STATUS_CRITICAL, dropout=True),
            ]

            sankey_edges = [
                ("mr", "coloc_support_both"), ("mr", "coloc_support_coloc_only"),
                ("mr", "coloc_support_pwcoco_only"), ("mr", "coloc_fail"),
                ("coloc_support_both", "finngen_pass"), ("coloc_support_both", "finngen_fail"),
                ("coloc_support_coloc_only", "finngen_pass"), ("coloc_support_coloc_only", "finngen_fail"),
                ("coloc_support_pwcoco_only", "finngen_pass"), ("coloc_support_pwcoco_only", "finngen_fail"),
                ("finngen_pass", "ukb_pass"), ("finngen_pass", "ukb_fail"),
                ("ukb_pass", "smr_both"), ("ukb_pass", "smr_bulk"),
                ("ukb_pass", "smr_sc"), ("ukb_pass", "smr_none"),
                ("smr_both", "three_trait_both"), ("smr_both", "three_trait_hyprcoloc_only"),
                ("smr_both", "three_trait_pwcoco_qtl_only"), ("smr_both", "three_trait_fail"),
                ("smr_sc", "three_trait_both"), ("smr_sc", "three_trait_hyprcoloc_only"),
                ("smr_sc", "three_trait_pwcoco_qtl_only"), ("smr_sc", "three_trait_fail"),
                ("smr_bulk", "three_trait_both"), ("smr_bulk", "three_trait_hyprcoloc_only"),
                ("smr_bulk", "three_trait_pwcoco_qtl_only"), ("smr_bulk", "three_trait_fail"),
            ]

            # an empty stage costs a label and a slot but carries no information,
            # so it is left out of the diagram (it stays in the selector below)
            drawn = [group for group in sankey_groups if group["proteins"]]
            node_index = {group["key"]: position for position, group in enumerate(drawn)}

            # Plotly keeps node labels inside the plot area, flipping the last
            # column's (the HyPrColoc/PWCoCo-QTL 3-trait support column's) to the
            # left of its nodes rather than letting them run into the margin - so
            # the right margin stays thin and the gap before the last column is the
            # widest, giving those flipped labels room to sit without overlapping
            # the SMR column's nodes. The early columns (mr -> coloc_support ->
            # finngen) carry the widest ribbons (the biggest drop-offs happen
            # here), which used to visually cover the trailing "(n)" on each
            # label - widened those gaps specifically rather than spacing every
            # column evenly.
            column_x = [0.02, 0.19, 0.37, 0.55, 0.68, 0.99]

            node_values = [len(group["proteins"]) for group in drawn]
            node_labels = [f"{group['label']} ({value})" for group, value in zip(drawn, node_values)]
            node_colors = [group["color"] for group in drawn]
            node_x = [column_x[group["column"]] for group in drawn]
            node_hover = [
                f"<b>{group['name']}: {value} target(s)</b><br><br>"
                f"{format_protein_list_html(group['proteins'])}"
                for group, value in zip(drawn, node_values)
            ]

            # `drawn` is already pass-lane-first within each column, so taking
            # each column in that order keeps the ribbons from crossing
            column_indices = [
                [node_index[group["key"]] for group in drawn if group["column"] == column]
                for column in range(len(column_x))
            ]
            node_y = layout_sankey_columns(column_indices, node_values, len(drawn))

            edges = [
                (node_index[source], node_index[target])
                for source, target in sankey_edges
                if source in node_index and target in node_index
            ]

            sankey_fig = go.Figure(go.Sankey(
                arrangement="fixed",
                textfont=dict(size=11, color="#2b2b33"),
                node=dict(
                    label=node_labels,
                    color=node_colors,
                    customdata=node_hover,
                    hovertemplate="%{customdata}<extra></extra>",
                    x=node_x,
                    y=node_y,
                    pad=18,
                    thickness=14,
                    line=dict(color="rgba(255,255,255,0.9)", width=0.8)
                ),
                link=dict(
                    source=[source for source, _ in edges],
                    target=[target for _, target in edges],
                    # a target can have more than 1 source (e.g. the 3 SMR support
                    # lanes all feed HyPrColoc), so each ribbon's value/hover must be
                    # the proteins actually shared between its own source and target -
                    # not the target's total, which would double- or triple-count them
                    value=[
                        len(drawn[source]["proteins"] & drawn[target]["proteins"])
                        for source, target in edges
                    ],
                    # ribbons take the colour of where they land, and drop-out
                    # ribbons sit fainter so the eye follows the surviving lane
                    color=[
                        hex_to_rgba(node_colors[target], 0.25 if drawn[target]["dropout"] else 0.4)
                        for _, target in edges
                    ],
                    customdata=[
                        f"<b>{drawn[source]['name']} → {drawn[target]['name']}</b><br>"
                        f"{len(drawn[source]['proteins'] & drawn[target]['proteins'])} target(s)<br><br>"
                        f"{format_protein_list_html(drawn[source]['proteins'] & drawn[target]['proteins'])}"
                        for source, target in edges
                    ],
                    hovertemplate="%{customdata}<extra></extra>"
                )
            ))

            sankey_fig.update_layout(
                template="plotly_white",
                height=420,
                margin=dict(l=14, r=20, t=18, b=18),
                hoverlabel=dict(align="left", bgcolor="white", bordercolor="rgba(0,0,0,0.15)", font=dict(size=12))
            )

            st.plotly_chart(sankey_fig, width="stretch")

        st.divider()
        st.subheader("Final target list")
        st.caption(
            "Two views of the target list, switched below, at 2 genuinely different "
            "depths - not just 2 presentations of the same targets. **Proteogenomic "
            "only** stops at cis-MR + COLOC/PWCoCo + FinnGen/UKB safety (pQTL + GWAS "
            "evidence only) and deliberately goes no further. **Multi-omics** requires "
            "SMR/HEIDI and HyPrColoc's 3-trait confirmation on top of that - since SMR "
            "already draws on eQTL data, anything that reaches SMR belongs to the "
            "Multi-omics side, not Proteogenomic."
        )

        final_targets_view = st.segmented_control(
            "Which target list to show",
            ["Multi-omics (pQTL + GWAS + eQTL)", "Proteogenomic only (pQTL + GWAS)"],
            default="Multi-omics (pQTL + GWAS + eQTL)",
            selection_mode="single",
            key="final_targets_view_selector"
        )

        # segmented_control returns None if the user clicks the selected option again to
        # deselect it - fall back to the stricter multi-omics (passed HyPrColoc) view
        # rather than an undefined state
        show_hyprcoloc_targets = final_targets_view != "Proteogenomic only (pQTL + GWAS)"

        if show_hyprcoloc_targets:
            st.success(
                "**Multi-omics targets** - on top of cis-MR, COLOC, FinnGen/UKB safety and "
                f"SMR/HEIDI, these also passed HyPrColoc (posterior probability ≥ {hyprcoloc_pp_threshold:.2f}), "
                "meaning all 3 omics layers - pQTL (proteomics), GWAS (genomics) and eQTL "
                "(transcriptomics) - share a single causal variant. **Top SNP** is HyPrColoc's "
                "own *candidate SNP* - the single variant it found shared across the pQTL, GWAS "
                "and eQTL signals - with alleles and betas aligned to the AD risk allele."
            )

            # smr_display carries every protein ever run through SMR, including ones that
            # never reached (or failed) cis-MR/COLOC/FinnGen/UKB safety/SMR-HEIDI - restrict
            # to smr_eligible_set (same set the Sankey gates on) so this view matches its
            # own "passed cis-MR + COLOC + safety + SMR/HEIDI" claim above
            base_targets = smr_display.copy()
            identity_cols = [
                col for col in ["topsnp", "topsnp_chr", "topsnp_bp", "a1", "a2", "b_gwas", "b_eqtl"]
                if col in base_targets.columns
            ]
            base_targets = base_targets.drop(columns=identity_cols)

            if "protein" in base_targets.columns:
                base_targets = base_targets[base_targets["protein"].astype(str).isin(smr_eligible_set)]

            snp_info = select_hyprcoloc_candidate_rows(hyprcoloc_display, hyprcoloc_pp_threshold)
            snp_info_cols = available_cols(
                snp_info,
                ["protein", "cell_type", "data_type", "candidate_snp", "a1", "a2",
                 "gwas_beta", "gwas_p", "pqtl_beta", "pqtl_p", "eqtl_beta", "eqtl_p", "posterior_prob"]
            )
            snp_info = snp_info[snp_info_cols].rename(columns={
                "candidate_snp": "topsnp",
                "gwas_beta": "b_gwas",
                "eqtl_beta": "b_eqtl",
                "posterior_prob": "hyprcoloc_posterior_prob"
            })

            if not snp_info.empty:
                snp_info["snp_source"] = "HyPrColoc candidate SNP"

            if not snp_info.empty and "a1" not in snp_info.columns:
                st.warning(
                    "This pQTL dataset's HyPrColoc results don't carry the candidate SNP's "
                    "aligned alleles/betas yet, so, rerun the HyPrColoc pipeline step "
                    "(bin/hyprcoloc_targets.py) to populate them. Showing the candidate SNP ID "
                    "only until then."
                )

            merge_cols = [col for col in ["protein", "cell_type", "data_type"] if col in base_targets.columns and col in snp_info.columns]
            final_targets = (
                base_targets.merge(snp_info, on=merge_cols, how="inner")
                if merge_cols and not snp_info.empty
                else base_targets.iloc[0:0]
            )
        else:
            st.info(
                "**Proteogenomic-only targets** - passed cis-MR, COLOC/PWCoCo and FinnGen/UKB "
                "safety on the pQTL + GWAS layers alone (proteomics + genomics). This view "
                "stops deliberately *before* SMR - SMR/HEIDI already draws on eQTL data, so a "
                "target that reaches SMR (whether or not it goes on to pass HyPrColoc) shows "
                "up in the **Multi-omics** view instead, not here. 1 row per target (there is "
                "no cell-type/tissue dimension without SMR/eQTL data). **Top SNP** is always "
                "the target's own top cis-pQTL SNP, aligned to the AD risk allele (p-values "
                "are only ever floored to 1e-300 when reported as exactly 0)."
            )

            # target-level only (no cell_type/eQTL dimension at all - this view never
            # touches smr_display/hyprcoloc_display) - straight from target_info, scoped to
            # the safety-cleared set (ukb_pass_set), independent of whether SMR was ever run
            proteogenomic_cols = available_cols(
                target_info, ["protein", "snp", "a1", "a2", "gwas_beta", "gwas_p", "pqtl_beta", "pqtl_p"]
            )
            final_targets = (
                target_info[proteogenomic_cols].rename(columns={"snp": "topsnp", "gwas_beta": "b_gwas"})
                if not target_info.empty and "protein" in proteogenomic_cols
                else pd.DataFrame()
            )

            if not final_targets.empty:
                final_targets = final_targets[final_targets["protein"].astype(str).isin(ukb_pass_set)].copy()
                final_targets["snp_source"] = "cis-pQTL top SNP"

        if final_targets.empty:
            if show_hyprcoloc_targets:
                st.info(
                    "No targets currently pass every gate including HyPrColoc (posterior "
                    f"probability ≥ {hyprcoloc_pp_threshold:.2f}) at the selected thresholds."
                )
            else:
                st.info(
                    "No targets currently pass cis-MR, COLOC/PWCoCo and FinnGen/UKB safety "
                    "at the selected thresholds."
                )
        else:
            sort_cols = [col for col in ["protein", "cell_type"] if col in final_targets.columns]
            final_targets = final_targets.sort_values(sort_cols, na_position="last") if sort_cols else final_targets

            with st.container(border=True):
                if "cell_type" in final_targets.columns:
                    metric1, metric2 = st.columns(2)
                    metric1.metric("Unique targets", safe_nunique(final_targets, "protein"))
                    metric2.metric("Target x cell/bulk rows", len(final_targets))
                else:
                    st.metric("Unique targets", safe_nunique(final_targets, "protein"))

            st.divider()

            final_cols = [
                "protein",
                "gene",
                "data_type",
                "cell_type",
                "snp_source",
                "topsnp",
                "a1",
                "a2",
                "b_gwas",
                "gwas_p",
                "pqtl_beta",
                "pqtl_p",
                "b_eqtl",
                "eqtl_p",
                "hyprcoloc_posterior_prob",
                "b_smr",
                "p_smr",
                "q_smr",
                "p_heidi"
            ]

            final_cols = available_cols(final_targets, final_cols)
            final_table = final_targets[final_cols].copy()

            final_column_names = {
                "protein": "Target",
                "gene": "Gene",
                "data_type": "Data type",
                "cell_type": "Cell type",
                "snp_source": "SNP source",
                "topsnp": "Top SNP",
                "a1": "Risk allele",
                "a2": "Other allele",
                "b_gwas": "GWAS beta",
                "gwas_p": "GWAS p-value",
                "pqtl_beta": "pQTL beta",
                "pqtl_p": "pQTL p-value",
                "b_eqtl": "eQTL beta",
                "eqtl_p": "eQTL p-value",
                "hyprcoloc_posterior_prob": "HyPrColoc posterior probability",
                "b_smr": "SMR beta",
                "p_smr": "SMR p-value",
                "q_smr": "SMR FDR",
                "p_heidi": "HEIDI p-value"
            }

            final_table = final_table.rename(columns=final_column_names)

            st.caption(
                "Betas are all aligned to the outcome (AD) risk allele shown in **Risk allele**. "
                "**SNP source** spells out, per row, which SNP that alignment (and the Top SNP / "
                "allele / beta columns) was computed at."
            )

            st.dataframe(
                final_table,
                hide_index=True,
                width="stretch",
                column_config={
                    "GWAS beta": st.column_config.NumberColumn(format="%.4f"),
                    "GWAS p-value": st.column_config.NumberColumn(format="%.3e"),
                    "pQTL beta": st.column_config.NumberColumn(format="%.4f"),
                    "pQTL p-value": st.column_config.NumberColumn(format="%.3e"),
                    "eQTL beta": st.column_config.NumberColumn(format="%.4f"),
                    "eQTL p-value": st.column_config.NumberColumn(format="%.3e"),
                    "HyPrColoc posterior probability": st.column_config.NumberColumn(format="%.4f"),
                    "SMR beta": st.column_config.NumberColumn(format="%.4f"),
                    "SMR p-value": st.column_config.NumberColumn(format="%.3e"),
                    "SMR FDR": st.column_config.NumberColumn(format="%.3e"),
                    "HEIDI p-value": st.column_config.NumberColumn(format="%.3e")
                }
            )

            st.download_button(
                label="Download final targets",
                data=final_table.to_csv(index=False, sep="\t"),
                file_name=f"{pqtl_dataset}_{phenotype}_final_targets.tsv",
                mime="text/tab-separated-values",
                key="download_final_targets",
                width="stretch"
            )

    with tab9:
        st.caption("PARALLEL METHOD · COMPLEMENTARY TO STAGE 2 (pQTL–GWAS COLOC)")
        st.subheader("PWCoCo (conditional coloc)")
        st.caption(
            "PWCoCo re-tests pQTL–GWAS colocalisation using a GCTA-COJO-style stepwise "
            "conditional analysis, which can separate multiple independent causal signals "
            "at a locus that standard pairwise COLOC (**Evidence by Stage → 2. pQTL–GWAS "
            "COLOC**) assumes is a single signal. It runs alongside standard COLOC on the "
            "same cis-MR-passing targets, not instead of it - a target that colocalises "
            "under either method is carried forward as a prioritised target (see the "
            "**Overview** tab), annotated with which method(s) supported it. It's kept as "
            "its own tab here, rather than nested under Evidence by Stage, since it's a "
            "parallel method rather than a downstream pipeline stage."
        )

        if pwcoco_display.empty:
            st.info("No PWCoCo results are available for this pQTL dataset yet.")
        else:
            pwcoco_filtered = pwcoco_display.copy()

            if "h4" in pwcoco_filtered.columns:
                pwcoco_filtered["h4"] = pd.to_numeric(pwcoco_filtered["h4"], errors="coerce")
                pwcoco_pass_rows = pwcoco_filtered[pwcoco_filtered["h4"].fillna(0) >= pp4]
            else:
                pwcoco_pass_rows = pwcoco_filtered.iloc[0:0]

            with st.container(border=True):
                metric1, metric2, metric3, metric4 = st.columns(4)
                metric1.metric("Unique targets tested", safe_nunique(pwcoco_filtered, "protein"))
                metric2.metric("Conditional signal rows", len(pwcoco_filtered))
                metric3.metric("Passing PP.H4 threshold", len(pwcoco_pass_rows))
                metric4.metric("Median PP.H4", safe_median(pwcoco_filtered, "h4"))
            st.caption("\"Unique targets tested\" above is every protein PWCoCo ran on - no cis-MR gate applied yet.")

            st.divider()
            st.subheader("Concordance with standard COLOC")
            st.caption("Scoped to targets that already passed cis-MR *and* cleared PP.H4 via COLOC and/or PWCoCo - not the broader \"tested\" count above.")

            if "coloc_support" in mr_coloc_pass.columns and not mr_coloc_pass.empty:
                support_counts = mr_coloc_pass.drop_duplicates("protein")["coloc_support"].value_counts()
                with st.container(border=True):
                    concordance_col1, concordance_col2, concordance_col3 = st.columns(3)
                    with concordance_col1:
                        st.metric("Both methods agree", int(support_counts.get("both", 0)))
                        st.badge("Both", color="green")
                    with concordance_col2:
                        st.metric("COLOC only", int(support_counts.get("coloc_only", 0)))
                        st.badge("COLOC only", color="blue")
                    with concordance_col3:
                        st.metric("PWCoCo only", int(support_counts.get("pwcoco_only", 0)))
                        st.badge("PWCoCo only", color="orange")
                st.caption(
                    "\"COLOC only\" / \"PWCoCo only\" targets are discordant between the two "
                    "methods but are still carried forward as prioritised targets, not dropped - "
                    "the disagreement usually reflects a genuine methodological difference "
                    "(single- vs multi-signal locus assumption), not necessarily absence of a "
                    "true target."
                )
            else:
                st.info("No cis-MR + COLOC/PWCoCo overlap to compare yet.")

            st.divider()
            st.subheader("PWCoCo results")
            st.caption(
                f"A row **passes** when its conditional PP.H4 ≥ {pp4:.2f} (same threshold as "
                "standard COLOC's PP.H4, set in the sidebar). A target can have more than 1 "
                "row when PWCoCo's stepwise selection found more than 1 conditionally "
                "independent signal at its locus."
            )

            # dataset1/dataset2 dropped entirely, not just deprioritised - they're the
            # internal scratch-file names PWCoCo.pwcoco() writes to a tempdir
            # ("sumstats1.txt"/"sumstats2.txt"), constant on every single row, so they
            # carry zero information here and are pure clutter, not an audit trail
            pwcoco_filtered = pwcoco_filtered.drop(columns=["dataset1", "dataset2"], errors="ignore")

            pwcoco_cols = [
                "protein",
                "snp1",
                "snp2",
                "nsnps",
                "h0",
                "h1",
                "h2",
                "h3",
                "h4",
                "log_abf_all"
            ]

            pwcoco_cols = available_cols(pwcoco_filtered, pwcoco_cols)
            remaining_pwcoco_cols = [col for col in pwcoco_filtered.columns if col not in pwcoco_cols]

            pwcoco_table = pwcoco_filtered[pwcoco_cols + remaining_pwcoco_cols].copy()

            if "h4" in pwcoco_table.columns:
                pwcoco_table = pwcoco_table.sort_values("h4", ascending=False, na_position="last")

            pwcoco_column_names = {
                "protein": "Target",
                "snp1": "Top SNP (dataset 1)",
                "snp2": "Top SNP (dataset 2)",
                "nsnps": "N SNPs",
                "h0": "PP.H0",
                "h1": "PP.H1",
                "h2": "PP.H2",
                "h3": "PP.H3",
                "h4": "PP.H4",
                "log_abf_all": "Log ABF (all)"
            }

            pwcoco_table = pwcoco_table.rename(columns=pwcoco_column_names)

            st.dataframe(
                pwcoco_table,
                width="stretch",
                hide_index=True,
                column_config={
                    "PP.H0": st.column_config.NumberColumn(format="%.4f"),
                    "PP.H1": st.column_config.NumberColumn(format="%.4f"),
                    "PP.H2": st.column_config.NumberColumn(format="%.4f"),
                    "PP.H3": st.column_config.NumberColumn(format="%.4f"),
                    "PP.H4": st.column_config.NumberColumn(format="%.4f")
                }
            )

            st.download_button(
                label="Download PWCoCo results",
                data=pwcoco_table.to_csv(index=False, sep="\t"),
                file_name=f"{pqtl_dataset}_{phenotype}_PWCoCo.tsv",
                mime="text/tab-separated-values",
                key="download_pwcoco_hits",
                width="stretch"
            )

    with tab10:
        st.caption("SAME QUESTION AS STAGE 6 (HYPRCOLOC) · ANSWERED A DIFFERENT WAY")
        st.subheader("PWCoCo-QTL: eQTL-level causal-variant triangulation")
        st.caption(
            "HyPrColoc and PWCoCo-QTL ask the identical question: across all 3 biological "
            "layers, pQTL (protein), eQTL (transcript) and GWAS (disease), is there a shared "
            "causal variant? They just answer it differently. **HyPrColoc** clusters all 3 "
            "traits at once, under the assumption that each trait has at most 1 causal variant "
            "in the region. **PWCoCo-QTL** instead runs PWCoCo pairwise, 3 times "
            "(pQTL-GWAS, eQTL-pQTL, eQTL-GWAS), each allowing more than 1 causal variant per "
            "trait via conditioning. A target is **triangulated** only when the exact same "
            "conditionally independent SNP clears the PP.H4 threshold in **all 3** of those "
            "pairwise analyses at once, not just some of them - anything less is not counted."
        )
        st.caption(
            "The 2 methods run side by side, not one after the other: a target reaches "
            "**Final Targets** if either one supports it. This is exactly how PWCoCo relates "
            "to standard pairwise COLOC on the pQTL-GWAS tab, one level up."
        )

        if hyprcoloc_testable_set:
            with st.container(border=True):
                st.markdown("**HyPrColoc vs. PWCoCo-QTL concordance**")
                concordance_col1, concordance_col2, concordance_col3, concordance_col4 = st.columns(4)
                with concordance_col1:
                    st.metric("Both agree", len(three_trait_both_set))
                    st.badge("Both", color="green")
                with concordance_col2:
                    st.metric("HyPrColoc only", len(three_trait_hyprcoloc_only_set))
                    st.badge("HyPrColoc only", color="blue")
                with concordance_col3:
                    st.metric("PWCoCo-QTL only", len(three_trait_pwcoco_qtl_only_set))
                    st.badge("PWCoCo-QTL only", color="orange")
                with concordance_col4:
                    st.metric("Neither", len(three_trait_fail_set))
                    st.badge("Neither", color="red")
            st.caption(
                "Scoped to every target with SMR support (both methods run on the same "
                "population). \"Neither\" is the only bucket that stops here; every other "
                "bucket reaches Final Targets, since only 1 of the 2 methods is required."
            )

        if pwcoco_eqtl_pqtl_display.empty and pwcoco_eqtl_gwas_display.empty:
            st.info("No PWCoCo-QTL results are available for this pQTL dataset yet.")
        else:
            st.divider()
            with st.container(border=True):
                metric1, metric2, metric3 = st.columns(3)
                metric1.metric("Targets tested (eQTL-pQTL)", safe_nunique(pwcoco_eqtl_pqtl_display, "protein"))
                metric2.metric("Targets tested (eQTL-GWAS)", safe_nunique(pwcoco_eqtl_gwas_display, "protein"))
                metric3.metric("Targets triangulated (all 3 combos)", len(triangulated_proteins))
            st.caption(
                "\"Targets tested\" is every protein PWCoCo was run on for that combo - no "
                "PP.H4 threshold applied yet. \"Triangulated\" requires all 3 combos to agree "
                "on 1 SNP, per the definition above."
            )

            st.divider()
            st.subheader("Shared-SNP triangulation")
            st.caption(
                "1 row per (target, SNP) pair where the same SNP clears the PP.H4 threshold "
                "(sidebar slider, currently "
                f"{pp4:.2f}) in pQTL-GWAS, eQTL-pQTL AND eQTL-GWAS simultaneously - every row "
                "here is a triangulated target, recomputed live as you move that slider."
            )

            if shared_snp_table.empty:
                st.info("No shared colocalising SNPs found across combos at this threshold yet.")
            else:
                shared_cols = available_cols(
                    shared_snp_table,
                    ["protein", "snp", "pqtl_gwas_h4", "eqtl_pqtl_h4", "eqtl_gwas_h4"]
                )
                shared_table = shared_snp_table[shared_cols].copy()

                shared_column_names = {
                    "protein": "Target",
                    "snp": "Shared SNP",
                    "pqtl_gwas_h4": "PP.H4 (pQTL-GWAS)",
                    "eqtl_pqtl_h4": "PP.H4 (eQTL-pQTL)",
                    "eqtl_gwas_h4": "PP.H4 (eQTL-GWAS)",
                }
                shared_table = shared_table.rename(columns=shared_column_names)

                st.dataframe(shared_table, width="stretch", hide_index=True)

                st.download_button(
                    label="Download shared-SNP triangulation",
                    data=shared_table.to_csv(index=False, sep="\t"),
                    file_name=f"{pqtl_dataset}_{phenotype}_PWCoCo_QTL_shared_snps.tsv",
                    mime="text/tab-separated-values",
                    key="download_pwcoco_qtl_shared",
                    width="stretch"
                )

            st.divider()
            st.subheader("Raw PWCoCo(eQTL-pQTL) and PWCoCo(eQTL-GWAS) results")
            st.caption(
                f"A row **passes** when its conditional PP.H4 is at least {pp4:.2f} (same "
                "threshold as the standard COLOC/PWCoCo tab). A target can have more than 1 row "
                "per combo when PWCoCo's stepwise selection found more than 1 conditionally "
                "independent signal, or when it has more than 1 SMR-passing eQTL source (e.g. "
                "multiple GTEx tissues)."
            )

            pwcoco_qtl_cols = [
                "protein", "eqtl_dataset", "cell_type",
                "snp1", "snp2", "nsnps", "h0", "h1", "h2", "h3", "h4", "log_abf_all"
            ]

            pwcoco_qtl_column_names = {
                "protein": "Target",
                "eqtl_dataset": "eQTL dataset",
                "cell_type": "Tissue / cell type",
                "snp1": "Top SNP (dataset 1)",
                "snp2": "Top SNP (dataset 2)",
                "nsnps": "N SNPs",
                "h0": "PP.H0",
                "h1": "PP.H1",
                "h2": "PP.H2",
                "h3": "PP.H3",
                "h4": "PP.H4",
                "log_abf_all": "Log ABF (all)"
            }

            eqtl_pqtl_col1, eqtl_gwas_col2 = st.columns(2)

            with eqtl_pqtl_col1:
                st.markdown("**eQTL-pQTL**")
                if pwcoco_eqtl_pqtl_display.empty:
                    st.info("No eQTL-pQTL PWCoCo results yet.")
                else:
                    cols = available_cols(pwcoco_eqtl_pqtl_display, pwcoco_qtl_cols)
                    table = pwcoco_eqtl_pqtl_display[cols].copy()
                    if "h4" in table.columns:
                        table = table.sort_values("h4", ascending=False, na_position="last")
                    table = table.rename(columns=pwcoco_qtl_column_names)
                    st.dataframe(table, width="stretch", hide_index=True)
                    st.download_button(
                        label="Download eQTL-pQTL PWCoCo results",
                        data=table.to_csv(index=False, sep="\t"),
                        file_name=f"{pqtl_dataset}_{phenotype}_PWCoCo_eQTL_pQTL.tsv",
                        mime="text/tab-separated-values",
                        key="download_pwcoco_eqtl_pqtl",
                        width="stretch"
                    )

            with eqtl_gwas_col2:
                st.markdown("**eQTL-GWAS**")
                if pwcoco_eqtl_gwas_display.empty:
                    st.info("No eQTL-GWAS PWCoCo results yet.")
                else:
                    cols = available_cols(pwcoco_eqtl_gwas_display, pwcoco_qtl_cols)
                    table = pwcoco_eqtl_gwas_display[cols].copy()
                    if "h4" in table.columns:
                        table = table.sort_values("h4", ascending=False, na_position="last")
                    table = table.rename(columns=pwcoco_qtl_column_names)
                    st.dataframe(table, width="stretch", hide_index=True)
                    st.download_button(
                        label="Download eQTL-GWAS PWCoCo results",
                        data=table.to_csv(index=False, sep="\t"),
                        file_name=f"{pqtl_dataset}_{phenotype}_PWCoCo_eQTL_GWAS.tsv",
                        mime="text/tab-separated-values",
                        key="download_pwcoco_eqtl_gwas",
                        width="stretch"
                    )

def main():
    p = argparse.ArgumentParser()
    p.add_argument("--db_name", required=True, type=str)
    p.add_argument("--port_number", required=True, type=str)
    p.add_argument("--phenotype", required=True, type=str)
    p.add_argument("--pqtl_dataset", required=True, type=str)
    args = p.parse_args()
    dashboard(
        db_name=args.db_name,
        port_number=args.port_number,
        phenotype=args.phenotype,
        pqtl_dataset=args.pqtl_dataset
    )


if __name__ == "__main__":
    main()