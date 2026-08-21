#!/usr/bin/env python3
import argparse
from pathlib import Path
import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
import streamlit as st

from drugmr import paths, registry

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

# single source of truth for "where am I in the pipeline" - the Overview tab's step
# map and every downstream tab's stage caption are both built from this list, so the
# two can never drift out of sync with each other or with the st.tabs() labels below
PIPELINE_STAGES = [
    dict(title="cis-MR", blurb="Mendelian randomisation of cis-instrumented protein abundance on the outcome."),
    dict(title="pQTL–GWAS COLOC", blurb="Pairwise colocalisation confirming the pQTL and GWAS signals share one causal variant."),
    dict(title="FinnGen PheWAS", blurb="Phenome-wide safety inference for adverse effects of the same risk allele in FinnGen."),
    dict(title="UKB PheWAS", blurb="Same phenome-wide safety scan, run in UK Biobank EHR-derived phenotypes."),
    dict(title="SMR (bulk/sc eQTL)", blurb="SMR + HEIDI test that the pQTL signal also acts through transcription."),
    dict(title="HyPrColoc (bulk/sc eQTL)", blurb="Joint colocalisation of pQTL + GWAS + eQTL signals onto one shared variant."),
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

def create_streamlit_ammenities(db_name: str, port_number: str):
    streamlit_dir = Path(".streamlit")
    streamlit_dir.mkdir(parents=True, exist_ok=True)

    secrets = f"""[connections.postgresql]
dialect = "postgresql"
host = "localhost"
port = "{port_number}"
database = "{db_name}"
username = ""
password = ""
"""

    # create local streamlit ammenities
    (streamlit_dir / "secrets.toml").write_text(secrets, encoding="utf-8")


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


def load_optional_tsv(file: Path, label: str):
    if not file.exists():
        st.warning(f"{label} result file not found: {file}")
        return pd.DataFrame()

    try:
        df = pd.read_csv(file, sep="\t", low_memory=False)
    except (OSError, pd.errors.ParserError, UnicodeDecodeError) as error:
        st.warning(f"{label} result file could not be read: {file}")
        st.exception(error)
        return pd.DataFrame()

    if df.empty:
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
    }


def resolve_dataset_files(project_dir: Path, phenotype: str, dataset_id: str, run_id: str = "latest"):
    """Resolve this dataset's 7 dashboard files via runs/registry.json first -
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

    if "bonferroni_significant" in df.columns:
        df["bonferroni_significant"] = (
            df["bonferroni_significant"]
            .astype(str)
            .str.lower()
            .isin(["true", "1", "yes"])
        )
    elif "p_bonferroni" in df.columns:
        df["bonferroni_significant"] = df["p_bonferroni"].fillna(np.inf) <= 0.05

    return df


def compute_phewas_safety_status(phewas_df: pd.DataFrame, proteins):
    """Per-protein FinnGen/UKB PheWAS safety flag for the prioritisation Sankey.

    A protein FAILS only if it has a Bonferroni-significant PheWAS association
    with beta_mr >= 0 (the same allele that raises AD risk also raises the
    safety-relevant phenotype). No significant association, a significant but
    protective (negative) beta, or no PheWAS coverage at all all count as PASS.
    """
    proteins = list(proteins)
    status = {protein: "pass" for protein in proteins}

    if phewas_df.empty or "protein" not in phewas_df.columns or "beta_mr" not in phewas_df.columns:
        return status

    if "bonferroni_significant" not in phewas_df.columns:
        return status

    df = phewas_df.copy()
    df["protein"] = df["protein"].astype(str)

    adverse = df["bonferroni_significant"].fillna(False) & (df["beta_mr"].fillna(0) >= 0)
    failing_proteins = set(df.loc[adverse, "protein"].unique())

    for protein in proteins:
        if protein in failing_proteins:
            status[protein] = "fail"

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


def subset_phewas_outcome(df: pd.DataFrame, outcome: str):
    if df.empty:
        return df

    if "outcome_trait" in df.columns:
        return df[df["outcome_trait"] == outcome].copy()

    if "pheno_id" in df.columns:
        return df[df["pheno_id"] == outcome].copy()

    return df.copy()


def probability_strength_color(value):
    """Tiering for an already-passing probability (COLOC PP.H4, HyPrColoc posterior).

    Every row shown has already cleared the sidebar threshold, so this isn't a
    pass/fail color - it's how comfortably a row cleared it, so the strongest
    colocalisation evidence is visually distinguishable from a borderline one.
    """
    if pd.isna(value):
        return "gray"

    if value >= 0.90:
        return "green"

    if value >= 0.80:
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


def render_prioritised_target_cards(df: pd.DataFrame, outcome: str, n_columns: int = 3):
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

                badge_col1, badge_col2 = st.columns(2)

                with badge_col1:
                    pp_label = f"PP.H4 {pp_h4:.2f}" if pd.notna(pp_h4) else "PP.H4 NA"
                    st.badge(pp_label, color=probability_strength_color(pp_h4))

                with badge_col2:
                    fdr_label = f"MR FDR {mr_fdr:.1e}" if pd.notna(mr_fdr) else "MR FDR NA"
                    st.badge(fdr_label, color=significance_strength_color(mr_fdr))

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
    source_name: str,
    source_description: str,
    n_endpoints: int,
    pqtl_dataset: str,
    outcome: str,
    key_prefix: str,
    stage_number: int
):
    stage_caption(stage_number)
    st.subheader(f"{source_name} PheWAS safety and repurposing profile")

    if phewas_outcome.empty:
        st.info(f"No local {source_name} PheWAS safety results are available for this outcome.")
        return

    if "protein" not in phewas_outcome.columns:
        st.error(f"The {source_name} PheWAS result file does not contain a protein column.")
        return

    phewas_targets = sorted(phewas_outcome["protein"].dropna().astype(str).unique())

    if len(phewas_targets) == 0:
        st.info(f"No proteins were found in the {source_name} PheWAS safety table.")
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

    target_phewas = phewas_outcome[
        phewas_outcome["protein"].astype(str) == selected_phewas_target
    ].copy()

    p_col = "p_mr" if "p_mr" in target_phewas.columns else None
    beta_col = "beta_mr" if "beta_mr" in target_phewas.columns else None
    bonferroni_col = "p_bonferroni" if "p_bonferroni" in target_phewas.columns else None

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

    if bonferroni_col is not None:
        target_phewas["bonferroni_significant"] = target_phewas[bonferroni_col].fillna(np.inf) <= 0.05
    elif "bonferroni_significant" not in target_phewas.columns:
        target_phewas["bonferroni_significant"] = False

    phenotype_col = "phenostring" if "phenostring" in target_phewas.columns else "phenocode"

    if phenotype_col not in target_phewas.columns:
        st.error(f"The {source_name} PheWAS result file does not contain a phenotype column.")
        return

    category_col = "category" if "category" in target_phewas.columns else None
    n_phenotypes = target_phewas[phenotype_col].nunique()
    n_nominal = int((target_phewas[p_col] < 0.05).sum())
    n_bonferroni = int(target_phewas["bonferroni_significant"].sum())

    with st.container(border=True):
        metric1, metric2, metric3 = st.columns(3)
        metric1.metric(f"{source_name} phenotypes tested", int(n_phenotypes))
        metric2.metric("Nominal associations", n_nominal)
        metric3.metric("Bonferroni-significant associations", n_bonferroni)

    st.caption(
        f"PheWAS MR estimates show the effect of genetically predicted protein abundance "
        f"on each {source_description}. Wald ratio is used for targets with one available "
        "cis-MR instrument and IVW is used for targets with > 1."
    )

    plot_kwargs = {
        "data_frame": target_phewas,
        "x": beta_col,
        "y": "minus_log10_p",
        "hover_name": phenotype_col,
        "symbol": "bonferroni_significant",
        "hover_data": {
            beta_col: ":.4f",
            p_col: ":.3e",
            "minus_log10_p": False,
            "bonferroni_significant": True
        },
        "labels": {
            beta_col: "PheWAS MR beta",
            "minus_log10_p": "-log10(PheWAS p-value)",
            "bonferroni_significant": "Bonferroni significant"
        },
        "title": f"{source_name} PheWAS profile: {selected_phewas_target}",
        "height": 600,
        "template": "plotly_white"
    }

    if "phenocode" in target_phewas.columns:
        plot_kwargs["hover_data"]["phenocode"] = True

    if category_col is not None:
        plot_kwargs["color"] = category_col
        plot_kwargs["labels"][category_col] = f"{source_name} category"

    st.divider()
    st.subheader("Phenome-wide association landscape")

    phewas_fig = px.scatter(**plot_kwargs)
    phewas_fig.add_hline(y=-np.log10(0.05 / n_endpoints), line_dash="dash", line_color="grey")
    phewas_fig.add_vline(x=0, line_dash="dash", line_color="grey")
    st.plotly_chart(phewas_fig, width="stretch")

    st.subheader("Bonferroni-significant PheWAS associations")
    top_phewas = target_phewas[target_phewas["bonferroni_significant"]].copy()
    top_phewas = top_phewas.sort_values(
        [beta_col, bonferroni_col if bonferroni_col is not None else p_col],
        ascending=[True, True]
    )

    if top_phewas.empty:
        st.info(
            f"No {source_name} phenotype associations survive Bonferroni correction across "
            f"{n_endpoints:,} endpoints for {selected_phewas_target}."
        )
    else:
        top_plot_kwargs = {
            "data_frame": top_phewas,
            "x": beta_col,
            "y": phenotype_col,
            "hover_data": {
                beta_col: ":.4f",
                p_col: ":.3e",
                "minus_log10_p": ":.3f"
            },
            "labels": {
                beta_col: "PheWAS MR beta",
                phenotype_col: ""
            },
            "title": "Bonferroni-significant PheWAS associations",
            "height": max(450, 45 * len(top_phewas)),
            "template": "plotly_white"
        }

        if "phenocode" in top_phewas.columns:
            top_plot_kwargs["hover_data"]["phenocode"] = True

        if category_col is not None:
            top_plot_kwargs["color"] = category_col
            top_plot_kwargs["labels"][category_col] = f"{source_name} category"

        top_phewas_fig = px.scatter(**top_plot_kwargs)
        top_phewas_fig.add_vline(x=0, line_dash="dash", line_color="grey")
        st.plotly_chart(top_phewas_fig, width="stretch")

    phewas_cols = [
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
        "p_bonferroni",
        "bonferroni_significant"
    ]

    phewas_cols = available_cols(target_phewas, phewas_cols)
    significant_phewas = target_phewas[target_phewas["bonferroni_significant"]].copy()
    significant_phewas = significant_phewas.sort_values(bonferroni_col if bonferroni_col is not None else p_col, ascending=True)

    if significant_phewas.empty:
        st.success(
            f"No {source_name} phenotype associations survive Bonferroni correction across "
            f"{n_endpoints:,} endpoints for this target."
        )
    else:
        st.dataframe(
            significant_phewas[phewas_cols],
            width="stretch",
            hide_index=True
        )

    with st.expander(f"View all {source_name} PheWAS associations"):
        remaining_cols = [col for col in target_phewas.columns if col not in phewas_cols]
        st.dataframe(
            target_phewas[phewas_cols + remaining_cols].sort_values(p_col, ascending=True),
            width="stretch",
            hide_index=True
        )

    st.download_button(
        label=f"Download {selected_phewas_target} {source_name} PheWAS results",
        data=target_phewas.to_csv(index=False, sep="\t"),
        file_name=f"{selected_phewas_target}_{outcome}_{source_name.replace(' ', '_')}_PheWAS.tsv",
        mime="text/tab-separated-values",
        key=f"{key_prefix}_download_phewas_{pqtl_dataset}_{outcome}_{selected_phewas_target}",
        width="stretch"
    )


def dashboard(db_name: str, port_number: str, phenotype: str, pqtl_dataset: str):
    mr_table = "cis_mr_results"
    coloc_table = "coloc_results"
    finngen_phewas_table = "finngen_phewas_safety"
    ukb_phewas_table = "ukb_phewas_safety"

    # main aesthetics
    # native Streamlit only - no HTML or CSS
    st.set_page_config(
        page_title=f"{db_name}",
        page_icon="",
        layout="wide",
        initial_sidebar_state="expanded"
    )

    try:
        conn = st.connection(
            "postgresql",
            type="sql",
            url=f"postgresql://localhost:{port_number}/{db_name}"
        )
    except Exception as error:
        st.error("The PostgreSQL connection could not be initialised.")
        st.exception(error)
        st.stop()

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
            _, dataset_result_files[pqtl_dataset] = resolve_dataset_files(
                project_dir, phenotype, pqtl_dataset, run_id=selected_run
            )

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

    # load local result files into PostgreSQL for the dashboard
    mr = load_required_tsv(mr_file, "cis-MR")
    coloc = load_required_tsv(coloc_file, "pQTL–GWAS COLOC")
    finngen_phewas = load_optional_tsv(finngen_phewas_file, "FinnGen PheWAS safety")
    ukb_phewas = load_optional_tsv(ukb_phewas_file, "UKB PheWAS safety")
    target_info = load_optional_tsv(target_info_file, "Harmonised target information")
    smr = load_optional_tsv(smr_file, "SMR (bulk/sc eQTL)")
    hyprcoloc = load_optional_tsv(hyprcoloc_file, "HyPrColoc (bulk/sc eQTL)")

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


    # refresh dashboard tables
    try:
        mr.to_sql(mr_table, conn.engine, if_exists="replace", index=False)
        coloc.to_sql(coloc_table, conn.engine, if_exists="replace", index=False)
    except Exception as error:
        st.error("The cis-MR or COLOC dashboard table could not be refreshed.")
        st.exception(error)
        st.stop()

    finngen_phewas_available = not finngen_phewas.empty
    ukb_phewas_available = not ukb_phewas.empty

    if finngen_phewas_available:
        try:
            finngen_phewas.to_sql(finngen_phewas_table, conn.engine, if_exists="replace", index=False)
        except Exception as error:
            st.warning("The FinnGen PheWAS dashboard table could not be refreshed.")
            st.exception(error)
            finngen_phewas_available = False

    if ukb_phewas_available:
        try:
            ukb_phewas.to_sql(ukb_phewas_table, conn.engine, if_exists="replace", index=False)
        except Exception as error:
            st.warning("The UKB PheWAS dashboard table could not be refreshed.")
            st.exception(error)
            ukb_phewas_available = False

    with st.sidebar.expander("Tracking", expanded=False):
        st.write(f"Loaded {len(mr)} rows into {mr_table}")
        st.write(f"Loaded {len(coloc)} rows into {coloc_table}")

        if finngen_phewas_available:
            st.write(f"Loaded {len(finngen_phewas)} rows into {finngen_phewas_table}")

        if ukb_phewas_available:
            st.write(f"Loaded {len(ukb_phewas)} rows into {ukb_phewas_table}")

        if not target_info.empty:
            st.write(f"Loaded {len(target_info)} harmonised top cis-hit rows")

        if not smr.empty:
            st.write(f"Loaded {len(smr)} SMR (bulk/sc eQTL) rows")

        if not hyprcoloc.empty:
            st.write(f"Loaded {len(hyprcoloc)} HyPrColoc (bulk/sc eQTL) rows")

    # load MR + COLOC results
    mr = conn.query(f"SELECT * FROM {mr_table};", ttl=0)
    coloc = conn.query(f"SELECT * FROM {coloc_table};", ttl=0)

    if finngen_phewas_available:
        finngen_phewas = prepare_phewas(conn.query(f"SELECT * FROM {finngen_phewas_table};", ttl=0))
    else:
        finngen_phewas = pd.DataFrame()

    if ukb_phewas_available:
        ukb_phewas = prepare_phewas(conn.query(f"SELECT * FROM {ukb_phewas_table};", ttl=0))
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

    with st.sidebar.expander("Target prioritisation thresholds", expanded=True):
        fdr = st.slider("MR FDR threshold", 0.0, 1.0, 0.05, 0.01)
        q_pval = st.slider("Minimum Cochran Q p-value", 0.0, 1.0, 0.05, 0.01)
        pp4 = st.slider("pQTL–GWAS COLOC PP.H4 threshold", 0.0, 1.0, 0.70, 0.01)
        smr_fdr_threshold = st.slider("SMR FDR (q_SMR) threshold", 0.0, 1.0, 0.05, 0.01)
        heidi_p_threshold = st.slider("Minimum HEIDI p-value", 0.0, 1.0, 0.01, 0.01)
        hyprcoloc_pp_threshold = st.slider("HyPrColoc posterior probability threshold", 0.0, 1.0, 0.50, 0.01)

    with st.sidebar.expander("Protein filter", expanded=True):
        protein = st.text_input(
            "Protein search",
            placeholder="e.g. PILRA, APOE or Q9UKJ1"
        )

    st.header(f"{dataset_name} → {outcome}")
    st.caption(
        f"N = {dataset_n:,} | MR FDR ≤ {fdr:.2f} | Q p ≥ {q_pval:.2f} | "
        f"pQTL–GWAS PP.H4 ≥ {pp4:.2f}"
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

    # proteins which pass both MR + COLOC thresholds
    mr_coloc_pass = mr_pass.merge(
        coloc_pass,
        on="protein",
        how="inner",
        suffixes=("_mr", "_pqtl_coloc")
    )

    # preserve assay-specific protein IDs
    # only remove fully duplicated merged rows
    mr_coloc_pass = mr_coloc_pass.drop_duplicates()

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

    # eQTLGen is no longer a configured bulk_eqtl_dataset (assets/config.yaml), but
    # stale rows from when it was configured persist in the combined SMR table since
    # compile_multi_omics_targets() only replaces rows for the dataset it's currently
    # processing - filter it out here so it never surfaces on the dashboard
    if "eqtl_dataset" in smr_display.columns:
        smr_display = smr_display[smr_display["eqtl_dataset"].str.lower() != "eqtlgen"]

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

    # round coloc posterior probs
    for col in coloc_numeric_cols:
        if col in coloc_pass.columns:
            coloc_pass[col] = coloc_pass[col].round(3)

        if col in mr_coloc_pass.columns:
            mr_coloc_pass[col] = mr_coloc_pass[col].round(3)

    # Overview-only PheWAS safety gate: a target that passes cis-MR + pQTL COLOC
    # is only "prioritised" on the front page if it also has no Bonferroni-
    # significant FinnGen or UKB PheWAS MR hit with beta_mr >= 0 (same allele
    # raising both the target and the adverse phenotype). Mirrors the
    # FinnGen/UKB stages of the Final Targets Sankey, but applied here so it
    # gates the Overview tab's card list too - mr_coloc_pass itself is left
    # untouched since the Sankey needs the un-gated COLOC-pass set to draw its
    # own FinnGen/UKB drop-off columns.
    mr_coloc_pass_proteins = (
        set(mr_coloc_pass["protein"].dropna().astype(str))
        if "protein" in mr_coloc_pass.columns else set()
    )

    overview_finngen_status = compute_phewas_safety_status(finngen_phewas_outcome, mr_coloc_pass_proteins)
    overview_finngen_safe = {
        protein for protein in mr_coloc_pass_proteins
        if overview_finngen_status.get(protein) != "fail"
    }

    overview_ukb_status = compute_phewas_safety_status(ukb_phewas_outcome, overview_finngen_safe)
    overview_phewas_safe_proteins = {
        protein for protein in overview_finngen_safe
        if overview_ukb_status.get(protein) != "fail"
    }

    mr_coloc_safe_pass = (
        mr_coloc_pass[mr_coloc_pass["protein"].astype(str).isin(overview_phewas_safe_proteins)].copy()
        if "protein" in mr_coloc_pass.columns else mr_coloc_pass.copy()
    )

    # main staged target counts
    n_tested = safe_nunique(mr_outcome, "protein")
    n_mr = safe_nunique(mr_pass, "protein")
    n_mr_coloc = safe_nunique(mr_coloc_pass, "protein")
    n_mr_coloc_safe = safe_nunique(mr_coloc_safe_pass, "protein")
    n_finngen_phewas = safe_nunique(finngen_phewas_outcome, "protein")
    n_ukb_phewas = safe_nunique(ukb_phewas_outcome, "protein")

    tab1, tab2, tab3, tab4, tab5, tab6, tab7, tab8 = st.tabs([
        "Overview",
        "1. cis-MR",
        "2. pQTL–GWAS COLOC",
        "3. FinnGen PheWAS",
        "4. UKB PheWAS",
        "5. SMR (bulk/sc eQTL)",
        "6. HyPrColoc (bulk/sc eQTL)",
        "7. Final Targets"
    ])

    with tab1:
        st.caption("OVERVIEW · START HERE")
        st.subheader("How this dashboard is organised")
        st.caption(
            "Every protein below moves through the same 7-stage pipeline, in order. A "
            "target only reaches a later stage once it has already passed everything "
            "before it, and each numbered tab shows that stage's own results in full."
        )

        pipeline_columns = st.columns(len(PIPELINE_STAGES))

        for position, (pipeline_column, stage) in enumerate(zip(pipeline_columns, PIPELINE_STAGES), start=1):
            with pipeline_column:
                with st.container(border=True):
                    st.markdown(f"**{position}. {stage['title']}**")
                    st.caption(stage["blurb"])

        st.divider()
        st.subheader("Target prioritisation")
        st.caption(
            "Targets move from cis-MR testing to MR support, then to shared pQTL–GWAS "
            "causal signal support through pairwise COLOC (stages 1-2 of the pipeline "
            "above), then through a FinnGen/UKB PheWAS safety check (stages 3-4). The "
            "complete flow through every stage, including SMR and HyPrColoc, is drawn as "
            "a Sankey diagram on the **7. Final Targets** tab."
        )

        with st.container(border=True):
            metric1, metric2, metric3, metric4 = st.columns(4)
            metric1.metric("Proteins tested by cis-MR", n_tested)
            metric2.metric("cis-MR supported", n_mr, f"{retention(n_mr, n_tested):.1f}% of tested", delta_color="off")
            metric3.metric("cis-MR + pQTL COLOC", n_mr_coloc, f"{retention(n_mr_coloc, n_mr):.1f}% retained", delta_color="off")
            metric4.metric(
                "+ PheWAS safe (FinnGen/UKB)",
                n_mr_coloc_safe,
                f"{retention(n_mr_coloc_safe, n_mr_coloc):.1f}% retained",
                delta_color="off"
            )

        st.divider()

        funnel_df = pd.DataFrame({
            "stage": [
                "Proteins tested by cis-MR",
                "cis-MR supported",
                "cis-MR + pQTL COLOC",
                "+ PheWAS safe (FinnGen/UKB)"
            ],
            "n_targets": [
                n_tested,
                n_mr,
                n_mr_coloc,
                n_mr_coloc_safe
            ]
        })

        funnel_fig = px.bar(
            funnel_df,
            x="n_targets",
            y="stage",
            orientation="h",
            text="n_targets",
            title="Progressive target prioritisation",
            labels={"n_targets": "Number of unique proteins", "stage": ""},
            height=420,
            template="plotly_white"
        )

        funnel_fig.update_yaxes(
            categoryorder="array",
            categoryarray=funnel_df["stage"][::-1]
        )
        funnel_fig.update_traces(textposition="outside")
        funnel_fig.update_layout(showlegend=False, margin=dict(l=20, r=40, t=60, b=20))
        st.plotly_chart(funnel_fig, width="stretch")

        st.divider()
        st.subheader("Prioritised targets")
        st.caption(
            "Every protein that has passed cis-MR and pairwise pQTL–GWAS COLOC at the "
            "thresholds set in the sidebar, and has no Bonferroni-significant FinnGen or "
            "UKB PheWAS MR hit with a beta in the same direction (i.e. no evidence the "
            "same risk allele also drives an adverse phenotype), 1 card each. Betas and "
            "alleles are harmonised to the outcome GWAS risk allele."
        )

        if not mr_coloc_safe_pass.empty:
            st.success(
                f"{n_mr_coloc_safe} unique target(s) passed the selected cis-MR and pairwise COLOC "
                "thresholds with no adverse Bonferroni-significant FinnGen/UKB PheWAS signal."
            )

            if n_mr_coloc_safe < n_mr_coloc:
                st.caption(
                    f"{n_mr_coloc - n_mr_coloc_safe} additional target(s) passed cis-MR + COLOC but were "
                    "excluded here for a Bonferroni-significant, same-direction FinnGen/UKB PheWAS hit - "
                    "see the Sankey diagram on the **7. Final Targets** tab for the full breakdown."
                )

            cards_df = mr_coloc_safe_pass.sort_values(
                ["pp_h4_abf", "mr_fdr_q"],
                ascending=[False, True],
                na_position="last"
            ) if "pp_h4_abf" in mr_coloc_safe_pass.columns else mr_coloc_safe_pass

            render_prioritised_target_cards(cards_df, outcome)

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
                "pp_h4_abf"
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
                "pp_h4_abf": "COLOC PP.H4"
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
                "excluded here for a Bonferroni-significant, same-direction FinnGen/UKB PheWAS hit - see the "
                "Sankey diagram on the **7. Final Targets** tab for the full breakdown."
            )
        else:
            st.info("No proteins currently pass both the selected cis-MR and pQTL COLOC thresholds.")

        if finngen_phewas_available or ukb_phewas_available:
            st.caption(
                f"FinnGen PheWAS results are available for {n_finngen_phewas} unique target(s); "
                f"UKB PheWAS results are available for {n_ukb_phewas} unique target(s)."
            )

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

        st.dataframe(
            mr_display[display_cols + remaining_cols],
            width="stretch",
            hide_index=True
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
        st.caption("Targets shown here pass both the selected cis-MR and pairwise pQTL–GWAS COLOC thresholds.")

        if not mr_coloc_pass.empty:
            with st.container(border=True):
                col1, col2, col3 = st.columns(3)
                col1.metric("Prioritised proteins", mr_coloc_pass["protein"].nunique())
                col2.metric("Median PP.H4", safe_median(mr_coloc_pass, "pp_h4_abf"))
                col3.metric("Median MR FDR", safe_median(mr_coloc_pass, "mr_fdr_q", scientific=True))

            prioritised_cols = [
                "protein",
                "mr_method",
                "n_instruments",
                "mr_beta",
                "mr_se",
                "mr_pval",
                "mr_fdr_q",
                "q_pval",
                "egger_intercept_pval",
                "top_snp",
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

            st.dataframe(
                mr_coloc_pass[prioritised_cols + remaining_cols],
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
            source_name="FinnGen",
            source_description="FinnGen phenotype (ICD coded)",
            n_endpoints=2511,
            pqtl_dataset=pqtl_dataset,
            outcome=outcome,
            key_prefix="finngen",
            stage_number=3
        )

    with tab5:
        render_phewas_section(
            phewas_outcome=ukb_phewas_outcome,
            mr_coloc_pass=mr_coloc_pass,
            source_name="UKB TOPMed",
            source_description="UK Biobank EHR-derived PheCode phenotype",
            n_endpoints=1419,
            pqtl_dataset=pqtl_dataset,
            outcome=outcome,
            key_prefix="ukb",
            stage_number=4
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
            smr_filtered = smr_display.copy()

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
            "aligned onto a common effect allele."
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
                hyprcoloc_filtered["all_3_traits_colocalise"] = (
                    traits_lower.str.contains("pqtl_") &
                    traits_lower.str.contains("gwas_") &
                    traits_lower.str.contains("eqtl_")
                )
            else:
                hyprcoloc_filtered["all_3_traits_colocalise"] = False

            if "posterior_prob" in hyprcoloc_filtered.columns:
                hyprcoloc_pass_rows = hyprcoloc_filtered[
                    hyprcoloc_filtered["all_3_traits_colocalise"] &
                    (pd.to_numeric(hyprcoloc_filtered["posterior_prob"], errors="coerce").fillna(0) >= hyprcoloc_pp_threshold)
                ]
            else:
                hyprcoloc_pass_rows = hyprcoloc_filtered.iloc[0:0]

            with st.container(border=True):
                metric1, metric2, metric3, metric4 = st.columns(4)
                metric1.metric("Unique targets tested", safe_nunique(hyprcoloc_filtered, "protein"))
                metric2.metric("Target x cell-type rows", len(hyprcoloc_filtered))
                metric3.metric("Passing HyPrColoc threshold", len(hyprcoloc_pass_rows))
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
                "cell_type",
                "traits",
                "all_3_traits_colocalise",
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
                "cell_type": "Cell type",
                "traits": "Clustered traits",
                "all_3_traits_colocalise": "pQTL + GWAS + eQTL cluster",
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
                "- **cis-MR** —> this flow starts from the proteins that already passed; the "
                "much larger screening drop-off across every tested protein is the funnel in "
                "the Overview tab.\n"
                "- **pQTL–GWAS COLOC** —> passes on the posterior-probability threshold set in "
                "the sidebar.\n"
                "- **FinnGen / UKB PheWAS safety** —> fails *only* when a target has a "
                "Bonferroni-significant association with `beta_mr >= 0`, i.e. the same allele "
                "that raises AD risk also raises the safety-relevant phenotype. No significant "
                "hit, a protective (negative) significant hit, or no PheWAS coverage at all "
                "each count as passing.\n"
                f"- **SMR support** —> requires SMR FDR (`q_SMR`) < {smr_fdr_threshold:.2f} and "
                f"HEIDI p-value > {heidi_p_threshold:.2f}, split by whether that support came "
                "from bulk/tissue eQTL data, single-cell eQTL data, or both.\n"
                "- **HyPrColoc** —> runs against whichever eQTL dataset(s) supported the target's "
                "SMR stage (bulk, single-cell, or both); no-SMR-support targets end at the SMR "
                "stage. Passes when a HyPrColoc cluster contains the pQTL, GWAS *and* eQTL trait "
                f"together (not just 2 of the 3) with posterior probability ≥ {hyprcoloc_pp_threshold:.2f}."
            )

        n_sankey_mr_pass = safe_nunique(mr_pass, "protein")

        # smr_eligible_set (passed everything up to and including SMR/HEIDI, HyPrColoc
        # notwithstanding) is populated below when the flow can be drawn - stays empty
        # otherwise so the target-list toggle further down always has something to check
        smr_eligible_set = set()

        if n_sankey_mr_pass == 0:
            st.info("No proteins currently pass the cis-MR threshold, so no flow can be drawn.")
        else:
            mr_pass_proteins = set(mr_pass["protein"].dropna().astype(str)) if "protein" in mr_pass.columns else set()
            coloc_pass_set = set(mr_coloc_pass["protein"].dropna().astype(str)) if "protein" in mr_coloc_pass.columns else set()

            finngen_status = compute_phewas_safety_status(finngen_phewas_outcome, coloc_pass_set)
            finngen_pass_set = {protein for protein in coloc_pass_set if finngen_status.get(protein) != "fail"}

            ukb_status = compute_phewas_safety_status(ukb_phewas_outcome, finngen_pass_set)
            ukb_pass_set = {protein for protein in finngen_pass_set if ukb_status.get(protein) != "fail"}

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

            # HyPrColoc now runs for any target with SMR support, whether that came
            # from bulk eQTL, single-cell eQTL, or both - no-SMR-support targets are
            # terminal at the SMR stage
            smr_eligible_set = both_set | sc_only_set | bulk_only_set
            hyprcoloc_status = compute_hyprcoloc_pass_status(hyprcoloc_display, smr_eligible_set, hyprcoloc_pp_threshold)
            hyprcoloc_pass_set = {protein for protein in smr_eligible_set if hyprcoloc_status.get(protein)}
            hyprcoloc_fail_set = smr_eligible_set - hyprcoloc_pass_set

            # single source of truth for the diagram, its hover text and the
            # drill-down selector below - every count is derived from these sets.
            # `label` is what the chart prints (kept short so five columns fit
            # without wrapping), `name` is the formal stage name used on hover
            # and in the selector.
            sankey_groups = [
                dict(key="mr", column=0, name="cis-MR passed", label="cis-MR passed",
                     proteins=mr_pass_proteins, color=STATUS_MUTED, dropout=False),
                dict(key="coloc_pass", column=1, name="COLOC passed", label="COLOC passed",
                     proteins=coloc_pass_set, color=STATUS_GOOD, dropout=False),
                dict(key="coloc_fail", column=1, name="COLOC failed", label="COLOC failed",
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
                dict(key="hyprcoloc_pass", column=5, name="HyPrColoc: pQTL+GWAS+eQTL colocalise", label="HyPrColoc passed",
                     proteins=hyprcoloc_pass_set, color=STATUS_GOOD, dropout=False),
                dict(key="hyprcoloc_fail", column=5, name="HyPrColoc: no shared 3-trait signal", label="HyPrColoc failed",
                     proteins=hyprcoloc_fail_set, color=STATUS_CRITICAL, dropout=True),
            ]

            sankey_edges = [
                ("mr", "coloc_pass"), ("mr", "coloc_fail"),
                ("coloc_pass", "finngen_pass"), ("coloc_pass", "finngen_fail"),
                ("finngen_pass", "ukb_pass"), ("finngen_pass", "ukb_fail"),
                ("ukb_pass", "smr_both"), ("ukb_pass", "smr_bulk"),
                ("ukb_pass", "smr_sc"), ("ukb_pass", "smr_none"),
                ("smr_both", "hyprcoloc_pass"), ("smr_both", "hyprcoloc_fail"),
                ("smr_sc", "hyprcoloc_pass"), ("smr_sc", "hyprcoloc_fail"),
                ("smr_bulk", "hyprcoloc_pass"), ("smr_bulk", "hyprcoloc_fail"),
            ]

            # an empty stage costs a label and a slot but carries no information,
            # so it is left out of the diagram (it stays in the selector below)
            drawn = [group for group in sankey_groups if group["proteins"]]
            node_index = {group["key"]: position for position, group in enumerate(drawn)}

            # Plotly keeps node labels inside the plot area, flipping the last
            # column's to the left of its nodes rather than letting them run into
            # the margin - so the right margin stays thin and the columns are
            # spaced apart instead. The gap before the last column is the widest
            # because that is the one place a right-hand label (the SMR stages')
            # meets a left-flipped one (the HyPrColoc stages').
            column_x = [0.02, 0.19, 0.35, 0.51, 0.67, 0.99]

            node_values = [len(group["proteins"]) for group in drawn]
            node_labels = [f"{group['label']} ({value})" for group, value in zip(drawn, node_values)]
            node_colors = [group["color"] for group in drawn]
            node_x = [column_x[group["column"]] for group in drawn]
            node_hover = [
                f"<b>{group['name']} — {value} target(s)</b><br><br>"
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
                textfont=dict(size=13, color="#2b2b33"),
                node=dict(
                    label=node_labels,
                    color=node_colors,
                    customdata=node_hover,
                    hovertemplate="%{customdata}<extra></extra>",
                    x=node_x,
                    y=node_y,
                    pad=18,
                    thickness=15,
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
            "Two views of the same list, switched below - both already need cis-MR + COLOC + "
            "FinnGen/UKB safety + SMR/HEIDI. They differ only in whether HyPrColoc is required "
            "too, and therefore in **which SNP** the Top SNP / allele / beta columns are "
            "reported at."
        )

        final_targets_view = st.segmented_control(
            "Which target list to show",
            ["Final targets (passed HyPrColoc)", "All targets (before HyPrColoc)"],
            default="Final targets (passed HyPrColoc)",
            selection_mode="single",
            key="final_targets_view_selector"
        )

        # segmented_control returns None if the user clicks the selected option again to
        # deselect it - fall back to the stricter "passed HyPrColoc" view rather than an
        # undefined state
        show_hyprcoloc_targets = final_targets_view != "All targets (before HyPrColoc)"

        base_targets = smr_display.copy()
        identity_cols = [
            col for col in ["topsnp", "topsnp_chr", "topsnp_bp", "a1", "a2", "b_gwas", "b_eqtl"]
            if col in base_targets.columns
        ]
        base_targets = base_targets.drop(columns=identity_cols)

        # smr_display carries every protein ever run through SMR, including ones that
        # never reached (or failed) cis-MR/COLOC/FinnGen/UKB safety/SMR-HEIDI - restrict
        # to the same smr_eligible_set the Sankey above already gated on, so both views
        # here actually match the "passed cis-MR + COLOC + safety + SMR/HEIDI" claim in
        # their captions instead of silently including ineligible proteins
        if "protein" in base_targets.columns:
            base_targets = base_targets[base_targets["protein"].astype(str).isin(smr_eligible_set)]

        if show_hyprcoloc_targets:
            st.success(
                "**Final targets** - on top of cis-MR, COLOC, FinnGen/UKB safety and SMR/HEIDI, "
                f"these also passed HyPrColoc (posterior probability ≥ {hyprcoloc_pp_threshold:.2f}). "
                "**Top SNP** is HyPrColoc's own *candidate SNP* - the single variant it found shared "
                "across the pQTL, GWAS and eQTL signals - with alleles and betas aligned to the AD "
                "risk allele."
            )

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
                "**All targets, before HyPrColoc** - passed cis-MR, COLOC, FinnGen/UKB safety "
                "and SMR/HEIDI, whether or not HyPrColoc has since confirmed a single shared "
                "causal variant. **Each row picks its own Top SNP**: HyPrColoc's candidate SNP "
                f"for rows where that target x cell-type/tissue actually passed HyPrColoc (posterior "
                f"probability ≥ {hyprcoloc_pp_threshold:.2f}), falling back to the target's top "
                "cis-pQTL SNP otherwise - **SNP source** says which applies to each row. Alleles "
                "and betas are aligned to the AD risk allele either way (p-values are only ever "
                "floored to 1e-300 when reported as exactly 0)."
            )

            # same passing rows as the "Final targets" branch above - a row lands here only if
            # this target x cell-type/tissue actually cleared the HyPrColoc threshold, so its
            # candidate SNP is trustworthy to use as this row's Top SNP
            hypr_pass = select_hyprcoloc_candidate_rows(hyprcoloc_display, hyprcoloc_pp_threshold)
            hypr_cols = available_cols(
                hypr_pass,
                ["protein", "cell_type", "data_type", "candidate_snp", "a1", "a2",
                 "gwas_beta", "gwas_p", "pqtl_beta", "pqtl_p", "eqtl_beta", "eqtl_p", "posterior_prob"]
            )
            hypr_pass = hypr_pass[hypr_cols].rename(columns={
                "candidate_snp": "hypr_topsnp",
                "a1": "hypr_a1",
                "a2": "hypr_a2",
                "gwas_beta": "hypr_b_gwas",
                "gwas_p": "hypr_gwas_p",
                "pqtl_beta": "hypr_pqtl_beta",
                "pqtl_p": "hypr_pqtl_p",
                "eqtl_beta": "b_eqtl",
                "eqtl_p": "eqtl_p",
                "posterior_prob": "hyprcoloc_posterior_prob"
            })

            if not hypr_pass.empty and "hypr_a1" not in hypr_pass.columns:
                st.warning(
                    "This pQTL dataset's HyPrColoc results don't carry the candidate SNP's "
                    "aligned alleles/betas yet - rerun the HyPrColoc pipeline step "
                    "(bin/hyprcoloc_targets.py) to populate them. Every row will fall back to "
                    "the top cis-pQTL SNP until then."
                )
                hypr_pass = hypr_pass.iloc[0:0]

            # target-level (not cell-type-specific) fallback for rows that didn't pass HyPrColoc
            pqtl_cols = available_cols(target_info, ["protein", "snp", "a1", "a2", "gwas_beta", "gwas_p", "pqtl_beta", "pqtl_p"])
            pqtl_fallback = (
                target_info[pqtl_cols].rename(columns={
                    "snp": "pqtl_topsnp",
                    "a1": "pqtl_a1",
                    "a2": "pqtl_a2",
                    "gwas_beta": "pqtl_b_gwas",
                    "gwas_p": "pqtl_gwas_p",
                    "pqtl_beta": "pqtl_pqtl_beta",
                    "pqtl_p": "pqtl_pqtl_p"
                })
                if not target_info.empty and "protein" in pqtl_cols
                else pd.DataFrame()
            )

            merge_cols = [col for col in ["protein", "cell_type", "data_type"] if col in base_targets.columns and col in hypr_pass.columns]
            final_targets = base_targets.merge(hypr_pass, on=merge_cols, how="left") if merge_cols else base_targets.copy()

            if not pqtl_fallback.empty:
                final_targets = final_targets.merge(pqtl_fallback, on="protein", how="left")

            hypr_topsnp_col = (
                final_targets["hypr_topsnp"] if "hypr_topsnp" in final_targets.columns
                else pd.Series(np.nan, index=final_targets.index)
            )
            passed_hyprcoloc_row = hypr_topsnp_col.notna()

            final_targets["snp_source"] = np.where(
                passed_hyprcoloc_row, "HyPrColoc candidate SNP", "cis-pQTL top SNP"
            )

            # collapse the hypr_*/pqtl_* pairs into the plain column names final_cols expects
            # below, picking the HyPrColoc value per row where it passed and the target-level
            # pQTL fallback everywhere else
            fallback_pairs = [
                ("topsnp", "hypr_topsnp", "pqtl_topsnp"),
                ("a1", "hypr_a1", "pqtl_a1"),
                ("a2", "hypr_a2", "pqtl_a2"),
                ("b_gwas", "hypr_b_gwas", "pqtl_b_gwas"),
                ("gwas_p", "hypr_gwas_p", "pqtl_gwas_p"),
                ("pqtl_beta", "hypr_pqtl_beta", "pqtl_pqtl_beta"),
                ("pqtl_p", "hypr_pqtl_p", "pqtl_pqtl_p"),
            ]

            for out_col, hypr_col, pqtl_col in fallback_pairs:
                hypr_series = final_targets[hypr_col] if hypr_col in final_targets.columns else pd.Series(np.nan, index=final_targets.index)
                pqtl_series = final_targets[pqtl_col] if pqtl_col in final_targets.columns else pd.Series(np.nan, index=final_targets.index)
                final_targets[out_col] = hypr_series.where(passed_hyprcoloc_row, pqtl_series)

            drop_cols = [col for _, hypr_col, pqtl_col in fallback_pairs for col in (hypr_col, pqtl_col) if col in final_targets.columns]
            final_targets = final_targets.drop(columns=drop_cols)

        if final_targets.empty:
            if show_hyprcoloc_targets:
                st.info(
                    "No targets currently pass every gate including HyPrColoc (posterior "
                    f"probability ≥ {hyprcoloc_pp_threshold:.2f}) at the selected thresholds."
                )
            else:
                st.info(
                    "No targets currently pass cis-MR, COLOC, FinnGen/UKB safety and SMR/HEIDI "
                    "at the selected thresholds."
                )
        else:
            final_targets = final_targets.sort_values(
                ["protein", "cell_type"],
                na_position="last"
            ) if "cell_type" in final_targets.columns else final_targets

            with st.container(border=True):
                metric1, metric2 = st.columns(2)
                metric1.metric("Unique targets", safe_nunique(final_targets, "protein"))
                metric2.metric("Target x cell/bulk rows", len(final_targets))

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


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--db_name", required=True, type=str)
    p.add_argument("--port_number", required=True, type=str)
    p.add_argument("--phenotype", required=True, type=str)
    p.add_argument("--pqtl_dataset", required=True, type=str)
    args = p.parse_args()
    create_streamlit_ammenities(args.db_name, args.port_number)
    dashboard(
        db_name=args.db_name,
        port_number=args.port_number,
        phenotype=args.phenotype,
        pqtl_dataset=args.pqtl_dataset
    )


if __name__ == "__main__":
    main()