#!/usr/bin/env python3
import argparse
from pathlib import Path
import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
import streamlit as st

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


def subset_phewas_outcome(df: pd.DataFrame, outcome: str):
    if df.empty:
        return df

    if "outcome_trait" in df.columns:
        return df[df["outcome_trait"] == outcome].copy()

    if "pheno_id" in df.columns:
        return df[df["pheno_id"] == outcome].copy()

    return df.copy()


def render_phewas_section(
    phewas_outcome: pd.DataFrame,
    mr_coloc_pass: pd.DataFrame,
    source_name: str,
    source_description: str,
    n_endpoints: int,
    pqtl_dataset: str,
    outcome: str,
    key_prefix: str
):
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

    # check which datasets have the required dashboard files
    dataset_result_files = {}
    available_datasets = []

    for dataset_id in dataset_names:
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

        finngen_phewas_file = project_dir / "results" / "PheWAS" / dataset_id / phenotype / f"{dataset_id}_{phenotype}_PheWAS.tsv"
        ukb_phewas_file = project_dir / "results" / "PheWAS_UKBB" / dataset_id / phenotype / f"{dataset_id}_{phenotype}_PheWAS.tsv"
        target_info_file = project_dir / "results" / "target_stats" / dataset_id / phenotype / f"{dataset_id}_{phenotype}_top_cis_hits.tsv"
        smr_file = project_dir / "results" / "SMR" / f"{dataset_id}_{phenotype}_final_multi_omics_targets.tsv"

        required_files = [
            mr_file,
            coloc_file
        ]

        if all(file is not None and file.exists() for file in required_files):
            available_datasets.append(dataset_id)
            dataset_result_files[dataset_id] = {
                "mr": mr_file,
                "coloc": coloc_file,
                "finngen_phewas": finngen_phewas_file,
                "ukb_phewas": ukb_phewas_file,
                "target_info": target_info_file,
                "smr": smr_file
            }

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
        dataset_col, dataset_info_col = st.columns([2, 1])

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

    # load local result files into PostgreSQL for the dashboard
    mr = load_required_tsv(mr_file, "cis-MR")
    coloc = load_required_tsv(coloc_file, "pQTL–GWAS COLOC")
    finngen_phewas = load_optional_tsv(finngen_phewas_file, "FinnGen PheWAS safety")
    ukb_phewas = load_optional_tsv(ukb_phewas_file, "UKB PheWAS safety")
    target_info = load_optional_tsv(target_info_file, "Harmonised target information")
    smr = load_optional_tsv(smr_file, "SMR (bulk/sc eQTL)")

    # standardise MR + pQTL COLOC columns before loading into PostgreSQL
    # avoids dataset-specific differences such as Wald_beta vs wald_beta
    mr = standardise_columns(mr)
    coloc = standardise_columns(coloc)

    if not target_info.empty:
        target_info = standardise_columns(target_info)

    if not smr.empty:
        smr = standardise_columns(smr)

    # make protein column consistent before loading into PostgreSQL
    if "protein_id" in mr.columns:
        mr = mr.rename(columns={"protein_id": "protein"})

    if "protein_id" in coloc.columns:
        coloc = coloc.rename(columns={"protein_id": "protein"})

    if not target_info.empty and "protein_id" in target_info.columns:
        target_info = target_info.rename(columns={"protein_id": "protein"})

    if not smr.empty and "protein_id" in smr.columns:
        smr = smr.rename(columns={"protein_id": "protein"})

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


    # SMR (bulk/sc eQTL) hits are not split by outcome_trait like mr/coloc are,
    # so they're carried through as their own table rather than an "_outcome" subset
    smr_display = smr.copy()

    # protein search
    if protein:
        mr_outcome = filter_protein(mr_outcome, protein)
        mr_pass = filter_protein(mr_pass, protein)
        coloc_pass = filter_protein(coloc_pass, protein)
        mr_coloc_pass = filter_protein(mr_coloc_pass, protein)
        finngen_phewas_outcome = filter_protein(finngen_phewas_outcome, protein)
        ukb_phewas_outcome = filter_protein(ukb_phewas_outcome, protein)
        smr_display = filter_protein(smr_display, protein)

    # round coloc posterior probs
    for col in coloc_numeric_cols:
        if col in coloc_pass.columns:
            coloc_pass[col] = coloc_pass[col].round(3)

        if col in mr_coloc_pass.columns:
            mr_coloc_pass[col] = mr_coloc_pass[col].round(3)

    # main staged target counts
    n_tested = safe_nunique(mr_outcome, "protein")
    n_mr = safe_nunique(mr_pass, "protein")
    n_mr_coloc = safe_nunique(mr_coloc_pass, "protein")
    n_finngen_phewas = safe_nunique(finngen_phewas_outcome, "protein")
    n_ukb_phewas = safe_nunique(ukb_phewas_outcome, "protein")

    tab1, tab2, tab3, tab4, tab5, tab6, tab7 = st.tabs([
        "Overview",
        "1. cis-MR",
        "2. pQTL–GWAS COLOC",
        "3. FinnGen PheWAS",
        "4. UKB PheWAS",
        "5. SMR (bulk/sc eQTL)",
        "6. Final Targets"
    ])

    with tab1:
        st.subheader("Target prioritisation")
        st.caption(
            "Targets move from cis-MR testing to MR support and then to shared pQTL–GWAS "
            "causal signal support through pairwise COLOC."
        )

        with st.container(border=True):
            metric1, metric2, metric3 = st.columns(3)
            metric1.metric("Proteins tested by cis-MR", n_tested)
            metric2.metric("cis-MR supported", n_mr, f"{retention(n_mr, n_tested):.1f}% of tested", delta_color="off")
            metric3.metric("cis-MR + pQTL COLOC", n_mr_coloc, f"{retention(n_mr_coloc, n_mr):.1f}% retained", delta_color="off")

        st.divider()

        funnel_df = pd.DataFrame({
            "stage": [
                "Proteins tested by cis-MR",
                "cis-MR supported",
                "cis-MR + pQTL COLOC"
            ],
            "n_targets": [
                n_tested,
                n_mr,
                n_mr_coloc
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

        if not mr_coloc_pass.empty:
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

            prioritised_cols = available_cols(mr_coloc_pass, prioritised_cols)

            if "pp_h4_abf" in mr_coloc_pass.columns:
                mr_coloc_pass = mr_coloc_pass.sort_values(
                    ["pp_h4_abf", "mr_fdr_q"],
                    ascending=[False, True],
                    na_position="last"
                )

            overview_table = mr_coloc_pass[prioritised_cols].copy()

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

            st.success(
                f"{n_mr_coloc} unique target(s) passed the selected cis-MR and pairwise COLOC thresholds."
            )

            st.caption(
                "All SNP effects are harmonised to the outcome GWAS risk allele. "
                "A positive pQTL beta means the risk allele increases protein abundance, "
                "whereas a negative pQTL beta means the risk allele decreases protein abundance."
            )

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

        else:
            st.info("No proteins currently pass both the selected cis-MR and pQTL COLOC thresholds.")

        if finngen_phewas_available or ukb_phewas_available:
            st.caption(
                f"FinnGen PheWAS results are available for {n_finngen_phewas} unique target(s); "
                f"UKB PheWAS results are available for {n_ukb_phewas} unique target(s)."
            )

    with tab2:
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
            key_prefix="finngen"
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
            key_prefix="ukb"
        )

    with tab6:
        st.subheader("SMR (bulk / single-cell eQTL)")
        st.caption(
            "Targets shown here already passed cis-MR + pQTL–GWAS COLOC, and additionally "
            "passed SMR (FDR-corrected) + HEIDI in the configured bulk and/or single-cell eQTL "
            "dataset(s). Alleles are aligned to the AD risk allele, same convention as the "
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
                "- **cis-MR** — this flow starts from the proteins that already passed; the "
                "much larger screening drop-off across every tested protein is the funnel in "
                "the Overview tab.\n"
                "- **pQTL–GWAS COLOC** — passes on the posterior-probability threshold set in "
                "the sidebar.\n"
                "- **FinnGen / UKB PheWAS safety** — fails *only* when a target has a "
                "Bonferroni-significant association with `beta_mr >= 0`, i.e. the same allele "
                "that raises AD risk also raises the safety-relevant phenotype. No significant "
                "hit, a protective (negative) significant hit, or no PheWAS coverage at all "
                "each count as passing.\n"
                f"- **SMR support** — requires SMR FDR (`q_SMR`) < {smr_fdr_threshold:.2f} and "
                f"HEIDI p-value > {heidi_p_threshold:.2f}, split by whether that support came "
                "from bulk/tissue eQTL data, single-cell eQTL data, or both."
            )

        n_sankey_mr_pass = safe_nunique(mr_pass, "protein")

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
            ]

            sankey_edges = [
                ("mr", "coloc_pass"), ("mr", "coloc_fail"),
                ("coloc_pass", "finngen_pass"), ("coloc_pass", "finngen_fail"),
                ("finngen_pass", "ukb_pass"), ("finngen_pass", "ukb_fail"),
                ("ukb_pass", "smr_both"), ("ukb_pass", "smr_bulk"),
                ("ukb_pass", "smr_sc"), ("ukb_pass", "smr_none"),
            ]

            # an empty stage costs a label and a slot but carries no information,
            # so it is left out of the diagram (it stays in the selector below)
            drawn = [group for group in sankey_groups if group["proteins"]]
            node_index = {group["key"]: position for position, group in enumerate(drawn)}

            # Plotly keeps node labels inside the plot area, flipping the last
            # column's to the left of its nodes rather than letting them run into
            # the margin - so the right margin stays thin and the columns are
            # spaced apart instead. The gap before the last column is the widest
            # because that is the one place a right-hand label (the UKB stage's)
            # meets a left-flipped one (the SMR stages').
            column_x = [0.02, 0.24, 0.44, 0.64, 0.99]

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
                    value=[node_values[target] for _, target in edges],
                    # ribbons take the colour of where they land, and drop-out
                    # ribbons sit fainter so the eye follows the surviving lane
                    color=[
                        hex_to_rgba(node_colors[target], 0.25 if drawn[target]["dropout"] else 0.4)
                        for _, target in edges
                    ],
                    customdata=[
                        f"<b>{drawn[source]['name']} → {drawn[target]['name']}</b><br>"
                        f"{node_values[target]} target(s)<br><br>"
                        f"{format_protein_list_html(drawn[target]['proteins'])}"
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

        if smr_display.empty:
            st.info("No SMR results are available for this pQTL dataset yet.")
        else:
            final_targets = smr_display.sort_values(
                ["protein", "cell_type"],
                na_position="last"
            ) if "cell_type" in smr_display.columns else smr_display.copy()

            # pQTL beta is per-protein (not per-cell-type) - bring it in from the
            # harmonised top cis-hit table, already aligned to the same risk allele
            # convention as the SMR file's own A1 (both scripts align A1 the same way)
            if not target_info.empty and "protein" in target_info.columns and "pqtl_beta" in target_info.columns:
                final_targets = final_targets.merge(
                    target_info[["protein", "pqtl_beta"]],
                    on="protein",
                    how="left"
                )

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
                "topsnp",
                "topsnp_chr",
                "topsnp_bp",
                "a1",
                "a2",
                "b_gwas",
                "pqtl_beta",
                "b_eqtl",
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
                "topsnp": "Top SNP",
                "topsnp_chr": "Chr",
                "topsnp_bp": "Position (bp)",
                "a1": "Risk allele",
                "a2": "Other allele",
                "b_gwas": "GWAS beta",
                "pqtl_beta": "pQTL beta",
                "b_eqtl": "eQTL beta",
                "b_smr": "SMR beta",
                "p_smr": "SMR p-value",
                "q_smr": "SMR FDR",
                "p_heidi": "HEIDI p-value"
            }

            final_table = final_table.rename(columns=final_column_names)

            st.caption("Betas are all aligned to the outcome (AD) risk allele shown in **Risk allele**.")

            st.dataframe(
                final_table,
                hide_index=True,
                width="stretch",
                column_config={
                    "Chr": st.column_config.NumberColumn(format="%d"),
                    "Position (bp)": st.column_config.NumberColumn(format="%d"),
                    "GWAS beta": st.column_config.NumberColumn(format="%.4f"),
                    "pQTL beta": st.column_config.NumberColumn(format="%.4f"),
                    "eQTL beta": st.column_config.NumberColumn(format="%.4f"),
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