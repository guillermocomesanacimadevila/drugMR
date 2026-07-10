#!/usr/bin/env python3
import argparse
import subprocess
import numpy as np
import plotly.express as px
import streamlit as st


def create_streamlit_ammenities(db_name: str, port_number: str):
    cmd = f"""
set -euo pipefail 
mkdir -p .streamlit
cat > .streamlit/secrets.toml <<EOF
[connections.postgresql]
dialect = "postgresql"
host = "localhost"
port = "{port_number}"
database = "{db_name}"
username = ""
password = ""
EOF
    """

    # run in terminal to create streamlit ammenities 
    subprocess.run(cmd, shell=True, check=True, executable="/bin/bash")


def dashboard(db_name: str, phenotype: str):
    mr_table = "cis_mr_results"
    coloc_table = "coloc_results"

    # main aesthetics
    st.set_page_config(page_title=f"{db_name}", layout="wide")
    conn = st.connection("postgresql", type="sql")

    # load MR + COLOC results
    mr = conn.query(f"SELECT * FROM {mr_table};", ttl="10m")
    coloc = conn.query(f"SELECT * FROM {coloc_table};", ttl="10m")

    # MR ammenities
    # if 1 instrument -> use Wald
    # otherwise -> use IVW
    mr["mr_method"] = np.where(mr["n_instruments"] == 1, "Wald ratio", "IVW")
    mr["mr_beta"] = np.where(mr["n_instruments"] == 1, mr["wald_beta"], mr["ivw_beta"])
    mr["mr_se"] = np.where(mr["n_instruments"] == 1, mr["wald_se"], mr["ivw_se"])
    mr["mr_pval"] = np.where(mr["n_instruments"] == 1, mr["wald_pval"], mr["ivw_pval"])
    mr["mr_fdr_q"] = np.where(mr["n_instruments"] == 1, mr["wald_fdr_q"], mr["ivw_fdr_q"])

    # make protein column consistent between MR and COLOC
    if "protein_id" in coloc.columns:
        coloc = coloc.rename(columns={"protein_id": "protein"})

    st.title(f"{db_name}: UKBB-PPP → {phenotype}")

    # sidebar filters
    outcome = st.sidebar.selectbox("Outcome", sorted(mr["outcome_trait"].dropna().unique()))
    fdr = st.sidebar.slider("MR FDR threshold", 0.0, 1.0, 1.0, 0.01)
    q_pval = st.sidebar.slider("Minimum Cochran Q p-value", 0.0, 1.0, 0.0, 0.01)
    egger_pval = st.sidebar.slider("Minimum Egger intercept p-value", 0.0, 1.0, 0.0, 0.01)
    pp4 = st.sidebar.slider("COLOC PP.H4 threshold", 0.0, 1.0, 0.0, 0.01)
    protein = st.sidebar.text_input("Protein search")

    # subset to selected outcome
    mr_x = mr[mr["outcome_trait"] == outcome].copy()
    coloc_x = coloc[coloc["outcome_trait"] == outcome].copy()

    # Wald FDR for 1 instrument
    # IVW FDR for multi-instrument proteins
    if "mr_fdr_q" in mr_x.columns:
        mr_x = mr_x[mr_x["mr_fdr_q"].fillna(np.inf) <= fdr]

    # apply Cochran Q only to IVW proteins
    # Wald proteins have no Cochran Q so keep them
    if "q_pval" in mr_x.columns:
        mr_x = mr_x[
            (
                (mr_x["mr_method"] == "IVW") &
                (mr_x["q_pval"].fillna(-np.inf) >= q_pval)
            )
            |
            (mr_x["mr_method"] == "Wald ratio")
        ]

    # apply Egger intercept only to IVW proteins
    # Wald proteins have no Egger intercept so keep them
    if "egger_intercept_pval" in mr_x.columns:
        mr_x = mr_x[
            (
                (mr_x["mr_method"] == "IVW") &
                (mr_x["egger_intercept_pval"].fillna(-np.inf) >= egger_pval)
            )
            |
            (mr_x["mr_method"] == "Wald ratio")
        ]

    # COLOC filter
    if "pp_h4_abf" in coloc_x.columns:
        coloc_x = coloc_x[coloc_x["pp_h4_abf"].fillna(0) >= pp4]

    # protein search
    if protein:
        mr_x = mr_x[mr_x["protein"].str.contains(protein, case=False, na=False)]
        coloc_x = coloc_x[coloc_x["protein"].str.contains(protein, case=False, na=False)]

    # proteins which pass both MR + COLOC thresholds
    prioritised = mr_x.merge(coloc_x, on="protein", how="inner", suffixes=("_mr", "_coloc"))

    # remove duplicate proteins if there are somehow repeated rows
    prioritised = prioritised.drop_duplicates(subset=["protein"])

    # round coloc posterior probs
    for col in ["pp_h0_abf", "pp_h1_abf", "pp_h2_abf", "pp_h3_abf", "pp_h4_abf"]:
        if col in coloc_x.columns:
            coloc_x[col] = coloc_x[col].round(3)

        if col in prioritised.columns:
            prioritised[col] = prioritised[col].round(3)

    # prioritised target metrics at the top
    st.subheader("Target prioritisation")

    metric1, metric2, metric3 = st.columns(3)

    with metric1:
        st.metric("MR proteins passing", len(mr_x))

    with metric2:
        st.metric("COLOC proteins passing", len(coloc_x))

    with metric3:
        st.metric("Targets passing MR + COLOC", len(prioritised))

    tab1, tab2, tab3 = st.tabs(["cis-MR results", "COLOC results", "Prioritised targets"])

    with tab1:
        n_ivw = (mr_x["mr_method"] == "IVW").sum()
        n_wald = (mr_x["mr_method"] == "Wald ratio").sum()

        col1, col2, col3 = st.columns(3)

        with col1:
            st.metric("MR proteins shown", len(mr_x))

        with col2:
            st.metric("IVW proteins", int(n_ivw))

        with col3:
            st.metric("Wald proteins", int(n_wald))

        # put main MR cols first
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

        display_cols = [col for col in display_cols if col in mr_x.columns]
        remaining_cols = [col for col in mr_x.columns if col not in display_cols]

        st.dataframe(mr_x[display_cols + remaining_cols], use_container_width=True)

        # primary MR volcano plot
        # uses IVW for multi-instrument proteins
        # uses Wald for single-instrument proteins
        plot_df = mr_x.copy()

        plot_df = plot_df[
            plot_df["mr_pval"].notna() &
            plot_df["mr_beta"].notna() &
            (plot_df["mr_pval"] > 0)
        ].copy()

        if not plot_df.empty:
            plot_df["minus_log10_mr_pval"] = -np.log10(plot_df["mr_pval"])
            plot_df["significant"] = plot_df["mr_fdr_q"] < 0.05

            fig = px.scatter(
                plot_df,
                x="mr_beta",
                y="minus_log10_mr_pval",
                hover_name="protein",
                color="significant",
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
                height=600
            )

            # raw p-value reference line
            # point significance itself is still defined using FDR < 0.05
            fig.add_hline(y=-np.log10(0.05), line_dash="dash", line_color="grey")
            fig.add_vline(x=0, line_dash="dash", line_color="grey")

            st.plotly_chart(fig, use_container_width=True)

        else:
            st.info("No MR results remain after applying the selected filters.")

    with tab2:
        st.metric("COLOC proteins shown", len(coloc_x))
        st.dataframe(coloc_x, use_container_width=True)

    with tab3:
        st.metric("Targets passing all selected thresholds", len(prioritised))

        if not prioritised.empty:
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
                "pp_h0_abf",
                "pp_h1_abf",
                "pp_h2_abf",
                "pp_h3_abf",
                "pp_h4_abf"
            ]

            prioritised_cols = [col for col in prioritised_cols if col in prioritised.columns]
            remaining_prioritised_cols = [col for col in prioritised.columns if col not in prioritised_cols]

            prioritised = prioritised.sort_values(
                by=["pp_h4_abf", "mr_fdr_q"],
                ascending=[False, True],
                na_position="last"
            )

            st.dataframe(
                prioritised[prioritised_cols + remaining_prioritised_cols],
                use_container_width=True
            )

            st.download_button(
                label="Download prioritised targets",
                data=prioritised.to_csv(index=False, sep="\t"),
                file_name=f"{outcome}_prioritised_targets.tsv",
                mime="text/tab-separated-values"
            )

            st.write("Proteins passing the current MR + COLOC thresholds:")

            for target in prioritised["protein"].tolist():
                st.code(target)

        else:
            st.info("No proteins currently pass both the selected MR and COLOC thresholds.")


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--db_name", required=True, type=str)
    p.add_argument("--port_number", required=True, type=str)
    p.add_argument("--phenotype", required=True, type=str)
    args = p.parse_args()
    create_streamlit_ammenities(args.db_name, args.port_number)
    dashboard(args.db_name, args.phenotype)


if __name__ == "__main__":
    main()