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
    st.set_page_config(page_title = f"{db_name}", layout = "wide")
    conn = st.connection("postgresql", type="sql")
    mr = conn.query(f"SELECT * FROM {mr_table};", ttl="10m")
    coloc = conn.query(f"SELECT * FROM {coloc_table};", ttl="10m")

    if "protein_id" in coloc.columns:
        coloc = coloc.rename(columns={"protein_id": "protein"})

    st.title(f"{db_name}: UKBB-PPP → {phenotype}")

    outcome = st.sidebar.selectbox("Outcome", sorted(mr["outcome_trait"].dropna().unique()))
    fdr = st.sidebar.slider("IVW FDR threshold", 0.0, 1.0, 1.0)
    pp4 = st.sidebar.slider("COLOC PP.H4 threshold", 0.0, 1.0, 0.0)
    protein = st.sidebar.text_input("Protein search")

    mr_x = mr[mr["outcome_trait"] == outcome].copy()
    coloc_x = coloc[coloc["outcome_trait"] == outcome].copy()

    if "ivw_fdr_q" in mr_x.columns:
        mr_x = mr_x[mr_x["ivw_fdr_q"] <= fdr]

    if "pp_h4_abf" in coloc_x.columns:
        coloc_x = coloc_x[coloc_x["pp_h4_abf"].fillna(0) >= pp4]

    if protein:
        mr_x = mr_x[mr_x["protein"].str.contains(protein, case=False, na=False)]
        coloc_x = coloc_x[coloc_x["protein"].str.contains(protein, case=False, na=False)]

    for col in ["pp_h0_abf", "pp_h1_abf", "pp_h2_abf", "pp_h3_abf", "pp_h4_abf"]:
        if col in coloc_x.columns:
            coloc_x[col] = coloc_x[col].round(3)

    tab1, tab2 = st.tabs(["cis-MR results", "COLOC results"])

    with tab1:
        st.metric("MR proteins shown", len(mr_x))
        st.dataframe(mr_x, use_container_width=True)

        plot_df = mr_x.copy()
        plot_df["minus_log10_ivw_pval"] = -np.log10(plot_df["ivw_pval"])
        plot_df["significant"] = plot_df["ivw_fdr_q"] < 0.05

        fig = px.scatter(
            plot_df,
            x="ivw_beta",
            y="minus_log10_ivw_pval",
            hover_name="protein",
            color="significant",
            title="IVW MR volcano plot",
            height=600
        )

        # add 0.05 thresh line 
        fig.add_hline(y=-np.log10(0.05), line_dash = "dash", line_color = "grey") # at FDR_q = 0.05   
        st.plotly_chart(fig, use_container_width=True)

    with tab2:
        st.metric("COLOC proteins shown", len(coloc_x))
        st.dataframe(coloc_x, use_container_width=True)

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