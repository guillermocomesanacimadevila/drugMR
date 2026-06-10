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
    table_id = "cis_mr_results"

    # main aesthetics
    st.set_page_config(page_title = f"{db_name}", layout = "wide")
    conn = st.connection("postgresql", type="sql")
    df = conn.query(f"SELECT * FROM {table_id};", ttl="10m")
    st.title(f"{db_name} cis-MR dashboard ({phenotype})")
    outcome = st.sidebar.selectbox("Outcome", sorted(df["outcome_trait"].dropna().unique()))
    fdr = st.sidebar.slider("IVW FDR threshold", 0.0, 1.0, 1.0)
    protein = st.sidebar.text_input("Protein search")
    x = df[df["outcome_trait"] == outcome].copy()

    if "ivw_fdr_q" in x.columns:
        x = x[x["ivw_fdr_q"] <= fdr]

    if protein:
        x = x[x["protein"].str.contains(protein, case=False, na=False)]

    st.metric("Proteins shown", len(x))
    st.dataframe(x, use_container_width=True)

    plot_df = df[df["outcome_trait"] == outcome].copy()
    plot_df["minus_log10_ivw_pval"] = -np.log10(plot_df["ivw_pval"])
    fig = px.scatter(
        plot_df,
        x="ivw_beta",
        y="minus_log10_ivw_pval",
        hover_name="protein",
        color="ivw_fdr_q",
        title="IVW MR volcano plot"
    )

    st.plotly_chart(fig, use_container_width=True)

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