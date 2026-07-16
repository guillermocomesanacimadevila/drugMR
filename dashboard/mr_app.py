#!/usr/bin/env python3
import argparse
import subprocess
from pathlib import Path
import numpy as np
import pandas as pd
import plotly.express as px
import streamlit as st
from sqlalchemy import inspect


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


def retention(current: int, previous: int):
    return 0.0 if previous == 0 else 100 * current / previous


def dashboard(db_name: str, phenotype: str):
    mr_table = "cis_mr_results"
    coloc_table = "coloc_results"
    smr_table = "single_cell_smr_results"
    eqtl_coloc_table = "eqtl_coloc_results"
    moloc_table = "moloc_results"
    overview_table = "multi_omics_overview"
    snp_table = "multi_omics_snp_evidence"

    # main aesthetics
    st.set_page_config(page_title=f"{db_name}", layout="wide")
    conn = st.connection("postgresql", type="sql")

    # corresponding multi-omics result files
    project_dir = Path(__file__).resolve().parent.parent
    smr_file = project_dir / f"results/SMR/SingleBrain/{phenotype}/ukb_ppp_{phenotype}_promising_targets_SMR.tsv"
    eqtl_coloc_file = project_dir / f"results/eQTL_coloc/ukb_ppp/SingleBrain/{phenotype}/ukb_ppp_{phenotype}_SingleBrain_all_eqtl_coloc.tsv"
    moloc_file = project_dir / f"results/QTL_moloc/ukb_ppp/SingleBrain/{phenotype}/ukb_ppp_{phenotype}_SingleBrain_moloc_summary.tsv"
    overview_file = project_dir / f"results/SMR/SingleBrain/{phenotype}/ukb_ppp_{phenotype}_multi_omics_overview.tsv"
    snp_file = project_dir / f"results/SMR/SingleBrain/{phenotype}/ukb_ppp_{phenotype}_multi_omics_snp_evidence.tsv"

    # create new dashboard tables if they do not exist yet
    tables = inspect(conn.engine).get_table_names()

    if smr_table not in tables:
        if not smr_file.exists():
            st.error(f"SMR result file not found: {smr_file}")
            st.stop()

        smr = pd.read_csv(smr_file, sep="\t")
        smr.to_sql(smr_table, conn.engine, if_exists="replace", index=False)
        st.write(f"[TRACKING] Loaded {len(smr)} rows into {smr_table}")

    if eqtl_coloc_table not in tables:
        if not eqtl_coloc_file.exists():
            st.error(f"eQTL COLOC result file not found: {eqtl_coloc_file}")
            st.stop()

        eqtl_coloc = pd.read_csv(eqtl_coloc_file, sep="\t")
        eqtl_coloc.to_sql(eqtl_coloc_table, conn.engine, if_exists="replace", index=False)
        st.write(f"[TRACKING] Loaded {len(eqtl_coloc)} rows into {eqtl_coloc_table}")

    if moloc_table not in tables:
        if not moloc_file.exists():
            st.error(f"MOLOC result file not found: {moloc_file}")
            st.stop()

        moloc = pd.read_csv(moloc_file, sep="\t")
        moloc.to_sql(moloc_table, conn.engine, if_exists="replace", index=False)
        st.write(f"[TRACKING] Loaded {len(moloc)} rows into {moloc_table}")

    if not overview_file.exists():
        st.error(f"Overview result file not found: {overview_file}")
        st.stop()

    overview = pd.read_csv(overview_file, sep="\t")
    overview.to_sql(overview_table, conn.engine, if_exists="replace", index=False)
    st.write(f"[TRACKING] Loaded {len(overview)} rows into {overview_table}")

    if not snp_file.exists():
        st.error(f"SNP evidence file not found: {snp_file}")
        st.stop()

    snp = pd.read_csv(snp_file, sep="\t")
    snp.to_sql(snp_table, conn.engine, if_exists="replace", index=False)
    st.write(f"[TRACKING] Loaded {len(snp)} rows into {snp_table}")

    # load MR + COLOC + multi-omics results
    mr = conn.query(f"SELECT * FROM {mr_table};", ttl=0)
    coloc = conn.query(f"SELECT * FROM {coloc_table};", ttl=0)
    smr = conn.query(f"SELECT * FROM {smr_table};", ttl=0)
    eqtl_coloc = conn.query(f"SELECT * FROM {eqtl_coloc_table};", ttl=0)
    moloc = conn.query(f"SELECT * FROM {moloc_table};", ttl=0)
    overview = conn.query(f"SELECT * FROM {overview_table};", ttl=0)
    snp = conn.query(f"SELECT * FROM {snp_table};", ttl=0)
    snp = snp[snp["phenotype"] == phenotype].copy()

    # standardise cols of QTL stuff
    # make SMR columns easier to use
    smr = smr.rename(columns={
        "Gene": "gene",
        "topSNP": "top_snp",
        "b_SMR": "b_smr",
        "se_SMR": "se_smr",
        "p_SMR": "p_smr",
        "p_HEIDI": "p_heidi",
        "q_SMR": "q_smr"
    })

    # make eQTL COLOC columns consistent
    eqtl_coloc = eqtl_coloc.rename(columns={
        "protein_id": "protein",
        "outcome_trait": "phenotype",
        "PP.H0.abf": "eqtl_pp_h0_abf",
        "PP.H1.abf": "eqtl_pp_h1_abf",
        "PP.H2.abf": "eqtl_pp_h2_abf",
        "PP.H3.abf": "eqtl_pp_h3_abf",
        "PP.H4.abf": "eqtl_pp_h4_abf"
    })

    # make MOLOC columns consistent
    moloc = moloc.rename(columns={"model": "moloc_model", "PPA": "moloc_ppa"})

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

    # available outcomes and default CLI phenotype
    outcomes = sorted(mr["outcome_trait"].dropna().unique())
    default_outcome = outcomes.index(phenotype) if phenotype in outcomes else 0

    # sidebar filters
    outcome = st.sidebar.selectbox("Outcome", outcomes, index=default_outcome)
    fdr = st.sidebar.slider("MR FDR threshold", 0.0, 1.0, 0.05, 0.01)
    q_pval = st.sidebar.slider("Minimum Cochran Q p-value", 0.0, 1.0, 0.05, 0.01)
    # egger_pval = st.sidebar.slider("Minimum Egger intercept p-value", 0.0, 1.0, 0.05, 0.01)
    pp4 = st.sidebar.slider("pQTL–GWAS COLOC PP.H4 threshold", 0.0, 1.0, 0.75, 0.01)
    smr_q = st.sidebar.slider("Single-cell SMR FDR threshold", 0.0, 1.0, 0.05, 0.01)
    heidi_p = st.sidebar.slider("Minimum HEIDI p-value", 0.0, 1.0, 0.01, 0.01)
    eqtl_pp4 = st.sidebar.slider("GWAS–eQTL COLOC PP.H4 threshold", 0.0, 1.0, 0.75, 0.01)
    moloc_ppa = st.sidebar.slider("Three-trait MOLOC PPA threshold", 0.0, 1.0, 0.70, 0.01)
    cell_types = sorted(smr["cell_type"].dropna().unique())
    selected_cell_types = st.sidebar.multiselect("Cell types", cell_types, default=cell_types)
    protein = st.sidebar.text_input("Protein search")

    st.title(f"{db_name}: UKBB-PPP → {outcome}")
    st.caption(
        # f"MR FDR ≤ {fdr:.2f} | Q p ≥ {q_pval:.2f} | Egger p ≥ {egger_pval:.2f} | "
        f"MR FDR ≤ {fdr:.2f} | Q p ≥ {q_pval:.2f} | "
        f"pQTL–GWAS PP.H4 ≥ {pp4:.2f} | SMR FDR ≤ {smr_q:.2f} | HEIDI p ≥ {heidi_p:.2f} | "
        f"GWAS–eQTL PP.H4 ≥ {eqtl_pp4:.2f} | MOLOC abc PPA ≥ {moloc_ppa:.2f}"
    )

    # subset everything to selected outcome
    mr_outcome = mr[mr["outcome_trait"] == outcome].copy()
    coloc_outcome = coloc[coloc["outcome_trait"] == outcome].copy()
    smr_outcome = smr[smr["phenotype"] == outcome].copy()
    eqtl_coloc_outcome = eqtl_coloc[eqtl_coloc["phenotype"] == outcome].copy()
    moloc_outcome = moloc[moloc["phenotype"] == outcome].copy()
    overview_outcome = overview[overview["phenotype"] == outcome].copy()
    snp_outcome = snp[snp["phenotype"] == outcome].copy()

    if "cell_type" in overview_outcome.columns:
        overview_outcome = overview_outcome[overview_outcome["cell_type"].isin(selected_cell_types)]

    if "cell_type" in snp_outcome.columns:
        snp_outcome = snp_outcome[snp_outcome["cell_type"].isin(selected_cell_types)]

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

    if "pp_h4_abf" in coloc_pass.columns:
        coloc_pass = coloc_pass[coloc_pass["pp_h4_abf"].fillna(0) >= pp4]

    # proteins which pass both MR + COLOC thresholds
    mr_coloc_pass = mr_pass.merge(coloc_pass, on="protein", how="inner", suffixes=("_mr", "_pqtl_coloc"))
    mr_coloc_pass = mr_coloc_pass.drop_duplicates(subset=["protein"])

    # STAGE 3
    # single-cell SMR + HEIDI for proteins already supported by MR + COLOC
    smr_pass = smr_outcome[
        (smr_outcome["q_smr"].fillna(np.inf) <= smr_q) &
        (smr_outcome["p_heidi"].fillna(-np.inf) >= heidi_p) &
        (smr_outcome["cell_type"].isin(selected_cell_types))
    ].copy()

    smr_stage = smr_pass.merge(mr_coloc_pass[["protein"]], on="protein", how="inner")
    smr_stage = smr_stage.drop_duplicates(subset=["protein", "cell_type"])

    # STAGE 4
    # GWAS - eQTL COLOC for SMR supported target x cell-type pairs
    eqtl_coloc_pass = eqtl_coloc_outcome[
        (eqtl_coloc_outcome["eqtl_pp_h4_abf"].fillna(0) >= eqtl_pp4) &
        (eqtl_coloc_outcome["cell_type"].isin(selected_cell_types))
    ].copy()

    smr_eqtl_coloc_stage = smr_stage.merge(
        eqtl_coloc_pass,
        on=["protein", "cell_type"],
        how="inner",
        suffixes=("_smr", "_eqtl_coloc")
    )

    smr_eqtl_coloc_stage = smr_eqtl_coloc_stage.drop_duplicates(subset=["protein", "cell_type"])

    # STAGE 5
    # pQTL - eQTL - GWAS MOLOC
    moloc_pass = moloc_outcome[
        (moloc_outcome["moloc_model"].astype(str).str.lower() == "abc") &
        (moloc_outcome["moloc_ppa"].fillna(0) >= moloc_ppa) &
        (moloc_outcome["cell_type"].isin(selected_cell_types))
    ].copy()

    final_multi_omics = smr_eqtl_coloc_stage.merge(
        moloc_pass,
        on=["protein", "cell_type"],
        how="inner",
        suffixes=("", "_moloc")
    )

    final_multi_omics = final_multi_omics.drop_duplicates(subset=["protein", "cell_type"])

    # protein search
    if protein:
        mr_outcome = mr_outcome[mr_outcome["protein"].str.contains(protein, case=False, na=False)]
        mr_pass = mr_pass[mr_pass["protein"].str.contains(protein, case=False, na=False)]
        coloc_pass = coloc_pass[coloc_pass["protein"].str.contains(protein, case=False, na=False)]
        mr_coloc_pass = mr_coloc_pass[mr_coloc_pass["protein"].str.contains(protein, case=False, na=False)]
        smr_stage = smr_stage[smr_stage["protein"].str.contains(protein, case=False, na=False)]
        smr_eqtl_coloc_stage = smr_eqtl_coloc_stage[smr_eqtl_coloc_stage["protein"].str.contains(protein, case=False, na=False)]
        final_multi_omics = final_multi_omics[final_multi_omics["protein"].str.contains(protein, case=False, na=False)]
        overview_outcome = overview_outcome[overview_outcome["protein"].str.contains(protein, case=False, na=False)]
        snp_outcome = snp_outcome[snp_outcome["protein"].str.contains(protein, case=False, na=False)]

    # round coloc posterior probs
    for col in ["pp_h0_abf", "pp_h1_abf", "pp_h2_abf", "pp_h3_abf", "pp_h4_abf"]:
        if col in coloc_pass.columns:
            coloc_pass[col] = coloc_pass[col].round(3)

        if col in mr_coloc_pass.columns:
            mr_coloc_pass[col] = mr_coloc_pass[col].round(3)

    for col in ["eqtl_pp_h0_abf", "eqtl_pp_h1_abf", "eqtl_pp_h2_abf", "eqtl_pp_h3_abf", "eqtl_pp_h4_abf"]:
        if col in eqtl_coloc_pass.columns:
            eqtl_coloc_pass[col] = eqtl_coloc_pass[col].round(3)

        if col in smr_eqtl_coloc_stage.columns:
            smr_eqtl_coloc_stage[col] = smr_eqtl_coloc_stage[col].round(3)

        if col in final_multi_omics.columns:
            final_multi_omics[col] = final_multi_omics[col].round(3)

    if "moloc_ppa" in final_multi_omics.columns:
        final_multi_omics["moloc_ppa"] = final_multi_omics["moloc_ppa"].round(3)

    # main staged target counts
    n_tested = mr_outcome["protein"].nunique()
    n_mr = mr_pass["protein"].nunique()
    n_mr_coloc = mr_coloc_pass["protein"].nunique()
    n_smr = smr_stage["protein"].nunique()
    n_eqtl_coloc = smr_eqtl_coloc_stage["protein"].nunique()
    n_final = final_multi_omics["protein"].nunique()

    tab1, tab2, tab3, tab4, tab5 = st.tabs([
        "Overview",
        "cis-MR results",
        "pQTL COLOC",
        "Single-cell evidence",
        "Final targets"
    ])

    with tab1:
        st.subheader("Target prioritisation")

        metric1, metric2, metric3 = st.columns(3)
        metric1.metric("Proteins tested by cis-MR", n_tested)
        metric2.metric("cis-MR supported", n_mr, f"{retention(n_mr, n_tested):.1f}% of tested", delta_color="off")
        metric3.metric("cis-MR + pQTL COLOC", n_mr_coloc, f"{retention(n_mr_coloc, n_mr):.1f}% retained", delta_color="off")

        metric4, metric5, metric6 = st.columns(3)
        metric4.metric("+ cell-type SMR/HEIDI", n_smr, f"{retention(n_smr, n_mr_coloc):.1f}% retained", delta_color="off")
        metric5.metric("+ GWAS–eQTL COLOC", n_eqtl_coloc, f"{retention(n_eqtl_coloc, n_smr):.1f}% retained", delta_color="off")
        metric6.metric("Final three-trait targets", n_final, f"{retention(n_final, n_eqtl_coloc):.1f}% retained", delta_color="off")

        funnel_df = pd.DataFrame({
            "stage": [
                "Proteins tested by cis-MR",
                "cis-MR supported",
                "cis-MR + pQTL COLOC",
                "Cell-type SMR + HEIDI",
                "GWAS–eQTL COLOC",
                "Three-trait MOLOC"
            ],
            "n_targets": [n_tested, n_mr, n_mr_coloc, n_smr, n_eqtl_coloc, n_final]
        })

        funnel_fig = px.bar(
            funnel_df,
            x="n_targets",
            y="stage",
            orientation="h",
            text="n_targets",
            title="Progressive target prioritisation",
            labels={"n_targets": "Number of unique proteins", "stage": ""},
            height=480
        )

        funnel_fig.update_yaxes(categoryorder="array", categoryarray=funnel_df["stage"][::-1])
        funnel_fig.update_traces(textposition="outside")
        st.plotly_chart(funnel_fig, use_container_width=True)

        if not final_multi_omics.empty:
            # we need to add some stuff here 
            overview_cols = [
                "protein",
                "gene_smr",
                "cell_type",
                "top_snp_smr",
                "q_smr",
                "p_heidi",
                "eqtl_pp_h4_abf",
                "moloc_model",
                "moloc_ppa"
            ]

            # we need another table with 
            # A1 (or risk allele) in GWAS
            # aligned to risk allele for pQTL and sc-eQTL 
            # beta for each 
            # pvalue on that specific thingy (of both top SNPs) (top pQTL SNP and top SMR SNP)

            # so (top pQTL SNP)
            # SNP (rsID)
            # CHR
            # POS (GRCh38)
            # protein_id
            # pheno_id
            # pqtl_dataset
            # A1 (risk allele) on GWAS
            # A2
            # Top pQTL SNP
            # GWAS beta
            # GWAS P
            # pQTL beta
            # pQTL P
            # sc-eQTL beta
            # sc-eQTL P
            # ---- NOW THE SAME FOR THE TOP SMR SNP 


            overview_cols = [col for col in overview_cols if col in final_multi_omics.columns]
            st.success(f"{n_final} unique target(s) passed all multi-omics evidence chain.")
            overview_cols = ["protein", "Gene", "cell_type", "pqtl_coloc_top_snp", "topSNP", "Wald_beta", "Wald_pval", "Wald_FDR_q", "pqtl_pp_h4", "q_SMR", "p_HEIDI", "eqtl_pp_h4", "moloc_model", "moloc_ppa"]
            overview_cols = [col for col in overview_cols if col in overview_outcome.columns]
            st.dataframe(overview_outcome[overview_cols], use_container_width=True, hide_index=True)
            st.subheader("SNP-level evidence")
            snp_cols = [
                "protein",
                "gene",
                "cell_type",
                "snp_type",
                "SNP",
                "CHR",
                "POS",
                "A1",
                "A2",
                "GWAS_beta",
                "GWAS_SE",
                "GWAS_P",
                "pQTL_beta",
                "pQTL_SE",
                "pQTL_P",
                "sc_eQTL_beta",
                "sc_eQTL_SE",
                "sc_eQTL_P"
            ]

            snp_cols = [col for col in snp_cols if col in snp_outcome.columns]
            st.dataframe(snp_outcome[snp_cols], use_container_width=True, hide_index=True)

        else:
            st.info("No targets currently pass the complete multi-omics evidence chain.")

    with tab2:
        show_all_mr = st.checkbox("Show all tested cis-MR proteins", value=False)
        mr_display = mr_outcome if show_all_mr else mr_pass
        n_ivw = (mr_display["mr_method"] == "IVW").sum()
        n_wald = (mr_display["mr_method"] == "Wald ratio").sum()

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

        display_cols = [col for col in display_cols if col in mr_display.columns]
        remaining_cols = [col for col in mr_display.columns if col not in display_cols]
        st.dataframe(mr_display[display_cols + remaining_cols], use_container_width=True, hide_index=True)

        # primary MR volcano plot
        plot_df = mr_display[
            mr_display["mr_pval"].notna() &
            mr_display["mr_beta"].notna() &
            (mr_display["mr_pval"] > 0)
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

            fig.add_hline(y=-np.log10(0.05), line_dash="dash", line_color="grey")
            fig.add_vline(x=0, line_dash="dash", line_color="grey")
            st.plotly_chart(fig, use_container_width=True)

        else:
            st.info("No MR results remain after applying the selected filters.")

    with tab3:
        st.subheader("pQTL–GWAS colocalisation")

        col1, col2 = st.columns(2)
        col1.metric("All COLOC-supported proteins", coloc_pass["protein"].nunique())
        col2.metric("Proteins supported by cis-MR + COLOC", mr_coloc_pass["protein"].nunique())

        coloc_tab1, coloc_tab2 = st.tabs(["All COLOC-supported proteins", "cis-MR + COLOC targets"])

        with coloc_tab1:
            st.dataframe(coloc_pass, use_container_width=True, hide_index=True)

        with coloc_tab2:
            if not mr_coloc_pass.empty:
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

                prioritised_cols = [col for col in prioritised_cols if col in mr_coloc_pass.columns]
                remaining_cols = [col for col in mr_coloc_pass.columns if col not in prioritised_cols]

                if "pp_h4_abf" in mr_coloc_pass.columns:
                    mr_coloc_pass = mr_coloc_pass.sort_values(["pp_h4_abf", "mr_fdr_q"], ascending=[False, True], na_position="last")

                st.dataframe(mr_coloc_pass[prioritised_cols + remaining_cols], use_container_width=True, hide_index=True)

                st.download_button(
                    label="Download cis-MR + COLOC targets",
                    data=mr_coloc_pass.to_csv(index=False, sep="\t"),
                    file_name=f"{outcome}_prioritised_targets.tsv",
                    mime="text/tab-separated-values"
                )

            else:
                st.info("No proteins currently pass both the selected cis-MR and pQTL COLOC thresholds.")

    with tab4:
        st.subheader("Single-cell multi-omics evidence")

        smr_tab, eqtl_tab, moloc_tab = st.tabs(["SMR + HEIDI", "GWAS–eQTL COLOC", "Three-trait MOLOC"])

        with smr_tab:
            col1, col2 = st.columns(2)
            col1.metric("Target × cell-type SMR hits", len(smr_stage))
            col2.metric("Unique proteins", smr_stage["protein"].nunique())

            if not smr_stage.empty:
                smr_cols = [
                    "protein",
                    "gene",
                    "cell_type",
                    "probeID",
                    "top_snp",
                    "b_smr",
                    "se_smr",
                    "p_smr",
                    "q_smr",
                    "p_heidi",
                    "nsnp_HEIDI",
                    "pqtl_dataset",
                    "eqtl_dataset"
                ]

                smr_cols = [col for col in smr_cols if col in smr_stage.columns]
                remaining_cols = [col for col in smr_stage.columns if col not in smr_cols]
                st.dataframe(smr_stage[smr_cols + remaining_cols], use_container_width=True, hide_index=True)

            else:
                st.info("No target × cell-type pairs pass the selected SMR and HEIDI thresholds.")

        with eqtl_tab:
            col1, col2 = st.columns(2)
            col1.metric("Target × cell-type eQTL COLOC hits", len(smr_eqtl_coloc_stage))
            col2.metric("Unique proteins", smr_eqtl_coloc_stage["protein"].nunique())

            if not smr_eqtl_coloc_stage.empty:
                eqtl_cols = [
                    "protein",
                    "gene_smr",
                    "cell_type",
                    "top_snp_smr",
                    "q_smr",
                    "p_heidi",
                    "top_snp_eqtl_coloc",
                    "eqtl_pp_h0_abf",
                    "eqtl_pp_h1_abf",
                    "eqtl_pp_h2_abf",
                    "eqtl_pp_h3_abf",
                    "eqtl_pp_h4_abf",
                    "n_eqtl_snps",
                    "n_gwas_snps"
                ]

                eqtl_cols = [col for col in eqtl_cols if col in smr_eqtl_coloc_stage.columns]
                remaining_cols = [col for col in smr_eqtl_coloc_stage.columns if col not in eqtl_cols]
                st.dataframe(smr_eqtl_coloc_stage[eqtl_cols + remaining_cols], use_container_width=True, hide_index=True)

            else:
                st.info("No SMR-supported target × cell-type pairs pass GWAS–eQTL COLOC.")

        with moloc_tab:
            col1, col2 = st.columns(2)
            col1.metric("Target × cell-type three-trait signals", len(final_multi_omics))
            col2.metric("Unique proteins", final_multi_omics["protein"].nunique())

            if not final_multi_omics.empty:
                moloc_cols = [
                    "protein",
                    "gene_smr",
                    "cell_type",
                    "top_snp_smr",
                    "q_smr",
                    "p_heidi",
                    "eqtl_pp_h4_abf",
                    "moloc_model",
                    "moloc_ppa",
                    "nsnps_moloc"
                ]

                moloc_cols = [col for col in moloc_cols if col in final_multi_omics.columns]
                remaining_cols = [col for col in final_multi_omics.columns if col not in moloc_cols]
                st.dataframe(final_multi_omics[moloc_cols + remaining_cols], use_container_width=True, hide_index=True)

            else:
                st.info("No target × cell-type pairs pass the selected three-trait MOLOC threshold.")

    with tab5:
        st.subheader("Final multi-omics drug targets")
        st.metric("Final unique targets", final_multi_omics["protein"].nunique())

        if not final_multi_omics.empty:
            selected_target = st.selectbox("Select final target", final_multi_omics["protein"].dropna().drop_duplicates().tolist())
            target_data = final_multi_omics[final_multi_omics["protein"] == selected_target].copy()
            selected_cell = st.selectbox("Cell type", target_data["cell_type"].dropna().drop_duplicates().tolist())
            target_row = target_data[target_data["cell_type"] == selected_cell].iloc[0]
            gene_col = "gene_smr" if "gene_smr" in target_row.index else "gene"

            if gene_col in target_row.index:
                st.markdown(f"### {target_row[gene_col]} `{selected_target}` in `{selected_cell}`")

            col1, col2, col3, col4 = st.columns(4)
            col1.metric("SMR FDR q", f"{target_row['q_smr']:.2e}")
            col2.metric("HEIDI p", f"{target_row['p_heidi']:.3f}")
            col3.metric("GWAS–eQTL PP.H4", f"{target_row['eqtl_pp_h4_abf']:.3f}")
            col4.metric("MOLOC abc PPA", f"{target_row['moloc_ppa']:.3f}")

            st.success(
                f"{selected_target} in {selected_cell} is supported by cis-pQTL MR, pQTL–GWAS colocalisation, "
                "cell-type SMR with a non-significant HEIDI test, GWAS–eQTL colocalisation and three-trait MOLOC."
            )

            final_cols = [
                "protein",
                "gene_smr",
                "cell_type",
                "probeID",
                "top_snp_smr",
                "b_smr",
                "se_smr",
                "p_smr",
                "q_smr",
                "p_heidi",
                "nsnp_HEIDI",
                "eqtl_pp_h4_abf",
                "moloc_model",
                "moloc_ppa",
                "nsnps_moloc"
            ]

            final_cols = [col for col in final_cols if col in final_multi_omics.columns]
            st.dataframe(final_multi_omics[final_cols], use_container_width=True, hide_index=True)

            with st.expander("View complete merged multi-omics evidence"):
                st.dataframe(final_multi_omics, use_container_width=True, hide_index=True)

            st.download_button(
                label="Download final multi-omics targets",
                data=final_multi_omics.to_csv(index=False, sep="\t"),
                file_name=f"{outcome}_final_multi_omics_targets.tsv",
                mime="text/tab-separated-values"
            )

            st.write("Final targets passing the complete evidence chain:")

            for target in final_multi_omics["protein"].drop_duplicates().tolist():
                cells = final_multi_omics[final_multi_omics["protein"] == target]["cell_type"].drop_duplicates().tolist()
                st.code(f"{target}: {', '.join(cells)}")

        else:
            st.info(
                "No targets currently pass cis-MR, pQTL COLOC, single-cell SMR/HEIDI, "
                "GWAS–eQTL COLOC and three-trait MOLOC using the selected thresholds."
            )


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
