#!/usr/bin/env python3
import argparse
import subprocess
from pathlib import Path
import numpy as np
import pandas as pd
import plotly.express as px
import streamlit as st
from sqlalchemy import inspect


# KEY CHANGES DOWN THE LINE WITH MORE PQTL DATASETS 
# -> CHANGE THE DASHBOARD FUNCT TO ADD EQTL AND PQL ARGS
# biomarker meta analysis: https://pmc.ncbi.nlm.nih.gov/articles/instance/12136742/pdf/nihpp-rs6597595v1.pdf

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


def dashboard(db_name: str, port_number: str, phenotype: str, pqtl_dataset: str):
    mr_table = "cis_mr_results"
    coloc_table = "coloc_results"
    smr_table = "single_cell_smr_results"
    eqtl_coloc_table = "eqtl_coloc_results"
    moloc_table = "moloc_results"
    overview_table = "multi_omics_overview"
    snp_table = "multi_omics_snp_evidence"
    phewas_table = "phewas_safety"

    # main aesthetics
    st.set_page_config(page_title=f"{db_name}", layout="wide")
    conn = st.connection("postgresql", type="sql", url=f"postgresql://localhost:{port_number}/{db_name}")

    # corresponding multi-omics result files
    project_dir = Path(__file__).resolve().parent.parent
    smr_file = project_dir / f"results/SMR/SingleBrain/{phenotype}/{pqtl_dataset}_{phenotype}_promising_targets_SMR.tsv"
    eqtl_coloc_file = project_dir / f"results/eQTL_coloc/{pqtl_dataset}/SingleBrain/{phenotype}/{pqtl_dataset}_{phenotype}_SingleBrain_all_eqtl_coloc.tsv"
    moloc_file = project_dir / f"results/QTL_moloc/{pqtl_dataset}/SingleBrain/{phenotype}/{pqtl_dataset}_{phenotype}_SingleBrain_moloc_summary.tsv"
    overview_file = project_dir / f"results/SMR/SingleBrain/{phenotype}/{pqtl_dataset}_{phenotype}_multi_omics_overview.tsv"
    snp_file = project_dir / f"results/SMR/SingleBrain/{phenotype}/{pqtl_dataset}_{phenotype}_multi_omics_snp_evidence.tsv"
    phewas_file = project_dir / f"results/PheWAS/{pqtl_dataset}/{phenotype}/{pqtl_dataset}_{phenotype}_PheWAS.tsv"
    # "{pqtl_dataset}_{pheno_id}_PheWAS.tsv"

    # create new dashboard tables if they do not exist yet
    tables = inspect(conn.engine).get_table_names()

    if not smr_file.exists():
        st.error(f"SMR result file not found: {smr_file}")
        st.stop()

    smr = pd.read_csv(smr_file, sep="\t")
    smr.to_sql(smr_table, conn.engine, if_exists="replace", index=False)
    st.write(
    f"[TRACKING] Loaded {len(smr)} rows into {smr_table} "
    f"for {pqtl_dataset}"
    )
    st.write("[TRACKING] Cell types found:",sorted(smr["cell_type"].dropna().unique()))

    
    if not eqtl_coloc_file.exists():
        st.error(f"eQTL COLOC result file not found: {eqtl_coloc_file}")
        st.stop()
    eqtl_coloc = pd.read_csv(eqtl_coloc_file, sep="\t")
    eqtl_coloc.to_sql(eqtl_coloc_table, conn.engine, if_exists="replace", index=False)

    if not moloc_file.exists():
        st.error(f"MOLOC result file not found: {moloc_file}")
        st.stop()

    moloc = pd.read_csv(moloc_file, sep="\t")
    moloc.to_sql(moloc_table, conn.engine, if_exists="replace", index=False)

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

    # phewas gist
    if not phewas_file.exists():
        st.warning(f"PheWAS safety result file not found: {phewas_file}")
        phewas_available = False
    else:
        phewas = pd.read_csv(phewas_file, sep="\t")
        if phewas.empty:
            st.warning(f"PheWAS safety result file is empty: {phewas_file}")
            phewas_available = False
        else:
            phewas.to_sql(phewas_table, conn.engine, if_exists="replace", index=False)
            st.write(f"[TRACKING] Loaded {len(phewas)} rows into {phewas_table}")
            phewas_available = True

    # load MR + COLOC + multi-omics results
    mr = conn.query(f"SELECT * FROM {mr_table};", ttl=0)
    coloc = conn.query(f"SELECT * FROM {coloc_table};", ttl=0)
    smr = conn.query(f"SELECT * FROM {smr_table};", ttl=0)
    eqtl_coloc = conn.query(f"SELECT * FROM {eqtl_coloc_table};", ttl=0)
    moloc = conn.query(f"SELECT * FROM {moloc_table};", ttl=0)
    overview = conn.query(f"SELECT * FROM {overview_table};", ttl=0)
    snp = conn.query(f"SELECT * FROM {snp_table};", ttl=0)
    snp = snp[snp["phenotype"] == phenotype].copy()
    
    # check whether phewas is avaulable cuz snp_targets.shape[0] might == 0
    if phewas_available:
        phewas = conn.query(f"SELECT * FROM {phewas_table};", ttl=0)
    else:
        phewas = pd.DataFrame()

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

    # make PheWAS columns consistent
    if not phewas.empty:
        phewas = phewas.rename(columns={
            "PROTEIN": "protein",
            "protein_id": "protein",
            "PHENO_ID": "pheno_id",
            "OUTCOME_TRAIT": "outcome_trait",
            "PHENOCODE": "phenocode",
            "PHENOSTRING": "phenostring",
            "CATEGORY": "category",
            "SNP": "snp",
            "RSID": "rsid",
            "METHOD": "method",
            "N_INSTRUMENTS": "n_instruments",
            "BETA_MR": "beta_mr",
            "SE_MR": "se_mr",
            "P_MR": "p_mr",
            "P_BONFERRONI": "p_bonferroni",
            "BONFERRONI_SIGNIFICANT": "bonferroni_significant"
        })

        # A1/A2 already come from the original outcome GWAS
        # do not overwrite them with FinnGen ALT/REF
        if "A1" in phewas.columns:
            phewas["A1"] = phewas["A1"].astype(str).str.upper()

        if "A2" in phewas.columns:
            phewas["A2"] = phewas["A2"].astype(str).str.upper()

        for col in [
            "n_instruments",
            "beta_mr",
            "se_mr",
            "p_mr",
            "p_bonferroni"
        ]:
            if col in phewas.columns:
                phewas[col] = pd.to_numeric(phewas[col], errors="coerce")

        if "bonferroni_significant" in phewas.columns:
            phewas["bonferroni_significant"] = (
                phewas["bonferroni_significant"]
                .astype(str)
                .str.lower()
                .isin(["true", "1", "yes"])
            )
        elif "p_bonferroni" in phewas.columns:
            phewas["bonferroni_significant"] = phewas["p_bonferroni"].fillna(np.inf) <= 0.05

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
    pp4 = st.sidebar.slider("pQTL–GWAS COLOC PP.H4 threshold", 0.0, 1.0, 0.70, 0.01)
    smr_q = st.sidebar.slider("Single-cell SMR FDR threshold", 0.0, 1.0, 0.05, 0.01)
    heidi_p = st.sidebar.slider("Minimum HEIDI p-value", 0.0, 1.0, 0.01, 0.01)
    eqtl_pp4 = st.sidebar.slider("GWAS–eQTL COLOC PP.H4 threshold", 0.0, 1.0, 0.70, 0.01)
    moloc_ppa = st.sidebar.slider("Three-trait MOLOC PPA threshold", 0.0, 1.0, 0.70, 0.01)
    cell_types = sorted(smr["cell_type"].dropna().unique())
    selected_cell_types = st.sidebar.multiselect("Cell types", cell_types, default=cell_types)
    protein = st.sidebar.text_input("Protein search")

    dataset_names = {"ukb_ppp": "UKB-PPP", "decode": "deCODE", "wu_csf": "WU-CSF"}
    dataset_name = dataset_names.get(pqtl_dataset, pqtl_dataset)
    st.title(f"{db_name}: {dataset_name} → {outcome}")
    
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

    if not phewas.empty:
        phewas_outcome = phewas.copy()

        if "outcome_trait" in phewas_outcome.columns:
            phewas_outcome = phewas_outcome[phewas_outcome["outcome_trait"] == outcome].copy()
        elif "pheno_id" in phewas_outcome.columns:
            phewas_outcome = phewas_outcome[phewas_outcome["pheno_id"] == outcome].copy()
    else:
        phewas_outcome = pd.DataFrame()

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

        if not phewas_outcome.empty and "protein" in phewas_outcome.columns:
            phewas_outcome = phewas_outcome[phewas_outcome["protein"].str.contains(protein, case=False, na=False)]

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

    tab1, tab2, tab3, tab4, tab5, tab6 = st.tabs([
        "Overview",
        "cis-MR results",
        "pQTL COLOC",
        "Single-cell evidence",
        "Final targets",
        "PheWAS safety"
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

            # so (top pQTL SNP)
            # SNP (rsID)
            # CHR
            # POS (GRCh38)
            # protein_id
            # pheno_id
            # pqtl_dataset
            # A1 (risk allele) on GWAS
            # A2
            # Top pQTL SNP
            # GWAS beta
            # GWAS P
            # pQTL beta
            # pQTL P
            # sc-eQTL beta
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


    with tab6:
        st.subheader("FinnGen PheWAS safety assessment")

        if phewas_outcome.empty:
            st.info("No local PheWAS safety results are available for this outcome.")

        elif "protein" not in phewas_outcome.columns:
            st.error("The PheWAS result file does not contain a protein column.")

        else:
            phewas_targets = sorted(phewas_outcome["protein"].dropna().astype(str).unique())

            if len(phewas_targets) == 0:
                st.info("No proteins were found in the PheWAS safety table.")

            else:
                default_phewas_target = 0
                final_target_names = final_multi_omics["protein"].dropna().astype(str).unique().tolist()

                for target in final_target_names:
                    if target in phewas_targets:
                        default_phewas_target = phewas_targets.index(target)
                        break

                selected_phewas_target = st.selectbox(
                    "Select target for PheWAS",
                    phewas_targets,
                    index=default_phewas_target,
                    key="selected_phewas_target"
                )

                target_phewas = phewas_outcome[
                    phewas_outcome["protein"].astype(str) == selected_phewas_target
                ].copy()

                p_col = None
                beta_col = None
                bonferroni_col = None

                for col in ["p_mr"]:
                    if col in target_phewas.columns:
                        p_col = col
                        break

                for col in ["beta_mr"]:
                    if col in target_phewas.columns:
                        beta_col = col
                        break

                for col in ["p_bonferroni"]:
                    if col in target_phewas.columns:
                        bonferroni_col = col
                        break

                if p_col is None or beta_col is None:
                    st.error(
                        "The PheWAS result file needs the MR effect column "
                        "(beta_mr) and the MR p-value column "
                        "(p_mr)."
                    )

                else:
                    target_phewas = target_phewas[
                        target_phewas[p_col].notna() &
                        target_phewas[beta_col].notna() &
                        (target_phewas[p_col] > 0)
                    ].copy()

                    if target_phewas.empty:
                        st.info(f"No valid PheWAS associations were found for {selected_phewas_target}.")

                    else:
                        target_phewas["minus_log10_p"] = -np.log10(target_phewas[p_col])

                        if bonferroni_col is not None:
                            target_phewas["bonferroni_significant"] = target_phewas[bonferroni_col].fillna(np.inf) <= 0.05
                        elif "bonferroni_significant" not in target_phewas.columns:
                            target_phewas["bonferroni_significant"] = False

                        phenotype_col = "phenostring" if "phenostring" in target_phewas.columns else "phenocode"
                        category_col = "category" if "category" in target_phewas.columns else None

                        n_phenotypes = target_phewas[phenotype_col].nunique()
                        n_nominal = int((target_phewas[p_col] < 0.05).sum())
                        n_bonferroni = int(target_phewas["bonferroni_significant"].sum())

                        metric1, metric2, metric3 = st.columns(3)
                        metric1.metric("FinnGen phenotypes tested", int(n_phenotypes))
                        metric2.metric("Nominal associations", n_nominal)
                        metric3.metric("Bonferroni-significant associations", n_bonferroni)
                        
                        st.caption(
                            "PheWAS MR estimates show the effect of genetically predicted protein levels "
                            "on each FinnGen phenotype (ICD-10 coded). Wald ratio is used for targets with "
                            "one COJO-selected cis-pQTL instrument and IVW is used for targets with more than one."
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
                            "title": f"FinnGen PheWAS profile: {selected_phewas_target}",
                            "height": 600
                        }

                        if "phenocode" in target_phewas.columns:
                            plot_kwargs["hover_data"]["phenocode"] = True

                        if category_col is not None:
                            plot_kwargs["color"] = category_col
                            plot_kwargs["labels"][category_col] = "FinnGen category"

                        phewas_fig = px.scatter(**plot_kwargs)
                        phewas_fig.add_hline(y=-np.log10(0.05 / 2511), line_dash="dash", line_color="grey")
                        phewas_fig.add_vline(x=0, line_dash="dash", line_color="grey")
                        st.plotly_chart(phewas_fig, use_container_width=True)

                        st.subheader("Bonferroni-significant PheWAS associations")

                        top_phewas = target_phewas[target_phewas["bonferroni_significant"]].copy()

                        if bonferroni_col is not None:
                            top_phewas = top_phewas.sort_values(bonferroni_col, ascending=True)
                        else:
                            top_phewas = top_phewas.sort_values(p_col, ascending=True)

                        top_phewas = top_phewas.sort_values(beta_col, ascending=True)

                        if top_phewas.empty:
                            st.info(
                                f"No FinnGen phenotype associations survive Bonferroni correction across "
                                f"2,511 ICD endpoints for {selected_phewas_target}."
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
                                "height": max(450, 45 * len(top_phewas))
                            }

                            if "phenocode" in top_phewas.columns:
                                top_plot_kwargs["hover_data"]["phenocode"] = True

                            if category_col is not None:
                                top_plot_kwargs["color"] = category_col
                                top_plot_kwargs["labels"][category_col] = "FinnGen category"

                            top_phewas_fig = px.scatter(**top_plot_kwargs)
                            top_phewas_fig.add_vline(x=0, line_dash="dash", line_color="grey")
                            st.plotly_chart(top_phewas_fig, use_container_width=True)

                        phewas_cols = [
                            "protein",
                            "method",
                            "n_instruments",
                            "rsid",
                            "A1",
                            "A2",
                            "phenocode",
                            "phenostring",
                            "category",
                            "beta_mr",
                            "se_mr",
                            "p_mr",
                            "p_bonferroni",
                            "bonferroni_significant"
                        ]

                        phewas_cols = [
                            col for col in phewas_cols
                            if col is not None and col in target_phewas.columns
                        ]

                        significant_phewas = target_phewas[
                            target_phewas["bonferroni_significant"]
                        ].copy()

                        if bonferroni_col is not None:
                            significant_phewas = significant_phewas.sort_values(bonferroni_col, ascending=True)
                        else:
                            significant_phewas = significant_phewas.sort_values(p_col, ascending=True)

                        if significant_phewas.empty:
                            st.success("No FinnGen phenotype associations survive Bonferroni correction across 2,511 ICD endpoints for this target.")
                        else:
                            st.dataframe(
                                significant_phewas[phewas_cols],
                                use_container_width=True,
                                hide_index=True
                            )

                        with st.expander("View all PheWAS associations"):
                            remaining_cols = [col for col in target_phewas.columns if col not in phewas_cols]
                            st.dataframe(
                                target_phewas[phewas_cols + remaining_cols].sort_values(p_col, ascending=True),
                                use_container_width=True,
                                hide_index=True
                            )

                        st.download_button(
                            label=f"Download {selected_phewas_target} PheWAS results",
                            data=target_phewas.to_csv(index=False, sep="\t"),
                            file_name=f"{selected_phewas_target}_{outcome}_FinnGen_PheWAS.tsv",
                            mime="text/tab-separated-values"
                        )


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--db_name", required=True, type=str)
    p.add_argument("--port_number", required=True, type=str)
    p.add_argument("--phenotype", required=True, type=str)
    p.add_argument("--pqtl_dataset", required=True, type=str)
    args = p.parse_args()
    create_streamlit_ammenities(args.db_name, args.port_number)
    dashboard(db_name=args.db_name, port_number=args.port_number, phenotype=args.phenotype, pqtl_dataset=args.pqtl_dataset)


if __name__ == "__main__":
    main()