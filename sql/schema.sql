/*
This is a collection of tables for the dashboard (see MIRO ERD)
- each table contains a specific entries respective to each run
- runs table -> keeps INSERT ... INTO runs -> so != fixed
- all other tables -> INSERT directly from TSVs

The multi-omics hits have to be compiled from a SQL query
stored within a PY script in bin/ -> inherited in dashboard

--

CREATE VIEW multi_omics_targets AS
SELECT
    mr.run_id, mr.protein, mr.pqtl_dataset,
    mr.primary_beta AS mr_beta, mr.primary_fdr_q,
    c.pp_h4_abf, c.coloc_pass,
    smr.eqtl_dataset, smr.cell_type,
    smr.b_gwas, smr.b_eqtl, smr.b_smr, smr.p_smr, smr.p_heidi,
    hc.posterior_prob AS hyprcoloc_pp, hc.candidate_snp
FROM cis_mr_results mr
JOIN coloc_results c ON c.run_id = mr.run_id AND c.protein = mr.protein
JOIN smr_results smr ON smr.run_id = mr.run_id AND smr.protein = mr.protein
LEFT JOIN hyprcoloc_results hc
    ON hc.run_id = mr.run_id AND hc.protein = mr.protein AND hc.eqtl_dataset = smr.eqtl_dataset;

*/


-- Table 1
CREATE TABLE runs (
    run_id VARCHAR NOT NULL,
    pheno_id VARCHAR NOT NULL,
    pqtl_dataset VARCHAR NOT NULL,
    created_at TIMESTAMPTZ,
    git_sha7 VARCHAR NOT NULL,
    mode VARCHAR NOT NULL,
    image_name VARCHAR NOT NULL,
    overwrite BOOLEAN,
    PRIMARY KEY (run_id)
);

CREATE INDEX idx_runs_pheno_dataset ON runs (pheno_id, pqtl_dataset);

-- cis-MR results
CREATE TABLE cis_mr_results (
    run_id VARCHAR NOT NULL REFERENCES runs(run_id),
    protein VARCHAR NOT NULL,
    pqtl_dataset VARCHAR NOT NULL,
    outcome_trait VARCHAR NOT NULL,
    n_instruments INTEGER,
    methods_run VARCHAR,
    primary_method VARCHAR,
    primary_beta DOUBLE PRECISION,
    primary_se DOUBLE PRECISION,
    primary_pval DOUBLE PRECISION,
    primary_fdr_q DOUBLE PRECISION,
    mde_or DOUBLE PRECISION,
    ivw_beta DOUBLE PRECISION,
    ivw_se DOUBLE PRECISION,
    ivw_pval DOUBLE PRECISION,
    ivw_fdr_q DOUBLE PRECISION,
    egger_beta DOUBLE PRECISION,
    egger_se DOUBLE PRECISION,
    egger_pval DOUBLE PRECISION,
    egger_fdr_q DOUBLE PRECISION,
    egger_intercept DOUBLE PRECISION,
    egger_intercept_pval DOUBLE PRECISION,
    wme_beta DOUBLE PRECISION,
    wme_se DOUBLE PRECISION,
    wme_pval DOUBLE PRECISION,
    wme_fdr_q DOUBLE PRECISION,
    wald_beta DOUBLE PRECISION,
    wald_se DOUBLE PRECISION,
    wald_pval DOUBLE PRECISION,
    wald_fdr_q DOUBLE PRECISION,
    q DOUBLE PRECISION,
    q_df INTEGER,
    q_pval DOUBLE PRECISION

);

CREATE INDEX idx_cis_mr_run ON cis_mr_results (run_id); -- to index cis_mr results with the run_id!!!!!

-- Table 3: pairwise Bayesian COLOC (pQTL-GWAS)
CREATE TABLE coloc_results (
    run_id VARCHAR NOT NULL REFERENCES runs(run_id),
    protein VARCHAR NOT NULL,
    pqtl_dataset VARCHAR NOT NULL,
    outcome_trait VARCHAR NOT NULL,
    top_snp VARCHAR,
    nsnps INTEGER,
    pp_h0_abf DOUBLE PRECISION,
    pp_h1_abf DOUBLE PRECISION,
    pp_h2_abf DOUBLE PRECISION,
    pp_h3_abf DOUBLE PRECISION,
    pp_h4_abf DOUBLE PRECISION,
    coloc_pass BOOLEAN,
    n_pqtl_snps INTEGER,
    n_gwas_snps INTEGER,
    n_cases INTEGER,
    n_controls INTEGER,
    n_gwas INTEGER,
    s_gwas DOUBLE PRECISION,
    pp4_threshold DOUBLE PRECISION
);

CREATE INDEX idx_coloc_run ON coloc_results (run_id);

-- Table 4: PWCoCo 
CREATE TABLE pwcoco_results (
    run_id VARCHAR NOT NULL REFERENCES runs(run_id),
    protein VARCHAR NOT NULL,
    pqtl_dataset VARCHAR NOT NULL,
    outcome_trait VARCHAR NOT NULL,
    snp1 VARCHAR,
    snp2 VARCHAR,
    nsnps INTEGER,
    h0 DOUBLE PRECISION,
    h1 DOUBLE PRECISION,
    h2 DOUBLE PRECISION,
    h3 DOUBLE PRECISION,
    h4 DOUBLE PRECISION,
    log_abf_all DOUBLE PRECISION
);

CREATE INDEX idx_pwcoco_run ON pwcoco_results (run_id);


/*
Now 2 PWCoCo tables 
* i.e. pQTL-eQTL
* i.e. eQTL-GWAS
*/

CREATE TABLE pwcoco_eqtl_pqtl_results (
    run_id VARCHAR NOT NULL REFERENCES runs(run_id),
    protein VARCHAR NOT NULL,
    pqtl_dataset VARCHAR NOT NULL,
    eqtl_dataset VARCHAR NOT NULL,
    cell_type VARCHAR,
    snp1 VARCHAR,
    snp2 VARCHAR,
    nsnps INTEGER,
    h0 DOUBLE PRECISION,
    h1 DOUBLE PRECISION,
    h2 DOUBLE PRECISION,
    h3 DOUBLE PRECISION,
    h4 DOUBLE PRECISION,
    log_abf_all DOUBLE PRECISION
);

CREATE INDEX idx_pwcoco_eqtl_pqtl_run ON pwcoco_eqtl_pqtl_results (run_id);

CREATE TABLE pwcoco_eqtl_gwas_results (
    run_id VARCHAR NOT NULL REFERENCES runs(run_id),
    protein VARCHAR NOT NULL,
    pqtl_dataset VARCHAR NOT NULL,
    outcome_trait VARCHAR NOT NULL,
    eqtl_dataset VARCHAR NOT NULL,
    cell_type VARCHAR,
    snp1 VARCHAR,
    snp2 VARCHAR,
    nsnps INTEGER,
    h0 DOUBLE PRECISION,
    h1 DOUBLE PRECISION,
    h2 DOUBLE PRECISION,
    h3 DOUBLE PRECISION,
    h4 DOUBLE PRECISION,
    log_abf_all DOUBLE PRECISION
);

CREATE INDEX idx_pwcoco_eqtl_gwas_run ON pwcoco_eqtl_gwas_results (run_id);

-- Table 7. SMR+HEIDI
CREATE TABLE smr_results (
    run_id VARCHAR NOT NULL REFERENCES runs(run_id),
    protein VARCHAR NOT NULL,
    pqtl_dataset VARCHAR NOT NULL,
    phenotype VARCHAR NOT NULL,
    eqtl_dataset VARCHAR NOT NULL,
    cell_type VARCHAR,
    data_type VARCHAR,
    gene VARCHAR,
    gene_id VARCHAR,
    qtl_name VARCHAR,
    probe_id VARCHAR,
    probe_chr VARCHAR,
    probe_bp BIGINT,
    top_snp VARCHAR,
    top_snp_chr VARCHAR,
    top_snp_bp BIGINT,
    a1 VARCHAR,
    a2 VARCHAR,
    freq DOUBLE PRECISION,
    b_gwas DOUBLE PRECISION,
    se_gwas DOUBLE PRECISION,
    p_gwas DOUBLE PRECISION,
    b_eqtl DOUBLE PRECISION,
    se_eqtl DOUBLE PRECISION,
    p_eqtl DOUBLE PRECISION,
    b_smr DOUBLE PRECISION,
    se_smr DOUBLE PRECISION,
    p_smr DOUBLE PRECISION,
    p_smr_multi DOUBLE PRECISION,
    p_heidi DOUBLE PRECISION,
    nsnp_heidi INTEGER,
    q_smr DOUBLE PRECISION,
    chr VARCHAR,
    start_bp BIGINT,
    end_bp BIGINT,
    strand VARCHAR,
    gwas_locus VARCHAR,
    lead_snp VARCHAR,
    lead_snp_bp BIGINT
);

CREATE INDEX idx_smr_run ON smr_results (run_id);

-- HyPrColoc res
CREATE TABLE hyprcoloc_results (
    run_id VARCHAR NOT NULL REFERENCES runs(run_id),
    protein VARCHAR NOT NULL,
    pqtl_dataset VARCHAR NOT NULL,
    outcome_trait VARCHAR NOT NULL,
    eqtl_dataset VARCHAR NOT NULL,
    cell_type VARCHAR,
    data_type VARCHAR,
    iteration INTEGER,
    traits VARCHAR,
    posterior_prob DOUBLE PRECISION,
    regional_prob DOUBLE PRECISION,
    candidate_snp VARCHAR,
    posterior_explained_by_snp DOUBLE PRECISION,
    dropped_trait VARCHAR,
    n_snps INTEGER,
    a1 VARCHAR,
    a2 VARCHAR,
    gwas_beta DOUBLE PRECISION,
    gwas_p DOUBLE PRECISION,
    pqtl_beta DOUBLE PRECISION,
    pqtl_p DOUBLE PRECISION,
    eqtl_beta DOUBLE PRECISION,
    eqtl_p DOUBLE PRECISION
);

CREATE INDEX idx_hyprcoloc_run ON hyprcoloc_results (run_id);

-- Table 9: FinnGen PheWAS safety screen
CREATE TABLE finngen_phewas_safety (
    run_id VARCHAR NOT NULL REFERENCES runs(run_id),
    protein VARCHAR NOT NULL,
    pqtl_dataset VARCHAR NOT NULL,
    method VARCHAR,
    n_instruments_original INTEGER,
    n_instruments_available INTEGER,
    n_instruments INTEGER,
    missing_instruments VARCHAR,
    instrument_completeness DOUBLE PRECISION,
    rsid VARCHAR,
    snp VARCHAR,
    a1 VARCHAR,
    a2 VARCHAR,
    ad_effect_allele_original VARCHAR,
    ad_other_allele_original VARCHAR,
    beta_ad_original VARCHAR,
    beta_ad VARCHAR,
    ad_a1_flipped VARCHAR,
    exposure_effect_allele VARCHAR,
    exposure_other_allele VARCHAR,
    finngen_ref VARCHAR,
    finngen_alt VARCHAR,
    beta_exposure_original VARCHAR,
    beta_exposure VARCHAR,
    exposure_a1_flipped VARCHAR,
    se_exposure VARCHAR,
    p_exposure VARCHAR,
    phewas_a1_flipped VARCHAR,
    beta_mr DOUBLE PRECISION,
    se_mr DOUBLE PRECISION,
    p_mr DOUBLE PRECISION,
    beta_phewas VARCHAR,
    se_phewas VARCHAR,
    p_phewas VARCHAR,
    phenocode VARCHAR,
    phenostring VARCHAR,
    category VARCHAR,
    n_endpoints_tested INTEGER,
    p_bonferroni DOUBLE PRECISION,
    bonferroni_significant BOOLEAN
);

CREATE INDEX idx_finngen_phewas_run ON finngen_phewas_safety (run_id);


-- Table 10: UKB PheWAS safety screen (fallback, gated on FinnGen's coverage manifest)
CREATE TABLE ukb_phewas_safety (
    run_id VARCHAR NOT NULL REFERENCES runs(run_id),
    protein VARCHAR NOT NULL,
    pqtl_dataset VARCHAR NOT NULL,
    method VARCHAR,
    n_instruments_original INTEGER,
    n_instruments_available INTEGER,
    n_instruments INTEGER,
    missing_instruments VARCHAR,
    instrument_completeness DOUBLE PRECISION,
    rsid VARCHAR,
    snp VARCHAR,
    a1 VARCHAR,
    a2 VARCHAR,
    ad_effect_allele_original VARCHAR,
    ad_other_allele_original VARCHAR,
    beta_ad_original VARCHAR,
    beta_ad VARCHAR,
    ad_a1_flipped VARCHAR,
    exposure_effect_allele VARCHAR,
    exposure_other_allele VARCHAR,
    ukb_ref VARCHAR,
    ukb_alt VARCHAR,
    beta_exposure_original VARCHAR,
    beta_exposure VARCHAR,
    exposure_a1_flipped VARCHAR,
    se_exposure VARCHAR,
    p_exposure VARCHAR,
    phewas_a1_flipped VARCHAR,
    beta_mr DOUBLE PRECISION,
    se_mr DOUBLE PRECISION,
    p_mr DOUBLE PRECISION,
    beta_phewas VARCHAR,
    se_phewas VARCHAR,
    p_phewas VARCHAR,
    phenocode VARCHAR,
    phenostring VARCHAR,
    category VARCHAR,
    n_endpoints_tested INTEGER,
    p_bonferroni DOUBLE PRECISION,
    bonferroni_significant BOOLEAN
);

CREATE INDEX idx_ukb_phewas_run ON ukb_phewas_safety (run_id);

-- Table 11: harmonised top cis-hit per protein (alleles aligned to the outcome
CREATE TABLE target_stats (
    run_id VARCHAR NOT NULL REFERENCES runs(run_id),
    protein VARCHAR NOT NULL,
    pqtl_dataset VARCHAR NOT NULL,
    outcome_trait VARCHAR NOT NULL,
    snp VARCHAR,
    a1 VARCHAR,
    a2 VARCHAR,
    frq DOUBLE PRECISION,
    gwas_beta DOUBLE PRECISION,
    gwas_p DOUBLE PRECISION,
    pqtl_beta DOUBLE PRECISION,
    pqtl_p DOUBLE PRECISION
);

CREATE INDEX idx_target_stats_run ON target_stats (run_id);