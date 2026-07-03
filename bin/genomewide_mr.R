#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(remotes)
  library(progress)
  library(TwoSampleMR)
  library(ieugwasr)
  library(genetics.binaRies)
  library(data.table)
})

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 7) {
  stop("Usage: Rscript genomewide_mr.R <M_id> <M_gwas> <M_N> <pheno_id> <pheno_gwas> <ref_bfile> <out_dir>")
}

M_id       <- args[1] # pTau
M_gwas     <- args[2] # dat/gwas/mediators/{M_id}.tsv
M_N        <- as.numeric(args[3]) # mediator GWAS N
pheno_id   <- args[4] # AD
pheno_gwas <- args[5] # dat/gwas/{pheno_id}.tsv
ref_bfile  <- args[6] # /Users/c.user/Desktop/neurobridge/ref/ldsc/1000G_EUR_Phase3_plink/1000G.EUR.QC.ALL
out_dir    <- args[7] # results/networkMR/M_Y/{pheno_id}

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

genomewide_mr <- function(M_id, M_gwas, M_N, pheno_id, pheno_gwas, ref_bfile, out_dir) {
  
  print(paste0("[TRACKING] Running genome-wide MR: ", M_id, " -> ", pheno_id))
  
  # read exposure mediator GWAS
  cat("> Reading exposure mediator GWAS...\n")
  df <- fread(M_gwas)
  df <- as.data.table(df)
  setorder(df, P)
  df <- df[!duplicated(SNP)]
  df[, N := M_N]
  df <- as.data.frame(df)
  print(dim(df))
  
  # read outcome GWAS
  cat("> Reading outcome GWAS...\n")
  df_pheno <- fread(pheno_gwas)
  df_pheno <- as.data.table(df_pheno)
  setorder(df_pheno, P)
  df_pheno <- df_pheno[!duplicated(SNP)]
  df_pheno <- as.data.frame(df_pheno)
  print(dim(df_pheno))
  
  exposure <- format_data(
    df,
    type              = "exposure",
    snp_col           = "SNP",
    beta_col          = "BETA",
    se_col            = "SE",
    effect_allele_col = "A1",
    other_allele_col  = "A2",
    eaf_col           = "FRQ",
    pval_col          = "P",
    samplesize_col    = "N",
    phenotype_col     = M_id
  )
  
  outcome <- format_data(
    df_pheno,
    type              = "outcome",
    snp_col           = "SNP",
    beta_col          = "BETA",
    se_col            = "SE",
    effect_allele_col = "A1",
    other_allele_col  = "A2",
    eaf_col           = "FRQ",
    pval_col          = "P",
    samplesize_col    = "N",
    phenotype_col     = pheno_id
  )
  
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Relevance assumption ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  exposure$pval.exposure <- as.numeric(exposure$pval.exposure)
  exposure <- exposure[exposure$pval.exposure < 5e-8, ]
  exposure <- exposure[
    exposure$eaf.exposure > 0.01 &
      exposure$eaf.exposure < 0.99,
  ]
  exposure$F <- (exposure$beta.exposure^2) / (exposure$se.exposure^2)
  exposure <- exposure[exposure$F >= 10, ]
  
  print(paste0("[TRACKING] Instruments after p/F filters: ", nrow(exposure)))
  
  if (nrow(exposure) == 0) {
    print(paste0("[CONCERN] No instruments after p/F filters for ", M_id))
    return(NULL)
  }
  
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Harmonise exposure and outcome ~~~~~~~~~~~~~~~~~~~~~~~
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  dat <- harmonise_data(exposure, outcome)
  
  if (nrow(dat) == 0) {
    print(paste0("[CONCERN] No harmonised SNPs for ", M_id))
    return(NULL)
  }
  
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # LD Clump -> Independent IVs ~~~~~~~~~~~~~~~~~~~~~~~~~~
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  dat.clump <- data.table::as.data.table(dat)
  dat.clump[, rsid := SNP]
  dat.clump[, pval := pval.exposure]
  
  dat.clump <- ld_clump(
    dat.clump,
    clump_kb = 10000,
    clump_r2 = 0.001,
    plink_bin = Sys.which("plink"),
    bfile = ref_bfile
  )
  
  print(paste0("[TRACKING] Instruments after clumping: ", nrow(dat.clump)))
  
  if (nrow(dat.clump) == 0) {
    print(paste0("[CONCERN] No instruments after clumping for ", M_id))
    return(NULL)
  }
  
  # ~~~~~~~~~~
  # run MR
  # ~~~~~~~~~~
  if (nrow(dat.clump) >= 3) {
    res.temp <- mr(
      dat.clump,
      method_list = c(
        "mr_ivw",
        "mr_egger_regression",
        "mr_weighted_median"
      )
    )
    
    res.pleio <- mr_pleiotropy_test(dat.clump)
    res.het <- mr_heterogeneity(dat.clump, method_list = c("mr_ivw"))
    res.temp <- data.table::as.data.table(res.temp)
    
  } else if (nrow(dat.clump) == 1) {
    res.temp <- mr(dat.clump, method_list = c("mr_wald_ratio"))
    res.pleio <- data.table(egger_intercept = NA_real_, pval = NA_real_)
    res.het <- data.table(Q = NA_real_, Q_df = NA_real_, Q_pval = NA_real_)
    res.temp <- data.table::as.data.table(res.temp)
    
  } else {
    print(paste0("[CONCERN] Not enough instruments for ", M_id))
    return(NULL)
  }
  
  res.temp <- dcast(
    res.temp,
    id.exposure + id.outcome ~ method,
    value.var = c("b", "se", "pval")
  )
  
  setnames(
    res.temp,
    old = grep("Inverse variance weighted", names(res.temp), value = TRUE),
    new = gsub("Inverse variance weighted", "IVW", grep("Inverse variance weighted", names(res.temp), value = TRUE))
  )
  
  setnames(
    res.temp,
    old = grep("MR Egger", names(res.temp), value = TRUE),
    new = gsub("MR Egger", "Egger", grep("MR Egger", names(res.temp), value = TRUE))
  )
  
  setnames(
    res.temp,
    old = grep("Weighted median", names(res.temp), value = TRUE),
    new = gsub("Weighted median", "WME", grep("Weighted median", names(res.temp), value = TRUE))
  )
  
  setnames(
    res.temp,
    old = c(
      "b_IVW", "se_IVW", "pval_IVW",
      "b_Egger", "se_Egger", "pval_Egger",
      "b_WME", "se_WME", "pval_WME",
      "b_Wald ratio", "se_Wald ratio", "pval_Wald ratio"
    ),
    new = c(
      "IVW_beta", "IVW_se", "IVW_pval",
      "Egger_beta", "Egger_se", "Egger_pval",
      "WME_beta", "WME_se", "WME_pval",
      "Wald_beta", "Wald_se", "Wald_pval"
    ),
    skip_absent = TRUE
  )
  
  res.temp[, mediator := M_id]
  res.temp[, mediator_N := M_N]
  res.temp[, outcome_trait := pheno_id]
  res.temp[, n_instruments := nrow(dat.clump)]
  res.temp[, egger_intercept := res.pleio$egger_intercept[1]]
  res.temp[, egger_intercept_pval := res.pleio$pval[1]]
  res.temp[, Q := res.het$Q[1]]
  res.temp[, Q_df := res.het$Q_df[1]]
  res.temp[, Q_pval := res.het$Q_pval[1]]
  
  keep_cols <- c(
    "mediator",
    "mediator_N",
    "outcome_trait",
    "n_instruments",
    "IVW_beta",
    "IVW_se",
    "IVW_pval",
    "Egger_beta",
    "Egger_se",
    "Egger_pval",
    "egger_intercept",
    "egger_intercept_pval",
    "WME_beta",
    "WME_se",
    "WME_pval",
    "Wald_beta",
    "Wald_se",
    "Wald_pval",
    "Q",
    "Q_df",
    "Q_pval"
  )
  
  keep_cols <- keep_cols[keep_cols %in% names(res.temp)]
  res.temp <- res.temp[, ..keep_cols]
  out_file <- file.path(out_dir, paste0(M_id, "_", pheno_id, "_genomewide_MR.tsv"))
  fwrite(res.temp, out_file, sep = "\t")
  print(paste0("[TRACKING] Saved genome-wide MR result: ", out_file))
  return(res.temp)
}

genomewide_mr(
  M_id       = M_id,
  M_gwas     = M_gwas,
  M_N        = M_N,
  pheno_id   = pheno_id,
  pheno_gwas = pheno_gwas,
  ref_bfile  = ref_bfile,
  out_dir    = out_dir
)