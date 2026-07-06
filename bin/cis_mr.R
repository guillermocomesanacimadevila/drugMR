#!/usr/bin/env Rscript

# slap onto docker image (env/Dockerfile)
# install.packages("remotes")
# install.packages("arrow")
# install.packages("progress")
# remotes::install_github("MRCIEU/TwoSampleMR")
# if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes")
# remotes::install_github("MRCIEU/genetics.binaRies")
# remotes::install_github("mrcieu/ieugwasr")


# TO DO'S
# MAKE A FUNCTION (WHICH WILL BE INHERITED WITHIN THE MR FUNCT) WHICH PER 1/N pQTL-MR -> adds an I or a symbol as kind of a progress bar

suppressPackageStartupMessages({
  library(remotes)
  library(progress)
  library(TwoSampleMR)
  library(ieugwasr)
  library(genetics.binaRies)
  library(data.table)
  library(arrow)
})

args <- commandArgs(trailingOnly = TRUE)
# args -> database (ukb-ppp) -> pheno1 -> out_dir

if (length(args) < 5) {
  stop("Usage: Rscript cis_mr.R <pQTL_dataset> <pqtl_dir> <pheno_id> <pheno_gwas> <ref_bfile>")
}

pqtl_dataset <- args[1] # ukb_ppp
pqtl_dir     <- args[2] # dat/cis_regions/{pqtl_dataset}
pheno_id     <- args[3] # AD
pheno_gwas   <- args[4] # dat/gwas/{pheno_id}
ref_bfile    <- args[5] # /Users/c.user/Desktop/neurobridge/ref/ldsc/1000G_EUR_Phase3_plink/1000G.EUR.QC.ALL"
local_results_dir = "../results/cis-MR"
# out_dir      <- args[6]

out_dir <- "./results/cis-MR"

# MR params
# clump_kb, clump_r2, clump_p1
# pval thresh, f_stat thresh

# create outdir 
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

mr_function <- function(pqtl_dataset, pqtl_dir, pheno_id, pheno_gwas, ref_bfile, out_dir) {
  
  # dataset specfic 
  if (pqtl_dataset == "ukb_ppp") {
    protein_dirs <- list.dirs(
      pqtl_dir,
      recursive = FALSE,
      full.names = TRUE
    )
    
    pb <- progress_bar$new(
      total = length(protein_dirs),
      format = "[:bar] :current/:total (:percent) ETA: :eta | :protein",
      clear = FALSE
    )
    
    # compile all res
    all_results <- list()
    
    for (i in protein_dirs) {
      protein <- basename(i)
      print(paste0("[TRACKING] Processing ", protein))
      pqtl_file <- file.path(i, "pqtl.parquet")
      gwas_file <- file.path(i, "gwas.parquet")
      
      if (!file.exists(pqtl_file)) {
        print(paste0("[CONCERN] Missing pqtl.parquet for ", protein))
        pb$tick(tokens = list(protein = protein))
        next
      }
      
      if (!file.exists(gwas_file)) {
        print(paste0("[CONCERN] Missing gwas.parquet for ", protein))
        pb$tick(tokens = list(protein = protein))
        next
      }
      
      # read exposure (i.e. pQTL)
      cat("> Reading exposure pQTLs from UKBB-PPP...\n")
      df <- arrow::read_parquet(pqtl_file)
      df <- as.data.table(df)
      setorder(df, P)
      df <- df[!duplicated(SNP)]
      df <- as.data.frame(df)
      print(dim(df))
      
      # read outcome data
      cat("> Reading outcome GWAS...\n")
      df_pheno <- arrow::read_parquet(gwas_file)
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
        phenotype_col     = protein
      )
      # check shape 
      dim(exposure)
      
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
      # check shape
      dim(outcome)
      
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
        print(paste0("No instruments after p/F filters for ", protein))
        pb$tick(tokens = list(protein = protein))
        next
      }
      
      # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      # LD Clump -> Ind IVs ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      dat <- harmonise_data(exposure, outcome)
      if (nrow(dat) > 0) {
        dat.clump <- data.table::as.data.table(dat)
        dat.clump[, rsid := SNP]
        dat.clump[, pval := pval.exposure]
        dat.clump <- ld_clump(
          dat.clump,
          clump_kb = 10000,
          clump_r2 = 0.001,
          plink_bin = Sys.which("plink"),
          # plink_bin = genetics.binaRies::get_plink_binary(),
          bfile = ref_bfile
        )
        
        print(paste0("[TRACKING] Instruments after clumping: ", nrow(dat.clump)))
        
        # HEREE ******* - SAVE as INSTRUMENTS.tsv (for that particular protein)
        # Then save onto results/IVs/protein-wide .parquet file with instruments
      } else {
        print(paste0("No harmonised SNPs for ", protein))
        pb$tick(tokens = list(protein = protein))
        next
      }
      
      # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      # Steiger filtering
      # Keep SNPs where R2_GX > R2_GY
      # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      
      dat.clump <- as.data.frame(dat.clump)
      dat.clump <- steiger_filtering(dat.clump)
      dat.clump <- dat.clump[dat.clump$steiger_dir == TRUE, ]
      
      print(paste0("[TRACKING] Instruments after Steiger filtering: ", nrow(dat.clump)))
      
      if (nrow(dat.clump) == 0) {
        print(paste0("No instruments after Steiger filtering for ", protein))
        pb$tick(tokens = list(protein = protein))
        next
      }
      
      dat.clump <- data.table::as.data.table(dat.clump)
      
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
        print(paste0("Not enough instruments for ", protein))
        pb$tick(tokens = list(protein = protein))
        next
      }
      
      res.temp <- dcast(res.temp, id.exposure + id.outcome ~ method, value.var = c("b", "se", "pval"))
      setnames(res.temp, old = grep("Inverse variance weighted", names(res.temp), value = TRUE), new = gsub("Inverse variance weighted", "IVW", grep("Inverse variance weighted", names(res.temp), value = TRUE)))
      setnames(res.temp, old = grep("MR Egger", names(res.temp), value = TRUE), new = gsub("MR Egger", "Egger", grep("MR Egger", names(res.temp), value = TRUE)))
      setnames(res.temp, old = grep("Weighted median", names(res.temp), value = TRUE), new = gsub("Weighted median", "WME", grep("Weighted median", names(res.temp), value = TRUE)))
      res.temp[, protein := protein]
      res.temp[, outcome_trait := pheno_id]
      res.temp[, n_instruments := nrow(dat.clump)]
      res.temp[, egger_intercept := res.pleio$egger_intercept[1]]
      res.temp[, egger_intercept_pval := res.pleio$pval[1]]
      res.temp[, Q := res.het$Q[1]]
      res.temp[, Q_df := res.het$Q_df[1]]
      res.temp[, Q_pval := res.het$Q_pval[1]]
      all_results[[protein]] <- res.temp
      pb$tick(tokens = list(protein = protein))
    }
    
    all_results <- rbindlist(all_results, fill = TRUE)
    
    if (nrow(all_results) == 0) {
      print("[CONCERN] No MR results generated.")
      return(NULL)
    }
    
    # reformat for shiny app / dashboard
    setnames(all_results,
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
    
    # check whether FDR correct or not
    if (length(protein_dirs) > 1) {
      if ("IVW_pval" %in% names(all_results)) {all_results[, IVW_FDR_q := p.adjust(IVW_pval, method = "fdr")]}
      if ("Egger_pval" %in% names(all_results)) {all_results[, Egger_FDR_q := p.adjust(Egger_pval, method = "fdr")]}
      if ("WME_pval" %in% names(all_results)) {all_results[, WME_FDR_q := p.adjust(WME_pval, method = "fdr")]}
      if ("Wald_pval" %in% names(all_results)) {all_results[, Wald_FDR_q := p.adjust(Wald_pval, method = "fdr")]}
    } else {
      if ("IVW_pval" %in% names(all_results)) {all_results[, IVW_FDR_q := IVW_pval]}
      if ("Egger_pval" %in% names(all_results)) {all_results[, Egger_FDR_q := Egger_pval]}
      if ("WME_pval" %in% names(all_results)) {all_results[, WME_FDR_q := WME_pval]}
      if ("Wald_pval" %in% names(all_results)) {all_results[, Wald_FDR_q := Wald_pval]}
    }
    
    keep_cols <- c(
      "protein",
      "outcome_trait",
      "n_instruments",
      "IVW_beta",
      "IVW_se",
      "IVW_pval",
      "IVW_FDR_q",
      "Egger_beta",
      "Egger_se",
      "Egger_pval",
      "Egger_FDR_q",
      "egger_intercept",
      "egger_intercept_pval",
      "WME_beta",
      "WME_se",
      "WME_pval",
      "WME_FDR_q",
      "Wald_beta",
      "Wald_se",
      "Wald_pval",
      "Wald_FDR_q",
      "Q",
      "Q_df",
      "Q_pval"
    )
    
    keep_cols <- keep_cols[keep_cols %in% names(all_results)]
    all_results <- all_results[, ..keep_cols]
    
    out_file <- file.path(out_dir, paste0(pqtl_dataset, "_", pheno_id, "_all_MR.tsv"))
    fwrite(all_results, out_file, sep = "\t")
    print(paste0("Saved all MR results: ", out_file))
  }
}

mr_function(
  pqtl_dataset = pqtl_dataset,
  pqtl_dir     = pqtl_dir,
  pheno_id     = pheno_id,
  pheno_gwas   = pheno_gwas,
  ref_bfile    = ref_bfile,
  out_dir      = out_dir
)

# rename output colnames 
# make ifs for if n proteins > 1 perform FDR if not !=


# first do an if pQTL_dataset == ukb_ppp
# read_parquet
# create progress bar for each protein (1/N)
# after results -> shiny app HTML report -> bring into local and open ./*html