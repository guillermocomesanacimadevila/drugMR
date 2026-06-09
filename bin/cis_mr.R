#!/usr/bin/env Rscript

# slap onto docker image (env/Dockerfile)
install.packages("remotes")
install.packages("progress")
remotes::install_github("MRCIEU/TwoSampleMR")

suppressPackageStartupMessages({
  library(remotes)
  library(progress)
  library(TwoSampleMR)
  library(data.table)
  library(arrow)
})

args <- commandArgs(trailingOnly = TRUE)
# args -> database (ukb-ppp) -> pheno1 -> out_dir

if (length(args) < 5) {
  stop("Usage: Rscript cis_mr.R <pQTL_dataset> <pqtl_dir> <pheno_id> <pheno_gwas> <out_dir>")
}

pqtl_dataset <- args[1]
pqtl_dir     <- args[2]
pheno_id     <- args[3]
pheno_gwas   <- args[4] 
out_dir      <- args[5]


# MR params
# clump_kb, clump_r2, clump_p1


# create outdir 
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

mr_function <- function(pqtl_dataset, pqtl_dir, pheno_id, pheno_gwas, out_dir) {
  # read outcome data
  df_pheno <- read.csv(pheno_gwas, sep = "\t")
  outcome <- read_outcome_data(
    df_pheno,
    type              = "outcome",
    sep               = "\t",
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
  dim(df_pheno)
  
  # dataset specfic 
  if (pqtl_dataset == "ukb_ppp") {
    proteins <- list.files(
      pqtl_dir,
      pattern = "\\.parquet$", 
      recursive = TRUE,
      full.names = TRUE
    )
    
    pb <- progress_bar$new(total = length(proteins), format = "[:bar] :current/:total (:percent) ETA: :eta")
    
    for (i in proteins) {
      pb$tick()
      protein <- tools::file_path_sans_ext(basename(i))
      print(paste0("[TRACKING] Processing ", protein))
      df <- arrow::read_parquet(i)
      print(dim(df))
      
      # read exposure (i.e. pQTL)
      exposure <- read_exposure_data(
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
      
      # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      # 1. Relevance assumption ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      exposure$pval_col <- as.numeric(exposure$pval_col)
      exposure <- filter(exposure, pval_col < 5e-8) # can set this as a parameter (CLI arg)
      exposure <- exposure[exposure$effect_allele_col > 0.01 & exposure$effect_allele_col < 0.99, ] # same here: can be set as params
      
      # F-stat
      exposure$F <- ((exposure$beta_col) ^ 2) / ((exposure$se_col) ^ 2)
      exposure < exposure[exposure$F >= 10, ]
      
      # harmonise 
      dat <- harmonise_data(exposure, outcome)
      
       
      
      
      
    }
    
  }
}

# first do an if pQTL_dataset == ukb_ppp

# read_parquet
# create progress bar for each protein (1/N)
# after results -> shiny app HTML report -> bring into local and open ./*html