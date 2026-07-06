#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(moloc)
  library(arrow)
  library(data.table)
  library(rjson)
})

# # - n_o of traits for protein X
# - pheno_id
# - pqtl_dataset
# - pqtl_dir (maybe)

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 3) {
  stop("Usage: Rscript moloc.R <pheno_id> <pqtl_dataset> <json_file>")
}

pheno_id     <- args[1]
pqtl_dataset <- args[2]
json_file    <- args[3]

out_dir <- file.path("results", "moloc", pqtl_dataset)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# df <- arrow::read_parquet(pqtl_file)
# df <- as.data.table(df)
# setnames(pheno_id,  c("BP", "FRQ"), c("POS", "MAF"))

moloc_runner <- function(pheno_id, json_file) {
  pqtls_dir <- file.path("dat", "cis_regions", pqtl_dataset)
  json_data <- fromJSON(file = json_file)
  cat("[TRACKING] Running MOLOC for", length(json_data), "proteins...\n")
  
  for (protein in names(json_data)) {
    cat("[TRACKING]", protein, "\n")
    traits <- json_data[[protein]]
    protein_dir <- file.path(pqtls_dir, protein)
    pqtl <- arrow::read_parquet(file.path(protein_dir, "pqtl.parquet"))
    pqtl <- as.data.table(pqtl)
    setnames(pqtl, c("BP", "FRQ"), c("POS", "MAF"))
    trait_dfs <- list(pqtl = pqtl)
    
    for (trait in traits) {
      cat("[TRACKING]   Loading:", trait, "\n")
      if (trait == pheno_id) {
        trait_file <- file.path(protein_dir, "gwas.parquet")
      } else {
        trait_file <- file.path(protein_dir, "mediators", paste0(trait, ".parquet"))
      }
      
      df <- arrow::read_parquet(trait_file)
      df <- as.data.table(df)
      setnames(df, c("BP", "FRQ"), c("POS", "MAF"))
      trait_dfs[[trait]] <- df
    }
    
    cat("[TRACKING]   Running MOLOC...\n")
    n_traits <- length(trait_dfs)
    
    # arjust priors in accordance with n_o of traits 
    priors <- c(1e-4, rep(1e-6, n_traits - 1))
    
    res <- moloc_test(
      trait_dfs,
      prior_var = c(0.01, 0.1, 0.5),
      priors = priors
    )
    
    # for 3 traits (default priors and prior_sigma^2)
    # res <- moloc_test(
    #   trait_dfs,
    #   prior_var = c(0.01, 0.1, 0.5),
    #   priors = c(1e-4, 1e-6, 1e-7)
    # )
    
    protein_out_dir <- file.path(out_dir, protein)
    dir.create(protein_out_dir, recursive = TRUE, showWarnings = FALSE)
    post <- as.data.table(res[[1]], keep.rownames = "model")
    best <- as.data.table(res[[3]], keep.rownames = "model")
    nsnps <- data.table(nsnps = res[[2]])
    fwrite(post, file.path(protein_out_dir, "posteriors.csv"))
    fwrite(best, file.path(protein_out_dir, "best_snps.csv"))
    fwrite(nsnps, file.path(protein_out_dir, "nsnps.csv"))
    cat("[TRACKING]   Saved:", protein_out_dir, "\n")
    
  }
  cat("[TRACKING] MOLOC complete.\n")
}

moloc_runner(
  pheno_id = pheno_id,
  json_file = json_file
)
