#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(coloc)
})

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 5) {
  stop("Usage: Rscript coloc.R <exposure_locus> <outcome_locus> <pheno_id> <protein_id> <pqtl_dataset>")
}


# we have defined that the exposure pQTL == quant ALWAYS
# The locusd definition schema is defined on previous script (we do NOT need to do that from coloc.R)

exposure_locus <- args[1] # straight tsv of pQTL with matching SNPs to GWAS 
outcome_locus  <- args[2] # straight tsv of GWAS with matching SNPs to pQTL
pheno_id       <- args[3]
protein_id     <- agrs[4] # loop overthis in backend py script for all prots within res dir where FDR_p < 0.05 
pqtl_dataset   <- args[5] 


exposure_def <- "quant"
outcome_def <- "cc" # might have to change this at some other stage
out_dir <- "../results"
dir.create(out_dir, showWarnings = FALSE, recursive = FALSE)
out_file <- file.path(out_dir, paste0(pheno_id, "_", protein_id, "_coloc.tsv"))

n <- c(
  exposure_locus$N[1],
  protein$N[1]
)




