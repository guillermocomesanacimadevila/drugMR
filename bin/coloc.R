#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(coloc)
  library(arrow)
})

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 7) {
  stop("Usage: Rscript coloc.R <pqtl_dataset> <protein_id> <pheno_id> <gwas_parquet> <pqtl_parquet> <n_cases> <n_controls>")
}

# we have defined that the exposure pQTL == quant ALWAYS
# The locus definition schema is defined on previous script (we do NOT need to do that from coloc.R)

pqtl_dataset <- args[1]
protein_id   <- args[2]
pheno_id     <- args[3]
gwas_file    <- args[4]
pqtl_file    <- args[5]
n_cases      <- as.numeric(args[6])
n_controls   <- as.numeric(args[7])

exposure_def <- "quant"
pp4_thresh   <- 0.75
outcome_def  <- "cc"   # might have to change this at some other stage
out_dir <- file.path("./results/coloc", pqtl_dataset)
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# out_file = just 1 file eventually, but for now 1 file per protein
# for each locus grab top pQTL SNP and add as a col (i.e. topSNP)
out_file <- file.path(out_dir, paste0(pheno_id, "_", protein_id, "_coloc.tsv"))
protein  <- as.data.table(read_parquet(pqtl_file))
gwas     <- as.data.table(read_parquet(gwas_file))

# grab top SNP based on pQTL for that given cis-region
protein <- protein[order(P), ]
protein <- protein[!duplicated(SNP)]
gwas <- gwas[order(P), ]
gwas <- gwas[!duplicated(SNP)]
top_snp <- protein$SNP[1]

# conform Ns
n_protein <- max(protein$N)
n_gwas <- n_cases + n_controls
s_gwas <- n_cases / n_gwas

# no of SNPs
n_pqtl_snps <- nrow(protein)
n_gwas_snps <- nrow(gwas)

print(paste0("[TRACKING] protein_id: ", protein_id))
print(paste0("[TRACKING] pQTL file: ", pqtl_file))
print(paste0("[TRACKING] GWAS file: ", gwas_file))
print(paste0("[TRACKING] pQTL SNPs: ", nrow(protein)))
print(paste0("[TRACKING] GWAS SNPs: ", nrow(gwas)))

# conform pQTL dataset
dataset1 <- list(
  snp     = protein$SNP,
  beta    = protein$BETA,
  varbeta = protein$SE^2,
  MAF     = protein$FRQ,
  N       = n_protein,
  type    = exposure_def
)

# conform GWAS dataset
dataset2 <- list(
  snp     = gwas$SNP,
  beta    = gwas$BETA,
  varbeta = gwas$SE^2,
  MAF     = gwas$FRQ,
  N       = n_gwas,
  s       = s_gwas,
  type    = outcome_def
)

# compile results with default priors
res <- coloc.abf(
  dataset1 = dataset1,
  dataset2 = dataset2,
  p1 = 1e-4,
  p2 = 1e-4,
  p12 = 1e-5
)

print(paste0("[RESULTS] COLOC results!: ", out_file))
print(res)

# make master results table 
summary <- as.data.table(as.list(res$summary))
summary[, protein_id := protein_id]
summary[, outcome_trait := pheno_id]
summary[, top_snp := top_snp]
summary[, n_pqtl_snps := n_pqtl_snps]
summary[, n_gwas_snps := n_gwas_snps]
summary[, n_cases := n_cases]
summary[, n_controls := n_controls]
summary[, n_gwas := n_gwas]
summary[, s_gwas := s_gwas]
summary[, pp4_threshold := pp4_thresh]
summary[, coloc_pass := PP.H4.abf >= pp4_thresh]

setcolorder(summary, c(
  "protein_id",
  "outcome_trait",
  "top_snp",
  "nsnps",
  "PP.H0.abf",
  "PP.H1.abf",
  "PP.H2.abf",
  "PP.H3.abf",
  "PP.H4.abf",
  "coloc_pass",
  "n_pqtl_snps",
  "n_gwas_snps",
  "n_cases",
  "n_controls",
  "n_gwas",
  "s_gwas",
  "pp4_threshold"
))

fwrite(summary, out_file, sep = "\t")
print(paste0("[DONE] Saved COLOC result: ", out_file))