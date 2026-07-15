#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(coloc)
  library(arrow)
})

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 5) {
  stop("Usage: Rscript coloc_eqtl.R <pheno_id> <pqtl_dataset> <eqtl_dataset> <n_cases> <n_controls>")
}

pheno_id <- args[1]
pqtl_dataset <- args[2]
eqtl_dataset <- args[3]
n_cases <- as.numeric(args[4])
n_controls <- as.numeric(args[5])

# we have defined that the exposure eQTL == quant ALWAYS
# GWAS and eQTL files have already been matched to the same SNPs from the Python script

exposure_def <- "quant"
pp4_thresh <- 0.75
outcome_def <- "cc" # might have to change this at some other stage
targets_file <- file.path("results", "SMR", eqtl_dataset, pheno_id, paste0(pqtl_dataset, "_", pheno_id, "_prepared_multi_omics_targets.tsv"))
out_dir <- file.path("./results/eQTL_coloc", pqtl_dataset, eqtl_dataset, pheno_id)
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
targets <- fread(targets_file)
all_results <- list()
print(paste0("[TRACKING] Running GWAS -> eQTL coloc for ", nrow(targets), " target x cell-type datasets"))

for (i in seq_len(nrow(targets))) {
  protein_id <- targets$protein[i]
  gene <- targets$gene[i]
  cell_type <- targets$cell_type[i]
  target_dir <- targets$target_dir[i]
  gwas_file <- file.path(target_dir, "gwas.parquet")
  eqtl_file <- file.path(target_dir, "sc_eqtl.parquet")
  out_file <- file.path(out_dir, paste0(pheno_id, "_", protein_id, "_", cell_type, "_eqtl_coloc.tsv"))
  eqtl <- as.data.table(read_parquet(eqtl_file))
  gwas <- as.data.table(read_parquet(gwas_file))

  # convert frequencies into proper MAF
  eqtl[, FRQ := as.numeric(FRQ)]
  gwas[, FRQ := as.numeric(FRQ)]
  eqtl[, MAF := pmin(FRQ, 1 - FRQ)]
  gwas[, MAF := pmin(FRQ, 1 - FRQ)]

  # remove invalid rows for coloc
  eqtl <- eqtl[!is.na(SNP) & !is.na(BETA) & !is.na(SE) & !is.na(MAF)]
  gwas <- gwas[!is.na(SNP) & !is.na(BETA) & !is.na(SE) & !is.na(MAF)]
  eqtl <- eqtl[is.finite(BETA) & is.finite(SE) & is.finite(MAF) & SE > 0 & MAF > 0 & MAF < 1]
  gwas <- gwas[is.finite(BETA) & is.finite(SE) & is.finite(MAF) & SE > 0 & MAF > 0 & MAF < 1]

  # remove duplicated SNPs
  eqtl <- eqtl[order(P), ]
  eqtl <- eqtl[!duplicated(SNP)]
  gwas <- gwas[order(P), ]
  gwas <- gwas[!duplicated(SNP)]

  # rematch GWAS and eQTL SNPs after frequency filtering
  shared_snps <- intersect(eqtl$SNP, gwas$SNP)
  eqtl <- eqtl[SNP %in% shared_snps]
  gwas <- gwas[SNP %in% shared_snps]
  setorder(eqtl, SNP)
  setorder(gwas, SNP)

  # grab top SNP based on eQTL for that given cis-region
  top_snp <- eqtl[order(P)]$SNP[1]

  # conform Ns
  n_eqtl <- eqtl$N[1]
  n_gwas <- n_cases + n_controls
  s_gwas <- n_cases / n_gwas

  # no of SNPs
  n_eqtl_snps <- nrow(eqtl)
  n_gwas_snps <- nrow(gwas)

  print("------------------------------------------------------------")
  print(paste0("[TRACKING] protein_id: ", protein_id))
  print(paste0("[TRACKING] gene: ", gene))
  print(paste0("[TRACKING] cell type: ", cell_type))
  print(paste0("[TRACKING] eQTL file: ", eqtl_file))
  print(paste0("[TRACKING] GWAS file: ", gwas_file))
  print(paste0("[TRACKING] eQTL N: ", n_eqtl))
  print(paste0("[TRACKING] eQTL SNPs: ", n_eqtl_snps))
  print(paste0("[TRACKING] GWAS SNPs: ", n_gwas_snps))
  print(paste0("[TRACKING] Top eQTL SNP: ", top_snp))

  # conform eQTL dataset
  dataset1 <- list(
    snp = eqtl$SNP,
    beta = eqtl$BETA,
    varbeta = eqtl$SE^2,
    MAF = eqtl$MAF,
    N = n_eqtl,
    type = exposure_def
  )

  # conform GWAS dataset
  dataset2 <- list(
    snp = gwas$SNP,
    beta = gwas$BETA,
    varbeta = gwas$SE^2,
    MAF = gwas$MAF,
    N = n_gwas,
    s = s_gwas,
    type = outcome_def
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
  summary[, gene := gene]
  summary[, outcome_trait := pheno_id]
  summary[, cell_type := cell_type]
  summary[, eqtl_dataset := eqtl_dataset]
  summary[, top_snp := top_snp]
  summary[, n_eqtl_snps := n_eqtl_snps]
  summary[, n_gwas_snps := n_gwas_snps]
  summary[, n_eqtl := n_eqtl]
  summary[, n_cases := n_cases]
  summary[, n_controls := n_controls]
  summary[, n_gwas := n_gwas]
  summary[, s_gwas := s_gwas]
  summary[, pp4_threshold := pp4_thresh]
  summary[, coloc_pass := PP.H4.abf >= pp4_thresh]

  setcolorder(summary, c(
    "protein_id",
    "gene",
    "outcome_trait",
    "cell_type",
    "eqtl_dataset",
    "top_snp",
    "nsnps",
    "PP.H0.abf",
    "PP.H1.abf",
    "PP.H2.abf",
    "PP.H3.abf",
    "PP.H4.abf",
    "coloc_pass",
    "n_eqtl_snps",
    "n_gwas_snps",
    "n_eqtl",
    "n_cases",
    "n_controls",
    "n_gwas",
    "s_gwas",
    "pp4_threshold"
  ))

  fwrite(summary, out_file, sep = "\t")
  all_results[[paste(protein_id, cell_type, sep = "_")]] <- summary
  print(paste0("[DONE] Saved COLOC result: ", out_file))
}

final_results <- rbindlist(all_results, fill = TRUE)
final_file <- file.path(out_dir, paste0(pqtl_dataset, "_", pheno_id, "_", eqtl_dataset, "_all_eqtl_coloc.tsv"))
fwrite(final_results, final_file, sep = "\t")

print("============================================================")
print(paste0("[DONE] Saved all GWAS -> eQTL coloc results: ", final_file))