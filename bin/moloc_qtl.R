#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(moloc)
  library(arrow)
  library(data.table)
})

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 3) {
  stop("Usage: Rscript moloc_qtl.R <pheno_id> <pqtl_dataset> <eqtl_dataset>")
}

pheno_id <- args[1]
pqtl_dataset <- args[2]
eqtl_dataset <- args[3]

targets_file <- file.path("results", "SMR", eqtl_dataset, pheno_id, paste0(pqtl_dataset, "_", pheno_id, "_prepared_multi_omics_targets.tsv"))
out_dir <- file.path("results", "QTL_moloc", pqtl_dataset, eqtl_dataset, pheno_id)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

moloc_runner <- function(pheno_id, pqtl_dataset, eqtl_dataset) {
  targets <- fread(targets_file)
  cat("[TRACKING] Running MOLOC for", nrow(targets), "target x cell-type datasets...\n")
  all_results <- list()
  for (i in seq_len(nrow(targets))) {
    protein <- targets$protein[i]
    gene <- targets$gene[i]
    cell_type <- targets$cell_type[i]
    target_dir <- targets$target_dir[i]

    cat("------------------------------------------------------------\n")
    cat("[TRACKING]", protein, "in", cell_type, "\n")

    gwas <- as.data.table(arrow::read_parquet(file.path(target_dir, "gwas.parquet")))
    pqtl <- as.data.table(arrow::read_parquet(file.path(target_dir, "pqtl.parquet")))
    eqtl <- as.data.table(arrow::read_parquet(file.path(target_dir, "sc_eqtl.parquet")))

    # moloc needs position and allele frequency named POS and MAF
    setnames(gwas, c("BP", "FRQ"), c("POS", "MAF"))
    setnames(pqtl, c("BP", "FRQ"), c("POS", "MAF"))
    setnames(eqtl, c("BP", "FRQ"), c("POS", "MAF"))

    # all datasets were already SNP matched in the Python script
    gwas <- gwas[order(SNP)]
    pqtl <- pqtl[order(SNP)]
    eqtl <- eqtl[order(SNP)]

    trait_dfs <- list(
      gwas = gwas,
      pqtl = pqtl,
      eqtl = eqtl
    )

    cat("[TRACKING] Running GWAS - pQTL - eQTL MOLOC with", nrow(gwas), "SNPs...\n")

    res <- moloc_test(
      trait_dfs,
      prior_var = c(0.01, 0.1, 0.5),
      priors = c(1e-4, 1e-6, 1e-7)
    )

    target_out_dir <- file.path(out_dir, protein, cell_type)
    dir.create(target_out_dir, recursive = TRUE, showWarnings = FALSE)
    post <- as.data.table(res[[1]], keep.rownames = "model")
    best <- as.data.table(res[[3]], keep.rownames = "model")
    nsnps <- data.table(nsnps = res[[2]])
    fwrite(post, file.path(target_out_dir, "posteriors.tsv"), sep = "\t")
    fwrite(best, file.path(target_out_dir, "best_snps.tsv"), sep = "\t")
    fwrite(nsnps, file.path(target_out_dir, "nsnps.tsv"), sep = "\t")
    best_model <- post[which.max(PPA)]
    best_model[, protein := protein]
    best_model[, gene := gene]
    best_model[, phenotype := pheno_id]
    best_model[, cell_type := cell_type]
    best_model[, pqtl_dataset := pqtl_dataset]
    best_model[, eqtl_dataset := eqtl_dataset]
    best_model[, nsnps := res[[2]]]
    all_results[[paste(protein, cell_type, sep = "_")]] <- best_model
    cat("[TRACKING] Best model:", best_model$model, "\n")
    cat("[TRACKING] Best model PPA:", best_model$PPA, "\n")
    cat("[TRACKING] Saved:", target_out_dir, "\n")
  }

  final_results <- rbindlist(all_results, fill = TRUE)
  final_file <- file.path(out_dir, paste0(pqtl_dataset, "_", pheno_id, "_", eqtl_dataset, "_moloc_summary.tsv"))
  fwrite(final_results, final_file, sep = "\t")

  cat("============================================================\n")
  cat("[TRACKING] MOLOC complete\n")
  cat("[TRACKING] Summary saved to", final_file, "\n")
}

moloc_runner(
  pheno_id = pheno_id,
  pqtl_dataset = pqtl_dataset,
  eqtl_dataset = eqtl_dataset
)