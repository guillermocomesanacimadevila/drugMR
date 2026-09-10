#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(hyprcoloc)
  library(arrow)
})

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 6) {
  stop("Usage: Rscript hyprcoloc.R <pqtl_dataset> <protein> <cell_type> <pheno_id> <trio_dir> <local_results_dir> [prior_1] [prior_c] [reg_thresh] [align_thresh] [equal_thresholds] [qtl_type]")
}

# trio_dir holds pqtl.parquet / gwas.parquet / eqtl.parquet for this protein x
# cell type - already matched to the same SNP set and aligned to a common
# effect allele (the GWAS A1) by bin/hyprcoloc_targets.py
pqtl_dataset      <- args[1]
protein           <- args[2]
cell_type         <- args[3]
pheno_id          <- args[4]
trio_dir          <- args[5]
local_results_dir <- args[6]

# HyPrColoc sensitivity plot tuning (params/*.yaml gates.hyprcoloc) - prior_c/reg_thresh/align_thresh
# are comma-separated vectors on the command line
HC_PRIOR_1          <- ifelse(length(args) >= 7, as.numeric(args[7]), 1e-4)
HC_PRIOR_C          <- if (length(args) >= 8) as.numeric(strsplit(args[8], ",")[[1]]) else c(0.05, 0.02, 0.01, 0.005)
HC_REG_THRESH       <- if (length(args) >= 9) as.numeric(strsplit(args[9], ",")[[1]]) else c(0.5, 0.6, 0.7)
HC_ALIGN_THRESH     <- if (length(args) >= 10) as.numeric(strsplit(args[10], ",")[[1]]) else c(0.5, 0.6, 0.7)
HC_EQUAL_THRESHOLDS <- ifelse(length(args) >= 11, as.logical(args[11]), TRUE)
QTL_TYPE            <- ifelse(length(args) >= 12, args[12], "eqtl")

# known qtl_type values get their real published abbreviation casing (eQTL,
# sQTL, mQTL, caQTL); an unregistered future type still gets a sane label
# instead of failing here
QTL_TYPE_LABELS <- c(eqtl = "eQTL", sqtl = "sQTL", mqtl = "mQTL", caqtl = "caQTL", pqtl = "pQTL")
QTL_LABEL <- if (QTL_TYPE %in% names(QTL_TYPE_LABELS)) QTL_TYPE_LABELS[[QTL_TYPE]] else toupper(QTL_TYPE)

# must match paths.hyprcoloc_dataset_out(...).parent.parent in
# bin/hyprcoloc_targets.py, which is where the caller looks for out_file
out_dir <- file.path(local_results_dir, "hyprcoloc", "by_eqtl_source")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
out_file <- file.path(out_dir, paste0(pheno_id, "_", protein, "_", cell_type, "_hyprcoloc.tsv"))

hyprcoloc_runner <- function(pqtl_dataset, protein, cell_type, pheno_id, trio_dir) {

  trait_files <- c(
    pqtl = file.path(trio_dir, "pqtl.parquet"),
    gwas = file.path(trio_dir, "gwas.parquet"),
    eqtl = file.path(trio_dir, "eqtl.parquet")
  )

  print(paste0("[TRACKING] Running HyPrColoc for ", protein, " x ", cell_type))

  dfs <- lapply(trait_files, function(f) as.data.table(read_parquet(f)))
  rsid <- dfs$gwas$SNP
  n_snps <- length(rsid)

  print(paste0("[TRACKING] Matched pQTL/GWAS/", QTL_LABEL, " SNPs: ", n_snps))

  if (n_snps < 2) {
    print(paste0("[CONCERN] Fewer than 2 matched SNPs for ", protein, " x ", cell_type, " - skipping"))
    return(NULL)
  }

  # trio_dir already carries the same SNP set (same order not guaranteed) so
  # re-align every trait onto the GWAS SNP order before building the matrices
  dfs <- lapply(dfs, function(df) df[match(rsid, df$SNP)])

  trait_names <- c(
    paste0("pQTL_", protein),
    paste0("GWAS_", pheno_id),
    paste0(QTL_LABEL, "_", cell_type)
  )

  betas <- cbind(dfs$pqtl$BETA, dfs$gwas$BETA, dfs$eqtl$BETA)
  ses   <- cbind(dfs$pqtl$SE, dfs$gwas$SE, dfs$eqtl$SE)
  colnames(betas) <- trait_names
  colnames(ses)   <- trait_names
  rownames(betas) <- rsid
  rownames(ses)   <- rsid

  result <- hyprcoloc::hyprcoloc(
    betas,
    ses,
    trait.names = trait_names,
    snp.id = rsid
  )

  res <- as.data.table(result$results)
  res[, protein := protein]
  res[, cell_type := cell_type]
  res[, pqtl_dataset := pqtl_dataset]
  res[, outcome_trait := pheno_id]
  res[, n_snps := n_snps]

  fwrite(res, out_file, sep = "\t")
  print(paste0("[DONE] Saved HyPrColoc result: ", out_file))

  sensitivity_file <- file.path(out_dir, paste0(pheno_id, "_", protein, "_", cell_type, "_sensitivity.pdf"))
  pdf(sensitivity_file)
  hyprcoloc::sensitivity.plot(
    effect.est       = betas,
    effect.se        = ses,
    trait.names      = trait_names,
    snp.id           = rsid,
    prior.1          = HC_PRIOR_1,
    prior.c          = HC_PRIOR_C,
    reg.thresh       = HC_REG_THRESH,
    align.thresh     = HC_ALIGN_THRESH,
    equal.thresholds = HC_EQUAL_THRESHOLDS
  )
  dev.off()
  print(paste0("[DONE] Saved HyPrColoc sensitivity plot: ", sensitivity_file))

  return(res)
}

result <- hyprcoloc_runner(
  pqtl_dataset = pqtl_dataset,
  protein      = protein,
  cell_type    = cell_type,
  pheno_id     = pheno_id,
  trio_dir     = trio_dir
)

print(result)