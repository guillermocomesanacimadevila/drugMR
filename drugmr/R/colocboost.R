#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(tools)
  library(colocboost)
})

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 4) {
  stop("Usage: Rscript colocboost.R <ld_file> <out_file> <sumstat_file_1> [<sumstat_file_2> ...]")
}

ld_file       <- args[1]
out_file      <- args[2]
sumstat_files <- args[3:length(args)]

trait_names <- file_path_sans_ext(basename(sumstat_files))
sumstat <- setNames(lapply(sumstat_files, fread), trait_names)
LD <- as.matrix(fread(ld_file, header = FALSE))
rownames(LD) <- colnames(LD) <- sumstat[[1]]$variant
res <- colocboost(sumstat = sumstat, LD = LD)
summary <- as.data.table(get_colocboost_summary(res)$cos_summary)

fwrite(summary, out_file, sep = "\t")
print(paste0("[DONE] Saved ColocBoost result (", length(sumstat), " traits: ", paste(trait_names, collapse = ", "), "): ", out_file))