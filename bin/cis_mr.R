#!/usr/bin/env Rscript

# slap onto docker image (env/Dockerfile)
# install.packages("remotes")
# install.packages("arrow")
# install.packages("progress")
# remotes::install_github("MRCIEU/TwoSampleMR")
# if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes")
# remotes::install_github("MRCIEU/genetics.binaRies")

suppressPackageStartupMessages({
  library(remotes)
  library(progress)
  library(TwoSampleMR)
  library(genetics.binaRies)
  library(data.table)
  library(arrow)
})

args <- commandArgs(trailingOnly = TRUE)
# args -> database (ukb-ppp) -> pheno1 -> out_dir

if (length(args) < 6) {
  stop("Usage: Rscript cis_mr.R <pQTL_dataset> <pqtl_dir> <pheno_id> <pheno_gwas> <ref_bfile> <out_dir>")
}

pqtl_dataset <- args[1]
pqtl_dir     <- args[2]
pheno_id     <- args[3]
pheno_gwas   <- args[4]
ref_bfile    <- args[5]
out_dir      <- args[6]


# MR params
# clump_kb, clump_r2, clump_p1
# pval thresh, f_stat thresh

# create outdir 
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

mr_function <- function(pqtl_dataset, pqtl_dir, pheno_id, pheno_gwas, ref_bfile, out_dir) {
  # read outcome data
  df_pheno <- read.csv(pheno_gwas, sep = "\t")
  cat("> Reading outcome GWAS...")
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
  dim(df_pheno)
  
  # dataset specfic 
  if (pqtl_dataset == "ukb_ppp") {
    proteins <- list.files(
      pqtl_dir,
      pattern = "\\.parquet$", 
      recursive = TRUE,
      full.names = TRUE
    )
    
    pb <- progress_bar$new(
      total = length(proteins),
      format = "[:bar] :current/:total (:percent) ETA: :eta | :protein",
      clear = FALSE
    )
    
    for (i in proteins) {
      pb$tick()
      protein <- tools::file_path_sans_ext(basename(i))
      print(paste0("[TRACKING] Processing ", protein))
      df <- arrow::read_parquet(i)
      print(dim(df))
      
      # read exposure (i.e. pQTL)
      cat("> Reading exposure pQTLs from UKBB-PPP...")
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
          plink_bin = genetics.binaRies::get_plink_binary(),
          bfile = ref_bfile
        )
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
        res.het   <- mr_heterogeneity(dat.clump, method_list = c("mr_ivw"))
        res.temp <- data.table::as.data.table(res.temp)
        res.temp[, egger_intercept := res.pleio$egger_intercept[1]]
        res.temp[, egger_intercept_pval := res.pleio$pval[1]]
        res.temp[, Q := res.het$Q[1]]
        res.temp[, Q_df := res.het$Q_df[1]]
        res.temp[, Q_pval := res.het$Q_pval[1]]
        
      } else if (nrow(dat.clump) == 1) {
        res.temp <- mr(dat.clump, method_list = c("mr_wald_ratio"))
        res.temp <- data.table::as.data.table(res.temp)
        res.temp[, egger_intercept := NA_real_]
        res.temp[, egger_intercept_pval := NA_real_]
        res.temp[, Q := NA_real_]
        res.temp[, Q_df := NA_real_]
        res.temp[, Q_pval := NA_real_]
        
      } else {
        print(paste0("Not enough instruments for ", protein))
        next
      }
      res.temp[, protein := protein]
      res.temp[, outcome_trait := pheno_id]
      res.temp[, n_instruments := nrow(dat.clump)]
      res.temp[, FDR_q := p.adjust(pval, method = "fdr")]
      out_file <- file.path(out_dir, paste0(protein, "_", pheno_id, "_MR.tsv"))
      fwrite(res.temp, out_file, sep = "\t")
      print(paste0("Saved MR results: ", out_file))
      pb$tick(tokens = list(protein = protein))
    } 
    
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

# first do an if pQTL_dataset == ukb_ppp
# read_parquet
# create progress bar for each protein (1/N)
# after results -> shiny app HTML report -> bring into local and open ./*html