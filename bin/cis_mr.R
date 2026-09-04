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

# FIXED IN THIS VERSION
# [1] n_instruments == 2 no longer silently dropped (was falling into the else -> next)
# [2] Egger / WME gated on IV count via MIN_IV_EGGER and MIN_IV_WME
#     (Egger needs residual df; InSIDE not defensible within one cis locus)
# [3] phenotype_col was being passed a VALUE not a COLUMN NAME -> exposure name silently defaulted
# [4] minimum detectable OR (MDE) reported per protein so nulls are interpretable

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
# results_dir mirrors drugmr/paths.py's out_dir (e.g. "runs/<run_id>/results") -
# must match what local.py's require_output() checks, hence not left hardcoded
results_dir <- ifelse(length(args) >= 6, args[6], "results")

out_dir <- file.path(results_dir, "cis-MR")

# let's just assume for now that the ldsc ref stuff is inside dat/ref
# ld <- ".dat/ref/ldsc/eur_w_ld_chr" -> for mediators
# hm3 <- ".dat/ref/ldsc/weights_hm3_no_hla" -> for mediators

# MR params
# clump_kb, clump_r2, clump_p1
# pval thresh, f_stat thresh
CLUMP_KB      <- 10000
CLUMP_R2      <- 0.001
PVAL_THRESH   <- 5e-8
F_THRESH      <- 10
MIN_IV_WME    <- 3
MIN_IV_EGGER  <- 3

# deCODE's own pQTL extraction has used both "FRQ" and "MAF" for the same allele-
# frequency column across different runs/scripts - pick whichever is actually present
# rather than hardcoding one, so a schema change upstream doesn't silently zero out
# every instrument via format_data()'s eaf.exposure being NA for a missing column
eaf_col_name <- function(df) if ("FRQ" %in% names(df)) "FRQ" else "MAF"

# create outdir 
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

mr_function <- function(pqtl_dataset, pqtl_dir, pheno_id, pheno_gwas, ref_bfile, out_dir) {
  
  # dataset specfic 
  supported_datasets <- c(
    "ukb_ppp",
    "decode",
    "wu_csf",
    "wingo_brain"
  )
  
  if (!pqtl_dataset %in% supported_datasets) {
    stop(
      paste0(
        "Unsupported pQTL dataset: ",
        pqtl_dataset,
        ". Supported datasets: ",
        paste(supported_datasets, collapse = ", ")
      )
    )
  }
  
  if (pqtl_dataset == "ukb_ppp") {
    dataset_label <- "UKBB-PPP"
  }
  
  if (pqtl_dataset == "decode") {
    dataset_label <- "deCODE"
  }
  
  if (pqtl_dataset == "wu_csf") {
    dataset_label <- "WU-CSF"
  }

  if (pqtl_dataset == "wingo_brain") {
    dataset_label <- "Wingo_Brain"
  }
  
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
  all_instruments <- list()

  
  # this is just so if one protein crashes later, we dont lose all the ones that worked
  out_file_running <- file.path(out_dir, paste0(pqtl_dataset, "_", pheno_id, "_all_MR.running.tsv"))
  out_file_final <- file.path(out_dir, paste0(pqtl_dataset, "_", pheno_id, "_all_MR.tsv"))

  # for saving instruments
  instruments_dir <- file.path(out_dir, "instruments")
  dir.create(instruments_dir, recursive = TRUE, showWarnings = FALSE)
  out_instruments_running <- file.path(instruments_dir, paste0(pqtl_dataset, "_", pheno_id, "_all_MR_instruments.running.tsv"))
  out_instruments_final <- file.path(instruments_dir, paste0(pqtl_dataset, "_", pheno_id, "_all_MR_instruments.tsv"))
 
  if (file.exists(out_file_running)) {
    file.remove(out_file_running)
  }

  if (file.exists(out_instruments_running)) {
    file.remove(out_instruments_running)
  }
  
  for (i in protein_dirs) {
    
    tryCatch({
      
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
      cat(paste0("> Reading exposure pQTLs from ", dataset_label, "...\n"))
      df <- arrow::read_parquet(pqtl_file)
      
      # skip empty pQTL datasets
      if (nrow(df) == 0) {
        print(paste0("[CONCERN] Empty pqtl.parquet for ", protein))
        pb$tick(tokens = list(protein = protein))
        next
      }
      
      df <- as.data.table(df)
      setorder(df, P)
      df <- df[!duplicated(SNP)]
      df <- as.data.frame(df)
      print(dim(df))
      
      # read outcome data
      cat("> Reading outcome GWAS...\n")
      df_pheno <- arrow::read_parquet(gwas_file)
      
      # skip empty GWAS datasets
      if (nrow(df_pheno) == 0) {
        print(paste0("[CONCERN] Empty gwas.parquet for ", protein))
        pb$tick(tokens = list(protein = protein))
        next
      }
      
      df_pheno <- as.data.table(df_pheno)
      setorder(df_pheno, P)
      df_pheno <- df_pheno[!duplicated(SNP)]
      df_pheno <- as.data.frame(df_pheno)
      print(dim(df_pheno))
      
      # FIX [3] format_data expects phenotype_col to be a COLUMN NAME, not the value itself.
      # Passing the protein string meant no such column existed -> silently defaulted to "exposure".
      df$phenotype <- protein
      df_pheno$phenotype <- pheno_id
      
      exposure <- format_data(
        df,
        type              = "exposure",
        snp_col           = "SNP",
        beta_col          = "BETA",
        se_col            = "SE",
        effect_allele_col = "A1",
        other_allele_col  = "A2",
        eaf_col           = eaf_col_name(df),
        pval_col          = "P",
        samplesize_col    = "N",
        phenotype_col     = "phenotype"
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
        eaf_col           = eaf_col_name(df_pheno),
        pval_col          = "P",
        samplesize_col    = "N",
        phenotype_col     = "phenotype"
      )
      # check shape
      dim(outcome)
      
      # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      # Relevance assumption ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      # NOTE F here is the squared z-statistic. With p < 5e-8 already imposed,
      # z^2 is ~30 by construction, so F >= 10 is non-binding. Kept for documentation.
      exposure$pval.exposure <- as.numeric(exposure$pval.exposure)
      exposure <- exposure[exposure$pval.exposure < PVAL_THRESH, ]
      exposure <- exposure[
        exposure$eaf.exposure > 0.01 &
          exposure$eaf.exposure < 0.99,
      ]
      exposure$F <- (exposure$beta.exposure^2) / (exposure$se.exposure^2)
      # NA F (e.g. missing SE upstream) must be dropped explicitly - exposure[NA, ]
      # keeps an all-NA row rather than excluding it, which silently corrupted
      # instruments and made harmonise_data() find "no harmonised SNPs" downstream
      exposure <- exposure[!is.na(exposure$F) & exposure$F >= F_THRESH, ]
      
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
        
        dat.clump <- tryCatch(
          {
            ld_clump(
              dat.clump,
              clump_kb = CLUMP_KB,
              clump_r2 = CLUMP_R2,
              plink_bin = Sys.which("plink"),
              # plink_bin = genetics.binaRies::get_plink_binary(),
              bfile = ref_bfile
            )
          },
          error = function(e) {
            print(paste0("[CONCERN] LD clumping failed for ", protein, " - ", e$message))
            return(NULL)
          }
        )
        
        if (is.null(dat.clump)) {
          print(paste0("[CONCERN] No clumped file / no LD clump results for ", protein))
          pb$tick(tokens = list(protein = protein))
          next
        }
        
        if (nrow(dat.clump) == 0) {
          print(paste0("[CONCERN] No instruments after clumping for ", protein))
          pb$tick(tokens = list(protein = protein))
          next
        }
        
        print(paste0("[TRACKING] Instruments after clumping: ", nrow(dat.clump)))
        
      } else {
        print(paste0("No harmonised SNPs for ", protein))
        pb$tick(tokens = list(protein = protein))
        next
      }
      
      # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      # Steiger filtering
      # Keep SNPs where R2_GX > R2_GY
      # NOTE near-vacuous for cis-pQTLs given the variance-explained asymmetry.
      # kept for documentation; reverse MR (AD -> protein) is the real directionality test.
      # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      
      dat.clump <- as.data.frame(dat.clump)
      
      dat.clump <- tryCatch(
        {
          steiger_filtering(dat.clump)
        },
        error = function(e) {
          print(paste0("[CONCERN] Steiger filtering failed for ", protein, " - ", e$message))
          return(NULL)
        }
      )
      
      if (is.null(dat.clump)) {
        pb$tick(tokens = list(protein = protein))
        next
      }
      
      dat.clump <- dat.clump[dat.clump$steiger_dir == TRUE, ]
      
      print(paste0("[TRACKING] Instruments after Steiger filtering: ", nrow(dat.clump)))
      
      if (nrow(dat.clump) == 0) {
        print(paste0("No instruments after Steiger filtering for ", protein))
        pb$tick(tokens = list(protein = protein))
        next
      }
      
      dat.clump <- data.table::as.data.table(dat.clump)

      # saving instruments
      instruments.temp <- data.table::copy(dat.clump)
      setorder(instruments.temp, pval.exposure)
      instruments.temp[, instrument_rank := seq_len(.N)]
      instruments.temp[, protein := protein]
      instruments.temp[, pqtl_dataset := pqtl_dataset]
      instruments.temp[, outcome_trait := pheno_id]
      instruments.temp[, selection_method := "LD_CLUMP"]
      instruments.temp[, clump_kb := CLUMP_KB]
      instruments.temp[, clump_r2 := CLUMP_R2]
      instruments.temp[, used_in_mr := TRUE]
      # z of the exposure effect -> sets the max r2 with a PAV that conditional
      # analysis can resolve downstream: r2_max = 1 - (5.45/z)^2
      instruments.temp[, z_exposure := beta.exposure / se.exposure]
      instruments.temp[, r2_max_resolvable := 1 - (5.45 / abs(z_exposure))^2]
      instrument_cols <- c(
        "protein",
        "pqtl_dataset",
        "outcome_trait",
        "instrument_rank",
        "SNP",
        "effect_allele.exposure",
        "other_allele.exposure",
        "eaf.exposure",
        "beta.exposure",
        "se.exposure",
        "pval.exposure",
        "samplesize.exposure",
        "F",
        "z_exposure",
        "r2_max_resolvable",
        "effect_allele.outcome",
        "other_allele.outcome",
        "eaf.outcome",
        "beta.outcome",
        "se.outcome",
        "pval.outcome",
        "samplesize.outcome",
        "rsq.exposure",
        "rsq.outcome",
        "steiger_dir",
        "steiger_pval",
        "mr_keep",
        "palindromic",
        "ambiguous",
        "selection_method",
        "clump_kb",
        "clump_r2",
        "used_in_mr"
      )
      
      instrument_cols <- instrument_cols[instrument_cols %in% names(instruments.temp)]
      instruments.temp <- instruments.temp[,
        ..instrument_cols
      ]
      
      all_instruments[[protein]] <- instruments.temp
      # save as it goes, because otherwise if one protein explodes
      # all previously selected instrument sets would be lost
      fwrite(instruments.temp, out_instruments_running, sep = "\t", append = file.exists(out_instruments_running), col.names = !file.exists(out_instruments_running))
      print(paste0("[TRACKING] Saved ", nrow(instruments.temp), " final cis-MR instruments for ", protein))
      
      # ~~~~~~~~~~
      # run MR
      # ~~~~~~~~~~
      # FIX [1] n == 2 used to fall through to the else and get dropped entirely.
      # FIX [2] method list now built from the IV count:
      #   IVW    >= 2
      #   WME    >= MIN_IV_WME    (breakdown point meaningless below this)
      #   Egger  >= MIN_IV_EGGER  (needs residual df; InSIDE indefensible in one cis locus)
      n_iv <- nrow(dat.clump)
      
      if (n_iv >= 2) {
        
        method_list <- c("mr_ivw")
        
        if (n_iv >= MIN_IV_WME) {
          method_list <- c(method_list, "mr_weighted_median")
        }
        
        if (n_iv >= MIN_IV_EGGER) {
          method_list <- c(method_list, "mr_egger_regression")
        }
        
        print(paste0("[TRACKING] n_IV = ", n_iv, " -> methods: ", paste(method_list, collapse = ", ")))
        
        res.temp <- tryCatch(
          {
            mr(
              dat.clump,
              method_list = method_list
            )
          },
          error = function(e) {
            print(paste0("[CONCERN] MR failed for ", protein, " - ", e$message))
            return(NULL)
          }
        )
        
        if (is.null(res.temp)) {
          pb$tick(tokens = list(protein = protein))
          next
        }
        
        # egger intercept only meaningful where egger itself was run
        if (n_iv >= MIN_IV_EGGER) {
          res.pleio <- tryCatch(
            {
              mr_pleiotropy_test(dat.clump)
            },
            error = function(e) {
              print(paste0("[CONCERN] Pleiotropy test failed for ", protein, " - ", e$message))
              data.table(egger_intercept = NA_real_, pval = NA_real_)
            }
          )
        } else {
          res.pleio <- data.table(egger_intercept = NA_real_, pval = NA_real_)
        }
        
        # Q is FLAGGED not used as a pass/fail gate - in cis, heterogeneity is
        # usually allelic heterogeneity (multiple independent regulatory variants),
        # not horizontal pleiotropy. and with 2-3 IVs it has no power anyway.
        res.het <- tryCatch(
          {
            mr_heterogeneity(dat.clump, method_list = c("mr_ivw"))
          },
          error = function(e) {
            print(paste0("[CONCERN] Heterogeneity test failed for ", protein, " - ", e$message))
            data.table(Q = NA_real_, Q_df = NA_real_, Q_pval = NA_real_)
          }
        )
        
        res.temp <- data.table::as.data.table(res.temp)
        
      } else if (n_iv == 1) {
        res.temp <- tryCatch(
          {
            mr(dat.clump, method_list = c("mr_wald_ratio"))
          },
          error = function(e) {
            print(paste0("[CONCERN] Wald ratio failed for ", protein, " - ", e$message))
            return(NULL)
          }
        )
        
        if (is.null(res.temp)) {
          pb$tick(tokens = list(protein = protein))
          next
        }
        
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
      res.temp[, pqtl_dataset := pqtl_dataset]
      res.temp[, outcome_trait := pheno_id]
      res.temp[, n_instruments := n_iv]
      res.temp[, methods_run := paste(unique(gsub("^(b|se|pval)_", "", grep("^(b|se|pval)_", names(res.temp), value = TRUE))), collapse = ",")]
      res.temp[, egger_intercept := res.pleio$egger_intercept[1]]
      res.temp[, egger_intercept_pval := res.pleio$pval[1]]
      res.temp[, Q := res.het$Q[1]]
      res.temp[, Q_df := res.het$Q_df[1]]
      res.temp[, Q_pval := res.het$Q_pval[1]]
      all_results[[protein]] <- res.temp
      
      # save as it goes, because otherwise if one protein explodes everything is gone
      fwrite(
        res.temp,
        out_file_running,
        sep = "\t",
        append = file.exists(out_file_running),
        col.names = !file.exists(out_file_running)
      )
      
      pb$tick(tokens = list(protein = protein))
      
    }, error = function(e) {
      protein <- basename(i)
      print(paste0("[CONCERN] Protein fully failed but moving on: ", protein, " - ", e$message))
      pb$tick(tokens = list(protein = protein))
    })
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
  
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # FIX [4] minimum detectable effect
  # OR detectable with 80% power at alpha = 0.05 given the primary SE.
  # without this a null is uninterpretable - "no effect" vs "no power".
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # skip_absent above means IVW_beta/Wald_beta (etc.) only exist as columns if
  # at least one protein in THIS batch actually used that method - true almost
  # always for the full ~2,900-protein X -> Y run (guaranteed some Wald-ratio
  # protein by sheer volume), but not guaranteed for a small pre-filtered batch
  # (e.g. NetworkMR's X -> M leg, just a handful of proteins) where every protein
  # could land on n_IV >= 2 and never trigger Wald ratio at all - backfill any
  # missing one as NA so the fifelse() calls below don't reference a genuinely
  # absent column
  for (col in c("IVW_beta", "IVW_se", "IVW_pval", "Wald_beta", "Wald_se", "Wald_pval")) {
    if (!col %in% names(all_results)) {
      all_results[, (col) := NA_real_]
    }
  }

  all_results[, primary_beta := fifelse(!is.na(IVW_beta), IVW_beta, Wald_beta)]
  all_results[, primary_se := fifelse(!is.na(IVW_se), IVW_se, Wald_se)]
  all_results[, primary_pval := fifelse(!is.na(IVW_pval), IVW_pval, Wald_pval)]
  all_results[, primary_method := fifelse(!is.na(IVW_beta), "IVW", "Wald ratio")]
  # all_results[, MDE_OR := exp((1.96 + 0.84) * primary_se)]
  
  # check whether FDR correct or not
  # FDR is applied WITHIN this pQTL dataset only. cross-dataset agreement is
  # treated as replication, not pooled into one testing family.
  if (length(protein_dirs) > 1) {
    if ("IVW_pval" %in% names(all_results)) {all_results[, IVW_FDR_q := p.adjust(IVW_pval, method = "fdr")]}
    if ("Egger_pval" %in% names(all_results)) {all_results[, Egger_FDR_q := p.adjust(Egger_pval, method = "fdr")]}
    if ("WME_pval" %in% names(all_results)) {all_results[, WME_FDR_q := p.adjust(WME_pval, method = "fdr")]}
    if ("Wald_pval" %in% names(all_results)) {all_results[, Wald_FDR_q := p.adjust(Wald_pval, method = "fdr")]}
    if ("primary_pval" %in% names(all_results)) {all_results[, primary_FDR_q := p.adjust(primary_pval, method = "fdr")]}
  } else {
    if ("IVW_pval" %in% names(all_results)) {all_results[, IVW_FDR_q := IVW_pval]}
    if ("Egger_pval" %in% names(all_results)) {all_results[, Egger_FDR_q := Egger_pval]}
    if ("WME_pval" %in% names(all_results)) {all_results[, WME_FDR_q := WME_pval]}
    if ("Wald_pval" %in% names(all_results)) {all_results[, Wald_FDR_q := Wald_pval]}
    if ("primary_pval" %in% names(all_results)) {all_results[, primary_FDR_q := primary_pval]}
  }
  
  keep_cols <- c(
    "protein",
    "pqtl_dataset",
    "outcome_trait",
    "n_instruments",
    "methods_run",
    "primary_method",
    "primary_beta",
    "primary_se",
    "primary_pval",
    "primary_FDR_q",
    # "MDE_OR",
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

  # save the final combined instrument set across all proteins
  all_instruments <- rbindlist(all_instruments, fill = TRUE, use.names = TRUE)
  if (nrow(all_instruments) > 0) {
    
    setorder(
      all_instruments,
      protein,
      instrument_rank
    )
    
    fwrite(
      all_instruments,
      out_instruments_final,
      sep = "\t"
    )
    
    print(
      paste0(
        "Saved all cis-MR instruments: ",
        out_instruments_final
      )
    )
    
    print(
      paste0(
        "Total instruments saved: ",
        nrow(all_instruments)
      )
    )
    
    print(
      paste0(
        "Proteins with instruments: ",
        uniqueN(all_instruments$protein)
      )
    )
    
  } else {
    
    print(
      "[CONCERN] No cis-MR instruments were generated."
    )
    
  }
  
  # instrument count distribution -> tells you immediately how often the
  # multi-instrument branch actually fires under r2 < 0.001
  print("[TRACKING] Instrument count distribution:")
  print(table(all_results$n_instruments))
  
  if (file.exists(out_file_running)) {
    print(
      paste0(
        "Saved running MR results too: ",
        out_file_running
      )
    )
  }
  
  if (file.exists(out_instruments_running)) {
    print(
      paste0(
        "Saved running instrument results too: ",
        out_instruments_running
      )
    )
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