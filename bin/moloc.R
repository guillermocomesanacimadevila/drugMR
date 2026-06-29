#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(remotes)
  library(moloc)
})

# Must have columns `SNP` or `CHR` and `POS`; `BETA`, `SE`; `N` and `MAF` (to estimate sdY);
# if a case control: `Ncases`.
# If eqtl/mqtl: `ProbeID`.
# If the regions are defined based on ProbeID, these must match 
# the ProbeID in the file.
# Optionally, "A1", "A2" if want to match alleles


moloc <- moloc_test(
  # data,
  prior_var=c(0.01, 0.1, 0.5),
  priors=c(1e-04, 1e-06, 1e-07)
)


