#!/usr/bin/env python3

# ETL pipeline workflow plan 
# 
# Download .xlxs from Dropbox dir/
# Extract all GC IDs from EBI GWAS catalogue alongside full names for GWAS phenotype
# Match with .xlsx to check overlap - rename the ones in EBI -> acronym
# Map each apptainer to its corresponding gene ID + Uniprot
# Preserve (ONLY cis-pQTLs)
# Map each gene to NCBI gene manifest
# Extract wg_csf manifest file with gene names, coords to hg38, etc... -> 100% overlap
# Download each file from the summary stats dir/
# Map apptainer
# ETL pipeline to -> .parquet (cis-region + biallelic + -INDELs) etc...
# Store as geneID_uniprotID.parquet (unlink() raw file) -> do the same for all proteins (cis-only)