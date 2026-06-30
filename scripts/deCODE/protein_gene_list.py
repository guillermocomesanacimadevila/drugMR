#!/usr/bin/env python3
import polars as pl
import argparse
from pathlib import Path

# PULL ALL FILE NAMES WITHIN pQTL DECODE PORTAL
# ONLY IN FORMAT (id1_id2_gene....)# EXTRACT GENE NAME
# BUT THIS SCRIPT IS JUST TO GRAB GENE NAMES AND OVERLAP THEM WITH NCBI REFERENCE
# % MISSINGNESS


# FOR TARGET PROTEIN WITH >1 APTAMER
# EXTRACT LOCUS COORDS FOR GENE THAT ENCODES PROTEIN Y
# GRAB LEAD SNP FOR BOTH AND COMPARE IT - LOWEST PVAL APTANER FOR SAME TARGET == KEEP
# OTHER == DISCARD 

# WORKFLOW STEPS (deCODE pQTLs -> GRCh38/hg38)
# id1_id2_gene - extract gene name - overlap with NCBI - % overlap
# for file in decode portal
# gunzip file.txt.gz
# extract gene name
# download file from decode portal 
# df.remove(df.isna())
# harmonise to LDSC-like colnames
# find gene name (and coordinates - chr:start-end on NCBI ref)
# grab cis-region (+/- 1mb from start and end coords)
# save as .parquet
# remove previous txt.gz file for protein X
# carry on with the next protein