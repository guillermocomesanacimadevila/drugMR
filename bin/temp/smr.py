#!/usr/bin/env python3
import argparse
import polars as pl
from pathlib import Path
import subprocess
import os 


# STEPS
# ETL pipeline for each single-cell eQTL dataset... (SingleBrain & OneK1K) ****
# For any given target with FDR_q < 0.05 and passes egger intercept and cochran Q & Coloc PP.H4 > 0.75 (or thresh...)
# Go to single-cell datasets (each cell type...) map same cis-region (and intersect with GWAS) - both files == parquet
# Save cis-region from sc-eQTL parquet onto .dat/cis_regions/{dataset}/target_id/...
# Run SMR with pre-imputed SMR ready files in datalake (for each cell type) only on cis-region
# FDR-correct for all tested targets within 1 dataset + P_HEIDI > 0.01 (FDR_pass; HEIDI_pass bools)
# Store within drugMR PostgreSQL database