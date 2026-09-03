#!/usr/bin/env python3
import os
import subprocess
from pathlib import Path

from drugmr import paths


def SMR(
    pheno_id: str,
    sumstats: str,
    ref_bfile: str,
    beqtl_summary: str,
    eqtl_dataset: str,
    peqtl_smr: float,
    peqtl_heidi: float,
    thread_num: int,
    maf: float,
    out_dir: str = "synthesis"
):
    ref_bfile = Path(ref_bfile)
    sumstats = Path(sumstats)
    beqtl_summary = Path(beqtl_summary)

    # eqtl_dataset can be stuff like SingleBrain/Ast
    # use full path for directory but only cell name for output prefix
    eqtl_dataset = Path(eqtl_dataset)
    # SMR(GWAS x eQTL) result for a given (pheno_id, eqtl_dataset) never depends on
    # pqtl_dataset, so this defaults to the shared synthesis/ tree rather than a
    # per-run out_dir - every pqtl_dataset run reuses the same computation instead
    # of re-running the smr binary from scratch
    raw_out_dir = paths.smr_raw_dir(eqtl_dataset, pheno_id, out_dir)
    os.makedirs(raw_out_dir, exist_ok=True)
    out_file = paths.smr_raw_prefix(eqtl_dataset, pheno_id, out_dir)
    print(f"[TRACKING] Running SMR on {pheno_id} using {eqtl_dataset}")

    cmd_smr = f"""
set -euo pipefail
smr \
  --bfile {ref_bfile} \
  --gwas-summary {sumstats} \
  --beqtl-summary {beqtl_summary} \
  --maf {maf} \
  --peqtl-smr {peqtl_smr} \
  --peqtl-heidi {peqtl_heidi} \
  --thread-num {thread_num} \
  --out {out_file}
"""

    subprocess.run(cmd_smr, shell=True, check=True, executable="/bin/bash")