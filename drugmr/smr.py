#!/usr/bin/env python3
import subprocess
from pathlib import Path
import os


def SMR(
    pheno_id: str,
    sumstats: str,
    ref_bfile: str,
    beqtl_summary: str,
    eqtl_dataset: str,
    peqtl_smr: float,
    peqtl_heidi: float,
    thread_num: int,
    maf: float
):
    ref_bfile = Path(ref_bfile)
    sumstats = Path(sumstats)
    beqtl_summary = Path(beqtl_summary)

    # eqtl_dataset can be stuff like SingleBrain/Ast
    # use full path for directory but only cell name for output prefix
    eqtl_dataset = Path(eqtl_dataset)
    eqtl_name = eqtl_dataset.name
    out_dir = Path(f"./results/SMR/{eqtl_dataset}/{pheno_id}")
    os.makedirs(out_dir, exist_ok=True)
    out_file = out_dir / f"{pheno_id}_{eqtl_name}"
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