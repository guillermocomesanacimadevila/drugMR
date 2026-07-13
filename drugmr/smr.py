#!/usr/bin/env python3
from pathlib import Path
import subprocess
import os


def SMR(pheno_id: str, sumstats: str, ref_bfile: str, beqtl_summary: str, eqtl_dataset: str, peqtl_smr, peqtl_heidi, thread_num: int, maf: float):
    ref_bfile = Path(ref_bfile)
    beqtl_summary = Path(beqtl_summary)
    out_dir = Path(f"./results/SMR/{eqtl_dataset}/{pheno_id}")
    os.makedirs(out_dir, exist_ok=True)
    out_prefix = out_dir / f"{pheno_id}_{eqtl_dataset}"

    if not Path(f"{ref_bfile}.bed").exists():
        raise FileNotFoundError(f"[CONCERN] Missing {ref_bfile}.bed")

    if not Path(f"{beqtl_summary}.besd").exists():
        raise FileNotFoundError(f"[CONCERN] Missing {beqtl_summary}.besd")

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
  --out {out_prefix}
"""

    print(f"[TRACKING] Running SMR on {pheno_id} using {eqtl_dataset}")

    subprocess.run(cmd_smr, shell=True, check=True, executable="/bin/bash")

    print(f"[TRACKING] Finished SMR for {pheno_id}")