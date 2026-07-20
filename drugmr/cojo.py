#!/usr/bin/env python3
import subprocess 
from pathlib import Path
import os 

# gcta64 \
#     --bfile dat/ref/1000G_EUR_Phase3_plink/1000G.EUR.QC.ALL \
#     --cojo-file BLNK.cojo.ma \
#     --cojo-slct \
#     --cojo-p 5e-8 \ 
#     --cojo-wind 10000 \
#     --cojo-collinear 0.9 \
#     --out results/cojo/BLNK

# wget GTCA into docker image

def COJO(ref_bfile, cojo_sumstats, p_thresh, wind_thresh, collinear_thresh, out_prefix):
    out_prefix = Path(out_prefix)
    os.makedirs(out_prefix.parent, exist_ok=True)
    cmd_cojo = f"""
set -euo pipefail
gcta \
    --bfile {ref_bfile} \
    --cojo-file {cojo_sumstats} \
    --cojo-slct \
    --cojo-p {p_thresh} \
    --cojo-wind {wind_thresh} \
    --cojo-collinear {collinear_thresh} \
    --out {out_prefix}
"""
    subprocess.run(cmd_cojo, shell=True, check=True, executable="/bin/bash")