#!/usr/bin/env python3
import subprocess

slurm_user = "c.c24102394"
host = "falconlogin.cf.ac.uk"

remote_cmd = """
set -euo pipefail

echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
echo "Connected to: $(hostname)"
echo "User: $(whoami)"
echo "Working directory: $(pwd)"
echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"

cd ~/drugMR
mkdir -p logs
mkdir -p dat/pQTL/wu_csf

sbatch scripts/wu_csf/wu_csf_preprocessing.sbatch
"""

print(f"Connecting to {slurm_user}@{host}...")

subprocess.run(["ssh", f"{slurm_user}@{host}", remote_cmd], check=True,)

print("Remote WU CSF preprocessing job submitted.")