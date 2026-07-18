#!/usr/bin/env python3
import subprocess

slurm_user = "c.user"
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
mkdir -p dat/pQTL/deCODE

sbatch scripts/deCODE/decode_preprocessing.sbatch
"""

print(f"Connecting to {slurm_user}@{host}...")

subprocess.run(["ssh", f"{slurm_user}@{host}", remote_cmd], check=True,)

print("Remote deCODE preprocessing job submitted.")