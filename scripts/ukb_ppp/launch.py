#!/usr/bin/env python3
import subprocess
from configparser import ConfigParser
from pathlib import Path

slurm_user = "c.c24102394"
host = "falconlogin.cf.ac.uk"
cfg = ConfigParser()
cfg.read(Path.home() / ".synapseConfig")
synapse_username = cfg["default"]["username"]
synapse_token = cfg["default"]["authtoken"]

remote_cmd = f"""
set -euo pipefail

echo ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
echo "Connected to: $(hostname)"
echo "User: $(whoami)"
echo "Working directory: $(pwd)"
echo ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

cd ~/drugMR
mkdir -p logs
export SYNAPSE_USERNAME='{synapse_username}'
export SYNAPSE_TOKEN='{synapse_token}'
sbatch scripts/ukb_ppp/ukb_ppp_download.sbatch
"""

print(f"Connecting to {slurm_user}@{host}...")
subprocess.run(["ssh", f"{slurm_user}@{host}", remote_cmd], check=True)
print("Remote sbatch job submitted.")