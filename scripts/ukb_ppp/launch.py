#!/usr/bin/env python3

import subprocess
from configparser import ConfigParser
from pathlib import Path

slurm_user = ""
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

apptainer exec \\
  --env SYNAPSE_USERNAME='{synapse_username}' \\
  --env SYNAPSE_TOKEN='{synapse_token}' \\
  env/drugmr.sif \\
  python scripts/ukb_ppp/list_proteins.py \\
    --synapse-username '{synapse_username}' \\
    --synapse-token '{synapse_token}'
"""

print(f"Connecting to {slurm_user}@{host}...")
subprocess.run(["ssh", f"{slurm_user}@{host}", remote_cmd], check=True)
print("Remote job finished.")