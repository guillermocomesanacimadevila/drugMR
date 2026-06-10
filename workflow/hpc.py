#!/usr/bin/env python3
import subprocess

def hpc(cmd: str, falcon_user: str):
    full_cmd = f"ssh {falcon_user}@falconlogin.cf.ac.uk '{cmd}'"
    print(full_cmd)
    subprocess.run(full_cmd, shell=True, check=True)

def get_remote_paths(falcon_user: str):
    remote = f"/shared/home1/{falcon_user}/drugMR"
    sif = f"{remote}/env/drugmr.sif"
    return remote, sif

def prepare_hpc(
    falcon_user: str,
    github_repo: str = "https://github.com/guillermocomesanacimadevila/drugMR.git",
    image_uri: str = "docker://ghcr.io/guillermocomesanacimadevila/drugmr:latest"
):
    remote, sif = get_remote_paths(falcon_user)

    hpc(f"""
set -euo pipefail

module load apptainer || true

if [ ! -d "{remote}" ]; then
    echo "[TRACKING] Cloning drugMR repo..."
    git clone "{github_repo}" "{remote}"
else
    echo "[TRACKING] drugMR repo already exists. Pulling latest changes..."
    cd "{remote}"
    git pull
fi

mkdir -p "{remote}/env"

if [ ! -f "{sif}" ]; then
    echo "[TRACKING] Pulling drugMR Apptainer image..."
    apptainer pull "{sif}" "{image_uri}"
else
    echo "[TRACKING] drugMR Apptainer image already exists."
fi
""", falcon_user)


# preprocessing pQTLs != part of the pipeline 