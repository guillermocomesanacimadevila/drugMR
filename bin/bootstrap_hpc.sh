#!/usr/bin/env bash
set -euo pipefail

# this will be executed from the main Jupyter notebook
# * log into SLURM cluster
# * module load apptainer python=3.12
# * if statement to check whether singularity container == there 
# * if so load - 
# * else - pull from my GHCR

REPO_URL="https://github.com/guillermocomesanacimadevila/drugMR.git"
REPO_DIR="$HOME/drugMR"
SIF_DIR="$REPO_DIR/env"
SIF_FILE="$SIF_DIR/drugmr.sif"
IMAGE_URI="docker://ghcr.io/guillermocomesanacimadevila/drugmr:latest"

if command -v module >/dev/null 2>&1; then
    module load apptainer || true
    module load python/3.12 || true
fi

command -v apptainer >/dev/null 2>&1 || { echo "ERROR: apptainer not found"; exit 1; }
command -v git >/dev/null 2>&1 || { echo "ERROR: git not found"; exit 1; }

if [ ! -d "$REPO_DIR" ]; then
    echo "[TRACKING] drugMR repo not found. Cloning..."
    git clone "$REPO_URL" "$REPO_DIR"
else
    echo "[TRACKING] drugMR repo already exists. Updating..."
    cd "$REPO_DIR"
    git pull
fi

mkdir -p "$SIF_DIR"

if [ -f "$SIF_FILE" ]; then
    echo "[TRACKING] Apptainer image already exists: $SIF_FILE"
else
    echo "[TRACKING] Pulling Apptainer image from GHCR..."
    apptainer pull "$SIF_FILE" "$IMAGE_URI"
fi

echo "[DONE] drugMR HPC environment ready."
echo "[INFO] Repo: $REPO_DIR"
echo "[INFO] SIF:  $SIF_FILE"