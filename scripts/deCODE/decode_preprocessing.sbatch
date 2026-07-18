#!/bin/bash --login

#SBATCH --job-name=decode
#SBATCH --account=SCWF00036_v_escott-price_34
#SBATCH --partition=compute
#SBATCH --output=/shared/home1/c.user/drugMR/logs/decode_%j.out
#SBATCH --error=/shared/home1/c.user/drugMR/logs/decode_%j.err
#SBATCH --time=2-00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G

set -euo pipefail

project_dir="/shared/home1/c.user/drugMR"
container="/usr/bin/apptainer"

cd "${project_dir}"

mkdir -p logs
mkdir -p dat/pQTL/deCODE

echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
echo "job started      : $(date)"
echo "host             : $(hostname)"
echo "user             : $(whoami)"
echo "working directory: $(pwd)"
echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"

if command -v module >/dev/null 2>&1; then
    module purge
fi

[ -x "${container}" ] || {
    echo "error: apptainer not found (${container})"
    exit 1
}

[ -f env/drugmr.sif ] || {
    echo "error: env/drugmr.sif not found"
    exit 1
}

[ -f scripts/deCODE/preprocesing_workflow.py ] || {
    echo "error: scripts/deCODE/preprocesing_workflow.py not found"
    exit 1
}

[ -f results/deCODE_manifest/deCODE_eur_pgwas_manifest.csv ] || {
    echo "error: results/deCODE_manifest/deCODE_eur_pgwas_manifest.csv not found"
    exit 1
}

echo "[TRACKING] Starting deCODE preprocessing..."

export POLARS_MAX_THREADS="${SLURM_CPUS_PER_TASK}"

"${container}" exec \
    --bind "${project_dir}:${project_dir}" \
    env/drugmr.sif \
    python scripts/deCODE/preprocesing_workflow.py

echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
echo "job finished: $(date)"
echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"