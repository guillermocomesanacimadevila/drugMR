#!/bin/bash --login

#SBATCH --job-name=wu_csf
#SBATCH --account=SCWF00036_v_escott-price_34
#SBATCH --partition=compute
#SBATCH --output=/shared/home1/c./drugMR/logs/wu_csf_%j.out
#SBATCH --error=/shared/home1/c./drugMR/logs/wu_csf_%j.err
#SBATCH --time=2-00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G

set -euo pipefail 

project_dir="/shared/home1/c./drugMR"
container="/usr/bin/apptainer"

cd "${project_dir}"

mkdir -p logs
mkdir -p dat/pQTL/wu_csf

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

[ -f scripts/wu_csf/csf_preprocessing_pipeline.py ] || {
    echo "error: scripts/wu_csf/csf_preprocessing_pipeline.py not found"
    exit 1
}

[ -f results/WS_CSF_manifest/wu_csf_manifest_autosomal.csv ] || {
    echo "error: results/WS_CSF_manifest/wu_csf_manifest_autosomal.csv not found"
    exit 1
}

echo "[TRACKING] Starting wu_csf preprocessing..."
export POLARS_MAX_THREADS="${SLURM_CPUS_PER_TASK}"

"${container}" exec \
    --bind "${project_dir}:${project_dir}" \
    env/drugmr.sif \
    python scripts/wu_csf/csf_preprocessing_pipeline.py

echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
echo "job finished: $(date)"
echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"