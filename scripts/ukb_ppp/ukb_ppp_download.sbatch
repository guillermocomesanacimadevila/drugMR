#!/bin/bash --login

#SBATCH --job-name=ukb_ppp
#SBATCH --account=SCWF00036_v_escott-price_34
#SBATCH --partition=compute
#SBATCH --output=/shared/home1/c.c24102394/drugMR/logs/ukb_ppp_%A_%a.out
#SBATCH --error=/shared/home1/c.c24102394/drugMR/logs/ukb_ppp_%A_%a.err
#SBATCH --time=0-05:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --array=1-10%2

set -euo pipefail

project_dir="/shared/home1/c.c24102394/drugMR"
container="/usr/bin/apptainer"
chunk_size=25

start=$(( (SLURM_ARRAY_TASK_ID - 1) * chunk_size + 1 ))
end=$(( SLURM_ARRAY_TASK_ID * chunk_size ))

cd "${project_dir}"

mkdir -p logs

echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
echo "job started      : $(date)"
echo "host             : $(hostname)"
echo "user             : $(whoami)"
echo "working directory: $(pwd)"
echo "array task       : ${SLURM_ARRAY_TASK_ID}"
echo "protein range    : ${start}-${end}"
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

[ -f scripts/ukb_ppp/list_proteins.py ] || {
    echo "error: scripts/ukb_ppp/list_proteins.py not found"
    exit 1
}

[ -n "${SYNAPSE_USERNAME:-}" ] || {
    echo "error: SYNAPSE_USERNAME not set"
    exit 1
}

[ -n "${SYNAPSE_TOKEN:-}" ] || {
    echo "error: SYNAPSE_TOKEN not set"
    exit 1
}

echo "starting protein processing..."

"${container}" exec \
    --env SYNAPSE_USERNAME="${SYNAPSE_USERNAME}" \
    --env SYNAPSE_TOKEN="${SYNAPSE_TOKEN}" \
    env/drugmr.sif \
    python scripts/ukb_ppp/list_proteins.py \
        --synapse-username "${SYNAPSE_USERNAME}" \
        --synapse-token "${SYNAPSE_TOKEN}" \
        --start-index "${start}" \
        --end-index "${end}"

echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
echo "job finished: $(date)"
echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"