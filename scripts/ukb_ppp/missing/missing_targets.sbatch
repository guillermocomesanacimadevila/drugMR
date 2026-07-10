#!/bin/bash --login

#SBATCH --job-name=ukb_missing
#SBATCH --account=SCWF00036_v_escott-price_34
#SBATCH --partition=compute
#SBATCH --output=/shared/home1/c.c24102394/drugMR/logs/ukb_missing_%j.out
#SBATCH --error=/shared/home1/c.c24102394/drugMR/logs/ukb_missing_%j.err
#SBATCH --time=0-10:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G

set -euo pipefail

project_dir="/shared/home1/c.c24102394/drugMR"
container="/usr/bin/apptainer"

cd "${project_dir}"
mkdir -p logs
mkdir -p dat/pQTL/ukb_ppp/missing

"${container}" exec \
    --env SYNAPSE_USERNAME="${SYNAPSE_USERNAME}" \
    --env SYNAPSE_TOKEN="${SYNAPSE_TOKEN}" \
    env/drugmr.sif \
    python scripts/ukb_ppp/missing/missing_targets.py \
        --synapse-username "${SYNAPSE_USERNAME}" \
        --synapse-token "${SYNAPSE_TOKEN}" \
        --gene-list scripts/ukb_ppp/missing/missing.txt