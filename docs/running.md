# Run the pipeline

## Local execution

Select the container profile installed on the machine:

```bash
source env/activate.sh

nextflow run main.nf \
  -profile docker \
  -params-file params/AD.test.yaml \
  --manifest_path assets/qtl_manifest.csv \
  -resume
```

Other local container profiles are `apptainer`, `singularity`, `podman`, `shifter`, and `charliecloud`.

## SLURM execution

There are two supported ways to run on SLURM. Do not mix them.

### Nextflow submits each pipeline task

Run Nextflow directly from a login or workflow node with the `falcon` profile. Nextflow remains a lightweight controller and submits each pipeline process to SLURM with its own resource request:

```bash
module load Java/17
cd /path/to/drugMR
source env/activate.sh

nextflow run main.nf \
  -profile falcon \
  -params-file params/AD.test.yaml \
  --manifest_path assets/qtl_manifest.csv \
  --slurm_account YOUR_ACCOUNT \
  --slurm_partition htc_genoa \
  --container_bind /path/to/external/data \
  -resume
```

This method creates child jobs named `nf-DRUGMR...`. Use it only where cluster policy permits a persistent Nextflow controller on the login or workflow node.

### Nextflow runs inside one submitted allocation

If Nextflow itself must be submitted with `sbatch`, use `-profile falcon,local`. The `falcon` profile enables the Falcon/Apptainer settings, and the trailing `local` profile makes pipeline processes run inside the allocation instead of submitting nested SLURM jobs:

```bash
#!/bin/bash
#SBATCH --job-name=drugmr_ad
#SBATCH --partition=htc_genoa
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --output=drugMR/%x_%j.out
#SBATCH --error=drugMR/%x_%j.err
#SBATCH --account YOUR_ACCOUNT

set -eo pipefail
module load Java/17
cd /path/to/drugMR
source env/activate.sh

nextflow run main.nf \
  -profile falcon,local \
  -params-file params/AD.test.yaml \
  --manifest_path assets/qtl_manifest.csv \
  --slurm_account YOUR_ACCOUNT \
  --slurm_partition htc_genoa \
  --container_bind /path/to/external/data \
  -resume
```

The profile order matters: `falcon,local` is correct because `local` must override the executor selected by `falcon`.

Do not use `sbatch` together with `-profile falcon` alone on a cluster with a one-node-per-user QOS. The controller allocation consumes the permitted node while its nested `nf-DRUGMR...` job remains pending with `QOSMaxNodePerUserLimit`.

Check jobs and accounting with:

```bash
squeue -u "$USER"
sacct -j JOB_ID --format=JobID,JobName,State,ExitCode,Elapsed,ReqMem,MaxRSS --units=G
```

Use `-resume` after a corrected or interrupted run. Nextflow reuses completed tasks whose inputs and commands have not changed.

## Find the run ID

Successful output is stored under `runs/<run_id>/`, where the ID includes the phenotype, pQTL dataset, date, and Git commit. List recent runs with:

```bash
ls -1dt runs/* | head
```

A successful run contains both its outputs and the effective parameters needed to reopen it later:

```text
runs/<run_id>/
├── manifest.json
├── params.lock.yaml
├── pipeline_info/
└── results/
```

The locked parameters are created only after the complete Nextflow workflow succeeds. Keep `params.lock.yaml` with the run when moving or archiving results.
