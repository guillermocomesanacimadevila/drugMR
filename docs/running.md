# Run drugMR

Run Nextflow from the root of the cloned repository. One command runs one params file. One params file contains one outcome trait and one pQTL dataset.

Use `-resume` when rerunning a corrected or interrupted analysis. Nextflow reuses completed tasks when their inputs and commands have not changed.

## Example: local

This example runs the pipeline and dashboard on the same computer.

First activate drugMR:

```bash
cd /path/to/drugMR
source env/activate.sh
```

Then run one configured outcome and pQTL pair:

```bash
nextflow run main.nf \
  -profile docker \
  -params-file params/AD.ukb_ppp.yaml \
  --manifest_path assets/qtl_manifest.csv \
  -resume
```

Replace `docker` with the installed container runtime when required:

```text
apptainer
singularity
podman
shifter
charliecloud
```

After a successful run, open that exact result in Python:

```python
import drugmr as dm

dm.results(
    config="runs/AD_ukb_ppp_YYYYMMDD_abcdef0/params.lock.yaml"
)
```

Docker must be running on the machine that calls `dm.results()` because it starts PostgreSQL with Docker Compose.

## Example: HPC (on a SLURM cluster!)

This example runs Nextflow on HPC and runs the dashboard on a local computer. Keep a complete clone of the repository on both machines.

SLURM and Nextflow control different levels of the run. SLURM allocates cluster resources. Nextflow decides where each pipeline process is executed. The `falcon` profile for example, selects the SLURM executor and Apptainer. The `local` profile selects the local Nextflow executor. Here, local means inside the machine or allocation where the Nextflow controller is already running. It does not mean your laptop.

Use `falcon` when Nextflow starts outside an allocation and should submit each process to SLURM. Use `falcon,local` when `sbatch` has already created an allocation for the whole workflow. In that second arrangement, the `falcon` part still supplies the Falcon and Apptainer settings, while the final `local` part stops Nextflow from submitting another SLURM job from inside the first SLURM job.

| How Nextflow is started | Profile | What happens |
| --- | --- | --- |
| Directly on a permitted login or workflow node | `falcon` | Each Nextflow process becomes a separate SLURM job with its own requested resources. |
| Through `sbatch run_drugmr.sbatch` | `falcon,local` | The controller and all processes run inside the resources requested by that one allocation. |

The second mode is useful on clusters with a one node per user policy because nested jobs can remain pending while the controller holds the permitted node. Its tradeoff is that the single `sbatch` request must provide enough CPUs, memory, and time for the largest process that will run inside it.

### Option 1: Nextflow submits each task to SLURM

Run this from an HPC login node or workflow node when cluster policy permits a persistent Nextflow controller:

```bash
module load Java/17
cd /shared/scratch/YOUR_PROJECT/pipelines/drugMR
source env/activate.sh

nextflow run main.nf \
  -profile falcon \
  -params-file params/SCZ.wingo.yaml \
  --manifest_path assets/qtl_manifest.csv \
  --slurm_account YOUR_ACCOUNT \
  --slurm_partition htc_genoa \
  --container_bind /shared/scratch/YOUR_PROJECT/data \
  -resume
```

Nextflow remains the controller and submits each process as its own SLURM job.

### Option 2: run Nextflow inside one SLURM allocation

Use this option when Nextflow itself must be submitted with `sbatch`:

```bash
#!/bin/bash
#SBATCH --job-name=drugmr_SCZ_wingo
#SBATCH --partition=htc_genoa
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --output=drugMR/%x_%j.out
#SBATCH --error=drugMR/%x_%j.err
#SBATCH --account=YOUR_ACCOUNT

set -eo pipefail

module load Java/17
cd /shared/scratch/YOUR_PROJECT/pipelines/drugMR
source env/activate.sh

nextflow run main.nf \
  -profile falcon,local \
  -params-file params/SCZ.wingo.yaml \
  --manifest_path assets/qtl_manifest.csv \
  --slurm_account YOUR_ACCOUNT \
  --slurm_partition htc_genoa \
  --container_bind /shared/scratch/YOUR_PROJECT/data \
  -resume
```

Save this as `run_drugmr.sbatch`, then submit it from the directory above the repository:

```bash
cd /shared/scratch/YOUR_PROJECT/pipelines
sbatch run_drugmr.sbatch
```

The profile order matters. `falcon,local` allows the Falcon settings and Apptainer configuration to load first, then makes each process run inside the existing allocation. Using `falcon` alone inside an `sbatch` job submits nested SLURM jobs. It could also be noted as `<profile>,local` instead of falcon if you're running drugMR on a non-SLURM cluster with different container software.

Monitor the controller job:

```bash
squeue -j JOB_ID
tail -f drugMR/drugmr_SCZ_wingo_JOB_ID.out
```

Inspect resource use after completion:

```bash
sacct -j JOB_ID \
  --format=JobID,JobName,State,ExitCode,Elapsed,ReqCPUS,ReqMem,MaxRSS \
  --units=G
```

## Output from one Nextflow run

A successful run ID contains the outcome, pQTL dataset, date, and Git revision:

```text
SCZ_wingo_brain_20260912_fe5675a
```

The complete portable run is written under `runs/<run_id>/`:

```text
runs/<run_id>/
├── manifest.json
├── params.lock.yaml
├── pipeline_info/
│   ├── execution_report_<timestamp>.html
│   ├── execution_timeline_<timestamp>.html
│   ├── execution_trace_<timestamp>.txt
│   └── pipeline_dag_<timestamp>.svg
└── results/
    ├── cis_mr/
    │   ├── mr.tsv
    │   └── instruments/mr_instruments.tsv
    ├── coloc/
    ├── pwcoco/summary/
    ├── target_stats/top_cis_hits.tsv
    ├── locus_data/
    │   ├── index.tsv
    │   └── <target>/
    │       ├── gwas.parquet
    │       ├── pqtl.parquet
    │       ├── ld.parquet when the reference contains the locus variants
    │       └── metadata.json
    ├── smr/
    ├── hyprcoloc/
    └── phewas/
```

Some downstream files are optional because they depend on passing targets, enabled stages, and available reference variants. `pipeline_info` records task status, CPU use, memory use, runtime, and the executed DAG.

`params.lock.yaml` records the effective parameters used for that exact run. `manifest.json` records its identity and provenance. These files are created only after the complete workflow succeeds, so failed analyses are not registered as completed runs.

List the newest successful runs with:

```bash
ls -1dt runs/* | head
```

## Fetch an HPC run to the local clone

Run this on the local computer, from its own drugMR clone:

```python
import drugmr as dm

config = dm.fetch_run(
    run_id="SCZ_wingo_brain_20260912_fe5675a",
    user="your_username",
    host="login.your-cluster.ac.uk", # without the @ (Do not include the @!)
    remote_root="/shared/scratch/YOUR_PROJECT/pipelines/drugMR/runs",
)

dm.results(config=config)
```

`dm.fetch_run()` copies the selected run into the local `runs/` directory and returns its local `params.lock.yaml`. The dashboard then uses the copied results and portable locus files. The original QTL and GWAS datasets can remain on HPC.
