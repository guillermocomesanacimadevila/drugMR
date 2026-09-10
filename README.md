![drugMR logo](docs/new_logo.png)

**drugMR: A multi-omics pipeline for genetically anchored drug target discovery**

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![Nextflow](https://img.shields.io/badge/nextflow-%E2%89%A526.04.0-23aa62?logo=nextflow&logoColor=white)](nextflow.config)
[![Python 3.12+](https://img.shields.io/badge/python-3.12%2B-blue?logo=python&logoColor=white)](pyproject.toml)
[![R](https://img.shields.io/badge/R-4%2B-blue?logo=r&logoColor=white)](env/Dockerfile)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://github.com/guillermocomesanacimadevila/drugMR/pkgs/container/drugmr)
[![run with apptainer](https://img.shields.io/badge/run%20with-apptainer-1d355c.svg?labelColor=000000)](https://apptainer.org/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/singularity/)
[![PostgreSQL](https://img.shields.io/badge/PostgreSQL-results%20store-blue?logo=postgresql&logoColor=white)](sql/schema.sql)
[![Streamlit](https://img.shields.io/badge/Streamlit-dashboard-blue?logo=streamlit&logoColor=white)](dashboard/mr_app.py)

## Contents

- [Introduction](#introduction)
- [Pipeline summary](#pipeline-summary)
- [Installation](#installation)
- [Usage](#usage)
- [Configuration](#configuration)
- [Runs, registry and synthesis](#runs-registry-and-synthesis)
- [Pipeline output](#pipeline-output)
- [Repository layout](#repository-layout)
- [Streamlit configuration](#streamlit-configuration)
- [Credits](#credits)
- [Citations](#citations)
- [License](#license)

## Introduction

**drugMR** takes an outcome GWAS and a panel of protein QTLs and returns a ranked, safety screened shortlist of druggable targets, run end to end with no manual intervention between stages. Each stage runs independently and passes forward only the targets that survived the stage before it: Mendelian randomisation, colocalisation, SMR, HyPrColoc, then a phenome wide PheWAS safety screen. It is demonstrated here on Alzheimer's disease but works against any outcome phenotype and any pQTL cohort registered in its dataset manifest. The pipeline ships with two independent, fully working ways to run it: a Nextflow pipeline and a Python orchestrator. Both write to the same `runs/<run_id>/` layout and share the same run registry, so results from either one are interchangeable.

## Pipeline summary

![drugMR pipeline DAG](docs/pipeline.png)

1. **GWAS QC**: harmonises and QCs the outcome GWAS.
2. **cis region prep**: matches pQTL cis regions to the outcome GWAS.
3. **cis-MR**: Wald ratio or inverse variance weighted MR per protein.
4. **Pairwise COLOC and PWCoCo**: colocalisation between the pQTL and the outcome GWAS, plus a conditional coloc check for loci with more than one causal signal.
5. **Top cis hit compilation**: aligns the top cis SNP per protein to the outcome risk allele.
6. **SMR and PWCoCo QTL**: eQTL to GWAS colocalisation via SMR and HEIDI, across bulk and single cell eQTL panels, plus SNP level pQTL to eQTL to GWAS triangulation.
7. **HyPrColoc**: clusters pQTL, GWAS and eQTL signals to confirm a shared causal variant.
8. **PheWAS**: FinnGen and UK Biobank phenome wide MR safety screen of the surviving targets.
9. **Results**: loads everything into PostgreSQL and serves it through a Streamlit dashboard.

Every gate threshold above (FDR cutoffs, `PP.H4.abf`, `q_SMR`, `p_HEIDI`) lives in the optional `gates` block of your params file, not hardcoded in the scripts. Completed stages are cached per run under `runs/<run_id>/results/` and reused unless `overwrite: true` is set.

## Installation

The pipeline itself never needs R, PostgreSQL, PLINK, GCTA, SMR, or PWCoCo installed on your machine directly. All of that lives inside the pipeline's own container image and is pulled automatically the first time you run it. What you do need on your own machine depends on which entry point you use.

### 1. Install Nextflow (version 26.04.0 or later)

```bash
curl -s https://get.nextflow.io | bash
chmod +x nextflow
sudo mv nextflow /usr/local/bin/
nextflow -version
```

Version 26.04.0 or later is required, not just recommended: `main.nf` uses the native `onComplete` section inside the entry workflow, which needs the strict syntax parser that Nextflow made its default behaviour from that version onward.

### 2. Install a container engine

At least one of the following, matching whichever profile you plan to use:

```bash
# Docker, for local runs
# see https://docs.docker.com/engine/install/ for your platform

# OR Apptainer, for HPC runs
# see https://apptainer.org/docs/admin/main/installation.html
```

### 3. Install Python 3.12 or later

```bash
python3 --version
```

Needed for the `drugmr` package itself and, if you plan to use it, the Python orchestrator (`dm.local()` / `dm.hpc()`).

### 4. PostgreSQL and R: no separate install needed

PostgreSQL 16 is started automatically via `docker compose up -d` (see `docker-compose.yml`) when you run `dm.results()`. R 4.4, along with PLINK, GCTA, SMR, and PWCoCo, is baked into the pipeline's own Docker image (`env/Dockerfile`) and only ever runs inside a container. Neither needs a manual install on your host.

### 5. Clone the repository and install the Python package

```bash
git clone --recurse-submodules https://github.com/guillermocomesanacimadevila/drugMR.git
cd drugMR
pip install -e .
```

Already cloned without `--recurse-submodules`? Run this from inside the repo instead of re-cloning:

```bash
git submodule update --init --recursive
```

### Quickstart (notebook)

```bash
chmod +x launch.sh
bash launch.sh
```

First run creates a `.venv`, installs `drugmr` plus Jupyter, and opens `notebooks/00_drugmr.ipynb` in Jupyter Lab. Every run after that just reactivates the same `.venv` and reopens the notebook.

## Usage

### Nextflow

```bash
nextflow run main.nf \
    -profile docker \
    -params-file params/AD.ukb_ppp.yaml \
    --manifest_path assets/qtl_manifest.csv
```

To resume a previous run instead of starting from scratch:

```bash
nextflow run main.nf \
    -profile docker \
    -params-file params/AD.ukb_ppp.yaml \
    --manifest_path assets/qtl_manifest.csv \
    -resume
```

`--manifest_path` points at the dataset registry (see `assets/qtl_manifest.csv`) and defaults to that same file, so it can be omitted for the standard registry or overridden to point at a different one entirely.

`-profile` picks the container engine. There is no default engine, you must pick one explicitly:

```
docker
apptainer
singularity
podman
shifter
charliecloud
```

An executor profile can be combined with it using a comma, for example:

```bash
nextflow run main.nf \
    -profile local,docker \
    -params-file params/AD.ukb_ppp.yaml \
    --manifest_path assets/qtl_manifest.csv
```

or, on Falcon, once the `falcon` profile is filled in with your own SLURM account (`process.queue`, `process.clusterOptions`, see `nextflow.config`):

```bash
nextflow run main.nf \
    -profile falcon,apptainer \
    -params-file params/AD.ukb_ppp.yaml \
    --manifest_path assets/qtl_manifest.csv
```

Postgres loading and dashboard serving are deliberately not a Nextflow stage. Run `dm.results(config=...)` afterwards, exactly as in the Python path below. If the run happened on a remote machine, pull `runs/<run_id>` back first:

```bash
rsync -avz <host>:<path>/runs/<run_id> ./runs/
```

Then run `dm.results()` locally as normal.

### Python

```python
import drugmr as dm

# run locally via Docker
dm.local(config="params/AD.ukb_ppp.yaml")

# or run on the Falcon HPC cluster via SLURM and Apptainer
dm.hpc(config="params/AD.ukb_ppp.yaml", falcon_user="your_falcon_username")

# load that run's results into PostgreSQL and launch the Streamlit dashboard
dm.results(config="params/AD.ukb_ppp.yaml")
```

See [`notebooks/00_drugmr.ipynb`](notebooks/00_drugmr.ipynb) for a worked example.

## Configuration

There is no single `config.yaml`. Every `(pheno_id, pqtl_dataset)` pair gets its own params file under `params/`, for example `params/AD.ukb_ppp.yaml` or `params/AD.wingo_brain.yaml`. `drugmr.config.Config` validates whatever you point it at against `params/schema.json`, so a mistyped column name fails immediately, not three hours into cis-MR. The Nextflow entry point reads the same file shape via `-params-file`, so one params file works for both entry points.

| Field | Purpose |
| --- | --- |
| `pheno_id`, `sumstats`, `n_cases`, `n_controls` | Outcome GWAS identity and sample size |
| `genome_build`, `target_build` | Source and target genome builds (liftover if they differ) |
| `snp_col` / `a1_col` / `a2_col` / `beta_col` / `se_col` / `p_col` / `pos_col` / `chr_col` / `af_col` | Column names in your outcome GWAS |
| `pqtl_dataset` | Which pQTL dataset to run, matching a dataset registered in `assets/qtl_manifest.csv` (add a row there to register a new one, no code changes needed) |
| `ref_bfile` | Reference panel (1000 Genomes) for cis-MR and SMR |
| `run_smr` | Master switch for the SMR step |
| `bulk_eqtl_datasets` | Pre-computed bulk eQTL datasets to ingest, for example `[MetaBrain, GTEx_v10]`. An empty list skips bulk SMR |
| `sc_eqtl_dataset` | Single cell eQTL dataset to run SMR against, for example `SingleBrain`. Empty skips single cell SMR |
| `maf`, `remove_mhc`, `remove_apoe` | QC filters applied to the GWAS and pQTLs |
| `overwrite` | Force every stage to rerun instead of reusing existing outputs |
| `gates` | Optional block of per stage statistical thresholds, see below |

```yaml
gates:
  cis_mr:
    wald_fdr_q: 0.05
    ivw_fdr_q: 0.05
    cochran_q_pval: 0.05
    egger_intercept_pval_min: 0
    min_instruments_for_ivw: 3
  coloc:
    pp4_threshold: 0.7
```

Omit the `gates` block entirely and the values above apply as defaults. HyPrColoc has no equivalent pipeline stage gate: its posterior probability threshold is applied interactively in the dashboard instead, so results can be explored at more than one threshold without rerunning the pipeline.

## Runs, registry and synthesis

Every run, from either entry point, gets its own stamped directory:

```
runs/AD_ukb_ppp_20260811_149fc55/
├── results/          # Everything cis-MR, COLOC, SMR and PheWAS wrote for this run
├── manifest.json     # What ran, when, against which commit (Python path only)
└── params.lock.yaml  # A frozen copy of the params file used (Python path only)
```

`runs/registry.json` maps `{pheno_id}__{pqtl_dataset}` to the latest successful run and its full history. It is only updated once every step of a run has actually succeeded, so it can never point at a half finished run. Preprocessing that is shared across every run for a given phenotype lives once in `dat/derived/<pheno_id>/`. Once you have run more than one pQTL dataset for the same phenotype, `synthesis/<pheno_id>/` holds the cross dataset target roll up.

## Pipeline output

`dm.results()` loads a run's results into PostgreSQL and launches the Streamlit dashboard (`dashboard/mr_app.py`), which includes a pQTL dataset selector, a run history picker, and sidebar filters for outcome, FDR/Q/PP.H4 thresholds, and protein search.

| Page | Contents |
| --- | --- |
| Overview | Target prioritisation funnel and pipeline stage guide |
| Target Profile | Full evidence trail for a single target, including a LocusZoom style regional plot |
| Evidence by Stage | Full per stage results tables, one sub tab per pipeline stage |
| Final Targets | Curated deliverable of targets passing every stage, with a Sankey diagram of the full branching |
| PWCoCo | Conditional colocalisation results for loci with more than one causal signal |
| PWCoCo QTL | SNP level triangulation across pQTL, eQTL and GWAS |

See also [`docs/RESULTS_SCHEMA.md`](docs/RESULTS_SCHEMA.md) for the exact schema of every results file.

## Repository layout

```
drugMR/
├── drugmr/              # Installable package: Config, paths, registry, SMR, PheWAS, PyTwoSampleMR, utils
├── bin/                 # Pipeline stage scripts (Python and R), invoked by both entry points
├── modules/local/       # Nextflow process definitions, one .nf file per stage
├── subworkflows/        # Nextflow subworkflows chaining modules together, one directory per stage
├── main.nf              # Nextflow entry point
├── nextflow.config      # Nextflow params, profiles, resource labels
├── conf/                # Nextflow resource config (base.config)
├── params/              # One params.yaml per (pheno_id, pqtl_dataset), plus schema.json
├── dat/                 # Input data: GWAS, pQTL, sc-eQTL, cis regions, reference panel
│   └── derived/         # Shared preprocessing reused across runs for a pheno_id
├── runs/                # One folder per run, plus registry.json
├── synthesis/           # Cross dataset roll ups per pheno_id
├── dashboard/           # Streamlit app (mr_app.py)
├── notebooks/           # Worked examples (00_drugmr.ipynb)
├── assets/              # qtl_manifest.csv and other run adjacent bits
├── env/                 # Dockerfile, requirements.txt
├── tools/               # Git submodules (pwcoco)
├── tests/               # Toy fixtures and per module Nextflow test configs
├── docs/                # Pipeline DAG, results schema
```

## Streamlit configuration

The dashboard reads results from PostgreSQL via `.streamlit/secrets.toml`, which is gitignored and regenerated automatically by `dm.results()`. The generated connection points at the Postgres instance started by `docker compose up -d` (user `drugmr_user`, database `drugmr`, port `5433`, see `docker-compose.yml`). Edit the constants at the top of `bin/write_streamlit_secrets.py` if your setup differs, or regenerate the file on its own:

```bash
python3 bin/write_streamlit_secrets.py
```

## Credits

drugMR was written by:

**Guillermo Comesaña Cimadevila**<sup>1,2</sup>, **Christian Pepler**<sup>2</sup>, **Marie-Joe Dib**<sup>3</sup>, **Dervis Salih**<sup>4</sup>, **Nicholas J. Bray**<sup>2</sup>, **Emily Simmonds**<sup>1</sup>, **Valentina Escott-Price**<sup>1,2</sup>

<sup>1</sup> UK Dementia Research Institute at Cardiff University, Cardiff, UK
<sup>2</sup> MRC Centre for Neuropsychiatric Genetics and Genomics, Cardiff University, Cardiff, UK
<sup>3</sup> Nascent Studio Ltd, London, UK
<sup>4</sup> UK Dementia Research Institute at University College London, London, UK

## Citations

If you use drugMR in your work, please cite it. A machine readable citation is provided in [`CITATION.cff`](CITATION.cff).

An extensive list of references for every statistical method and tool this pipeline depends on (TwoSampleMR, coloc, HyPrColoc, SMR, GCTA, PLINK, Nextflow, and the container engines) is in [`CITATIONS.md`](CITATIONS.md).

## License

Released under the [MIT License](LICENSE).
