![drugMR logo](docs/new_logo.png)

**drugMR: a multi omics pipeline for genetically anchored drug target discovery**

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![Nextflow](https://img.shields.io/badge/nextflow-%E2%89%A526.04.0-23aa62?logo=nextflow&logoColor=white)](nextflow.config)
[![Python 3.12+](https://img.shields.io/badge/python-3.12%2B-blue?logo=python&logoColor=white)](pyproject.toml)
[![R](https://img.shields.io/badge/R-4%2B-blue?logo=r&logoColor=white)](env/Dockerfile)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://github.com/guillermocomesanacimadevila/drugMR/pkgs/container/drugmr)
[![run with apptainer](https://img.shields.io/badge/run%20with-apptainer-1d355c.svg?labelColor=000000)](https://apptainer.org/)
[![PostgreSQL](https://img.shields.io/badge/PostgreSQL-results%20store-blue?logo=postgresql&logoColor=white)](sql/schema.sql)
[![Streamlit](https://img.shields.io/badge/Streamlit-dashboard-blue?logo=streamlit&logoColor=white)](dashboard/mr_app.py)

## Introduction

drugMR takes an outcome GWAS and a panel of protein QTLs and returns a ranked, safety screened shortlist of druggable targets, end to end with no manual steps between stages. Each stage passes forward only the targets that survived the one before it: Mendelian randomisation, colocalisation, SMR, HyPrColoc, then a phenome wide PheWAS safety screen.

drugMR supports:

* Any outcome GWAS and any pQTL cohort registered in its dataset manifest, no code changes needed to add a new one
* Bulk and single cell QTL panels (eQTL today, the manifest and every downstream stage are generic to any QTL type)
* Two independent, fully working entry points: a Nextflow pipeline and a Python orchestrator, writing to the same run layout and sharing the same registry
* Local execution via Docker, or HPC execution via SLURM and Apptainer over SSH
* Per stage statistical thresholds set entirely from your params file, nothing hardcoded in the scripts

drugMR is not tied to Alzheimer's disease, Cardiff/UKDRI own compute, or any single pQTL provider. Every threshold, every dataset, and every compute target is a parameter, not an assumption baked into the code. Developed by Guillermo Comesaña Cimadevila under the leadership of Valentina Escott-Price at the UK Dementia Research Institute, this pipeline was built for a genome wide Alzheimer's disease drug target search but is general enough for any GWAS backed target discovery project.

## Quick start

1. Install [Nextflow](https://www.nextflow.io/) (`>=26.04.0`)
2. Install any of [Docker](https://docs.docker.com/engine/install/), [Apptainer](https://apptainer.org/), [Singularity](https://sylabs.io/singularity/), Podman, Shifter or Charliecloud for full pipeline reproducibility. R, PostgreSQL, PLINK, GCTA, SMR and PWCoCo never need installing directly, they all live inside the pipeline's own container image, pulled automatically on first run.

Download the pipeline and run it against the bundled example dataset with a single command:

```bash
git clone --recurse-submodules https://github.com/guillermocomesanacimadevila/drugMR.git

nextflow run /path/to/cloned/drugMR/main.nf -profile docker -params-file params/AD.ukb_ppp.yaml
```

## Run on your own data

Every `(pheno_id, pqtl_dataset)` pair gets its own params file under `params/`. `drugmr.config.Config` validates it against `params/schema.json`, so a mistyped column name fails immediately rather than hours into cis-MR.

```yaml
pheno_id: AD
sumstats: dat/gwas/AD_bellenguez.tsv
n_cases: 39106
n_controls: 401577
pqtl_dataset: ukb_ppp
ref_bfile: dat/ref/1000G_EUR_Phase3_plink/1000G.EUR.QC.ALL
snp_col: SNP
a1_col: A1
a2_col: A2
beta_col: BETA
se_col: SE
p_col: P
pos_col: BP
chr_col: CHR
af_col: FRQ
genome_build: GRCh38
target_build: GRCh38
run_smr: true
bulk_qtl_datasets: [MetaBrain, GTEx_v10]
sc_qtl_dataset: SingleBrain
```

`pqtl_dataset` and every entry in `bulk_qtl_datasets` / `sc_qtl_dataset` must match a dataset registered in `assets/qtl_manifest.csv`. Registering a new one is a new row in that manifest (parquet, csv, tsv or txt all work), no code changes.

Run it via Nextflow:

```bash
nextflow run main.nf \
  -profile <docker,apptainer,singularity,podman,shifter,charliecloud> \
  -params-file <path/to/params.yaml>
```

for example:

```bash
nextflow run main.nf \
  -profile docker \
  -params-file params/AD.ukb_ppp.yaml
```

`-profile` can be combined with an executor profile, for example `-profile local,docker` on your own machine or, on a SLURM cluster with Apptainer:

```bash
nextflow run main.nf \
  -profile falcon,apptainer \
  -params-file params/AD.ukb_ppp.yaml
```

or resume a previous run instead of starting fresh:

```bash
nextflow run main.nf \
  -profile docker \
  -params-file params/AD.ukb_ppp.yaml \
  -resume
```

Or run it via the Python orchestrator:

```python
import drugmr as dm

dm.local(config="params/AD.ukb_ppp.yaml")

dm.hpc(
    config="params/AD.ukb_ppp.yaml",
    falcon_user="your_username",
    host="falconlogin.cf.ac.uk",       # any SLURM + Apptainer cluster works
    remote_repo_root="/shared/home1/{falcon_user}/drugMR",
)

dm.results(config="params/AD.ukb_ppp.yaml")
```

Postgres loading and dashboard serving are deliberately not a Nextflow stage. Run `dm.results(config=...)` afterwards regardless of which entry point produced the run. If it ran on a machine other than the one you want the dashboard on, pull the run across first:

```python
dm.fetch_run(run_id, host="your_cluster.ac.uk", remote_root="/path/to/drugMR/runs")
```

See [`notebooks/00_drugmr.ipynb`](notebooks/00_drugmr.ipynb) for a worked example.

## Pipeline summary

* GWAS QC and cis region extraction against the pQTL panel
* Cis-MR: Wald ratio or inverse variance weighted, per protein
* Pairwise colocalisation (coloc.abf) and PWCoCo for loci with more than one causal signal
* SMR and HEIDI against bulk and single cell QTL panels, plus SNP level pQTL to QTL to GWAS triangulation
* HyPrColoc, clustering pQTL, GWAS and QTL signals onto a shared causal variant
* PheWAS: FinnGen and UK Biobank phenome wide safety screen of the surviving targets
* Results loaded into PostgreSQL and served through an interactive Streamlit dashboard

Gate thresholds (FDR cutoffs, `PP.H4.abf`, `q_SMR`, `p_HEIDI`, and the rest) live in the optional `gates` block of your params file. Completed stages are cached per run under `runs/<run_id>/results/` and reused unless `overwrite: true` is set. See [`docs/RESULTS_SCHEMA.md`](docs/RESULTS_SCHEMA.md) for the exact schema of every results file.

## Credits

drugMR was written by:

**Guillermo Comesaña Cimadevila**<sup>1,2</sup>, **Marie-Joe Dib**<sup>3</sup>, **Matthew Bracher-Smith**<sup>1</sup>, **Christian Pepler**<sup>2</sup>, **Dervis Salih**<sup>4</sup>, **Nicholas Bray**<sup>2</sup>, **Emily Simmonds**<sup>1</sup>, **Valentina Escott-Price**<sup>1,2</sup>

<sup>1</sup> UK Dementia Research Institute at Cardiff University, Cardiff, UK
<sup>2</sup> MRC Centre for Neuropsychiatric Genetics and Genomics, Cardiff University, Cardiff, UK
<sup>3</sup> Nascent Studio Ltd, London, UK
<sup>4</sup> UK Dementia Research Institute at University College London, London, UK

The pipeline is built on Nextflow, runs one container per process for reproducibility, and follows nf-core structural conventions throughout.

## Contributions and support

Contributions are welcome. Fork the repository, open a pull request, and flag it to Guillermo (ComesanaCimadevilaG@cardiff.ac.uk).

## Citations

If you use drugMR in your work, please cite it. A machine readable citation is provided in [`CITATION.cff`](CITATION.cff). An extensive list of references for every statistical method and tool this pipeline depends on (TwoSampleMR, coloc, HyPrColoc, SMR, GCTA, PLINK, Nextflow, and the container engines) is in [`CITATIONS.md`](CITATIONS.md).

## License

Released under the [MIT License](LICENSE).
