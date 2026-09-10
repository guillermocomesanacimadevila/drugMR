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

drugMR takes an outcome GWAS and a panel of protein QTLs and returns a ranked, safety screened shortlist of druggable targets, end to end with no manual steps between stages: Mendelian randomisation, colocalisation, SMR, HyPrColoc, then a phenome wide PheWAS safety screen. It runs via a Nextflow pipeline or a Python orchestrator, against any outcome GWAS and any pQTL cohort registered in its dataset manifest, bulk or single cell. The downstream SMR, PWCoCo and HyPrColoc stages are generic to any QTL type, eQTL, sQTL, mQTL or otherwise, for genuine multi-omics triangulation.

## Quick start

1. Install [`Nextflow`](https://www.nextflow.io/) (`>=26.04.0`)
2. Install any of [`Docker`](https://docs.docker.com/engine/install/), [`Apptainer`](https://apptainer.org/) or [`Singularity`](https://sylabs.io/singularity/)

```bash
git clone --recurse-submodules https://github.com/guillermocomesanacimadevila/drugMR.git

nextflow run /path/to/cloned/drugMR/main.nf -profile docker -params-file params/AD.ukb_ppp.yaml --manifest_path assets/qtl_manifest.csv
```

## Run on your own data

Every `(pheno_id, pqtl_dataset)` pair gets its own params file under `params/`, validated against `params/schema.json`.

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

`pqtl_dataset`, `bulk_qtl_datasets` and `sc_qtl_dataset` must match a dataset registered in `assets/qtl_manifest.csv`. Registering a new one is a new row there, parquet, csv, tsv or txt all work, no code changes.

```bash
nextflow run main.nf \
  -profile <docker,apptainer,singularity,podman,shifter,charliecloud> \
  -params-file <path/to/params.yaml> \
  --manifest_path <path/to/qtl_manifest.csv>
```

or via Python:

```python
import drugmr as dm

dm.local(config="params/AD.ukb_ppp.yaml")
dm.hpc(config="params/AD.ukb_ppp.yaml", falcon_user="your_username")
dm.results(config="params/AD.ukb_ppp.yaml")
```

Postgres loading and dashboard serving are deliberately not a Nextflow stage. Run `dm.results()` afterwards regardless of which entry point produced the run.

## Pipeline summary

* GWAS QC and cis region extraction
* Cis-MR (Wald ratio, inverse variance weighted)
* Pairwise colocalisation and PWCoCo
* SMR and HEIDI, bulk and single cell QTL panels
* HyPrColoc
* PheWAS (FinnGen, UK Biobank)
* PostgreSQL and an interactive Streamlit dashboard

Gate thresholds live in the optional `gates` block of your params file. See [`docs/RESULTS_SCHEMA.md`](docs/RESULTS_SCHEMA.md) for the results schema.

## Credits

drugMR was written by:

**Guillermo Comesaña Cimadevila**<sup>1,2</sup>, **Marie-Joe Dib**<sup>3</sup>, **Matthew Bracher-Smith**<sup>1</sup>, **Christian Pepler**<sup>2</sup>, **Dervis Salih**<sup>4</sup>, **Nicholas Bray**<sup>2</sup>, **Emily Simmonds**<sup>1</sup>, **Valentina Escott-Price**<sup>1,2</sup>

<sup>1</sup> UK Dementia Research Institute at Cardiff University, Cardiff, UK
<sup>2</sup> MRC Centre for Neuropsychiatric Genetics and Genomics, Cardiff University, Cardiff, UK
<sup>3</sup> Nascent Studio Ltd, London, UK
<sup>4</sup> UK Dementia Research Institute at University College London, London, UK

## Contributions and support

Contributions are welcome. Fork the repository, open a pull request, and flag it to Guillermo (ComesanaCimadevilaG@cardiff.ac.uk).

## Citations

If you use drugMR in your work, please cite it. A machine readable citation is provided in [`CITATION.cff`](CITATION.cff). An extensive list of references for every tool this pipeline depends on is in [`CITATIONS.md`](CITATIONS.md).

## License

Released under the [MIT License](LICENSE).
