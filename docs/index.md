# drugMR

## Overview

**drugMR** is a multi-omics pipeline for genetically anchored drug target discovery. It takes an outcome GWAS and a panel of protein QTLs and produces a ranked, safety-screened shortlist of druggable targets.

The pipeline is built using [Nextflow](https://www.nextflow.io), a workflow tool to run tasks across multiple compute infrastructures in a portable manner. It uses Docker/Apptainer/Singularity containers for reproducible execution, and connects Mendelian randomisation, colocalisation, SMR and HEIDI, PWCoCo, HyPrColoc, and phenome-wide safety screening into a single reproducible run.

## Documentation

- [Install](getting-started.md): set up the pinned environment and choose a container runtime.
- [Configure a run](configuration.md): define an outcome GWAS, register QTL datasets, and set analysis thresholds.
- [Run the pipeline](running.md): launch and monitor an analysis.
- [Results and dashboard](results-dashboard.md): fetch completed runs, load PostgreSQL, and explore the Streamlit dashboard.
- [Output schema](RESULTS_SCHEMA.md): column reference for every results TSV.
- [Troubleshooting](troubleshooting.md): diagnose SLURM, memory, paths, containers, and interrupted runs.

## Workflow

![drugMR analysis pipeline](pipeline.png)

1. Quality-control the outcome GWAS and extract cis regions.
2. Run cis-MR and pairwise colocalisation against the selected pQTL panel.
3. Triangulate promising targets with registered bulk and single-cell QTL datasets.
4. Run PheWAS safety screens using FinnGen and UK Biobank data.
5. Load a completed run into PostgreSQL and explore it in the Streamlit dashboard.

## Choose where each part runs

The computational pipeline can run locally or through SLURM. The results dashboard runs on whichever machine calls `dm.results()`.

| Situation | What to do |
| --- | --- |
| Pipeline and dashboard on one machine | Run Nextflow, then pass the run's `params.lock.yaml` to `dm.results()`. |
| Pipeline on HPC, dashboard on your computer | Run Nextflow on HPC, assign the value returned by `dm.fetch_run()` to `config`, then call `dm.results(config=config)` locally. |
| Pipeline and dashboard on HPC | Pass the run's `params.lock.yaml` to `dm.results()` on HPC; use SSH port forwarding to view Streamlit. |

```{note}
`dm.fetch_run()` is only needed when the Nextflow run is on a remote HPC or cloud machine and the dashboard will run elsewhere. It transfers the run's locked parameter snapshot along with its results.

`dm.results()` always needs a configuration. For an exact local run, use `dm.results(config="runs/<run_id>/params.lock.yaml")`. To select the latest successful run matching a regular parameter file, use `dm.results(config="params/<file>.yaml")`.
```

```{toctree}
:maxdepth: 2
:caption: Getting started
:hidden:

getting-started
configuration
running
```

```{toctree}
:maxdepth: 2
:caption: Reference
:hidden:

results-dashboard
RESULTS_SCHEMA
troubleshooting
```
