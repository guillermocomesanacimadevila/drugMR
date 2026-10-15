# Changelog

All notable changes to drugMR are listed here. Versions follow [semantic versioning](https://semver.org).

## [Unreleased]

## [1.0.0] (2026-10-15)

First release. Container image `ghcr.io/guillermocomesanacimadevila/drugmr:v1.0.0`, pinned by digest in `nextflow.config`.

### Added

- Nextflow DSL2 pipeline: GWAS QC, cis regions, cis-MR, coloc, PWCoCo, PheWAS (FinnGen with UK Biobank fallback), bulk and single cell SMR and HEIDI, QTL PWCoCo and HyPrColoc.
- Profiles for local, SLURM, AWS Batch and GCP Batch, with Docker, Apptainer, Singularity, Podman, Shifter and Charliecloud.
- QTL manifest for registering pQTL, eQTL, sQTL and single cell QTL datasets.
- Run registry under `runs/<run_id>/` with locked parameters and provenance.
- PostgreSQL loader and Streamlit dashboard.
- Synthetic demo dataset (`-profile demo`).
- `ResultsExplorer` (`drugmr.interpret`): read a completed run's results tables directly in Python, without PostgreSQL or the dashboard, from the run overview to an evidence summary per target, with forest and regional plots.
- `notebooks/01_interpreting_results.ipynb`: a guided tour of every `ResultsExplorer` method on the synthetic demo run, with outputs.
- Documentation site, `CONTRIBUTING.md`, and issue templates for bug reports and feature requests.

[Unreleased]: https://github.com/guillermocomesanacimadevila/drugMR/compare/v1.0.0...HEAD
[1.0.0]: https://github.com/guillermocomesanacimadevila/drugMR/releases/tag/v1.0.0
