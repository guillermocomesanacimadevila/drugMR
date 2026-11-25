# Changelog

All notable changes to drugMR are listed here. Versions follow [semantic versioning](https://semver.org).

## [Unreleased]

## [1.1.0] (2026-11-25)

Container image `ghcr.io/guillermocomesanacimadevila/drugmr:v1.1.0`, pinned by digest in `nextflow.config`. It is the v1.0.0 image with two Python packages added (`pymysql`, `chembl_webresource_client`); every other package version is unchanged.

### Added

- Locus conditioning: for each target passing cis-MR and colocalisation, the outcome GWAS is conditioned on other independent signals in the cis window with GCTA-COJO, and cis-MR and coloc are rerun. Thresholds sit in a new `gates.locus_conditioning` block.
- Druggability as the final stage: every target passing cis-MR and colocalisation is looked up in Pharos, ChEMBL, Open Targets and B3DB, and each drug's direction is compared with the MR. The potency threshold sits in a new `gates.druggability` block, and `--skip_druggability true` turns the stage off on offline compute nodes.

### Changed

- The `drugmr` package is split into `methods/`, `io/` and `runners/` subpackages. `import drugmr as dm` works as before.

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

[Unreleased]: https://github.com/guillermocomesanacimadevila/drugMR/compare/v1.1.0...HEAD
[1.1.0]: https://github.com/guillermocomesanacimadevila/drugMR/releases/tag/v1.1.0
[1.0.0]: https://github.com/guillermocomesanacimadevila/drugMR/releases/tag/v1.0.0
