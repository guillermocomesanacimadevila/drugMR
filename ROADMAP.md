# drugMR engineering roadmap

Order top to bottom is priority. Sections split by when to tackle them, not by topic.

---

## Now

### 1. Commit today's cleanup
- [x] removed `tools/ukbppp_dl` submodule (not used anywhere in the actual pipeline, and it was broken: tracked as a gitlink but missing from `.gitmodules`, which is why `git status` kept failing)
- [x] fixed `launch.sh` (`python3.12 -v venv .venv` was wrong, should be `-m`)
- [x] updated `README.md` (`git clone --recurse-submodules`, notebook quickstart section)
- [x] staged untracking of `analysis/colocboost/*` so it actually respects the existing `analysis/*` gitignore rule going forward
- [ ] commit all of the above as its own commit, separate from the colocboost/sample-overlap work already sitting in the tree

### 2. docker-compose.yml for postgres
- wraps the native postgres setup `bin/load_db_into_postgres.py` already assumes, nothing more
- no design decisions needed here, it's just containerizing what already exists
- do this before the schema work below so migrations get tested against a throwaway container instead of the native install

---

## Next

### 3. SQL schema + ERD
- `load_db_into_postgres.py` right now just does `df.to_sql(table, engine, if_exists="replace")` off whatever columns happen to be in that day's tsv. no primary keys, no foreign keys, no constraints, table gets wiped and rebuilt from scratch every load
- draw the ERD first (miro is fine), the join key across `cis_mr_results` / `coloc_results` / hyprcoloc / phewas tables needs to be nailed down properly: `protein` x `pqtl_dataset` x `pheno_id`
- then write an actual migration, stop letting pandas guess column types every single load

### 4. R consolidation
- `cis_mr.R` and `genomewide_mr.R` duplicate `format_data`, `ld_clump`, `steiger_filtering`, and the IVW/Egger/WME dcast+rename block almost line for line
- pull the shared bits into one sourced helper `.R` file
- orchestration logic stays in `local.py` / `hpc.py`, not `drugmr/__init__.py`

### 5. HPC: submit one job per run, not subprocess calls scattered across hpc.py
- right now `hpc.py` fires off `subprocess.run(cmd, shell=True, ...)` stage by stage, straight from python, talking to SLURM directly inline
- instead: a single baseline submission script in `bin/` (something like `bin/submit_job.sh`) that `hpc.py` calls once per run. one run = one job, not a scatter of ad hoc subprocess calls
- take a `project_id` as a required param on submission, so every job on Falcon is attributable and traceable back to a run, instead of anonymous sbatch calls
- worth doing at the same time as item 6 below (shell=True cleanup), since this is exactly the code that needs it

### 6. generic QTL ingestion: stop hardcoding pQTL/eQTL datasets
- right now the pipeline only knows about 4 named pQTL cohorts (ukb_ppp, decode, wu_csf, wingo_brain) and a handful of named eQTL panels (MetaBrain, GTEx_v10, SingleBrain), each with its own ingestion code under `scripts/<cohort>/`
- goal: let a user bring any QTL dataset, as long as it's structured as parquet
  - pQTL dataset = one or more parquet files
  - eQTL dataset = one or a series of parquet files (one per tissue / cell type)
- user declares their own column names (snp/beta/se/p/a1/a2/chr/pos, same idea as the GWAS column mapping already in `params/schema.json`)
- build this as a **qtl manifest csv** under `assets/`, one row per dataset, that the user fills in themselves: dataset name, file path(s), column name mapping
- pipeline reads the manifest and handles any dataset generically from there, no more writing a new per-cohort script under `scripts/` every time someone wants to bring their own QTL data

---

## Parked (come back to this later)

### 7. CI, properly this time
- **git workflow**: still pushing straight to `main`, no PR flow, no branch protection. this is the actual root cause behind commits getting silently clobbered on a force-push, not something `ci.yml` alone fixes
- **ci.yml**: only runs `pytest`. none of the R scripts get touched at all (`cis_mr.R`, `coloc.R`, `hyprcoloc.R`, `moloc.R`), despite that being where the FRQ/MAF bug and the GRCh37/38 build mismatches actually lived
- two tiers to do this properly:
  - cheap: parse-check every `bin/*.R` script and confirm required packages load, run inside the existing `ghcr.io/guillermocomesanacimadevila/drugmr` image, no fixture data needed
  - real: toy fixture datasets (small per-protein pQTL slice, small GWAS slice, small plink bfile subset) so the R scripts actually run end to end in CI, not just parse

### 8. remove shell=True from subprocess.run()
- 13 call sites total: mostly `drugmr/local.py` and `drugmr/hpc.py`, plus `bin/coloc_targets.py`, `bin/load_db_into_postgres.py`, `bin/assort_network_mr.py`, `drugmr/smr.py`
- nothing here takes untrusted input right now so it's not an active exploit, but it's still silently trusting string-built shell commands and should get cleaned up
- do this together with item 5, since the hpc.py rewrite touches the exact same lines

### 9. output format
- runs/ + registry.json + params/ restructure is done (2026-08-11): `dat/derived/` for shared preprocessing, `runs/<run_id>/` per run with `manifest.json` + `registry.json`, `synthesis/` for cross-dataset rollups, `params/` schema-validated with a `gates:` block replacing the old hardcoded thresholds. dashboard reads through the registry now, not path-guessing
- what never got done in that pass: outputs are still `.tsv`, not `.parquet`/`.gz`

### 10. upload dat/ref onto Dropbox or Zenodo
- reference data doesn't belong committed to the repo and right now it only exists locally

---

## Skip, or revisit only if it actually becomes painful

### 11. workflow engine (Nextflow or Snakemake)
- original plan here was Nextflow, still nothing Nextflow-shaped anywhere in the repo
- `local.py` / `hpc.py` already give you a chunk of what a workflow engine buys: output-existence checks before rerunning a stage, and the `.running.tsv` incremental-save pattern in `cis_mr.R` is basically a poor man's DAG cache
- if this ever does get picked up: snakemake over nextflow. it's python native, and its wildcard system maps almost directly onto the `results/{stage}/{pqtl_dataset}/...` naming already in use. rules can just `shell:` call the existing `bin/*.R` / `bin/*.py` scripts unmodified
- only worth doing if SLURM submission through `hpc.py` actually turns into a real bottleneck, not before
