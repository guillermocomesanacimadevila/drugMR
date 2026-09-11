# Results and dashboard

Use Python for the post-pipeline steps:

```python
import drugmr as dm
```

## Fetch a run from HPC or cloud

Skip this step when the completed run already exists on the machine where the dashboard will run.

```python
config = dm.fetch_run(
    run_id="AD_test_20260911_abcdef0",
    user="your_username",
    host="login.your-cluster.ac.uk",
    remote_root="/path/to/drugMR/runs",
)
```

`fetch_run()` copies `remote_root/<run_id>/` into the local `runs/` directory with `rsync`. This includes the results, pipeline reports, `manifest.json`, and `params.lock.yaml`. It validates those metadata files, registers the run locally, and returns the local path to the locked parameter snapshot for `dm.results()`.

The `remote_root` value ends at the remote `runs` directory. `fetch_run()` appends the supplied `run_id` itself.

### Prefer the command line for SSH prompts

If SSH requests a password or private-key passphrase, run the fetch in a local terminal. Notebook interfaces do not always expose an interactive SSH prompt reliably.

```bash
python -c 'import drugmr as dm; config = dm.fetch_run(
    run_id="AD_test_20260911_abcdef0",
    user="your_username",
    host="login.your-cluster.ac.uk",
    remote_root="/path/to/drugMR/runs"
); print(config)'
```

SSH keys avoid repeated password prompts and work well with `rsync`.

## Load results and open the dashboard

`dm.results()` always requires a `config` argument. Choose one of the following forms according to where the pipeline ran.

### After fetching a remote run

```python
dm.results(config=config)
```

Here, `config` is the path returned by `dm.fetch_run()`.

### After a local Nextflow run

To open one exact local run, point directly to its saved snapshot:

```python
dm.results(config="runs/AD_test_20260911_abcdef0/params.lock.yaml")
```

To open the latest successful run matching a phenotype and pQTL dataset, use the original parameter file:

```python
dm.results(config="params/AD.test.yaml")
```

Do not call bare `dm.results()`: there is no reliable way to know which phenotype, pQTL dataset, or historical run the user intends to open.

This command:

1. Reads the run identity and settings from the supplied configuration.
2. Selects the exact run when given `params.lock.yaml`. A regular `params/*.yaml` instead selects the latest successful run matching its `pheno_id` and `pqtl_dataset`.
3. Starts PostgreSQL with `docker compose up -d`.
4. Writes the Streamlit database connection settings.
5. Loads the run's MR and colocalisation results into PostgreSQL.
6. Starts the Streamlit dashboard.

Run it locally after `fetch_run()` when you want the dashboard on your computer. You may instead run it directly on the HPC if Docker and the dashboard environment are available there. To view an HPC-hosted dashboard, forward its Streamlit port through SSH according to your cluster's access rules.

!!! warning
    Fetch only completed runs. The successful Nextflow controller creates `manifest.json` and `params.lock.yaml` and records the run only after all required stages finish. Older runs created before parameter snapshots were introduced still require their original file under `params/`.
