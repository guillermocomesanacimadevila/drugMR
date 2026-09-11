# Results and dashboard

Use Python for the post-pipeline steps:

```python
import drugmr as dm
```

## Fetch a run from HPC or cloud

Skip this step when the completed run already exists on the machine where the dashboard will run.

```python
dm.fetch_run(
    run_id="AD_test_20260911_abcdef0",
    user="your_username",
    host="login.your-cluster.ac.uk",
    remote_root="/path/to/drugMR/runs",
)
```

`fetch_run()` copies `remote_root/<run_id>/` into the local `runs/` directory with `rsync`, validates the run manifest, and registers the run locally.

### Prefer the command line for SSH prompts

If SSH requests a password or private-key passphrase, run the fetch in a local terminal. Notebook interfaces do not always expose an interactive SSH prompt reliably.

```bash
python -c 'import drugmr as dm; dm.fetch_run(
    run_id="AD_test_20260911_abcdef0",
    user="your_username",
    host="login.your-cluster.ac.uk",
    remote_root="/path/to/drugMR/runs"
)'
```

SSH keys avoid repeated password prompts and work well with `rsync`.

## Load results and open the dashboard

```python
dm.results(config="params/AD.test.yaml")
```

This command:

1. Finds the latest successful registered run matching `pheno_id` and `pqtl_dataset` in the configuration.
2. Starts PostgreSQL with `docker compose up -d`.
3. Writes the Streamlit database connection settings.
4. Loads the run's MR and colocalisation results into PostgreSQL.
5. Starts the Streamlit dashboard.

Run it locally after `fetch_run()` when you want the dashboard on your computer. You may instead run it directly on the HPC if Docker and the dashboard environment are available there. To view an HPC-hosted dashboard, forward its Streamlit port through SSH according to your cluster's access rules.

!!! warning
    Fetch only completed runs. The successful Nextflow controller records the run in its registry after all required stages finish.
