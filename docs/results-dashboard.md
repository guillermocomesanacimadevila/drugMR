# Results and dashboard

Use Python for the post-pipeline steps:

```python
import drugmr as dm
```

## Fetch a run from HPC or cloud

Skip this step when the completed run already exists on the machine where the dashboard will run.

```python
config = dm.fetch_run(
    run_id="AD_test_20260911_143022_abcdef0",
    user="your_username",
    host="login.your-cluster.ac.uk",
    remote_root="/path/to/drugMR/runs",
)
```

`fetch_run()` copies `remote_root/<run_id>/` into the local `runs/` directory with `rsync`. This includes the results, pipeline reports, `manifest.json`, and `params.lock.yaml`. It validates those metadata files, registers the run locally, and returns the local path to the locked parameter snapshot for `dm.results()`.

The `remote_root` value ends at the remote `runs` directory. `fetch_run()` appends the supplied `run_id` itself.

A worked example of this fetch and results flow is in `notebooks/00_drugmr.ipynb`. Use the notebook for the fetch step only when SSH keys are already configured for passwordless login. If SSH requests a password or private-key passphrase, run the fetch in a local terminal instead, since the notebook kernel does not reliably expose an interactive SSH prompt and the fetch can hang or fail silently. The notebook is fine for `dm.results()` once the run is already local.

```bash
python -c 'import drugmr as dm; config = dm.fetch_run(
    run_id="AD_test_20260911_143022_abcdef0",
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
dm.results(config="runs/AD_test_20260911_143022_abcdef0/params.lock.yaml")
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

Run `dm.results(config=config)` on the local computer after `fetch_run()`. The machine running the dashboard must support Docker Compose because drugMR starts PostgreSQL in Docker. Standard HPC environments use Apptainer for pipeline jobs and do not provide a user Docker daemon, so the dashboard should not be launched there.

```{warning}
Fetch only completed runs. The successful Nextflow controller creates `manifest.json` and `params.lock.yaml` and records the run only after all required stages finish. Older runs created before parameter snapshots were introduced still require their original file under `params/`.
```

## Explore results in Python

`ResultsExplorer` reads a completed run's results tables directly, without PostgreSQL or Docker. It works on any run under `runs/`, fetched or local. Every method returns a table indexed by gene symbol, reuses the pipeline's own gates, and accepts an optional protein (symbol in any case, or full protein ID).

```python
from drugmr.interpret import ResultsExplorer

explorer = ResultsExplorer(
    run_id="AD_test_20260911_143022_abcdef0",
    root="runs",
)

explorer.funnel()
explorer.targets()
explorer.target("TREM2")
```

`root` is the directory holding the runs. Pass `root="../runs"` when working from inside `notebooks/`.

Tables show the key columns by default. Pass `detailed=True` to any table method to get every column, for example `explorer.cis_mr(detailed=True)`.

| Method | Returns |
| --- | --- |
| `run_info()` | Run identity, provenance, outcome GWAS inputs and every gate used |
| `funnel()` | Proteins remaining after each stage, with bulk and single cell xQTL support kept apart |
| `cis_mr()` | One row per protein: OR with 95% CI, P, FDR, effect and drug direction |
| `instruments()` | One row per instrument SNP: alleles, pQTL and outcome effects, F statistic, Steiger |
| `sensitivity()` | IVW, MR Egger and weighted median side by side for multi instrument proteins |
| `colocalisation()` | COLOC and PWCoCo per protein, the best PWCoCo pair and which method supports it |
| `coloc_priors()` | COLOC PP4 under each prior set |
| `lead_variants()` | Top cis variant with outcome and pQTL effects on the same allele |
| `smr()` | SMR and HEIDI per protein and tissue or cell type, with direction against cis-MR |
| `hyprcoloc()` | Clusters sharing one causal variant across pQTL, GWAS and xQTL |
| `pwcoco_qtl()` | PWCoCo of each xQTL against the pQTL and against the outcome GWAS |
| `directions()` | One row per target: whether every line of evidence points the same way |
| `phewas()` | Bonferroni significant PheWAS traits, classed as possible additional indications or adverse effects |
| `targets()` | One row per prioritised target with every stage summarised |
| `target(protein)` | Every table above for one protein, in reading order |
| `forest_plot()` | Odds ratios with 95% CI for the prioritised targets |
| `regional_plot(protein)` | Outcome GWAS and pQTL across the cis window, coloured by LD with the lead SNP |

A guided walkthrough of every method on a completed run is in `notebooks/01_interpreting_results.ipynb`.
