<div class="drugmr-hero" markdown>

<span class="drugmr-eyebrow">GENETICALLY ANCHORED DRUG DISCOVERY</span>

# From association to actionable target

drugMR is a reproducible multi-omics pipeline that turns outcome GWAS and protein QTL panels into a ranked, safety-screened shortlist of druggable targets.

[Get started](getting-started.md){ .md-button .md-button--primary }
[Explore the workflow](#workflow){ .md-button }

<div class="drugmr-pills">
  <span>Nextflow DSL2</span><span>Cis-MR</span><span>Colocalisation</span><span>SMR + HEIDI</span><span>PheWAS</span>
</div>

</div>

<div class="drugmr-stat-grid" markdown>

<div class="drugmr-stat" markdown>
**End to end**

One workflow from GWAS QC to dashboard-ready evidence.
</div>

<div class="drugmr-stat" markdown>
**Multi-omics**

Triangulate protein, bulk-tissue and single-cell QTL evidence.
</div>

<div class="drugmr-stat" markdown>
**Reproducible**

Containerised execution, locked parameters and registered runs.
</div>

</div>

## Workflow

![drugMR analysis pipeline](pipeline.png){ .drugmr-pipeline }

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

!!! note
    `dm.fetch_run()` is only needed when the Nextflow run is on a remote HPC or cloud machine and the dashboard will run elsewhere. It transfers the run's locked parameter snapshot along with its results.

    `dm.results()` always needs a configuration. For an exact local run, use `dm.results(config="runs/<run_id>/params.lock.yaml")`. To select the latest successful run matching a regular parameter file, use `dm.results(config="params/<file>.yaml")`.
