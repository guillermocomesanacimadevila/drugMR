# drugMR

**A multi-omics pipeline for genetically anchored drug target discovery.**

![drugMR analysis pipeline](pipeline.png)

drugMR takes an outcome GWAS and a panel of protein QTLs and produces a ranked, safety-screened shortlist of druggable targets. Its Nextflow workflow connects Mendelian randomisation, colocalisation, SMR and HEIDI, PWCoCo, HyPrColoc, and phenome-wide safety screening.

[Install drugMR](getting-started.md){ .md-button .md-button--primary }
[View on GitHub](https://github.com/guillermocomesanacimadevila/drugMR){ .md-button }

## Workflow

1. Quality-control the outcome GWAS and extract cis regions.
2. Run cis-MR and pairwise colocalisation against the selected pQTL panel.
3. Triangulate promising targets with registered bulk and single-cell QTL datasets.
4. Run PheWAS safety screens using FinnGen and UK Biobank data.
5. Load a completed run into PostgreSQL and explore it in the Streamlit dashboard.

## Choose where each part runs

The computational pipeline can run locally or through SLURM. The results dashboard runs on whichever machine calls `dm.results()`.

| Situation | What to do |
| --- | --- |
| Pipeline and dashboard on one machine | Run Nextflow, then call `dm.results()` there. |
| Pipeline on HPC, dashboard on your computer | Run Nextflow on HPC, call `dm.fetch_run()` on your computer, then call `dm.results()` locally. |
| Pipeline and dashboard on HPC | Run Nextflow and `dm.results()` on HPC; use SSH port forwarding to view Streamlit. |

!!! note
    `dm.fetch_run()` is only a transfer step. Use it when the Nextflow run is on a remote HPC or cloud machine and the dashboard will run elsewhere.
