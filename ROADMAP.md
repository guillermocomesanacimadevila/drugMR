# drugMR DevOps roadmap

## Add CI workflow (GitHub actions + pytest + toy datasets in test/) -> accomodate .gitignore

## Add a Nextflow DSL2 wrapper (modules/ and workflows/) - after the "science" aspect is perfected

## Re-structure out/ directory properly (outs in .parquet or .gz)
### - DONE (2026-08-11): dat/derived/ (shared preprocessing) + runs/<run_id>/ (per-run,
###   with manifest.json + registry.json) + synthesis/ (cross-dataset) + params/
###   (schema-validated, one file per pheno_id x pqtl_dataset, with a gates: block
###   replacing several hardcoded/debug-broken thresholds). Dashboard reads via the
###   registry now, not path-guessing.
### - STILL OPEN: outputs are still .tsv, not .parquet/.gz - format conversion was not
###   part of this pass.

## Upload reference data (./dat/ref) onto Dropbox/Zenodo

## Remove shell=True from subprocess.run()