# Troubleshooting

## Exit status 137 or `OUT_OF_MEMORY`

Exit status 137 commonly means SLURM killed the process after it reached its memory allocation. Confirm with:

```bash
sacct -j JOB_ID --format=JobID,State,ExitCode,ReqMem,MaxRSS,Elapsed --units=G
```

On SLURM, use `-profile slurm` when Nextflow is launched from a login or workflow node and should submit every process as a separate SLURM job. If Nextflow itself is launched by an `sbatch` script, use `-profile slurm,local` so processes run inside that allocation. Using `sbatch` with `-profile slurm` alone creates nested jobs and can deadlock under a one-node-per-user QOS with `QOSMaxNodePerUserLimit`.

## The controller exits while a child job is still running

Inspect the child directly with `squeue` and `sacct`. Do not submit another copy while the child is active. After it finishes, rerun the controller with `-resume`. A complete conversion found under `synthesis/qtl_esd/` is recovered and moved beside its manifest-declared QTL source; later runs reuse it there.

## Follow a process log

Nextflow prints each failed process work directory. Its unbuffered task output is in `.command.out`:

```bash
tail -20 tmp/nf-work/XX/HASH/.command.out
```

Use `tail -f` to follow it and press ++ctrl+c++ to stop following. This does not cancel the SLURM job.

## An SMR results file has a header but no rows

This is expected when no target passed the upstream cis-MR and colocalisation gates for that pQTL and QTL dataset pair. drugMR still writes a schema valid, header-only TSV so the declared Nextflow output exists, rather than failing the process. Check the upstream `cis_mr` and `coloc` results and gate thresholds for that run if targets were expected. This is not a pipeline failure.

## Graphviz warnings

Graphviz is only needed to render the execution DAG. A missing Graphviz installation does not explain an analysis process failure. The execution report and timeline warnings should be investigated separately if those artifacts are required.

## A reference file exists but is reported missing

Relative paths are resolved from the repository, while a Nextflow task executes inside its work directory. Keep current code and pass reference paths through the configuration. If the data is outside the repository, bind the parent directory into the container with `--container_bind`.

## GWAS and pQTL positions do not match

**Symptom.** Regional plots show the GWAS peak shifted from the pQTL peak, or cis regions keep far fewer GWAS variants than expected. For one SNP, the GWAS and pQTL positions differ by the same amount across a region.

**Cause.** `genome_build` in the params file does not match the GWAS file. The most common case is a GWAS already in GRCh38 declared as GRCh37, so it is lifted a second time.

**Fix.** Check a few known rsIDs against dbSNP, set `genome_build` to the real build, and rerun from GWAS QC.

## A rerun reuses results you expected to change

Nextflow `-resume` reruns a task when its inputs, script or parameters change. Two caches sit outside Nextflow and are reused across runs:

| Cache | Reused when | To force a rebuild |
| --- | --- | --- |
| BESD files built from QTL parquet files (`synthesis/qtl_esd/` and next to the manifest path) | The BESD files for that dataset already exist. | Delete that dataset's `.besd`, `.esi` and `.epi` files. |
| SMR results (`synthesis/SMR/`) | A non empty `.smr` file exists for the same phenotype and QTL dataset. | Delete the matching `.smr` files. |

If you change the QTL data, `gates.smr.p_qtl_smr`, `gates.smr.p_qtl_heidi` or `gates.smr.diff_freq_prop`, delete the files above first, or the old results are reused. The SMR pass thresholds (`p_smr_threshold`, `p_heidi_threshold`) are applied after SMR and take effect without deleting anything.

Editing an unrelated row of `assets/qtl_manifest.csv` does not invalidate cached cis regions. Each run stages only the manifest row for its own pQTL dataset.
