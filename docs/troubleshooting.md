# Troubleshooting

## Exit status 137 or `OUT_OF_MEMORY`

Exit status 137 commonly means SLURM killed the process after it reached its memory allocation. Confirm with:

```bash
sacct -j JOB_ID --format=JobID,State,ExitCode,ReqMem,MaxRSS,Elapsed --units=G
```

On SLURM, use `-profile falcon`. Using only `-profile apptainer` runs tasks inside the controller allocation and prevents the per-process SLURM resource rules from taking effect.

## The controller exits while a child job is still running

Inspect the child directly with `squeue` and `sacct`. Do not submit another copy while the child is active. After it finishes, rerun the controller with `-resume`; generated SMR files under `synthesis/qtl_esd/` are reusable.

## Follow a process log

Nextflow prints each failed process work directory. Its unbuffered task output is in `.command.out`:

```bash
tail -20 tmp/nf-work/XX/HASH/.command.out
```

Use `tail -f` to follow it and press ++ctrl+c++ to stop following. This does not cancel the SLURM job.

## Graphviz warnings

Graphviz is only needed to render the execution DAG. A missing Graphviz installation does not explain an analysis process failure. The execution report and timeline warnings should be investigated separately if those artifacts are required.

## A reference file exists but is reported missing

Relative paths are resolved from the repository, while a Nextflow task executes inside its work directory. Keep current code and pass reference paths through the configuration. If the data is outside the repository, bind the parent directory into the container with `--container_bind`.
