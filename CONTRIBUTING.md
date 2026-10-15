# Contributing to drugMR

Thank you for your interest in drugMR. Bug reports, questions and pull requests are all welcome.

## Reporting a problem

Open an [issue](https://github.com/guillermocomesanacimadevila/drugMR/issues/new/choose) and choose the template that fits. For a bug, please include:

- the exact `nextflow run` command and the params file you used
- the profile (`docker`, `apptainer`, `singularity`, ...) and where it ran (laptop, HPC, cloud)
- the `run_id` and, if the run finished, its `params.lock.yaml`
- the end of `.nextflow.log` and the error message of the failing process

Before opening an issue, check the [troubleshooting page](https://guillermocomesanacimadevila.github.io/drugMR/troubleshooting.html).

## Setting up a development copy

```bash
git clone --recurse-submodules https://github.com/guillermocomesanacimadevila/drugMR.git
cd drugMR
./env/bootstrap.sh
source env/activate.sh
```

`env/bootstrap.sh` creates `.venv` with the pinned Python packages, installs the `drugmr` package in editable mode and installs a matching Nextflow.

## Checking a change

Every change should pass the same end to end test that CI runs on each push:

```bash
nextflow run main.nf \
    -profile test,docker \
    -params-file tests/data/params.yaml

python tests/check_results.py
```

For changes that affect results, also run the synthetic demo (see *Try the demo* in the README) and compare its outputs with a run from the previous commit.

## Pull requests

1. Create a branch from `main`.
2. Keep each pull request to one change, and describe what it changes and why.
3. Run the test above and fill in the pull request checklist.
4. Add a line to `CHANGELOG.md` under *Unreleased*.
5. Update the documentation in `docs/` when a parameter, output column or behaviour changes.

## Code style

- Plain, readable code: explicit loops and `if` blocks, one step per line, descriptive names, short comments.
- File and column names shared between scripts are defined in `drugmr/paths.py`; change them there and update every reader.
- British English in documentation, comments and user facing text.

## Questions

For anything else, open a [discussion or issue](https://github.com/guillermocomesanacimadevila/drugMR/issues) or see the [Reach out page](https://guillermocomesanacimadevila.github.io/drugMR/reach-out.html).
