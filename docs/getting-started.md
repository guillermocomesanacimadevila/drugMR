# Install

## Requirements

- Python 3.12 or newer
- Java for Nextflow
- Docker, Apptainer, Singularity, Podman, Shifter, or Charliecloud for pipeline containers
- Docker Compose and PostgreSQL 16 or newer for the results dashboard
- `rsync` on both machines when fetching a remote run

Download the reference data from the [drugMR Zenodo record](https://doi.org/10.5281/zenodo.22706705), then clone and bootstrap the repository:

```bash
git clone --recurse-submodules https://github.com/guillermocomesanacimadevila/drugMR.git
cd drugMR
./env/bootstrap.sh
source env/activate.sh
```

Run `source env/activate.sh` in every new shell. On an HPC system, load an available Java module first if Java is not already on `PATH`:

```bash
module avail java
module load Java/17
source env/activate.sh
```

Confirm the installation:

```bash
python -c 'import drugmr; print("drugMR is ready")'
nextflow -version
```

## Reference and input locations

Paths in the parameter file and QTL manifest may be repository-relative or absolute. If input data lives outside the repository, expose its parent directory to the container with `--container_bind`.

```bash
--container_bind /shared/scratch/PROJECT/data
```
