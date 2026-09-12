# Install drugMR

drugMR runs from a complete clone of this repository. The source repository is approximately 70 MB before reference data (1.51 Gb unzipped), QTL data, containers, work files, and results are added.

## What you need

1. Python 3.12 or newer.
2. R 4.4 or newer.
3. Java 17 or newer for Nextflow.
4. Nextflow 26.04.0 or newer.
5. One supported container runtime. The available profiles are Docker, Apptainer, Singularity, Podman, Shifter, and Charliecloud.
6. Docker Compose and PostgreSQL 16 or newer on the machine that will run the dashboard.
7. `rsync` on both machines when a completed run will be fetched from HPC.

## Example: local installation

Clone the whole repository on the machine that will run Nextflow and the dashboard:

```bash
git clone --recurse-submodules https://github.com/guillermocomesanacimadevila/drugMR.git
cd drugMR
./env/bootstrap.sh
source env/activate.sh
```

Download the reference bundle:

```bash
wget -c -O ref.zip "https://zenodo.org/records/22706705/files/ref.zip?download=1"
unzip ref.zip -d ref
rm ref.zip
```

Activate the environment in every new terminal:

```bash
cd /path/to/drugMR
source env/activate.sh
```

Confirm the installation:

```bash
python -c 'import drugmr; print("drugMR is ready")'
nextflow -version
```

## Example: HPC installation

For an HPC analysis with a local dashboard, clone the whole repository twice:

1. Clone it on the HPC system. This copy runs Nextflow and stores the full analysis.
2. Clone it on your local computer. This copy receives completed runs and starts PostgreSQL and Streamlit.

On the HPC system:

```bash
cd /shared/scratch/YOUR_PROJECT/pipelines
git clone --recurse-submodules https://github.com/guillermocomesanacimadevila/drugMR.git
cd drugMR
module load Java/17
./env/bootstrap.sh
source env/activate.sh
```

On the local computer:

```bash
cd ~/bioinformatics
git clone --recurse-submodules https://github.com/guillermocomesanacimadevila/drugMR.git
cd drugMR
./env/bootstrap.sh
source env/activate.sh
```

Keep both clones on the same Git revision where possible:

```bash
git rev-parse --short HEAD
```

Reference and input data do not need to be stored inside the repository. Paths in params files and the QTL manifest may be relative to the repository or absolute. When data is outside the repository, expose its parent directory to the container:

```bash
--container_bind /shared/scratch/YOUR_PROJECT/data
```

Do not copy the large HPC datasets to the local clone just to use the dashboard. `dm.fetch_run()` transfers the completed run, including the portable locus data needed by the dashboard.
