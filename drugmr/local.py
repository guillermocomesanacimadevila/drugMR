#!/usr/bin/env python3
import sys
import subprocess
from pathlib import Path

# LOCAL RESULTS / DASHBOARD STUFF
# To do's (After Greece)
# Need to create local running functions including a pulling docker from container function
# So then still the QC+MR runs in micromamba Docker env
# PostgreSQL db pulling and dashboard == jupyter (with .toml in ./)
# Goal == have a flagging variable (local/hpc)

def cmd_base(cmd):
    """
    Baseline function for parsing and running CLI-based .py scripts
    """
    return subprocess.run(cmd, shell=True, check=True, executable="/bin/bash")

def results(
    pqtl_dataset: str,
    pheno_id: str,
    db_id: str = "drugmr",
    dashboard_script: str = "dashboard/mr_app.py",
    db_script: str = "bin/load_db_into_postgres.py",
    port_number: int = 5432,
):
    project_root = Path(__file__).resolve().parents[1]
    mr_res = project_root / "results" / "cis-MR" / f"{pqtl_dataset}_{pheno_id}_all_MR.tsv"
    coloc_res = project_root / "results" / "coloc" / pqtl_dataset / f"{pqtl_dataset}_{pheno_id}_all_coloc.tsv"
    db_script = project_root / db_script
    dashboard_script = project_root / dashboard_script

    print("[TRACKING] Loading MR results into PostgreSQL...")

    subprocess.run(
        [
            sys.executable,
            str(db_script),
            "--results_file",
            str(mr_res),
            "--db_id",
            db_id,
            "--pqtl_dataset",
            pqtl_dataset,
            "--pheno_id",
            pheno_id,
            "--table",
            "cis_mr_results",
        ],
        check=True,
    )

    print("[TRACKING] Loading COLOC results into PostgreSQL...")

    subprocess.run(
        [
            sys.executable,
            str(db_script),
            "--results_file",
            str(coloc_res),
            "--db_id",
            db_id,
            "--pqtl_dataset",
            pqtl_dataset,
            "--pheno_id",
            pheno_id,
            "--table",
            "coloc_results",
        ],
        check=True,
    )

    print("[TRACKING] Launching Streamlit dashboard...")

    subprocess.run(
        [
            sys.executable,
            "-m",
            "streamlit",
            "run",
            str(dashboard_script),
            "--",
            "--db_name",
            db_id,
            "--port_number",
            str(port_number),
            "--phenotype",
            pheno_id,
        ],
        check=True,
    )

# secondment functions
# if run with local -> load up docker container
# i.e. check whether the container exists within current local machine 
# if so -> load docker container -> and run within that
# produce the same output as in the cloud -> which then run local scripts to load into postgres db and then dashboard
# to check docker container - create a function within local()

def local(
    pheno_id: str,
    sumstats: str,
    n_cases: int,
    n_controls: int,
    pqtl_dataset: str,
    pqtl_dir: str,
    ref_bfile: str,
    snp_col: str,
    a1_col: str,
    a2_col: str,
    beta_col: str,
    se_col: str,
    p_col: str,
    pos_col: str,
    chr_col: str,
    af_col: str,
    genome_build: str,
    target_build: str,
    out_dir: str = "results", # default dir in ./results within drugmR/
    maf: float = 0.01, # set default at 0.01 
    info_threshold: float | None = None,
    info_col: str | None = None,
    mediators: bool = False,
    mediator_manifest: str = "",
    remove_mhc: bool = True,
    remove_apoe: bool = False,
    image_uri: str = "ghcr.io/guillermocomesanacimadevila/drugmr:latest",
    image_name: str = "ghcr.io/guillermocomesanacimadevila/drugmr:latest"
):
    
    # set projectDir()
    project_root = Path(__file__).resolve().parents[1] # i.e. "Users/.../drugMR"

    def check_docker():
        print("[TRACKING] Checking Docker...")
        cmd = f"""
set -euo pipefail

if ! command -v docker >/dev/null 2>&1; then
    echo "[ERROR] Mate, install Docker before you run this locally."
    exit 1
fi

if ! docker info >/dev/null 2>&1; then
    echo "[ERROR] Mate, Docker is installed but it is not running."
    echo "[HINT] Open Docker Desktop and try again."
    exit 1
fi

echo "[TRACKING] Docker is installed and running."

if docker image inspect "{image_name}" >/dev/null 2>&1; then
    echo "[TRACKING] drugMR Docker image already exists locally."
else
    echo "[TRACKING] drugMR Docker image not found locally."
    echo "[TRACKING] Pulling drugMR image from GHCR..."
    docker pull "{image_uri}"
fi
        """
        cmd_base(cmd)

    # call check docker function
    check_docker()

    info_args = ""
    if info_col is not None:
        info_args += f"--info-col {info_col}"
    if info_threshold is not None:
        info_args += f" --info-threshold {info_threshold}"

    flag_args = ""
    if remove_mhc:
        flag_args += " --remove_mhc"
    if remove_apoe:
        flag_args += " --remove_apoe"

    # running individual modules
    print("[TRACKING] Running GWAS QC locally via Docker...")

    cmd_qc = f"""
set -euo pipefail 
docker run --rm \\
  -v "{project_root}:/work" \\
  -w /work \\
  "{image_name}" \\
  python bin/qc_gwas.py \\
    --pheno-id {pheno_id} \\
    --sumstats {sumstats} \\
    --out-dir {out_dir} \\
    --maf {maf} \\
    --snp-col {snp_col} \\
    --a1-col {a1_col} \\
    --a2-col {a2_col} \\
    --beta-col {beta_col} \\
    --se-col {se_col} \\
    --p-col {p_col} \\
    --pos-col {pos_col} \\
    --chr-col {chr_col} \\
    --af_col {af_col} \\
    --genome_build {genome_build} \\
    --target_build {target_build} \\
    --n_cases {n_cases} \\
    --n_controls {n_controls} \\
    --falcon-user local \\
    {info_args} \\
    {flag_args}
"""
    cmd_base(cmd_qc)

    # mediators stuff
    if mediators:
        print("[TRACKING] Qceing mediators locally via Docker...")

        mediator_args = f"--mediator-manifest {mediator_manifest}"
        if remove_mhc:
            mediator_args += " --remove_mhc"
        if remove_apoe:
            mediator_args += " --remove_apoe"

        cmd_m_qc = f"""
set -euo pipefail 
docker run --rm \\
  -v "{project_root}:/work" \\
  -w /work \\
  "{image_name}" \\
  python bin/arrange_mediators.py \\
    --mediators \\
    {mediator_args} \\
    --maf {maf}
"""
        cmd_base(cmd_m_qc)

    else:
        print("[TRACKING] No mediators specificed, running drugMR without them then!")

    # cis-region module
    print("[TRACKING] Preparing cis-regions locally...")
    
    cmd_cis = f"""
set -euo pipefail
docker run --rm \\
  -v "{project_root}:/work" \\
  -w /work \\
  "{image_name}" \\
  python bin/prep_cis_regions.py \\
    --pqtl_dataset {pqtl_dataset} \\
    --pheno_id {pheno_id} \\
    --pqtl_dir {pqtl_dir}
"""
    cmd_base(cmd_cis)

    cis_dir = project_root / "dat" / "cis_regions" / pqtl_dataset
    print(f"[TRACKING] Checking cis-region output: {cis_dir}")

    if not cis_dir.exists():
        raise FileNotFoundError(f"cis-region directory not created: {cis_dir}")

    n_cis = len(list(cis_dir.glob("*/pqtl.parquet")))
    print(f"[TRACKING] cis-region loci generated: {n_cis}")

    if n_cis == 0:
        raise RuntimeError("No cis-region files generated. Check pqtl_dir path.")

    # cis-MR module 
    print("[TRACKING] Running cis-MR locally via Docker...")

    cmd_mr = f"""
set -euo pipefail
docker run --rm \\
  -v "{project_root}:/work" \\
  -w /work \\
  "{image_name}" \\
  Rscript bin/cis_mr.R \\
    {pqtl_dataset} \\
    dat/cis_regions/{pqtl_dataset} \\
    {pheno_id} \\
    {out_dir}/QC/{pheno_id}/{pheno_id}.tsv \\
    {ref_bfile}
"""
    cmd_base(cmd_mr)

    mr_out = project_root / "results" / "cis-MR" / f"{pqtl_dataset}_{pheno_id}_all_MR.tsv"

    if not mr_out.exists():
        print(f"[CONCERN] MR results file not found: {mr_out}")
        print("[CONCERN] Skipping COLOC because no MR results were generated.")
        print("[DONE] Local Docker run completed.")
        return

    # CMD COLOC TARGETS
    # CMD RUN COLOC (Pairwise)
    # Need to test


    # networkMR (HERE)
        # networkMR
    if mediators:
        print("[TRACKING] Running NetworkMR with mediators!")

        cmd_network_mr = f"""
set -euo pipefail
docker run --rm \\
  -v "{project_root}:/work" \\
  -w /work \\
  -e PYTHONPATH=. \\
  "{image_name}" \\
  python bin/assort_network_mr.py \\
    --pheno_id {pheno_id} \\
    --pheno_gwas {out_dir}/QC/{pheno_id}/{pheno_id}.tsv \\
    --ref_bfile {ref_bfile} \\
    --pqtl_dataset {pqtl_dataset} \\
    --pqtl_dir {pqtl_dir} \\
    --run_genomewide_mr \\
    --run_cis_mr_X_M \\
    --run_network_mr
"""
        cmd_base(cmd_network_mr)
    else:
        print("[TRACKING] No mediators specified, skipping NetworkMR.")


    # coloc target module
    # ******** RE-DO -> AD if mediators:
    print("[TRACKING] Running COLOC locally...")

    # cmd_coloc with mediators
    cmd_coloc_with_mediators = f"""
set -euo pipefail
docker run --rm \\
  -v "{project_root}:/work" \\
  -w /work \\
  "{image_name}" \\
  python bin/coloc_targets.py \\
    --pqtl_dataset {pqtl_dataset} \\
    --local_results_dir results/cis-MR \\
    --pqtl_dir dat/cis_regions \\
    --pheno_id {pheno_id} \\
    --n_cases {n_cases} \\
    --n_controls {n_controls} \\
    --mediators \\
    --mediator_manifest {mediator_manifest}
"""
    
    # without mediators
    cmd_coloc_without_mediators = f"""
set -euo pipefail
docker run --rm \\
  -v "{project_root}:/work" \\
  -w /work \\
  "{image_name}" \\
  python bin/coloc_targets.py \\
    --pqtl_dataset {pqtl_dataset} \\
    --local_results_dir results/cis-MR \\
    --pqtl_dir dat/cis_regions \\
    --pheno_id {pheno_id} \\
    --n_cases {n_cases} \\
    --n_controls {n_controls}
"""
    
    if mediators:
        cmd_base(cmd_coloc_with_mediators)
    else:
        cmd_base(cmd_coloc_without_mediators)

    print("[DONE] Local Docker run completed.")