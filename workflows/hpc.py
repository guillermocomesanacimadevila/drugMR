#!/usr/bin/env python3
import subprocess
from pathlib import Path

# * Notes for myself before going to Greece
# the git clone thingy
# remember QC run for GWAS as well
# singularity script (just pull and run)
# NO prep data 
# dashboard pull into local and do it there rather than in HPC
# pull TSV output also into local 
# then just script running stuff - for each part as a sequence with an main() in sequence as well (with appropaite ifs as checks and prints)

def hpc(cmd: str, falcon_user: str):
    full_cmd = f"ssh {falcon_user}@falconlogin.cf.ac.uk '{cmd}'"
    print(full_cmd)
    subprocess.run(full_cmd, shell=True, check=True)

def get_remote_paths(falcon_user: str):
    remote = f"/shared/home1/{falcon_user}/drugMR"
    sif = f"{remote}/env/drugmr.sif"
    return remote, sif

def clone_repo(falcon_user: str):
    hpc("""
set -euo pipefail

if [ -d "$HOME/drugMR" ]; then
    echo "[TRACKING] I found the directory!"
    cd "$HOME/drugMR"
    git pull
else
    echo "[CONCERN] Yowza! I cannot see the drugMR directory..."
    echo "[TRACKING] Cloning from GitHub..."
    git clone https://github.com/guillermocomesanacimadevila/drugMR.git "$HOME/drugMR"
fi
""", falcon_user)

def container_checks(falcon_user: str):
    hpc("""
set -euo pipefail

if [ ! -d "$HOME/drugMR" ]; then
    git clone https://github.com/guillermocomesanacimadevila/drugMR.git "$HOME/drugMR"
fi

cd "$HOME/drugMR"
git pull

chmod +x bin/bootstrap_hpc.sh
bash bin/bootstrap_hpc.sh
""", falcon_user)

# NOW -> FUNCTIONS TO RUN EACH SCRIPT FROM THE PIPELINE 

# QC GWAS

def run_gwas_qc(
    falcon_user: str,
    pheno_id: str,
    sumstats: str,
    out_dir: str,
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
    n_cases: int,
    n_controls: int,
    maf: float = 0.01,
    info_threshold: float | None = None,
    info_col: str | None = None,
    remove_mhc: bool = True,
    remove_apoe: bool = False
):
    remote, sif = get_remote_paths(falcon_user)

    info_args = ""
    if info_col is not None:
        info_args += f" --info-col {info_col}"
    if info_threshold is not None:
        info_args += f" --info-threshold {info_threshold}"

    flag_args = ""
    if remove_mhc:
        flag_args += " --remove_mhc"
    if remove_apoe:
        flag_args += " --remove_apoe"

    hpc(f"""
set -euo pipefail
cd "{remote}"

apptainer exec --bind "{remote}:/work" "{sif}" \\
bash -c "cd /work && python bin/qc_gwas.py \\
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
  --n_cases {n_cases} \\
  --n_controls {n_controls} \\
  --falcon-user {falcon_user} \\
  {info_args} \\
  {flag_args}"
""", falcon_user)

# RUN MR 

def run_cis_mr(
    falcon_user: str,
    pqtl_dataset: str,
    pqtl_dir: str,
    pheno_id: str,
    pheno_gwas: str,
    ref_bfile: str
):
    remote, sif = get_remote_paths(falcon_user)

    hpc(f"""
set -euo pipefail
cd "{remote}"

apptainer exec --bind "{remote}:/work" "{sif}" \\
bash -c "cd /work && Rscript bin/cis_mr.R \\
  {pqtl_dataset} \\
  {pqtl_dir} \\
  {pheno_id} \\
  {pheno_gwas} \\
  {ref_bfile}"
""", falcon_user)

# SLAP ONTO POSTGRESQL DB

def load_postgres(
    falcon_user: str,
    pqtl_dataset: str,
    pheno_id: str,
    db_id: str = "drugmr"
):
    remote, sif = get_remote_paths(falcon_user)
    mr_res = f"results/cis-MR/{pqtl_dataset}_{pheno_id}_all_MR.tsv"

    hpc(f"""
set -euo pipefail
cd "{remote}"

apptainer exec --bind "{remote}:/work" "{sif}" \\
bash -c "cd /work && python bin/load_db_into_postgres.py \\
  --mr_res {mr_res} \\
  --db_id {db_id} \\
  --pqtl_dataset {pqtl_dataset} \\
  --pheno_id {pheno_id}"
""", falcon_user)

# PULL RESULTS INTO LOCAL

def pull_results_local(
    falcon_user: str,
    pqtl_dataset: str,
    pheno_id: str,
    local_results_dir: str = "results/cis-MR"
):
    remote, _ = get_remote_paths(falcon_user)
    remote_file = f"{remote}/results/cis-MR/{pqtl_dataset}_{pheno_id}_all_MR.tsv"
    local_results_dir = Path(local_results_dir)
    local_results_dir.mkdir(parents=True, exist_ok=True)

    cmd = f"scp {falcon_user}@falconlogin.cf.ac.uk:{remote_file} {local_results_dir}/"
    print(cmd)
    subprocess.run(cmd, shell=True, check=True)

# STREAMLIT DASHBOARD

def run_dashboard_local(
    db_name: str,
    phenotype: str,
    port_number: int = 5432
):
    cmd = f"""
python -m streamlit run dashboard/mr_app.py -- \\
  --db_name {db_name} \\
  --port_number {port_number} \\
  --phenotype {phenotype}
"""
    print(cmd)
    subprocess.run(cmd, shell=True, check=True)

# CHECK OUTPUTS

def check_outputs(
    falcon_user: str,
    pqtl_dataset: str,
    pheno_id: str
):
    remote, _ = get_remote_paths(falcon_user)
    mr_res = f"results/cis-MR/{pqtl_dataset}_{pheno_id}_all_MR.tsv"

    hpc(f"""
set -euo pipefail
cd "{remote}"

echo "[TRACKING] Checking MR output..."
ls -lh results/cis-MR/
head -5 "{mr_res}"
""", falcon_user)

