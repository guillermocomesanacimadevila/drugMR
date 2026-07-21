#!/usr/bin/env python3
import subprocess
from pathlib import Path
from drugmr.config import Config

# * Notes for myself before going to Greece
# the git clone thingy
# remember QC run for GWAS as well
# singularity script (just pull and run)
# NO prep data 
# dashboard pull into local and do it there rather than in HPC
# pull TSV output also into local 
# then just script running stuff - for each part as a sequence with an main() in sequence as well (with appropaite ifs as checks and prints)

def ssh(cmd: str, falcon_user: str, allowed_returncodes: tuple = (0,)):
    full_cmd = f"ssh {falcon_user}@falconlogin.cf.ac.uk '{cmd}'"
    result = subprocess.run(full_cmd, shell=True, executable="/bin/bash", capture_output=True, text=True)
    if result.stdout:
        print(result.stdout)
    if result.returncode not in allowed_returncodes:
        print("[ERROR] Falcon command failed.")
        print(result.stderr)
        raise subprocess.CalledProcessError(result.returncode, full_cmd)
    return result

def get_remote_paths(falcon_user: str):
    remote = f"/shared/home1/{falcon_user}/drugMR"
    sif = f"{remote}/env/drugmr.sif"
    return remote, sif

def check_remote_output(
    falcon_user: str,
    path: str,
    step: str,
    overwrite: bool = False
):
    # run step if overwrite == True
    if overwrite:
        print(f"[TRACKING] Overwrite enabled - rerunning {step}...")
        return False

    remote, _ = get_remote_paths(falcon_user)

    result = ssh(f"""
set -euo pipefail
cd "{remote}"

if [ -s "{path}" ]; then
    echo "[TRACKING] {step} already completed: {path}"
    exit 0
fi

exit 3
""", falcon_user, allowed_returncodes=(0, 3))

    if result.returncode == 0:
        print(f"[TRACKING] Skipping {step}...")
        return True

    print(f"[TRACKING] No existing {step} output found - running step...")
    return False

def check_remote_cis_regions(
    falcon_user: str,
    pqtl_dataset: str,
    overwrite: bool = False
):
    # run step if overwrite == True
    if overwrite:
        print("[TRACKING] Overwrite enabled - rerunning cis-region preparation...")
        return False

    remote, _ = get_remote_paths(falcon_user)

    result = ssh(f"""
set -euo pipefail
cd "{remote}"

n_cis=$(find "dat/cis_regions/{pqtl_dataset}" -mindepth 2 -maxdepth 2 -name "pqtl.parquet" 2>/dev/null | wc -l)

if [ "$n_cis" -gt 0 ]; then
    echo "[TRACKING] cis-regions already completed: $n_cis loci found"
    exit 0
fi

exit 3
""", falcon_user, allowed_returncodes=(0, 3))

    if result.returncode == 0:
        print("[TRACKING] Skipping cis-region preparation...")
        return True

    print("[TRACKING] No complete cis-region output found - running step...")
    return False

def check_remote_cojo(
    falcon_user: str,
    pqtl_dataset: str,
    pheno_id: str,
    overwrite: bool = False
):
    # run step if overwrite == True
    if overwrite:
        print("[TRACKING] Overwrite enabled - rerunning GCTA-COJO...")
        return False

    remote, _ = get_remote_paths(falcon_user)

    result = ssh(f"""
set -euo pipefail
cd "{remote}"

n_cojo=$(find "results/COJO/{pqtl_dataset}/{pheno_id}" -mindepth 2 -maxdepth 2 -name "*.jma.cojo" -size +0c 2>/dev/null | wc -l)

if [ "$n_cojo" -gt 0 ]; then
    echo "[TRACKING] GCTA-COJO already completed: $n_cojo loci found"
    exit 0
fi

exit 3
""", falcon_user, allowed_returncodes=(0, 3))

    if result.returncode == 0:
        print("[TRACKING] Skipping GCTA-COJO...")
        return True

    print("[TRACKING] No complete GCTA-COJO outputs found - running step...")
    return False

def require_remote_output(
    falcon_user: str,
    path: str,
    step: str,
    required_for: str
):
    remote, _ = get_remote_paths(falcon_user)

    ssh(f"""
set -euo pipefail
cd "{remote}"

if [ ! -f "{path}" ]; then
    echo "[ERROR] {required_for} cannot run because {step} output was not found: {path}"
    exit 1
fi

if [ ! -s "{path}" ]; then
    echo "[ERROR] {required_for} cannot run because {step} output is empty: {path}"
    exit 1
fi

echo "[TRACKING] Required {step} output found for {required_for}"
""", falcon_user)

def clone_repo(falcon_user: str):
    ssh("""
set -euo pipefail

echo 'Hello Falcon HPC!'
if [ -d "$HOME/drugMR" ]; then
    echo "[TRACKING] I found the directory!"
    cd "$HOME/drugMR"

    echo "[TRACKING] Resetting Falcon repo to GitHub main..."
    git fetch origin main
    git reset --hard origin/main
    git clean -fd \
      -e dat/ \
      -e results/ \
      -e work/ \
      -e assets/config.yaml
else
    echo "[CONCERN] Yowza! I cannot see the drugMR directory..."
    echo "[TRACKING] Cloning from GitHub..."
    git clone https://github.com/guillermocomesanacimadevila/drugMR.git "$HOME/drugMR"
fi
""", falcon_user)

def container_checks(falcon_user: str):
    ssh("""
set -euo pipefail

if [ ! -d "$HOME/drugMR" ]; then
    git clone https://github.com/guillermocomesanacimadevila/drugMR.git "$HOME/drugMR"
fi

cd "$HOME/drugMR"
# git pull

chmod +x bin/bootstrap_hpc.sh
bash bin/bootstrap_hpc.sh
""", falcon_user)


# NOW -> FUNCTIONS TO RUN EACH SCRIPT FROM THE PIPELINE 

# **************************
# **************************
# ANALYTICS PIPELINE - START
# **************************
# **************************

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
    target_build: str,
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

    ssh(f"""
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
  --target_build {target_build} \\
  --n_cases {n_cases} \\
  --n_controls {n_controls} \\
  --falcon-user {falcon_user} \\
  {info_args} \\
  {flag_args}"
""", falcon_user)


# mediators
def run_mediator_qc(
    falcon_user: str,
    mediator_manifest: str,
    maf: float = 0.01,
    remove_mhc: bool = True,
    remove_apoe: bool = False,
    overwrite: bool = True,
):
    remote, sif = get_remote_paths(falcon_user)

    flag_args = ""
    if remove_mhc:
        flag_args += " --remove_mhc"
    if remove_apoe:
        flag_args += " --remove_apoe"
    if overwrite:
        flag_args += " --overwrite"

    ssh(f"""
set -euo pipefail
cd "{remote}"

apptainer exec --bind "{remote}:/work" "{sif}" \\
bash -c "cd /work && python bin/arrange_mediators.py \\
  --mediators \\
  --mediator-manifest {mediator_manifest} \\
  --maf {maf} \\
  {flag_args}"
""", falcon_user)



# *********** Extract cis-regions from pQTLs
def prep_cis_regions(
    falcon_user: str,
    pheno_id: str,
    pqtl_dataset: str,
    pqtl_dir: str
):
    remote, sif = get_remote_paths(falcon_user)

    ssh(f"""
set -euo pipefail 
cd "{remote}"

apptainer exec --bind "{remote}:/work" "{sif}" \\
bash -c "cd /work && python bin/prep_cis_regions.py \\
  --pqtl_dataset {pqtl_dataset} \\
  --pheno_id {pheno_id} \\
  --pqtl_dir {pqtl_dir}"
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

    ssh(f"""
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

def run_network_mr(
    falcon_user: str,
    pheno_id: str,
    pheno_gwas: str,
    ref_bfile: str,
    pqtl_dataset: str,
    pqtl_dir: str,
):
    remote, sif = get_remote_paths(falcon_user)

    ssh(f"""
set -euo pipefail
cd "{remote}"

apptainer exec --bind "{remote}:/work" \\
  --env PYTHONPATH=. \\
  "{sif}" \\
  python bin/assort_network_mr.py \\
    --pheno_id {pheno_id} \\
    --pheno_gwas {pheno_gwas} \\
    --ref_bfile {ref_bfile} \\
    --pqtl_dataset {pqtl_dataset} \\
    --pqtl_dir {pqtl_dir} \\
    --run_genomewide_mr \\
    --run_cis_mr_X_M \\
    --run_network_mr
""", falcon_user)


# RUN COLOC
def run_coloc_without_mediators(
    falcon_user: str,
    pqtl_dataset: str,
    pheno_id: str,
    n_cases: int,
    n_controls: int
):
    remote, sif = get_remote_paths(falcon_user)

    ssh(f"""
set -euo pipefail
cd "{remote}"

apptainer exec --bind "{remote}:/work" \\
  --env PYTHONPATH=. \\
  "{sif}" \\
  bash -c "cd /work && python bin/coloc_targets.py \\
    --pqtl_dataset {pqtl_dataset} \\
    --local_results_dir results/cis-MR \\
    --pqtl_dir dat/cis_regions/{pqtl_dataset} \\
    --pheno_id {pheno_id} \\
    --n_cases {n_cases} \\
    --n_controls {n_controls}"
""", falcon_user)


def run_coloc_with_mediators(
    falcon_user: str,
    pqtl_dataset: str,
    pheno_id: str,
    n_cases: int,
    n_controls: int,
    mediators: bool = False,
    mediator_manifest: str = ""
):

    remote, sif = get_remote_paths(falcon_user)

    ssh(f"""
set -euo pipefail
cd "{remote}"

apptainer exec --bind "{remote}:/work" \\
  --env PYTHONPATH=. \\
  "{sif}" \\
  bash -c "cd /work && python bin/coloc_targets.py \\
    --pqtl_dataset {pqtl_dataset} \\
    --local_results_dir results/cis-MR \\
    --pqtl_dir dat/cis_regions/{pqtl_dataset} \\
    --pheno_id {pheno_id} \\
    --n_cases {n_cases} \\
    --n_controls {n_controls} \\
    --mediators \\
    --mediator_manifest {mediator_manifest}"
""", falcon_user)


# def run_smr(
#         falcon_user: str, 
#         pqtl_dataset: str,
#         pheno_id: str,
# ):

def run_smr(
    falcon_user: str,
    pqtl_dataset: str, 
    pheno_id: str,
    eqtl_dataset: str,
    sumstats: str,
    ref_bfile: str,
    maf: float
):
    remote, sif = get_remote_paths(falcon_user)

    ssh(f"""
set -euo pipefail 
cd "{remote}"

apptainer exec --bind "{remote}:/work" \\
  --env PYTHONPATH=. \\
  "{sif}" \\
  bash -c "cd /work && python bin/sort_single_cell_smr.py \\
    --pqtl_dataset {pqtl_dataset} \\
    --pheno_id {pheno_id} \\
    --eqtl_dataset {eqtl_dataset} \\
    --sumstats {sumstats} \\
    --maf {maf} \\
    --ref_bfile {ref_bfile}"
""", falcon_user)


def run_multi_omics_coloc(
    falcon_user: str,
    pheno_id: str,
    pqtl_dataset: str,
    eqtl_dataset: str,
    n_cases: int,
    n_controls: int
):
    remote, sif = get_remote_paths(falcon_user)

    ssh(f"""
set -euo pipefail
cd "{remote}"

apptainer exec --bind "{remote}:/work" \\
  --env PYTHONPATH=. \\
  "{sif}" \\
  bash -c "cd /work && python bin/assort_moloc_for_sc_hits.py \\
    --pheno_id {pheno_id} \\
    --pqtl_dataset {pqtl_dataset} \\
    --eqtl_dataset {eqtl_dataset} \\
    --n_cases {n_cases} \\
    --n_controls {n_controls}"
""", falcon_user)

def run_multi_omics_summary(
    falcon_user: str,
    pheno_id: str,
    pqtl_dataset: str,
    eqtl_dataset: str
):
    remote, sif = get_remote_paths(falcon_user)

    ssh(f"""
set -euo pipefail
cd "{remote}"
apptainer exec --bind "{remote}:/work" \\
  --env PYTHONPATH=. \\
  "{sif}" \\
  bash -c "cd /work && python bin/summarise_multi_omics.py \\
    --pheno_id {pheno_id} \\
    --pqtl_dataset {pqtl_dataset} \\
    --eqtl_dataset {eqtl_dataset}"
""", falcon_user)

def run_cojo(
    falcon_user: str,
    pheno_id: str,
    pqtl_dataset: str,
    eqtl_dataset: str,
    ref_bfile: str,
):
    remote, sif = get_remote_paths(falcon_user)

    ssh(f"""
set -euo pipefail
cd "{remote}"

apptainer exec --bind "{remote}:/work" \\
  --env PYTHONPATH=/work \\
  "{sif}" \\
  bash -c "cd /work && python bin/cojo_on_pqtls.py \\
    --pheno_id {pheno_id} \\
    --pqtl_dataset {pqtl_dataset} \\
    --eqtl_dataset {eqtl_dataset} \\
    --ref_bfile {ref_bfile}"
""", falcon_user)

# RUN PHEWAS CHECKS FOR SAFETY (LOCALLY) -> API != WORK IN SLURM HPC
# ******************************************************************

def phewas_safety(
    pheno_id: str,
    pqtl_dataset: str,
    eqtl_dataset: str,
    local_results_dir: str = "results",
    overwrite: bool = False
):
    project_root = Path(__file__).resolve().parents[1]
    local_results_dir = Path(local_results_dir)

    if not local_results_dir.is_absolute():
        local_results_dir = project_root / local_results_dir

    top_snp_file = (
        local_results_dir
        / "SMR"
        / eqtl_dataset
        / pheno_id
        / f"{pqtl_dataset}_{pheno_id}_multi_omics_snp_evidence.tsv"
    )

    phewas_out = (
        local_results_dir
        / "PheWAS"
        / pqtl_dataset
        / pheno_id
        / f"{pqtl_dataset}_{pheno_id}_PheWAS.tsv"
    )

    if phewas_out.exists() and phewas_out.stat().st_size > 0 and not overwrite:
        print(f"[TRACKING] PheWAS safety analysis already completed: {phewas_out}")
        print("[TRACKING] Skipping PheWAS safety analysis...")
        return

    if overwrite:
        print("[TRACKING] Overwrite enabled - rerunning PheWAS safety analysis...")
    else:
        print("[TRACKING] No existing PheWAS safety output found - running step...")

    phewas_out.parent.mkdir(parents=True, exist_ok=True)

    cmd = f"""
set -euo pipefail
cd "{project_root}"
python bin/phewas_cis_pqtls.py \\
  --pheno_id {pheno_id} \\
  --pqtl_dataset {pqtl_dataset} \\
  --eqtl_dataset {eqtl_dataset}
"""

    print(f"[TRACKING] PheWAS SNP evidence input found: {top_snp_file}")
    print("[TRACKING] Running PheWAS safety analysis locally...")
    subprocess.run(cmd, shell=True, check=True, executable="/bin/bash")
    print(f"[TRACKING] PheWAS safety results found: {phewas_out}")


# ******************************************************************
# ******************************************************************


# **************************
# **************************
# ANALYTICS PIPELINE - END
# **************************
# **************************


# Database functs and dashboard assortments
# SLAP ONTO POSTGRESQL DB

def load_postgres(
    falcon_user: str,
    pqtl_dataset: str,
    pheno_id: str,
    db_id: str = "drugmr"
):
    remote, sif = get_remote_paths(falcon_user)
    mr_res = f"results/cis-MR/{pqtl_dataset}_{pheno_id}_all_MR.tsv"
    coloc_res = f"results/coloc/{pqtl_dataset}/{pqtl_dataset}_{pheno_id}_all_coloc.tsv"

    ssh(f"""
set -euo pipefail
cd "{remote}"

apptainer exec --bind "{remote}:/work" "{sif}" \\
bash -c "cd /work && python bin/load_db_into_postgres.py \\
  --results_file {mr_res} \\
  --db_id {db_id} \\
  --pqtl_dataset {pqtl_dataset} \\
  --pheno_id {pheno_id} \\
  --table cis_mr_results"

apptainer exec --bind "{remote}:/work" "{sif}" \\
bash -c "cd /work && python bin/load_db_into_postgres.py \\
  --results_file {coloc_res} \\
  --db_id {db_id} \\
  --pqtl_dataset {pqtl_dataset} \\
  --pheno_id {pheno_id} \\
  --table coloc_results"
""", falcon_user)


# PULL RESULTS INTO LOCAL
def pull_results_local(
    falcon_user: str,
    pqtl_dataset: str,
    pheno_id: str,
    eqtl_dataset: str,
    local_results_dir: str = "results",
    overwrite: bool = True
):
    remote, _ = get_remote_paths(falcon_user)
    remote_mr = f"{remote}/results/cis-MR/{pqtl_dataset}_{pheno_id}_all_MR.tsv"
    remote_coloc = f"{remote}/results/coloc/{pqtl_dataset}/{pqtl_dataset}_{pheno_id}_all_coloc.tsv"
    remote_smr = f"{remote}/results/SMR/{eqtl_dataset}/{pheno_id}"
    remote_eqtl_coloc = f"{remote}/results/eQTL_coloc/{pqtl_dataset}/{eqtl_dataset}/{pheno_id}"
    remote_moloc = f"{remote}/results/QTL_moloc/{pqtl_dataset}/{eqtl_dataset}/{pheno_id}"
    remote_cojo = f"{remote}/results/COJO/{pqtl_dataset}/{pheno_id}"
    local_results_dir = Path(local_results_dir)
    local_mr_dir = local_results_dir / "cis-MR"
    local_coloc_dir = local_results_dir / "coloc" / pqtl_dataset
    local_smr_dir = local_results_dir / "SMR" / eqtl_dataset / pheno_id
    local_eqtl_coloc_dir = local_results_dir / "eQTL_coloc" / pqtl_dataset / eqtl_dataset / pheno_id
    local_moloc_dir = local_results_dir / "QTL_moloc" / pqtl_dataset / eqtl_dataset / pheno_id
    local_cojo_dir = local_results_dir / "COJO" / pqtl_dataset / pheno_id
    local_mr_dir.mkdir(parents=True, exist_ok=True)
    local_coloc_dir.mkdir(parents=True, exist_ok=True)
    local_smr_dir.mkdir(parents=True, exist_ok=True)
    local_eqtl_coloc_dir.mkdir(parents=True, exist_ok=True)
    local_moloc_dir.mkdir(parents=True, exist_ok=True)
    local_cojo_dir.mkdir(parents=True, exist_ok=True)
    local_mr = local_mr_dir / f"{pqtl_dataset}_{pheno_id}_all_MR.tsv"
    local_coloc = local_coloc_dir / f"{pqtl_dataset}_{pheno_id}_all_coloc.tsv"
    for remote_file, local_file in [
        (remote_mr, local_mr),
        (remote_coloc, local_coloc),
    ]:
        if local_file.exists() and not overwrite:
            print(f"[TRACKING] {local_file} already exists locally. Skipping pull.")
            continue

        if local_file.exists() and overwrite:
            print(f"[TRACKING] {local_file} already exists locally. Overwriting...")

        cmd = f"scp {falcon_user}@falconlogin.cf.ac.uk:{remote_file} {local_file}"
        print(cmd)
        subprocess.run(cmd, shell=True, check=True)
        print(f"[DONE] Pulled results into {local_file}")

    # pull all compiled SMR results
    if any(local_smr_dir.iterdir()) and not overwrite:
        print(f"[TRACKING] {local_smr_dir} already contains SMR results. Skipping pull.")
    else:
        cmd = (
            f"scp -r "
            f"{falcon_user}@falconlogin.cf.ac.uk:{remote_smr}/. "
            f"{local_smr_dir}/"
        )
        print(cmd)
        subprocess.run(cmd, shell=True, check=True)
        print(f"[DONE] Pulled SMR results into {local_smr_dir}")

    # pull all GWAS - sc-eQTL COLOC results
    if any(local_eqtl_coloc_dir.iterdir()) and not overwrite:
        print(f"[TRACKING] {local_eqtl_coloc_dir} already contains eQTL COLOC results. Skipping pull.")
    else:
        cmd = f"scp -r {falcon_user}@falconlogin.cf.ac.uk:{remote_eqtl_coloc}/. {local_eqtl_coloc_dir}/"
        print(cmd)
        subprocess.run(cmd, shell=True, check=True)
        print(f"[DONE] Pulled eQTL COLOC results into {local_eqtl_coloc_dir}")

    # pull all GWAS - pQTL - sc-eQTL MOLOC results
    if any(local_moloc_dir.iterdir()) and not overwrite:
        print(f"[TRACKING] {local_moloc_dir} already contains MOLOC results. Skipping pull.")
    else:
        cmd = f"scp -r {falcon_user}@falconlogin.cf.ac.uk:{remote_moloc}/. {local_moloc_dir}/"
        print(cmd)
        subprocess.run(cmd, shell=True, check=True)
        print(f"[DONE] Pulled MOLOC results into {local_moloc_dir}")

    # pull all GCTA-COJO results
    if any(local_cojo_dir.iterdir()) and not overwrite:
        print(f"[TRACKING] {local_cojo_dir} already contains COJO results. Skipping pull.")
    else:
        cmd = f"scp -r {falcon_user}@falconlogin.cf.ac.uk:{remote_cojo}/. {local_cojo_dir}/"
        print(cmd)
        subprocess.run(cmd, shell=True, check=True)
        print(f"[DONE] Pulled COJO results into {local_cojo_dir}")


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
    pheno_id: str,
    eqtl_dataset: str
):
    remote, _ = get_remote_paths(falcon_user)
    mr_res = f"results/cis-MR/{pqtl_dataset}_{pheno_id}_all_MR.tsv"
    coloc_res = f"results/coloc/{pqtl_dataset}/{pqtl_dataset}_{pheno_id}_all_coloc.tsv"
    smr_res = f"results/SMR/{eqtl_dataset}/{pheno_id}/{pqtl_dataset}_{pheno_id}_promising_targets_SMR.tsv"
    final_targets = f"results/SMR/{eqtl_dataset}/{pheno_id}/{pqtl_dataset}_{pheno_id}_final_multi_omics_targets.tsv"
    prepared_multi_omics = f"results/SMR/{eqtl_dataset}/{pheno_id}/{pqtl_dataset}_{pheno_id}_prepared_multi_omics_targets.tsv"
    eqtl_coloc = f"results/eQTL_coloc/{pqtl_dataset}/{eqtl_dataset}/{pheno_id}/{pqtl_dataset}_{pheno_id}_{eqtl_dataset}_all_eqtl_coloc.tsv"
    moloc = f"results/QTL_moloc/{pqtl_dataset}/{eqtl_dataset}/{pheno_id}/{pqtl_dataset}_{pheno_id}_{eqtl_dataset}_moloc_summary.tsv"
    cojo_dir = f"results/COJO/{pqtl_dataset}/{pheno_id}"

    ssh(f"""
set -euo pipefail
cd "{remote}"

echo "[TRACKING] Checking MR output..."
if [ -s "{mr_res}" ]; then
    ls -lh "{mr_res}"
    head -5 "{mr_res}"
else
    echo "[CONCERN] MR output not found or empty"
fi

echo "[TRACKING] Checking COLOC output..."
if [ -s "{coloc_res}" ]; then
    ls -lh "{coloc_res}"
    head -5 "{coloc_res}"
else
    echo "[CONCERN] COLOC output not found or empty"
fi

echo "[TRACKING] Checking SMR output..."
if [ -d "results/SMR/{eqtl_dataset}/{pheno_id}" ]; then
    ls -lh "results/SMR/{eqtl_dataset}/{pheno_id}/"
else
    echo "[CONCERN] SMR output directory not found"
fi

if [ -s "{smr_res}" ]; then
    head -5 "{smr_res}"
else
    echo "[CONCERN] Promising target SMR file not found or empty"
fi

if [ -s "{final_targets}" ]; then
    echo "[TRACKING] Final multi-omics targets found!"
    head -5 "{final_targets}"
else
    echo "[CONCERN] No final multi-omics target TSV found or file is empty"
fi

echo "[TRACKING] Checking prepared multi-omics targets..."
if [ -s "{prepared_multi_omics}" ]; then
    head -5 "{prepared_multi_omics}"
else
    echo "[CONCERN] Prepared multi-omics target manifest not found or empty"
fi

echo "[TRACKING] Checking GWAS - sc-eQTL COLOC output..."
if [ -s "{eqtl_coloc}" ]; then
    head -5 "{eqtl_coloc}"
else
    echo "[CONCERN] GWAS - sc-eQTL COLOC output not found or empty"
fi

echo "[TRACKING] Checking GWAS - pQTL - sc-eQTL MOLOC output..."
if [ -s "{moloc}" ]; then
    head -5 "{moloc}"
else
    echo "[CONCERN] GWAS - pQTL - sc-eQTL MOLOC output not found or empty"
fi

echo "[TRACKING] Checking GCTA-COJO output..."
n_cojo=$(find "{cojo_dir}" -mindepth 2 -maxdepth 2 -name "*.jma.cojo" -size +0c 2>/dev/null | wc -l)
if [ "$n_cojo" -gt 0 ]; then
    echo "[TRACKING] GCTA-COJO outputs found for $n_cojo loci"
    find "{cojo_dir}" -mindepth 2 -maxdepth 2 -name "*.jma.cojo" -size +0c -print
else
    echo "[CONCERN] No GCTA-COJO .jma.cojo outputs found"
fi
""", falcon_user)


# Function to run all the HPC gist
def hpc(config: str = "assets/config.yaml"):
    cfg = Config(config)
    falcon_user = cfg.falcon_user
    pheno_id = cfg.pheno_id
    sumstats = cfg.sumstats
    n_cases = cfg.n_cases
    n_controls = cfg.n_controls
    pqtl_dataset = cfg.pqtl_dataset
    pqtl_dir = cfg.pqtl_dir
    ref_bfile = cfg.ref_bfile
    snp_col = cfg.snp_col
    eqtl_dataset = cfg.eqtl_dataset
    a1_col = cfg.a1_col
    a2_col = cfg.a2_col
    beta_col = cfg.beta_col
    se_col = cfg.se_col
    p_col = cfg.p_col
    pos_col = cfg.pos_col
    chr_col = cfg.chr_col
    af_col = cfg.af_col
    genome_build = cfg.genome_build
    target_build = cfg.target_build
    out_dir = getattr(cfg, "out_dir", "results")
    maf = getattr(cfg, "maf", 0.01)
    mediators = getattr(cfg, "mediators", False)
    mediator_manifest = getattr(cfg, "mediator_manifest", "")
    info_threshold = getattr(cfg, "info_threshold", None)
    info_col = getattr(cfg, "info_col", None)
    remove_mhc = getattr(cfg, "remove_mhc", True)
    remove_apoe = getattr(cfg, "remove_apoe", False)
    local_results_dir = getattr(cfg, "local_results_dir", "results")
    overwrite = getattr(cfg, "overwrite", False)

    # define all outputs first so pipeline knows what has already been ran
    qc_out = f"{out_dir}/QC/{pheno_id}/{pheno_id}.tsv"
    mr_out = f"results/cis-MR/{pqtl_dataset}_{pheno_id}_all_MR.tsv"
    coloc_out = f"results/coloc/{pqtl_dataset}/{pqtl_dataset}_{pheno_id}_all_coloc.tsv"
    promising_smr_out = f"results/SMR/{eqtl_dataset}/{pheno_id}/{pqtl_dataset}_{pheno_id}_promising_targets_SMR.tsv"
    final_smr_out = f"results/SMR/{eqtl_dataset}/{pheno_id}/{pqtl_dataset}_{pheno_id}_final_multi_omics_targets.tsv"
    prepared_multi_omics_out = f"results/SMR/{eqtl_dataset}/{pheno_id}/{pqtl_dataset}_{pheno_id}_prepared_multi_omics_targets.tsv"
    eqtl_coloc_out = f"results/eQTL_coloc/{pqtl_dataset}/{eqtl_dataset}/{pheno_id}/{pqtl_dataset}_{pheno_id}_{eqtl_dataset}_all_eqtl_coloc.tsv"
    moloc_out = f"results/QTL_moloc/{pqtl_dataset}/{eqtl_dataset}/{pheno_id}/{pqtl_dataset}_{pheno_id}_{eqtl_dataset}_moloc_summary.tsv"
    summary_out = f"results/SMR/{eqtl_dataset}/{pheno_id}/{pqtl_dataset}_{pheno_id}_multi_omics_overview.tsv"
    snp_evidence_out = f"results/SMR/{eqtl_dataset}/{pheno_id}/{pqtl_dataset}_{pheno_id}_multi_omics_snp_evidence.tsv"
    cojo_dir = f"results/COJO/{pqtl_dataset}/{pheno_id}"

    # change this where NetworkMR saves its final compiled output
    network_mr_out = f"results/network-MR/{pqtl_dataset}/{pqtl_dataset}_{pheno_id}_network_MR.tsv"

    print("[TRACKING] Preparing Falcon repo...")
    clone_repo(falcon_user)

    print("[TRACKING] Preparing Falcon env...")
    container_checks(falcon_user)

    if not check_remote_output(
        falcon_user=falcon_user,
        path=qc_out,
        step="GWAS QC",
        overwrite=overwrite
    ):
        print("[TRACKING] Running GWAS QC...")
        run_gwas_qc(
            falcon_user=falcon_user,
            pheno_id=pheno_id,
            sumstats=sumstats,
            out_dir=out_dir,
            snp_col=snp_col,
            a1_col=a1_col,
            a2_col=a2_col,
            beta_col=beta_col,
            se_col=se_col,
            p_col=p_col,
            pos_col=pos_col,
            chr_col=chr_col,
            af_col=af_col,
            genome_build=genome_build,
            target_build=target_build,
            n_cases=n_cases,
            n_controls=n_controls,
            maf=maf,
            info_threshold=info_threshold,
            info_col=info_col,
            remove_mhc=remove_mhc,
            remove_apoe=remove_apoe,
        )

    require_remote_output(
        falcon_user=falcon_user,
        path=qc_out,
        step="GWAS QC",
        required_for="cis-region preparation"
    )

    if mediators:
        print("[TRACKING] Running mediator QC...")
        run_mediator_qc(
            falcon_user=falcon_user,
            mediator_manifest=mediator_manifest,
            maf=maf,
            remove_mhc=remove_mhc,
            remove_apoe=remove_apoe,
            overwrite=overwrite,
        )
    else:
        print("[TRACKING] No mediators specificed, running drugMR without them then!")

    if not check_remote_cis_regions(
        falcon_user=falcon_user,
        pqtl_dataset=pqtl_dataset,
        overwrite=overwrite
    ):
        print("[TRACKING] Preparing cis-regions...")
        prep_cis_regions(
            falcon_user=falcon_user,
            pheno_id=pheno_id,
            pqtl_dataset=pqtl_dataset,
            pqtl_dir=pqtl_dir,
        )

    if not check_remote_output(
        falcon_user=falcon_user,
        path=mr_out,
        step="cis-MR",
        overwrite=overwrite
    ):
        print("[TRACKING] Running cis-MR...")
        run_cis_mr(
            falcon_user=falcon_user,
            pqtl_dataset=pqtl_dataset,
            pqtl_dir=f"dat/cis_regions/{pqtl_dataset}",
            pheno_id=pheno_id,
            pheno_gwas=qc_out,
            ref_bfile=ref_bfile,
        )

    require_remote_output(
        falcon_user=falcon_user,
        path=mr_out,
        step="cis-MR",
        required_for="COLOC"
    )

    if mediators:
        require_remote_output(
            falcon_user=falcon_user,
            path=mr_out,
            step="cis-MR",
            required_for="NetworkMR"
        )

        if not check_remote_output(
            falcon_user=falcon_user,
            path=network_mr_out,
            step="NetworkMR",
            overwrite=overwrite
        ):
            print("[TRACKING] Running NetworkMR with mediators...")
            run_network_mr(
                falcon_user=falcon_user,
                pheno_id=pheno_id,
                pheno_gwas=qc_out,
                ref_bfile=ref_bfile,
                pqtl_dataset=pqtl_dataset,
                pqtl_dir=f"dat/cis_regions/{pqtl_dataset}",
            )
    else:
        print("[TRACKING] No mediators specified, skipping NetworkMR.")

    if not check_remote_output(
        falcon_user=falcon_user,
        path=coloc_out,
        step="COLOC",
        overwrite=overwrite
    ):
        print("[TRACKING] Running COLOC...")

        if mediators:
            run_coloc_with_mediators(
                falcon_user=falcon_user,
                pqtl_dataset=pqtl_dataset,
                pheno_id=pheno_id,
                n_cases=n_cases,
                n_controls=n_controls,
                mediators=mediators,
                mediator_manifest=mediator_manifest    
            )
        else:
            run_coloc_without_mediators(
                falcon_user=falcon_user,
                pqtl_dataset=pqtl_dataset,
                pheno_id=pheno_id,
                n_cases=n_cases,
                n_controls=n_controls,
            )

    require_remote_output(
        falcon_user=falcon_user,
        path=mr_out,
        step="cis-MR",
        required_for="single-cell SMR"
    )

    require_remote_output(
        falcon_user=falcon_user,
        path=coloc_out,
        step="COLOC",
        required_for="single-cell SMR"
    )

    # SMR
    if not check_remote_output(
        falcon_user=falcon_user,
        path=final_smr_out,
        step="single-cell SMR",
        overwrite=overwrite
    ):
        print(f"[TRACKING] Running single-cell SMR for {eqtl_dataset}...")
        run_smr(
            falcon_user=falcon_user,
            pqtl_dataset=pqtl_dataset,
            eqtl_dataset=eqtl_dataset,
            maf=maf,
            ref_bfile=ref_bfile,
            sumstats=qc_out,
            pheno_id=pheno_id
        )

    require_remote_output(
        falcon_user=falcon_user,
        path=final_smr_out,
        step="single-cell SMR",
        required_for="multi-omics QTL colocalisation"
    )

    # multi-omics QTL colocalisation
    # GWAS - sc-eQTL pairwise coloc
    # GWAS - pQTL - sc-eQTL MOLOC
    if not check_remote_output(
        falcon_user=falcon_user,
        path=moloc_out,
        step="multi-omics QTL colocalisation",
        overwrite=overwrite
    ):
        print("[TRACKING] Running multi-omics QTL colocalisation...")
        run_multi_omics_coloc(
            falcon_user=falcon_user,
            pheno_id=pheno_id,
            pqtl_dataset=pqtl_dataset,
            eqtl_dataset=eqtl_dataset,
            n_cases=n_cases,
            n_controls=n_controls
        )

    require_remote_output(
        falcon_user=falcon_user,
        path=prepared_multi_omics_out,
        step="multi-omics target preparation",
        required_for="pipeline completion"
    )

    require_remote_output(
        falcon_user=falcon_user,
        path=eqtl_coloc_out,
        step="GWAS - sc-eQTL COLOC",
        required_for="pipeline completion"
    )

    require_remote_output(
        falcon_user=falcon_user,
        path=moloc_out,
        step="GWAS - pQTL - sc-eQTL MOLOC",
        required_for="pipeline completion"
    )

    # grab onto master df for 2 key SNPs
    overview_done = check_remote_output(
        falcon_user=falcon_user,
        path=summary_out,
        step="multi-omics overview",
        overwrite=overwrite
    )

    snp_done = check_remote_output(
        falcon_user=falcon_user,
        path=snp_evidence_out,
        step="multi-omics SNP evidence",
        overwrite=overwrite
    )

    if not (overview_done and snp_done):
        print("[TRACKING] Building dashboard-ready multi-omics tables...")
        run_multi_omics_summary(
            falcon_user=falcon_user,
            pheno_id=pheno_id,
            pqtl_dataset=pqtl_dataset,
            eqtl_dataset=eqtl_dataset,
        )

    require_remote_output(
        falcon_user=falcon_user,
        path=summary_out,
        step="multi-omics overview",
        required_for="pipeline completion"
    )

    require_remote_output(
        falcon_user=falcon_user,
        path=snp_evidence_out,
        step="multi-omics SNP evidence",
        required_for="pipeline completion"
    )

    # GCTA-COJO
    if not check_remote_cojo(
        falcon_user=falcon_user,
        pqtl_dataset=pqtl_dataset,
        pheno_id=pheno_id,
        overwrite=overwrite
    ):
        print("[TRACKING] Running GCTA-COJO...")
        run_cojo(
            falcon_user=falcon_user,
            pheno_id=pheno_id,
            pqtl_dataset=pqtl_dataset,
            eqtl_dataset=eqtl_dataset,
            ref_bfile=ref_bfile,
        )

    if not check_remote_cojo(
        falcon_user=falcon_user,
        pqtl_dataset=pqtl_dataset,
        pheno_id=pheno_id,
        overwrite=False
    ):
        raise RuntimeError(
            f"No GCTA-COJO outputs were produced in {cojo_dir}"
        )

    print("[TRACKING] Checking outputs...")
    check_outputs(
        falcon_user=falcon_user,
        pqtl_dataset=pqtl_dataset,
        pheno_id=pheno_id,
        eqtl_dataset=eqtl_dataset,
    )

    print("[TRACKING] Pulling results locally...")
    pull_results_local(
        falcon_user=falcon_user,
        pqtl_dataset=pqtl_dataset,
        pheno_id=pheno_id,
        eqtl_dataset=eqtl_dataset,
        local_results_dir=local_results_dir,
        overwrite=overwrite,
    )

    print("[TRACKING] Running PheWAS safety analysis locally...")
    phewas_safety(
        pheno_id=pheno_id,
        pqtl_dataset=pqtl_dataset,
        eqtl_dataset=eqtl_dataset,
        local_results_dir=local_results_dir,
        overwrite=overwrite,
    )

    print(f"[TRACKING] Expected promising target SMR output: {promising_smr_out}")
    print(f"[TRACKING] Expected final multi-omics target output: {final_smr_out}")
    print(f"[TRACKING] Expected prepared multi-omics target output: {prepared_multi_omics_out}")
    print(f"[TRACKING] Expected GWAS - sc-eQTL COLOC output: {eqtl_coloc_out}")
    print(f"[TRACKING] Expected GWAS - pQTL - sc-eQTL MOLOC output: {moloc_out}")
    print(f"[TRACKING] Expected GCTA-COJO output directory: {cojo_dir}")
    print("[DONE] drugMR pipeline completed successfully.")