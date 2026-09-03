#!/usr/bin/env python3
import subprocess
from datetime import datetime
from pathlib import Path

from drugmr import paths, registry
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
      -e runs/ \
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

apptainer exec --bind "{remote}:/work" \\
  --env PYTHONPATH=. \\
  "{sif}" \\
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

apptainer exec --bind "{remote}:/work" \\
  --env PYTHONPATH=. \\
  "{sif}" \\
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
    ref_bfile: str,
    out_dir: str = "results"
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
  {ref_bfile} \\
  {out_dir}"
""", falcon_user)

def run_network_mr(
    falcon_user: str,
    pheno_id: str,
    pheno_gwas: str,
    ref_bfile: str,
    pqtl_dataset: str,
    pqtl_dir: str,
    local_results_dir: str = "results",
    ivw_fdr_q: float = 0.05,
    egger_intercept_pval_min: float = 0,
    cochran_q_pval: float = 0.05,
    m_y_pval_threshold: float = 0.05,
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
    --run_network_mr \\
    --local_results_dir {local_results_dir} \\
    --ivw_fdr_q {ivw_fdr_q} \\
    --egger_intercept_pval_min {egger_intercept_pval_min} \\
    --cochran_q_pval {cochran_q_pval} \\
    --m_y_pval_threshold {m_y_pval_threshold}
""", falcon_user)


# RUN COLOC
def run_coloc_without_mediators(
    falcon_user: str,
    pqtl_dataset: str,
    pheno_id: str,
    n_cases: int,
    n_controls: int,
    local_results_dir: str = "results",
    wald_fdr_q: float = 0.05,
    ivw_fdr_q: float = 0.05,
    cochran_q_pval: float = 0.05,
    egger_intercept_pval_min: float = 0,
    min_instruments_for_ivw: int = 3,
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
    --local_results_dir {local_results_dir} \\
    --pqtl_dir dat/cis_regions/{pqtl_dataset} \\
    --pheno_id {pheno_id} \\
    --n_cases {n_cases} \\
    --n_controls {n_controls} \\
    --wald_fdr_q {wald_fdr_q} \\
    --ivw_fdr_q {ivw_fdr_q} \\
    --cochran_q_pval {cochran_q_pval} \\
    --egger_intercept_pval_min {egger_intercept_pval_min} \\
    --min_instruments_for_ivw {min_instruments_for_ivw}"
""", falcon_user)


# RUN PWCoCo
def run_pwcoco(
    falcon_user: str,
    pqtl_dataset: str,
    pheno_id: str,
    ref_bfile: str,
    n_cases: int,
    n_controls: int,
    local_results_dir: str = "results",
    cochran_q_pval: float = 0.05,
    wald_fdr_q: float = 0.05,
):
    remote, sif = get_remote_paths(falcon_user)

    ssh(f"""
set -euo pipefail
cd "{remote}"

apptainer exec --bind "{remote}:/work" \\
  --env PYTHONPATH=. \\
  "{sif}" \\
  bash -c "cd /work && python bin/pwcoco_wrapper.py \\
    --pqtl_dataset {pqtl_dataset} \\
    --pheno_id {pheno_id} \\
    --ref_bfile {ref_bfile} \\
    --n_cases {n_cases} \\
    --n_controls {n_controls} \\
    --local_results_dir {local_results_dir} \\
    --cochran_q_pval {cochran_q_pval} \\
    --wald_fdr_q {wald_fdr_q}"
""", falcon_user)


# RUN PWCoCo (eQTL-informed) - eQTL-pQTL / eQTL-GWAS PWCoCo on every SMR-passing
# target, then compared for shared colocalising SNPs against the pQTL-GWAS PWCoCo
# above (see project_pwcoco_wiring memory / bin/pwcoco_qtl_wrapper.py)
def run_pwcoco_qtl(
    falcon_user: str,
    pqtl_dataset: str,
    pheno_id: str,
    ref_bfile: str,
    n_cases: int,
    n_controls: int,
    local_results_dir: str = "results",
):
    remote, sif = get_remote_paths(falcon_user)

    ssh(f"""
set -euo pipefail
cd "{remote}"

apptainer exec --bind "{remote}:/work" \\
  --env PYTHONPATH=. \\
  "{sif}" \\
  bash -c "cd /work && python bin/pwcoco_qtl_wrapper.py \\
    --pqtl_dataset {pqtl_dataset} \\
    --pheno_id {pheno_id} \\
    --ref_bfile {ref_bfile} \\
    --n_cases {n_cases} \\
    --n_controls {n_controls} \\
    --local_results_dir {local_results_dir}"
""", falcon_user)


def run_coloc_with_mediators(
    falcon_user: str,
    pqtl_dataset: str,
    pheno_id: str,
    n_cases: int,
    n_controls: int,
    mediators: bool = False,
    mediator_manifest: str = "",
    local_results_dir: str = "results",
    ivw_fdr_q: float = 0.05,
    pp4_threshold: float = 0.7,
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
    --local_results_dir {local_results_dir} \\
    --pqtl_dir dat/cis_regions/{pqtl_dataset} \\
    --pheno_id {pheno_id} \\
    --n_cases {n_cases} \\
    --n_controls {n_controls} \\
    --mediators \\
    --mediator_manifest {mediator_manifest} \\
    --ivw_fdr_q {ivw_fdr_q} \\
    --pp4_threshold {pp4_threshold}"
""", falcon_user)

# RUN SMR (bulk or single-cell, depending on eqtl_mode)
# named run_smr_step (not run_smr) to avoid clashing with the run_smr config flag in hpc()
def run_smr_step(
    falcon_user: str,
    pqtl_dataset: str,
    eqtl_dataset: str,
    eqtl_mode: str,
    pheno_id: str,
    sumstats: str,
    ref_bfile: str,
    maf: float,
    local_results_dir: str = "results"
):
    remote, sif = get_remote_paths(falcon_user)

    ssh(f"""
set -euo pipefail
cd "{remote}"

apptainer exec --bind "{remote}:/work" \\
  --env PYTHONPATH=. \\
  "{sif}" \\
  bash -c "cd /work && python bin/sort_smr.py \\
    --pheno_id {pheno_id} \\
    --sumstats {sumstats} \\
    --pqtl_dataset {pqtl_dataset} \\
    --eqtl_dataset {eqtl_dataset} \\
    --eqtl_mode {eqtl_mode} \\
    --ref_bfile {ref_bfile} \\
    --maf {maf} \\
    --local_results_dir {local_results_dir}"
""", falcon_user)


# RUN HyPrColoc (bulk and/or single-cell eQTL) - for every target x cell-type/tissue
# hit in the combined final multi-omics target table for the given eqtl_dataset, runs
# a 3-trait (pQTL / GWAS / eQTL) HyPrColoc restricted to that target's cis-region
def run_hyprcoloc_step(
    falcon_user: str,
    pqtl_dataset: str,
    pheno_id: str,
    eqtl_dataset: str,
    local_results_dir: str = "results"
):
    remote, sif = get_remote_paths(falcon_user)

    ssh(f"""
set -euo pipefail
cd "{remote}"

apptainer exec --bind "{remote}:/work" \\
  --env PYTHONPATH=. \\
  "{sif}" \\
  bash -c "cd /work && python bin/hyprcoloc_targets.py \\
    --pqtl_dataset {pqtl_dataset} \\
    --pheno_id {pheno_id} \\
    --eqtl_dataset {eqtl_dataset} \\
    --local_results_dir {local_results_dir}"
""", falcon_user)


# get final snp-wide hits
def compile_top_hits(
    falcon_user: str,
    pheno_id: str,
    pqtl_dataset: str,
    local_results_dir: str = "results"
):
    remote, sif = get_remote_paths(falcon_user)

    ssh(f"""
set -euo pipefail
cd "{remote}"

apptainer exec --bind "{remote}:/work" \\
  --env PYTHONPATH=/work \\
  "{sif}" \\
  bash -c "cd /work && python bin/compile_cis_hit_info.py \\
    --pheno_id {pheno_id} \\
    --pqtl_dataset {pqtl_dataset} \\
    --local_results_dir {local_results_dir}"
""", falcon_user)

# RUN PHEWAS CHECKS FOR SAFETY (LOCALLY) -> API != WORK IN SLURM HPC
# ******************************************************************

def phewas_safety_finngen(
    pheno_id: str,
    pqtl_dataset: str,
    local_results_dir: str = "results",
    overwrite: bool = False
):
    project_root = Path(__file__).resolve().parents[1]
    local_results_dir = Path(local_results_dir)

    if not local_results_dir.is_absolute():
        local_results_dir = project_root / local_results_dir

    top_snp_file = paths.coloc_out(pqtl_dataset, pheno_id, out_dir=str(local_results_dir))
    phewas_out = paths.phewas_out(pqtl_dataset, pheno_id, out_dir=str(local_results_dir))

    if phewas_out.exists() and phewas_out.stat().st_size > 0 and not overwrite:
        print(f"[TRACKING] FinnGen PheWAS safety analysis already completed: {phewas_out}")
        print("[TRACKING] Skipping FinnGen PheWAS safety analysis...")
        return

    if overwrite:
        print("[TRACKING] Overwrite enabled - rerunning FinnGen PheWAS safety analysis...")
    else:
        print("[TRACKING] No existing FinnGen PheWAS safety output found - running step...")

    if not top_snp_file.exists():
        raise FileNotFoundError(
            f"FinnGen PheWAS cannot run because pairwise COLOC output was not found: {top_snp_file}"
        )

    if top_snp_file.stat().st_size == 0:
        raise RuntimeError(
            f"FinnGen PheWAS cannot run because pairwise COLOC output is empty: {top_snp_file}"
        )

    phewas_out.parent.mkdir(parents=True, exist_ok=True)

    cmd = f"""
set -euo pipefail
cd "{project_root}"
python bin/phewas_cis_pqtls.py \\
  --pheno_id {pheno_id} \\
  --pqtl_dataset {pqtl_dataset} \\
  --local_results_dir {local_results_dir}
"""

    print(f"[TRACKING] FinnGen PheWAS pairwise COLOC input found: {top_snp_file}")
    print("[TRACKING] Running FinnGen PheWAS safety analysis locally...")
    subprocess.run(cmd, shell=True, check=True, executable="/bin/bash")
    print(f"[TRACKING] FinnGen PheWAS safety results found: {phewas_out}")


def phewas_safety_ukbb(
    pheno_id: str,
    pqtl_dataset: str,
    local_results_dir: str = "results",
    overwrite: bool = False
):
    project_root = Path(__file__).resolve().parents[1]
    local_results_dir = Path(local_results_dir)

    if not local_results_dir.is_absolute():
        local_results_dir = project_root / local_results_dir

    top_snp_file = paths.coloc_out(pqtl_dataset, pheno_id, out_dir=str(local_results_dir))
    phewas_out = paths.phewas_ukbb_out(pqtl_dataset, pheno_id, out_dir=str(local_results_dir))
    # UKB PheWAS is a fallback - only run for targets with zero FinnGen instrument
    # coverage - and reads this manifest internally to build that fallback set
    finngen_coverage_file = paths.phewas_finngen_coverage_out(pqtl_dataset, pheno_id, out_dir=str(local_results_dir))

    if phewas_out.exists() and phewas_out.stat().st_size > 0 and not overwrite:
        print(f"[TRACKING] UKBB PheWAS safety analysis already completed: {phewas_out}")
        print("[TRACKING] Skipping UKBB PheWAS safety analysis...")
        return

    if overwrite:
        print("[TRACKING] Overwrite enabled - rerunning UKBB PheWAS safety analysis...")
    else:
        print("[TRACKING] No existing UKBB PheWAS safety output found - running step...")

    if not top_snp_file.exists():
        raise FileNotFoundError(
            f"UKBB PheWAS cannot run because pairwise COLOC output was not found: {top_snp_file}"
        )

    if top_snp_file.stat().st_size == 0:
        raise RuntimeError(
            f"UKBB PheWAS cannot run because pairwise COLOC output is empty: {top_snp_file}"
        )

    if not finngen_coverage_file.exists():
        raise FileNotFoundError(
            f"UKBB PheWAS cannot run because the FinnGen PheWAS coverage manifest was not found: "
            f"{finngen_coverage_file}. Run FinnGen PheWAS first."
        )

    phewas_out.parent.mkdir(parents=True, exist_ok=True)

    cmd = f"""
set -euo pipefail
cd "{project_root}"
python bin/ukb_phewas.py \\
  --pheno_id {pheno_id} \\
  --pqtl_dataset {pqtl_dataset} \\
  --local_results_dir {local_results_dir}
"""

    print(f"[TRACKING] UKBB PheWAS pairwise COLOC input found: {top_snp_file}")
    print("[TRACKING] Running UKBB PheWAS safety analysis locally...")
    subprocess.run(cmd, shell=True, check=True, executable="/bin/bash")
    print(f"[TRACKING] UKBB PheWAS safety results found: {phewas_out}")


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
    db_id: str = "drugmr",
    local_results_dir: str = "results"
):
    remote, sif = get_remote_paths(falcon_user)
    mr_res = str(paths.mr_out(pqtl_dataset, pheno_id, local_results_dir))
    coloc_res = str(paths.coloc_out(pqtl_dataset, pheno_id, local_results_dir))

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
    local_results_dir: str = "results",
    overwrite: bool = True
):
    remote, _ = get_remote_paths(falcon_user)
    remote_mr = f"{remote}/{paths.mr_out(pqtl_dataset, pheno_id)}"
    remote_coloc = f"{remote}/{paths.coloc_out(pqtl_dataset, pheno_id)}"
    remote_target_stats = f"{remote}/{paths.target_stats_out(pqtl_dataset, pheno_id)}"
    remote_smr = f"{remote}/{paths.smr_final_targets_out(pqtl_dataset, pheno_id)}"
    remote_hyprcoloc = f"{remote}/{paths.hyprcoloc_out(pqtl_dataset, pheno_id)}"
    local_results_dir = Path(local_results_dir)
    local_mr = paths.mr_out(pqtl_dataset, pheno_id, out_dir=str(local_results_dir))
    local_coloc = paths.coloc_out(pqtl_dataset, pheno_id, out_dir=str(local_results_dir))
    local_target_stats = paths.target_stats_out(pqtl_dataset, pheno_id, out_dir=str(local_results_dir))
    local_smr = paths.smr_final_targets_out(pqtl_dataset, pheno_id, out_dir=str(local_results_dir))
    local_hyprcoloc = paths.hyprcoloc_out(pqtl_dataset, pheno_id, out_dir=str(local_results_dir))
    local_mr.parent.mkdir(parents=True, exist_ok=True)
    local_coloc.parent.mkdir(parents=True, exist_ok=True)
    local_target_stats.parent.mkdir(parents=True, exist_ok=True)
    local_smr.parent.mkdir(parents=True, exist_ok=True)
    local_hyprcoloc.parent.mkdir(parents=True, exist_ok=True)
    for remote_file, local_file in [
        (remote_mr, local_mr),
        (remote_coloc, local_coloc),
        (remote_target_stats, local_target_stats),
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

    # SMR is optional (bulk and/or single-cell, gated by run_smr) so only pull it
    # down if it was actually produced remotely
    if local_smr.exists() and not overwrite:
        print(f"[TRACKING] {local_smr} already exists locally. Skipping pull.")
    else:
        remote_smr_check = check_remote_output(
            falcon_user=falcon_user,
            path=str(paths.smr_final_targets_out(pqtl_dataset, pheno_id)),
            step="SMR",
            overwrite=False
        )

        if remote_smr_check:
            cmd = f"scp {falcon_user}@falconlogin.cf.ac.uk:{remote_smr} {local_smr}"
            print(cmd)
            subprocess.run(cmd, shell=True, check=True)
            print(f"[DONE] Pulled results into {local_smr}")
        else:
            print("[TRACKING] No remote SMR output found - skipping SMR pull.")

    # HyPrColoc is also optional (gated by bulk_eqtl_datasets / sc_eqtl_dataset)
    # so only pull it down if it was actually produced remotely
    if local_hyprcoloc.exists() and not overwrite:
        print(f"[TRACKING] {local_hyprcoloc} already exists locally. Skipping pull.")
    else:
        remote_hyprcoloc_check = check_remote_output(
            falcon_user=falcon_user,
            path=str(paths.hyprcoloc_out(pqtl_dataset, pheno_id)),
            step="HyPrColoc",
            overwrite=False
        )

        if remote_hyprcoloc_check:
            cmd = f"scp {falcon_user}@falconlogin.cf.ac.uk:{remote_hyprcoloc} {local_hyprcoloc}"
            print(cmd)
            subprocess.run(cmd, shell=True, check=True)
            print(f"[DONE] Pulled results into {local_hyprcoloc}")
        else:
            print("[TRACKING] No remote HyPrColoc output found - skipping HyPrColoc pull.")


# STREAMLIT DASHBOARD
def run_dashboard_local(
    db_name: str,
    phenotype: str,
    pqtl_dataset: str,
    port_number: int = 5432
):
    # cwd is pinned to project_root (same fix as drugmr/local.py's results()) so
    # Streamlit reliably finds <project_root>/.streamlit/config.toml (the custom
    # theme) regardless of the caller's own working directory, instead of
    # silently falling back to Streamlit defaults.
    project_root = Path(__file__).resolve().parents[1]
    cmd = f"""
python -m streamlit run dashboard/mr_app.py -- \\
  --db_name {db_name} \\
  --port_number {port_number} \\
  --phenotype {phenotype} \\
  --pqtl_dataset {pqtl_dataset}
"""
    print(cmd)
    subprocess.run(cmd, shell=True, check=True, cwd=str(project_root))


# CHECK OUTPUTS
def check_outputs(
    falcon_user: str,
    pqtl_dataset: str,
    pheno_id: str,
    local_results_dir: str = "results"
):
    remote, _ = get_remote_paths(falcon_user)
    mr_res = str(paths.mr_out(pqtl_dataset, pheno_id, local_results_dir))
    coloc_res = str(paths.coloc_out(pqtl_dataset, pheno_id, local_results_dir))
    target_stats_res = str(paths.target_stats_out(pqtl_dataset, pheno_id, local_results_dir))
    smr_res = str(paths.smr_final_targets_out(pqtl_dataset, pheno_id, local_results_dir))
    hyprcoloc_res = str(paths.hyprcoloc_out(pqtl_dataset, pheno_id, local_results_dir))

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

echo "[TRACKING] Checking top cis-hit compilation output..."
if [ -s "{target_stats_res}" ]; then
    ls -lh "{target_stats_res}"
    head -5 "{target_stats_res}"
else
    echo "[CONCERN] Top cis-hit compilation output not found or empty"
fi

echo "[TRACKING] Checking SMR output..."
if [ -s "{smr_res}" ]; then
    ls -lh "{smr_res}"
    head -5 "{smr_res}"
else
    echo "[CONCERN] SMR output not found or empty (SMR may not be configured for this run)"
fi

echo "[TRACKING] Checking HyPrColoc output..."
if [ -s "{hyprcoloc_res}" ]; then
    ls -lh "{hyprcoloc_res}"
    head -5 "{hyprcoloc_res}"
else
    echo "[CONCERN] HyPrColoc output not found or empty (HyPrColoc may not be configured for this run)"
fi

""", falcon_user)


# Function to run all the HPC gist
def hpc(config: str, run_id: str = None):
    # config has no default on purpose - there's no single correct params file
    # anymore now that each (pheno_id, pqtl_dataset) pair has its own under
    # params/ (e.g. params/AD.wingo_brain.yaml) - pass one explicitly.
    #
    # run_id defaults to None, which keeps the deterministic
    # (pheno_id, pqtl_dataset, day, remote commit) behaviour below. Pass an
    # existing runs/<run_id> value explicitly to resume/retry into that same
    # run dir instead of starting a fresh one.
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
    maf = getattr(cfg, "maf", 0.01)
    mediators = getattr(cfg, "mediators", False)
    mediator_manifest = getattr(cfg, "mediator_manifest", "")
    info_threshold = getattr(cfg, "info_threshold", None)
    info_col = getattr(cfg, "info_col", None)
    remove_mhc = getattr(cfg, "remove_mhc", True)
    remove_apoe = getattr(cfg, "remove_apoe", False)
    overwrite = getattr(cfg, "overwrite", False)
    run_smr = getattr(cfg, "run_smr", True)
    bulk_eqtl_datasets = getattr(cfg, "bulk_eqtl_datasets", [])
    sc_eqtl_dataset = getattr(cfg, "sc_eqtl_dataset", "")

    # cis-MR / coloc gate thresholds - see params/schema.json's gates block;
    # defaults match what bin/coloc_targets.py used to hardcode
    wald_fdr_q = cfg.gate("cis_mr", "wald_fdr_q", 0.05)
    ivw_fdr_q = cfg.gate("cis_mr", "ivw_fdr_q", 0.05)
    cochran_q_pval = cfg.gate("cis_mr", "cochran_q_pval", 0.05)
    egger_intercept_pval_min = cfg.gate("cis_mr", "egger_intercept_pval_min", 0)
    min_instruments_for_ivw = cfg.gate("cis_mr", "min_instruments_for_ivw", 3)
    pp4_threshold = cfg.gate("coloc", "pp4_threshold", 0.7)
    m_y_pval_threshold = cfg.gate("network_mr", "m_y_pval_threshold", 0.05)

    print("[TRACKING] Preparing Falcon repo...")
    clone_repo(falcon_user)

    print("[TRACKING] Preparing Falcon env...")
    container_checks(falcon_user)

    # run_id uses the REMOTE repo's HEAD (post clone_repo() reset), since that's the
    # code version that actually executes the pipeline - not this local machine's HEAD.
    # Deterministic for a given (pheno_id, pqtl_dataset, day, remote commit): rerunning
    # today against the same remote commit reuses the same runs/<run_id>/ dir (and its
    # check_remote_output() skip behavior) both on Falcon and in the pulled-down local
    # copy - unless run_id is passed explicitly, in which case that existing run dir is
    # reused as-is.
    remote, _ = get_remote_paths(falcon_user)
    git_sha_result = ssh(f'cd "{remote}" && git rev-parse --short=7 HEAD', falcon_user)
    git_sha7 = git_sha_result.stdout.strip()
    date_str = datetime.now().strftime("%Y%m%d")
    if run_id is None:
        run_id = paths.make_run_id(pheno_id, pqtl_dataset, date_str, git_sha7)
    elif not run_id.startswith(f"{pheno_id}_{pqtl_dataset}_"):
        raise ValueError(
            f"run_id {run_id!r} does not match config's pheno_id={pheno_id!r}, "
            f"pqtl_dataset={pqtl_dataset!r} - refusing to write into a run dir "
            "for a different (pheno_id, pqtl_dataset) pair."
        )
    out_dir = str(paths.run_results_dir(run_id))
    local_results_dir = out_dir
    print(f"[TRACKING] run_id: {run_id}")

    # define all outputs first so pipeline knows what has already been ran
    qc_out = str(paths.qc_out(pheno_id))
    mr_out = str(paths.mr_out(pqtl_dataset, pheno_id, out_dir))
    coloc_out = str(paths.coloc_out(pqtl_dataset, pheno_id, out_dir))

    # PWCoCo (conditional coloc) - runs alongside coloc_out above, not instead of it
    # (see project_pwcoco_wiring memory); a failure here is logged and does not halt
    # the run, since standard COLOC is the required path and PWCoCo is a complementary
    # annotation on top of it
    pwcoco_out = str(paths.pwcoco_out(pqtl_dataset, pheno_id, out_dir))

    # PWCoCo (eQTL-informed) - eQTL-pQTL / eQTL-GWAS on SMR-passing targets, gated
    # on the same complementary-not-required basis as pwcoco_out above
    pwcoco_qtl_out = str(paths.pwcoco_eqtl_pqtl_out(pqtl_dataset, pheno_id, out_dir))

    # change this where NetworkMR saves its final compiled output
    # NetworkMR is gated directly on its actual final output (the mediation
    # estimates file coloc_with_mediators() also reads) - no separate gate literal
    network_mr_out = str(paths.network_mr_mediation_estimates_out(pqtl_dataset, pheno_id, out_dir))
    target_stats_out = str(paths.target_stats_out(pqtl_dataset, pheno_id, out_dir))

    # SMR (bulk and/or single-cell) - promising target output per eQTL mode
    # bulk_eqtl_datasets is a list (MetaBrain / GTEx_v10 etc. are pre-computed
    # separately under results/SMR/bulk/{dataset}/) so its per-dataset outputs are built
    # inside the SMR step below rather than up front here
    smr_sc_out = str(paths.smr_sc_out(pqtl_dataset, pheno_id, sc_eqtl_dataset, out_dir))

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
            out_dir=str(paths.qc_out(pheno_id).parent),
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
            out_dir=out_dir,
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
                local_results_dir=out_dir,
                ivw_fdr_q=ivw_fdr_q,
                egger_intercept_pval_min=egger_intercept_pval_min,
                cochran_q_pval=cochran_q_pval,
                m_y_pval_threshold=m_y_pval_threshold,
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
                mediator_manifest=mediator_manifest,
                local_results_dir=out_dir,
                ivw_fdr_q=ivw_fdr_q,
                pp4_threshold=pp4_threshold
            )
        else:
            run_coloc_without_mediators(
                falcon_user=falcon_user,
                pqtl_dataset=pqtl_dataset,
                pheno_id=pheno_id,
                n_cases=n_cases,
                n_controls=n_controls,
                local_results_dir=out_dir,
                wald_fdr_q=wald_fdr_q,
                ivw_fdr_q=ivw_fdr_q,
                cochran_q_pval=cochran_q_pval,
                egger_intercept_pval_min=egger_intercept_pval_min,
                min_instruments_for_ivw=min_instruments_for_ivw
            )

    require_remote_output(
        falcon_user=falcon_user,
        path=coloc_out,
        step="COLOC",
        required_for="Top cis-hit compilation"
    )

    # PWCoCo (conditional coloc) - complementary to standard COLOC above, not a
    # replacement (see project_pwcoco_wiring memory): runs on the same cis-MR-passing
    # targets and its results are joined against coloc_out downstream (dashboard
    # coloc_support annotation), not used to gate anything in this orchestration.
    if not check_remote_output(
        falcon_user=falcon_user,
        path=pwcoco_out,
        step="PWCoCo",
        overwrite=overwrite
    ):
        print("[TRACKING] Running PWCoCo...")
        try:
            run_pwcoco(
                falcon_user=falcon_user,
                pqtl_dataset=pqtl_dataset,
                pheno_id=pheno_id,
                ref_bfile=ref_bfile,
                n_cases=n_cases,
                n_controls=n_controls,
                local_results_dir=out_dir,
                cochran_q_pval=cochran_q_pval,
                wald_fdr_q=wald_fdr_q,
            )
        except subprocess.CalledProcessError as error:
            print(f"[CONCERN] PWCoCo run failed - continuing without it: {error}")

    # compile final hits
    if not check_remote_output(
        falcon_user=falcon_user,
        path=target_stats_out,
        step="Top cis-hit compilation",
        overwrite=overwrite
    ):
        print("[TRACKING] Compiling harmonised top cis-hit table...")
        compile_top_hits(
            falcon_user=falcon_user,
            pheno_id=pheno_id,
            pqtl_dataset=pqtl_dataset,
            local_results_dir=out_dir
        )

    require_remote_output(
        falcon_user=falcon_user,
        path=target_stats_out,
        step="Top cis-hit compilation",
        required_for="Dashboard target information"
    )

    # SMR module (bulk and/or single-cell eQTL, run right after coloc + top-cis-hit compilation)
    # -> targets which survive cis-MR + COLOC are checked against SMR + HEIDI in the
    #    configured eQTL dataset(s), alleles aligned to the AD risk allele
    if run_smr:
        if bulk_eqtl_datasets:
            # bulk eQTL SMR (MetaBrain / GTEx_v10) is pre-computed elsewhere -
            # bin/sort_smr.py ingests results/SMR/bulk/{dataset}/ rather than re-running SMR
            for bulk_dataset in bulk_eqtl_datasets:
                smr_bulk_out = str(paths.smr_bulk_out(pqtl_dataset, pheno_id, bulk_dataset, out_dir))

                if not check_remote_output(
                    falcon_user=falcon_user,
                    path=smr_bulk_out,
                    step=f"Bulk SMR ({bulk_dataset})",
                    overwrite=overwrite
                ):
                    print(f"[TRACKING] Ingesting pre-computed bulk eQTL SMR for {bulk_dataset}...")
                    run_smr_step(
                        falcon_user=falcon_user,
                        pqtl_dataset=pqtl_dataset,
                        eqtl_dataset=bulk_dataset,
                        eqtl_mode="bulk",
                        pheno_id=pheno_id,
                        sumstats=qc_out,
                        ref_bfile=ref_bfile,
                        maf=maf,
                        local_results_dir=out_dir
                    )
        else:
            print("[TRACKING] No bulk_eqtl_datasets specified, skipping bulk SMR.")

        if sc_eqtl_dataset:
            if not check_remote_output(
                falcon_user=falcon_user,
                path=smr_sc_out,
                step="Single-cell SMR",
                overwrite=overwrite
            ):
                print("[TRACKING] Running single-cell eQTL SMR...")
                run_smr_step(
                    falcon_user=falcon_user,
                    pqtl_dataset=pqtl_dataset,
                    eqtl_dataset=sc_eqtl_dataset,
                    eqtl_mode="single_cell",
                    pheno_id=pheno_id,
                    sumstats=qc_out,
                    ref_bfile=ref_bfile,
                    maf=maf,
                    local_results_dir=out_dir
                )
        else:
            print("[TRACKING] No sc_eqtl_dataset specified, skipping single-cell SMR.")
    else:
        print("[TRACKING] run_smr is False, skipping SMR entirely.")

    # PWCoCo (eQTL-informed) - eQTL-pQTL / eQTL-GWAS PWCoCo on every SMR-passing
    # target, then compared for shared colocalising SNPs against the pQTL-GWAS
    # PWCoCo above (see project_pwcoco_wiring memory) - runs only when SMR did,
    # since it depends on smr_final_targets_out; non-fatal like PWCoCo above.
    if run_smr:
        if not check_remote_output(
            falcon_user=falcon_user,
            path=pwcoco_qtl_out,
            step="PWCoCo (eQTL)",
            overwrite=overwrite
        ):
            print("[TRACKING] Running PWCoCo (eQTL)...")
            try:
                run_pwcoco_qtl(
                    falcon_user=falcon_user,
                    pqtl_dataset=pqtl_dataset,
                    pheno_id=pheno_id,
                    ref_bfile=ref_bfile,
                    n_cases=n_cases,
                    n_controls=n_controls,
                    local_results_dir=out_dir,
                )
            except subprocess.CalledProcessError as error:
                print(f"[CONCERN] PWCoCo (eQTL) run failed - continuing without it: {error}")

    # HyPrColoc (bulk and/or single-cell eQTL) - run right after SMR so the
    # combined final multi-omics target table (bulk + single-cell) is complete.
    # Each dataset is run (and gated) independently so bulk and single-cell compose.
    hyprcoloc_eqtl_datasets = list(bulk_eqtl_datasets) + ([sc_eqtl_dataset] if sc_eqtl_dataset else [])

    if run_smr and hyprcoloc_eqtl_datasets:
        for hc_dataset in hyprcoloc_eqtl_datasets:
            hc_dataset_out = str(paths.hyprcoloc_dataset_out(pqtl_dataset, hc_dataset, pheno_id, out_dir))

            if not check_remote_output(
                falcon_user=falcon_user,
                path=hc_dataset_out,
                step=f"HyPrColoc ({hc_dataset})",
                overwrite=overwrite
            ):
                print(f"[TRACKING] Running HyPrColoc for {hc_dataset}...")
                run_hyprcoloc_step(
                    falcon_user=falcon_user,
                    pqtl_dataset=pqtl_dataset,
                    pheno_id=pheno_id,
                    eqtl_dataset=hc_dataset,
                    local_results_dir=out_dir
                )
    else:
        print("[TRACKING] No bulk_eqtl_datasets or sc_eqtl_dataset specified (or run_smr is False), skipping HyPrColoc.")

    print("[TRACKING] Checking outputs...")
    check_outputs(
        falcon_user=falcon_user,
        pqtl_dataset=pqtl_dataset,
        pheno_id=pheno_id,
        local_results_dir=out_dir
    )

    print("[TRACKING] Pulling results locally...")
    pull_results_local(
        falcon_user=falcon_user,
        pqtl_dataset=pqtl_dataset,
        pheno_id=pheno_id,
        local_results_dir=local_results_dir,
        overwrite=overwrite,
    )

    print("[TRACKING] Running FinnGen PheWAS safety analysis locally...")
    phewas_safety_finngen(
        pheno_id=pheno_id,
        pqtl_dataset=pqtl_dataset,
        local_results_dir=local_results_dir,
        overwrite=overwrite,
    )

    print("[TRACKING] Running UKBB PheWAS safety analysis locally...")
    phewas_safety_ukbb(
        pheno_id=pheno_id,
        pqtl_dataset=pqtl_dataset,
        local_results_dir=local_results_dir,
        overwrite=overwrite,
    )

    print(f"[TRACKING] Expected cis-MR output: {mr_out}")
    print(f"[TRACKING] Expected pairwise COLOC output: {coloc_out}")
    print(f"[TRACKING] Expected top cis-hit output: {target_stats_out}")

    # Reached only if every check_remote_output()/require_remote_output() gate above
    # passed - an ssh/apptainer failure earlier raises and this is never reached, so
    # the registry can never point at a partial/failed run.
    project_root = Path(__file__).resolve().parents[1]
    registry.write_manifest(
        run_id,
        {
            "pheno_id": pheno_id,
            "pqtl_dataset": pqtl_dataset,
            "git_sha7": git_sha7,
            "date": date_str,
            "created_at": datetime.now().isoformat(),
            "mode": "hpc",
            "falcon_user": falcon_user,
            "overwrite": overwrite,
        },
        root=str(project_root / "runs"),
    )
    registry.record_successful_run(pheno_id, pqtl_dataset, run_id, root=str(project_root / "runs"))
    print(f"[TRACKING] Recorded successful run in registry: {run_id}")
    print("[DONE] drugMR pipeline completed successfully.")
