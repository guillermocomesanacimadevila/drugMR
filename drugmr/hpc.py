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

def ssh(cmd: str, user: str, host: str, allowed_returncodes: tuple = (0,)):
    # cmd is passed as its own argv element, not embedded in a locally-shell-parsed
    # string - no local shell involved, so nothing in cmd (built from interpolated
    # pheno_id/pqtl_dataset/paths elsewhere in this module) can break out of local
    # quoting. ssh itself still hands cmd to the REMOTE host's shell to interpret,
    # which is inherent to how ssh runs a command - that's expected, not a shell=True
    # concern (cmd is our own multi-line bash script, not untrusted external input).
    result = subprocess.run(
        ["ssh", f"{user}@{host}", cmd],
        capture_output=True,
        text=True,
    )
    if result.stdout:
        print(result.stdout)
    if result.returncode not in allowed_returncodes:
        print("[ERROR] Remote command failed.")
        print(result.stderr)
        raise subprocess.CalledProcessError(result.returncode, cmd)
    return result

def get_remote_paths(user: str, remote_repo_root: str):
    remote = remote_repo_root.format(user=user)
    sif = f"{remote}/env/drugmr.sif"
    return remote, sif

def check_remote_output(
    user: str,
    host: str,
    remote_repo_root: str,
    path: str,
    step: str,
    overwrite: bool = False
):
    # run step if overwrite == True
    if overwrite:
        print(f"[TRACKING] Overwrite enabled - rerunning {step}...")
        return False

    remote, _ = get_remote_paths(user, remote_repo_root)

    result = ssh(f"""
set -euo pipefail
cd "{remote}"

if [ -s "{path}" ]; then
    echo "[TRACKING] {step} already completed: {path}"
    exit 0
fi

exit 3
""", user, host, allowed_returncodes=(0, 3))

    if result.returncode == 0:
        print(f"[TRACKING] Skipping {step}...")
        return True

    print(f"[TRACKING] No existing {step} output found - running step...")
    return False

def check_remote_cis_regions(
    user: str,
    host: str,
    remote_repo_root: str,
    pqtl_dataset: str,
    overwrite: bool = False
):
    # run step if overwrite == True
    if overwrite:
        print("[TRACKING] Overwrite enabled - rerunning cis-region preparation...")
        return False

    remote, _ = get_remote_paths(user, remote_repo_root)

    result = ssh(f"""
set -euo pipefail
cd "{remote}"

n_cis=$(find "dat/cis_regions/{pqtl_dataset}" -mindepth 2 -maxdepth 2 -name "pqtl.parquet" 2>/dev/null | wc -l)

if [ "$n_cis" -gt 0 ]; then
    echo "[TRACKING] cis-regions already completed: $n_cis loci found"
    exit 0
fi

exit 3
""", user, host, allowed_returncodes=(0, 3))

    if result.returncode == 0:
        print("[TRACKING] Skipping cis-region preparation...")
        return True

    print("[TRACKING] No complete cis-region output found - running step...")
    return False

def require_remote_output(
    user: str,
    host: str,
    remote_repo_root: str,
    path: str,
    step: str,
    required_for: str
):
    remote, _ = get_remote_paths(user, remote_repo_root)

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
""", user, host)

def clone_repo(user: str, host: str, remote_repo_root: str):
    remote, _ = get_remote_paths(user, remote_repo_root)

    ssh(f"""
set -euo pipefail

echo 'Hello HPC!'
if [ -d "{remote}" ]; then
    echo "[TRACKING] I found the directory!"
    cd "{remote}"

    echo "[TRACKING] Resetting remote repo to GitHub main..."
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
    git clone https://github.com/guillermocomesanacimadevila/drugMR.git "{remote}"
fi
""", user, host)

def container_checks(user: str, host: str, remote_repo_root: str):
    remote, _ = get_remote_paths(user, remote_repo_root)

    ssh(f"""
set -euo pipefail

if [ ! -d "{remote}" ]; then
    git clone https://github.com/guillermocomesanacimadevila/drugMR.git "{remote}"
fi

cd "{remote}"
# git pull

chmod +x bin/bootstrap_hpc.sh
bash bin/bootstrap_hpc.sh "{remote}"
""", user, host)


# NOW -> FUNCTIONS TO RUN EACH SCRIPT FROM THE PIPELINE 

# **************************
# **************************
# ANALYTICS PIPELINE - START
# **************************
# **************************

# QC GWAS
def run_gwas_qc(
    user: str,
    host: str,
    remote_repo_root: str,
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
    remote, sif = get_remote_paths(user, remote_repo_root)

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
  --user {user} \\
  {info_args} \\
  {flag_args}"
""", user, host)


# *********** Extract cis-regions from pQTLs
def prep_cis_regions(
    user: str,
    host: str,
    remote_repo_root: str,
    pheno_id: str,
    pqtl_dataset: str,
):
    remote, sif = get_remote_paths(user, remote_repo_root)

    ssh(f"""
set -euo pipefail
cd "{remote}"

apptainer exec --bind "{remote}:/work" \\
  --env PYTHONPATH=. \\
  "{sif}" \\
bash -c "cd /work && python bin/prep_cis_regions.py \\
  --pqtl_dataset {pqtl_dataset} \\
  --pheno_id {pheno_id}"
""", user, host)


# RUN MR
def run_cis_mr(
    user: str,
    host: str,
    remote_repo_root: str,
    pqtl_dataset: str,
    pqtl_dir: str,
    pheno_id: str,
    pheno_gwas: str,
    ref_bfile: str,
    out_dir: str = "results",
    clump_kb: int = 10000,
    clump_r2: float = 0.001,
    instrument_pval_threshold: float = 5.0e-8,
    min_f_stat: float = 10,
    apply_steiger_filter: bool = False,
    maf: float = 0.01,
):
    remote, sif = get_remote_paths(user, remote_repo_root)

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
  {out_dir} \\
  {clump_kb} \\
  {clump_r2} \\
  {instrument_pval_threshold} \\
  {min_f_stat} \\
  {apply_steiger_filter} \\
  {maf}"
""", user, host)

# RUN COLOC
def run_coloc(
    user: str,
    host: str,
    remote_repo_root: str,
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
    pp4_threshold: float = 0.7,
    p1: float = 1e-4,
    p2: float = 1e-4,
    p12: float = 1e-5,
):
    remote, sif = get_remote_paths(user, remote_repo_root)

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
    --min_instruments_for_ivw {min_instruments_for_ivw} \\
    --pp4_threshold {pp4_threshold} \\
    --p1 {p1} \\
    --p2 {p2} \\
    --p12 {p12}"
""", user, host)


# RUN PWCoCo
def run_pwcoco(
    user: str,
    host: str,
    remote_repo_root: str,
    pqtl_dataset: str,
    pheno_id: str,
    ref_bfile: str,
    n_cases: int,
    n_controls: int,
    local_results_dir: str = "results",
    wald_fdr_q: float = 0.05,
    ivw_fdr_q: float = 0.05,
    cochran_q_pval: float = 0.05,
    egger_intercept_pval_min: float = 0,
    min_instruments_for_ivw: int = 3,
):
    remote, sif = get_remote_paths(user, remote_repo_root)

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
    --wald_fdr_q {wald_fdr_q} \\
    --ivw_fdr_q {ivw_fdr_q} \\
    --cochran_q_pval {cochran_q_pval} \\
    --egger_intercept_pval_min {egger_intercept_pval_min} \\
    --min_instruments_for_ivw {min_instruments_for_ivw}"
""", user, host)


# RUN PWCoCo (QTL-informed) - QTL-pQTL / QTL-GWAS PWCoCo on every SMR-passing
# target, then compared for shared colocalising SNPs against the pQTL-GWAS PWCoCo
# above (see project_pwcoco_wiring memory / bin/pwcoco_qtl_wrapper.py)
def run_pwcoco_qtl(
    user: str,
    host: str,
    remote_repo_root: str,
    pqtl_dataset: str,
    pheno_id: str,
    ref_bfile: str,
    n_cases: int,
    n_controls: int,
    local_results_dir: str = "results",
    pp4_threshold: float = 0.7,
):
    remote, sif = get_remote_paths(user, remote_repo_root)

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
    --local_results_dir {local_results_dir} \\
    --pp4_threshold {pp4_threshold}"
""", user, host)


# RUN SMR (bulk or single-cell, depending on qtl_mode)
# named run_smr_step (not run_smr) to avoid clashing with the run_smr config flag in hpc()
def run_smr_step(
    user: str,
    host: str,
    remote_repo_root: str,
    pqtl_dataset: str,
    qtl_dataset: str,
    qtl_mode: str,
    pheno_id: str,
    sumstats: str,
    ref_bfile: str,
    maf: float,
    local_results_dir: str = "results",
    wald_fdr_q: float = 0.05,
    ivw_fdr_q: float = 0.05,
    cochran_q_pval: float = 0.05,
    p_qtl_smr: float = 5.0e-8,
    p_qtl_heidi: float = 1.57e-3,
    p_smr_threshold: float = 0.05,
    p_heidi_threshold: float = 0.01,
):
    remote, sif = get_remote_paths(user, remote_repo_root)

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
    --qtl_dataset {qtl_dataset} \\
    --qtl_mode {qtl_mode} \\
    --ref_bfile {ref_bfile} \\
    --maf {maf} \\
    --local_results_dir {local_results_dir} \\
    --wald_fdr_q {wald_fdr_q} \\
    --ivw_fdr_q {ivw_fdr_q} \\
    --cochran_q_pval {cochran_q_pval} \\
    --p_qtl_smr {p_qtl_smr} \\
    --p_qtl_heidi {p_qtl_heidi} \\
    --p_smr_threshold {p_smr_threshold} \\
    --p_heidi_threshold {p_heidi_threshold}"
""", user, host)


# RUN HyPrColoc (bulk and/or single-cell QTL) - for every target x cell-type/tissue
# hit in the combined final multi-omics target table for the given qtl_dataset, runs
# a 3-trait (pQTL / GWAS / QTL) HyPrColoc restricted to that target's cis-region
def run_hyprcoloc_step(
    user: str,
    host: str,
    remote_repo_root: str,
    pqtl_dataset: str,
    pheno_id: str,
    qtl_dataset: str,
    local_results_dir: str = "results",
    prior_1: float = 1e-4,
    prior_c: list[float] = (0.05, 0.02, 0.01, 0.005),
    reg_thresh: list[float] = (0.5, 0.6, 0.7),
    align_thresh: list[float] = (0.5, 0.6, 0.7),
    equal_thresholds: bool = True,
):
    remote, sif = get_remote_paths(user, remote_repo_root)
    prior_c_arg = ",".join(str(v) for v in prior_c)
    reg_thresh_arg = ",".join(str(v) for v in reg_thresh)
    align_thresh_arg = ",".join(str(v) for v in align_thresh)

    ssh(f"""
set -euo pipefail
cd "{remote}"

apptainer exec --bind "{remote}:/work" \\
  --env PYTHONPATH=. \\
  "{sif}" \\
  bash -c "cd /work && python bin/hyprcoloc_targets.py \\
    --pqtl_dataset {pqtl_dataset} \\
    --pheno_id {pheno_id} \\
    --qtl_dataset {qtl_dataset} \\
    --local_results_dir {local_results_dir} \\
    --prior_1 {prior_1} \\
    --prior_c {prior_c_arg} \\
    --reg_thresh {reg_thresh_arg} \\
    --align_thresh {align_thresh_arg} \\
    --equal_thresholds {equal_thresholds}"
""", user, host)


# get final snp-wide hits
def compile_top_hits(
    user: str,
    host: str,
    remote_repo_root: str,
    pheno_id: str,
    pqtl_dataset: str,
    local_results_dir: str = "results"
):
    remote, sif = get_remote_paths(user, remote_repo_root)

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
""", user, host)

# RUN PHEWAS CHECKS FOR SAFETY (LOCALLY) -> API != WORK IN SLURM HPC
# ******************************************************************

def phewas_safety_finngen(
    pheno_id: str,
    pqtl_dataset: str,
    local_results_dir: str = "results",
    overwrite: bool = False,
    coloc_threshold: float = 0,
    bonferroni_alpha: float = 0.05,
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

    cmd = [
        "python", "bin/phewas_cis_pqtls.py",
        "--pheno_id", pheno_id,
        "--pqtl_dataset", pqtl_dataset,
        "--local_results_dir", str(local_results_dir),
        "--coloc_threshold", str(coloc_threshold),
        "--bonferroni_alpha", str(bonferroni_alpha),
    ]

    print(f"[TRACKING] FinnGen PheWAS pairwise COLOC input found: {top_snp_file}")
    print("[TRACKING] Running FinnGen PheWAS safety analysis locally...")
    subprocess.run(cmd, check=True, cwd=str(project_root))
    print(f"[TRACKING] FinnGen PheWAS safety results found: {phewas_out}")


def phewas_safety_ukbb(
    pheno_id: str,
    pqtl_dataset: str,
    local_results_dir: str = "results",
    overwrite: bool = False,
    coloc_threshold: float = 0,
    bonferroni_alpha: float = 0.05,
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

    cmd = [
        "python", "bin/ukb_phewas.py",
        "--pheno_id", pheno_id,
        "--pqtl_dataset", pqtl_dataset,
        "--local_results_dir", str(local_results_dir),
        "--coloc_threshold", str(coloc_threshold),
        "--bonferroni_alpha", str(bonferroni_alpha),
    ]

    print(f"[TRACKING] UKBB PheWAS pairwise COLOC input found: {top_snp_file}")
    print("[TRACKING] Running UKBB PheWAS safety analysis locally...")
    subprocess.run(cmd, check=True, cwd=str(project_root))
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
    user: str,
    host: str,
    remote_repo_root: str,
    run_id: str,
    pqtl_dataset: str,
    pheno_id: str,
    db_id: str = "drugmr",
    local_results_dir: str = "results"
):
    remote, sif = get_remote_paths(user, remote_repo_root)
    mr_res = str(paths.mr_out(pqtl_dataset, pheno_id, local_results_dir))
    coloc_res = str(paths.coloc_out(pqtl_dataset, pheno_id, local_results_dir))

    ssh(f"""
set -euo pipefail
cd "{remote}"

apptainer exec --bind "{remote}:/work" "{sif}" \\
bash -c "cd /work && python bin/load_db_into_postgres.py \\
  --results_file {mr_res} \\
  --db_id {db_id} \\
  --run_id {run_id} \\
  --pqtl_dataset {pqtl_dataset} \\
  --table cis_mr_results"

apptainer exec --bind "{remote}:/work" "{sif}" \\
bash -c "cd /work && python bin/load_db_into_postgres.py \\
  --results_file {coloc_res} \\
  --db_id {db_id} \\
  --run_id {run_id} \\
  --pqtl_dataset {pqtl_dataset} \\
  --table coloc_results"
""", user, host)


# PULL RESULTS INTO LOCAL
def pull_results_local(
    user: str,
    host: str,
    remote_repo_root: str,
    run_id: str,
    pqtl_dataset: str,
    pheno_id: str,
    overwrite: bool = True,
):
    remote, _ = get_remote_paths(user, remote_repo_root)
    relative_results_dir = paths.run_results_dir(run_id)
    remote_results_dir = f"{remote}/{relative_results_dir}"

    project_root = Path(__file__).resolve().parents[1]
    local_results_dir = project_root / relative_results_dir

    required_outputs = (
        ("cis-MR", paths.mr_out),
        ("colocalisation", paths.coloc_out),
        ("target statistics", paths.target_stats_out),
    )

    # Check the current remote run before copying. Existing local files
    # must not conceal missing outputs on the cluster.
    for step, output_path in required_outputs:
        require_remote_output(
            user=user,
            host=host,
            remote_repo_root=remote_repo_root,
            path=str(
                output_path(
                    pqtl_dataset,
                    pheno_id,
                    out_dir=str(relative_results_dir),
                )
            ),
            step=step,
            required_for="Result retrieval",
        )

    local_results_dir.mkdir(parents=True, exist_ok=True)

    cmd = ["rsync", "-avz"]
    if not overwrite:
        cmd.append("--ignore-existing")

    cmd.extend([
        f"{user}@{host}:{remote_results_dir}/",
        f"{local_results_dir}/",
    ])

    print(
        f"[TRACKING] Fetching {user}@{host}:{remote_results_dir}/ "
        f"-> {local_results_dir}"
    )
    subprocess.run(cmd, check=True)

    for step, output_path in required_outputs:
        local_file = output_path(
            pqtl_dataset,
            pheno_id,
            out_dir=str(local_results_dir),
        )
        if not local_file.is_file() or local_file.stat().st_size == 0:
            raise RuntimeError(
                f"{step} output missing or empty after retrieval: {local_file}"
            )

    print(f"[DONE] Retrieved results for run: {run_id}")


# STREAMLIT DASHBOARD
def run_dashboard_local(
    db_name: str,
    phenotype: str,
    pqtl_dataset: str,
    port_number: int = 5433
):
    # cwd is pinned to project_root (same fix as drugmr/local.py's results()) so
    # Streamlit reliably finds <project_root>/.streamlit/config.toml (the custom
    # theme) regardless of the caller's own working directory, instead of
    # silently falling back to Streamlit defaults.
    project_root = Path(__file__).resolve().parents[1]
    cmd = [
        "python", "-m", "streamlit", "run", "dashboard/mr_app.py", "--",
        "--db_name", str(db_name),
        "--port_number", str(port_number),
        "--phenotype", phenotype,
        "--pqtl_dataset", pqtl_dataset,
    ]
    print(cmd)
    subprocess.run(cmd, check=True, cwd=str(project_root))


# CHECK OUTPUTS
def check_outputs(
    user: str,
    host: str,
    remote_repo_root: str,
    pqtl_dataset: str,
    pheno_id: str,
    local_results_dir: str = "results"
):
    remote, _ = get_remote_paths(user, remote_repo_root)
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

""", user, host)


# Function to run all the HPC gist
def hpc(
    config: str,
    user: str,
    host: str,
    remote_repo_root: str,
    run_id: str = None,
):
    # config has no default on purpose - there's no single correct params file
    # anymore now that each (pheno_id, pqtl_dataset) pair has its own under
    # params/ (e.g. params/AD.wingo_brain.yaml) - pass one explicitly.
    #
    # user is a per-invocation credential, not an analysis parameter -
    # it doesn't belong in a (pheno_id, pqtl_dataset) params file, so it's a
    # real argument here instead of cfg.user.
    #
    # run_id defaults to None, which keeps the deterministic
    # (pheno_id, pqtl_dataset, day, remote commit) behaviour below. Pass an
    # existing runs/<run_id> value explicitly to resume/retry into that same
    # run dir instead of starting a fresh one.
    #
    # Connection details are required per invocation and passed to each helper.
    # remote_repo_root can use {user} as a placeholder for the username.
    cfg = Config(config)
    pheno_id = cfg.pheno_id
    sumstats = cfg.sumstats
    n_cases = cfg.n_cases
    n_controls = cfg.n_controls
    pqtl_dataset = cfg.pqtl_dataset
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
    info_threshold = getattr(cfg, "info_threshold", None)
    info_col = getattr(cfg, "info_col", None)
    remove_mhc = getattr(cfg, "remove_mhc", True)
    remove_apoe = getattr(cfg, "remove_apoe", False)
    overwrite = getattr(cfg, "overwrite", False)
    run_smr = getattr(cfg, "run_smr", True)
    bulk_qtl_datasets = getattr(cfg, "bulk_qtl_datasets", [])
    sc_qtl_dataset = getattr(cfg, "sc_qtl_dataset", "")

    # cis-MR / coloc gate thresholds - see params/schema.json's gates block;
    # defaults match what bin/coloc_targets.py used to hardcode
    wald_fdr_q = cfg.gate("cis_mr", "wald_fdr_q", 0.05)
    ivw_fdr_q = cfg.gate("cis_mr", "ivw_fdr_q", 0.05)
    cochran_q_pval = cfg.gate("cis_mr", "cochran_q_pval", 0.05)
    egger_intercept_pval_min = cfg.gate("cis_mr", "egger_intercept_pval_min", 0)
    min_instruments_for_ivw = cfg.gate("cis_mr", "min_instruments_for_ivw", 3)
    apply_steiger_filter = cfg.gate("cis_mr", "apply_steiger_filter", False)
    clump_kb = cfg.gate("cis_mr", "clump_kb", 10000)
    clump_r2 = cfg.gate("cis_mr", "clump_r2", 0.001)
    instrument_pval_threshold = cfg.gate("cis_mr", "instrument_pval_threshold", 5.0e-8)
    min_f_stat = cfg.gate("cis_mr", "min_f_stat", 10)
    pp4_threshold = cfg.gate("coloc", "pp4_threshold", 0.7)
    p1 = cfg.gate("coloc", "p1", 1e-4)
    p2 = cfg.gate("coloc", "p2", 1e-4)
    p12 = cfg.gate("coloc", "p12", 1e-5)
    p_qtl_smr = cfg.gate("smr", "p_qtl_smr", 5.0e-8)
    p_qtl_heidi = cfg.gate("smr", "p_qtl_heidi", 1.57e-3)
    p_smr_threshold = cfg.gate("smr", "p_smr_threshold", 0.05)
    p_heidi_threshold = cfg.gate("smr", "p_heidi_threshold", 0.01)
    hc_prior_1 = cfg.gate("hyprcoloc", "prior_1", 1e-4)
    hc_prior_c = cfg.gate("hyprcoloc", "prior_c", [0.05, 0.02, 0.01, 0.005])
    hc_reg_thresh = cfg.gate("hyprcoloc", "reg_thresh", [0.5, 0.6, 0.7])
    hc_align_thresh = cfg.gate("hyprcoloc", "align_thresh", [0.5, 0.6, 0.7])
    hc_equal_thresholds = cfg.gate("hyprcoloc", "equal_thresholds", True)
    pwcoco_pp4_threshold = cfg.gate("pwcoco", "pp4_threshold", 0.7)
    bonferroni_alpha = cfg.gate("phewas", "bonferroni_alpha", 0.05)
    phewas_coloc_threshold = cfg.gate("phewas", "coloc_threshold", 0)

    print("[TRACKING] Preparing remote repo...")
    clone_repo(user, host, remote_repo_root)

    print("[TRACKING] Preparing remote env...")
    container_checks(user, host, remote_repo_root)

    # run_id uses the REMOTE repo's HEAD (post clone_repo() reset), since that's the
    # code version that actually executes the pipeline - not this local machine's HEAD.
    # Deterministic for a given (pheno_id, pqtl_dataset, day, remote commit): rerunning
    # today against the same remote commit reuses the same runs/<run_id>/ dir (and its
    # check_remote_output() skip behavior) both on the cluster and in the pulled-down local
    # copy - unless run_id is passed explicitly, in which case that existing run dir is
    # reused as-is.
    remote, _ = get_remote_paths(user, remote_repo_root)
    git_sha_result = ssh(f'cd "{remote}" && git rev-parse --short=7 HEAD', user, host)
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

    # PWCoCo (QTL-informed) - QTL-pQTL / QTL-GWAS on SMR-passing targets, gated
    # on the same complementary-not-required basis as pwcoco_out above
    pwcoco_qtl_out = str(paths.pwcoco_eqtl_pqtl_out(pqtl_dataset, pheno_id, out_dir))

    target_stats_out = str(paths.target_stats_out(pqtl_dataset, pheno_id, out_dir))

    # SMR (bulk and/or single-cell) - promising target output per QTL mode
    # bulk_qtl_datasets is a list (MetaBrain / GTEx_v10 etc. are pre-computed
    # separately under results/SMR/bulk/{dataset}/) so its per-dataset outputs are built
    # inside the SMR step below rather than up front here
    smr_sc_out = str(paths.smr_sc_out(pqtl_dataset, pheno_id, sc_qtl_dataset, out_dir))

    if not check_remote_output(
        user=user,
        host=host,
        remote_repo_root=remote_repo_root,
        path=qc_out,
        step="GWAS QC",
        overwrite=overwrite
    ):
        print("[TRACKING] Running GWAS QC...")
        run_gwas_qc(
            user=user,
            host=host,
            remote_repo_root=remote_repo_root,
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
        user=user,
        host=host,
        remote_repo_root=remote_repo_root,
        path=qc_out,
        step="GWAS QC",
        required_for="cis-region preparation"
    )

    if not check_remote_cis_regions(
        user=user,
        host=host,
        remote_repo_root=remote_repo_root,
        pqtl_dataset=pqtl_dataset,
        overwrite=overwrite
    ):
        print("[TRACKING] Preparing cis-regions...")
        prep_cis_regions(
            user=user,
            host=host,
            remote_repo_root=remote_repo_root,
            pheno_id=pheno_id,
            pqtl_dataset=pqtl_dataset,
        )

    if not check_remote_output(
        user=user,
        host=host,
        remote_repo_root=remote_repo_root,
        path=mr_out,
        step="cis-MR",
        overwrite=overwrite
    ):
        print("[TRACKING] Running cis-MR...")
        run_cis_mr(
            user=user,
            host=host,
            remote_repo_root=remote_repo_root,
            pqtl_dataset=pqtl_dataset,
            pqtl_dir=f"dat/cis_regions/{pqtl_dataset}",
            pheno_id=pheno_id,
            pheno_gwas=qc_out,
            ref_bfile=ref_bfile,
            out_dir=out_dir,
            clump_kb=clump_kb,
            clump_r2=clump_r2,
            instrument_pval_threshold=instrument_pval_threshold,
            min_f_stat=min_f_stat,
            apply_steiger_filter=apply_steiger_filter,
            maf=maf,
        )

    require_remote_output(
        user=user,
        host=host,
        remote_repo_root=remote_repo_root,
        path=mr_out,
        step="cis-MR",
        required_for="COLOC"
    )

    if not check_remote_output(
        user=user,
        host=host,
        remote_repo_root=remote_repo_root,
        path=coloc_out,
        step="COLOC",
        overwrite=overwrite
    ):
        print("[TRACKING] Running COLOC...")
        run_coloc(
            user=user,
            host=host,
            remote_repo_root=remote_repo_root,
            pqtl_dataset=pqtl_dataset,
            pheno_id=pheno_id,
            n_cases=n_cases,
            n_controls=n_controls,
            local_results_dir=out_dir,
            wald_fdr_q=wald_fdr_q,
            ivw_fdr_q=ivw_fdr_q,
            cochran_q_pval=cochran_q_pval,
            egger_intercept_pval_min=egger_intercept_pval_min,
            min_instruments_for_ivw=min_instruments_for_ivw,
            pp4_threshold=pp4_threshold,
            p1=p1,
            p2=p2,
            p12=p12,
        )

    require_remote_output(
        user=user,
        host=host,
        remote_repo_root=remote_repo_root,
        path=coloc_out,
        step="COLOC",
        required_for="Top cis-hit compilation"
    )

    # PWCoCo (conditional coloc) - complementary to standard COLOC above, not a
    # replacement (see project_pwcoco_wiring memory): runs on the same cis-MR-passing
    # targets and its results are joined against coloc_out downstream (dashboard
    # coloc_support annotation), not used to gate anything in this orchestration.
    if not check_remote_output(
        user=user,
        host=host,
        remote_repo_root=remote_repo_root,
        path=pwcoco_out,
        step="PWCoCo",
        overwrite=overwrite
    ):
        print("[TRACKING] Running PWCoCo...")
        try:
            run_pwcoco(
                user=user,
                host=host,
                remote_repo_root=remote_repo_root,
                pqtl_dataset=pqtl_dataset,
                pheno_id=pheno_id,
                ref_bfile=ref_bfile,
                n_cases=n_cases,
                n_controls=n_controls,
                local_results_dir=out_dir,
                wald_fdr_q=wald_fdr_q,
                ivw_fdr_q=ivw_fdr_q,
                cochran_q_pval=cochran_q_pval,
                egger_intercept_pval_min=egger_intercept_pval_min,
                min_instruments_for_ivw=min_instruments_for_ivw,
            )
        except subprocess.CalledProcessError as error:
            print(f"[CONCERN] PWCoCo run failed - continuing without it: {error}")

    # compile final hits
    if not check_remote_output(
        user=user,
        host=host,
        remote_repo_root=remote_repo_root,
        path=target_stats_out,
        step="Top cis-hit compilation",
        overwrite=overwrite
    ):
        print("[TRACKING] Compiling harmonised top cis-hit table...")
        compile_top_hits(
            user=user,
            host=host,
            remote_repo_root=remote_repo_root,
            pheno_id=pheno_id,
            pqtl_dataset=pqtl_dataset,
            local_results_dir=out_dir
        )

    require_remote_output(
        user=user,
        host=host,
        remote_repo_root=remote_repo_root,
        path=target_stats_out,
        step="Top cis-hit compilation",
        required_for="Dashboard target information"
    )

    # SMR module (bulk and/or single-cell QTL, run right after coloc + top-cis-hit compilation)
    # -> targets which survive cis-MR + COLOC are checked against SMR + HEIDI in the
    #    configured QTL dataset(s), alleles aligned to the AD risk allele
    if run_smr:
        if bulk_qtl_datasets:
            # bulk QTL SMR (MetaBrain / GTEx_v10) is pre-computed elsewhere -
            # bin/sort_smr.py ingests results/SMR/bulk/{dataset}/ rather than re-running SMR
            for bulk_dataset in bulk_qtl_datasets:
                smr_bulk_out = str(paths.smr_bulk_out(pqtl_dataset, pheno_id, bulk_dataset, out_dir))

                if not check_remote_output(
                    user=user,
                    host=host,
                    remote_repo_root=remote_repo_root,
                    path=smr_bulk_out,
                    step=f"Bulk SMR ({bulk_dataset})",
                    overwrite=overwrite
                ):
                    print(f"[TRACKING] Ingesting pre-computed bulk QTL SMR for {bulk_dataset}...")
                    run_smr_step(
                        user=user,
                        host=host,
                        remote_repo_root=remote_repo_root,
                        pqtl_dataset=pqtl_dataset,
                        qtl_dataset=bulk_dataset,
                        qtl_mode="bulk",
                        pheno_id=pheno_id,
                        sumstats=qc_out,
                        ref_bfile=ref_bfile,
                        maf=maf,
                        local_results_dir=out_dir,
                        wald_fdr_q=wald_fdr_q,
                        ivw_fdr_q=ivw_fdr_q,
                        cochran_q_pval=cochran_q_pval,
                        p_qtl_smr=p_qtl_smr,
                        p_qtl_heidi=p_qtl_heidi,
                        p_smr_threshold=p_smr_threshold,
                        p_heidi_threshold=p_heidi_threshold,
                    )
        else:
            print("[TRACKING] No bulk_qtl_datasets specified, skipping bulk SMR.")

        if sc_qtl_dataset:
            if not check_remote_output(
                user=user,
                host=host,
                remote_repo_root=remote_repo_root,
                path=smr_sc_out,
                step="Single-cell SMR",
                overwrite=overwrite
            ):
                print("[TRACKING] Running single-cell QTL SMR...")
                run_smr_step(
                    user=user,
                    host=host,
                    remote_repo_root=remote_repo_root,
                    pqtl_dataset=pqtl_dataset,
                    qtl_dataset=sc_qtl_dataset,
                    qtl_mode="single_cell",
                    pheno_id=pheno_id,
                    sumstats=qc_out,
                    ref_bfile=ref_bfile,
                    maf=maf,
                    local_results_dir=out_dir,
                    wald_fdr_q=wald_fdr_q,
                    ivw_fdr_q=ivw_fdr_q,
                    cochran_q_pval=cochran_q_pval,
                    p_qtl_smr=p_qtl_smr,
                    p_qtl_heidi=p_qtl_heidi,
                    p_smr_threshold=p_smr_threshold,
                    p_heidi_threshold=p_heidi_threshold,
                )
        else:
            print("[TRACKING] No sc_qtl_dataset specified, skipping single-cell SMR.")
    else:
        print("[TRACKING] run_smr is False, skipping SMR entirely.")

    # PWCoCo (QTL-informed) - QTL-pQTL / QTL-GWAS PWCoCo on every SMR-passing
    # target, then compared for shared colocalising SNPs against the pQTL-GWAS
    # PWCoCo above (see project_pwcoco_wiring memory) - runs only when SMR did,
    # since it depends on smr_final_targets_out; non-fatal like PWCoCo above.
    if run_smr:
        if not check_remote_output(
            user=user,
            host=host,
            remote_repo_root=remote_repo_root,
            path=pwcoco_qtl_out,
            step="PWCoCo (QTL)",
            overwrite=overwrite
        ):
            print("[TRACKING] Running PWCoCo (QTL)...")
            try:
                run_pwcoco_qtl(
                    user=user,
                    host=host,
                    remote_repo_root=remote_repo_root,
                    pqtl_dataset=pqtl_dataset,
                    pheno_id=pheno_id,
                    ref_bfile=ref_bfile,
                    n_cases=n_cases,
                    n_controls=n_controls,
                    local_results_dir=out_dir,
                    pp4_threshold=pwcoco_pp4_threshold,
                )
            except subprocess.CalledProcessError as error:
                print(f"[CONCERN] PWCoCo (QTL) run failed - continuing without it: {error}")

    # HyPrColoc (bulk and/or single-cell QTL) - run right after SMR so the
    # combined final multi-omics target table (bulk + single-cell) is complete.
    # Each dataset is run (and gated) independently so bulk and single-cell compose.
    hyprcoloc_qtl_datasets = list(bulk_qtl_datasets) + ([sc_qtl_dataset] if sc_qtl_dataset else [])

    if run_smr and hyprcoloc_qtl_datasets:
        for hc_dataset in hyprcoloc_qtl_datasets:
            hc_dataset_out = str(paths.hyprcoloc_dataset_out(pqtl_dataset, hc_dataset, pheno_id, out_dir))

            if not check_remote_output(
                user=user,
                host=host,
                remote_repo_root=remote_repo_root,
                path=hc_dataset_out,
                step=f"HyPrColoc ({hc_dataset})",
                overwrite=overwrite
            ):
                print(f"[TRACKING] Running HyPrColoc for {hc_dataset}...")
                run_hyprcoloc_step(
                    user=user,
                    host=host,
                    remote_repo_root=remote_repo_root,
                    pqtl_dataset=pqtl_dataset,
                    pheno_id=pheno_id,
                    qtl_dataset=hc_dataset,
                    local_results_dir=out_dir,
                    prior_1=hc_prior_1,
                    prior_c=hc_prior_c,
                    reg_thresh=hc_reg_thresh,
                    align_thresh=hc_align_thresh,
                    equal_thresholds=hc_equal_thresholds,
                )
    else:
        print("[TRACKING] No bulk_qtl_datasets or sc_qtl_dataset specified (or run_smr is False), skipping HyPrColoc.")

    print("[TRACKING] Checking outputs...")
    check_outputs(
        user=user,
        host=host,
        remote_repo_root=remote_repo_root,
        pqtl_dataset=pqtl_dataset,
        pheno_id=pheno_id,
        local_results_dir=out_dir
    )

    print("[TRACKING] Pulling results locally...")
    pull_results_local(
        user=user,
        host=host,
        remote_repo_root=remote_repo_root,
        run_id=run_id,
        pqtl_dataset=pqtl_dataset,
        pheno_id=pheno_id,
        overwrite=overwrite,
    )

    print("[TRACKING] Running FinnGen PheWAS safety analysis locally...")
    phewas_safety_finngen(
        pheno_id=pheno_id,
        pqtl_dataset=pqtl_dataset,
        local_results_dir=local_results_dir,
        overwrite=overwrite,
        coloc_threshold=phewas_coloc_threshold,
        bonferroni_alpha=bonferroni_alpha,
    )

    print("[TRACKING] Running UKBB PheWAS safety analysis locally...")
    phewas_safety_ukbb(
        pheno_id=pheno_id,
        pqtl_dataset=pqtl_dataset,
        local_results_dir=local_results_dir,
        overwrite=overwrite,
        coloc_threshold=phewas_coloc_threshold,
        bonferroni_alpha=bonferroni_alpha,
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
            "user": user,
            "host": host,
            "remote_repo_root": remote,
            "overwrite": overwrite,
        },
        root=str(project_root / "runs"),
    )
    registry.record_successful_run(pheno_id, pqtl_dataset, run_id, root=str(project_root / "runs"))
    print(f"[TRACKING] Recorded successful run in registry: {run_id}")
    print("[DONE] drugMR pipeline completed successfully.")
