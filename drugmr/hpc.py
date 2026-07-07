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

def ssh(cmd: str, falcon_user: str):
    full_cmd = f"ssh {falcon_user}@falconlogin.cf.ac.uk '{cmd}'"
    result = subprocess.run(full_cmd, shell=True, executable="/bin/bash", capture_output=True, text=True)
    if result.stdout:
        print(result.stdout)
    if result.returncode != 0:
        print("[ERROR] Falcon command failed.")
        print(result.stderr)
        raise subprocess.CalledProcessError(result.returncode, full_cmd)

def get_remote_paths(falcon_user: str):
    remote = f"/shared/home1/{falcon_user}/drugMR"
    sif = f"{remote}/env/drugmr.sif"
    return remote, sif

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
    git clean -fd
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

apptainer exec --bind "{remote}:/work" "{sif}" \\
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

apptainer exec --bind "{remote}:/work" "{sif}" \\
bash -c "cd /work && python bin/coloc_targets.py \\
  --pqtl_dataset {pqtl_dataset} \\
  --local_results_dir results/cis-MR \\
  --pqtl_dir dat/cis_regions \\
  --pheno_id {pheno_id} \\
  --n_cases {n_cases} \\
  --n_controls {n_controls} \\
  --mediators \\
  --mediator_manifest {mediator_manifest}"
""", falcon_user)



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
    local_results_dir: str = "results",
    overwrite: bool = True
):
    remote, _ = get_remote_paths(falcon_user)
    remote_mr = f"{remote}/results/cis-MR/{pqtl_dataset}_{pheno_id}_all_MR.tsv"
    remote_coloc = f"{remote}/results/coloc/{pqtl_dataset}/{pqtl_dataset}_{pheno_id}_all_coloc.tsv"
    local_results_dir = Path(local_results_dir)
    local_mr_dir = local_results_dir / "cis-MR"
    local_coloc_dir = local_results_dir / "coloc" / pqtl_dataset
    local_mr_dir.mkdir(parents=True, exist_ok=True)
    local_coloc_dir.mkdir(parents=True, exist_ok=True)
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
    coloc_res = f"results/coloc/{pqtl_dataset}/{pqtl_dataset}_{pheno_id}_all_coloc.tsv"

    ssh(f"""
set -euo pipefail
cd "{remote}"

echo "[TRACKING] Checking MR output..."
ls -lh results/cis-MR/
head -5 "{mr_res}"

echo "[TRACKING] Checking COLOC output..."
ls -lh results/coloc/{pqtl_dataset}/
head -5 "{coloc_res}"
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
    overwrite = getattr(cfg, "overwrite", True)

    print("[TRACKING] Preparing Falcon repo...")
    clone_repo(falcon_user)

    print("[TRACKING] Preparing Falcon env...")
    container_checks(falcon_user)

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

    print("[TRACKING] Preparing cis-regions...")
    prep_cis_regions(
        falcon_user=falcon_user,
        pheno_id=pheno_id,
        pqtl_dataset=pqtl_dataset,
        pqtl_dir=pqtl_dir,
    )

    print("[TRACKING] Running cis-MR...")
    run_cis_mr(
        falcon_user=falcon_user,
        pqtl_dataset=pqtl_dataset,
        pqtl_dir=f"dat/cis_regions/{pqtl_dataset}",
        pheno_id=pheno_id,
        pheno_gwas=f"results/QC/{pheno_id}/{pheno_id}.tsv",
        ref_bfile=ref_bfile,
    )

    if mediators:
        print("[TRACKING] Running NetworkMR with mediators...")
        run_network_mr(
            falcon_user=falcon_user,
            pheno_id=pheno_id,
            pheno_gwas=f"results/QC/{pheno_id}/{pheno_id}.tsv",
            ref_bfile=ref_bfile,
            pqtl_dataset=pqtl_dataset,
            pqtl_dir=pqtl_dir,
        )
    else:
        print("[TRACKING] No mediators specified, skipping NetworkMR.")

    # ******** RE-DO -> ADD if mediators:
    # ******** RE-DO -> ADD if mediators:
    # ******** RE-DO -> ADD if mediators:
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

    print("[TRACKING] Checking outputs...")
    check_outputs(
        falcon_user=falcon_user,
        pqtl_dataset=pqtl_dataset,
        pheno_id=pheno_id,
    )

    print("[TRACKING] Pulling results locally...")
    pull_results_local(
        falcon_user=falcon_user,
        pqtl_dataset=pqtl_dataset,
        pheno_id=pheno_id,
        local_results_dir=local_results_dir,
        overwrite=overwrite,
    )

    print("[DONE] drugMR pipeline completed successfully.")