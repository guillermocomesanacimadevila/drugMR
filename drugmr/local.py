#!/usr/bin/env python3
import sys
import subprocess
from pathlib import Path
from drugmr.config import Config

# LOCAL RESULTS / DASHBOARD STUFF
# Need to create local running functions including a pulling docker from container function
# So then still the QC+MR runs in micromamba Docker env
# PostgreSQL db pulling and dashboard == jupyter (with .toml in ./)
# Goal == have a flagging variable (local/hpc)

def cmd_base(cmd):
    """
    Baseline function for parsing and running CLI-based .py scripts
    """
    return subprocess.run(cmd, shell=True, check=True, executable="/bin/bash")

def check_output(path: Path, step: str, overwrite: bool = False):
    # run step if overwrite == True
    if overwrite:
        print(f"[TRACKING] Overwrite enabled - rerunning {step}...")
        return False

    # run step if output does not exist
    if not path.exists():
        print(f"[TRACKING] No existing {step} output found - running step...")
        return False

    # run step if output exists but is empty
    if path.stat().st_size == 0:
        print(f"[CONCERN] {step} output exists but is empty - rerunning step...")
        return False

    print(f"[TRACKING] {step} already completed: {path}")
    print(f"[TRACKING] Skipping {step}...")
    return True

def check_cis_regions(cis_dir: Path, overwrite: bool = False):
    # run step if overwrite == True
    if overwrite:
        print("[TRACKING] Overwrite enabled - rerunning cis-region preparation...")
        return False

    # run step if cis-region directory does not exist
    if not cis_dir.exists():
        print("[TRACKING] No existing cis-region directory found - running step...")
        return False

    # check whether any pQTL cis-regions actually exist
    n_cis = len(list(cis_dir.glob("*/pqtl.parquet")))

    if n_cis == 0:
        print("[CONCERN] cis-region directory exists but no pqtl.parquet files were found - rerunning step...")
        return False

    print(f"[TRACKING] cis-regions already completed: {n_cis} loci found")
    print("[TRACKING] Skipping cis-region preparation...")
    return True

def require_output(path: Path, step: str, required_for: str):
    # do not run downstream step where required upstream output does not exist
    if not path.exists():
        raise FileNotFoundError(
            f"{required_for} cannot run because {step} output was not found: {path}"
        )

    # do not run downstream step where required upstream output is empty
    if path.stat().st_size == 0:
        raise RuntimeError(
            f"{required_for} cannot run because {step} output is empty: {path}"
        )

def results(
    config: str = "assets/config.yaml",
    db_id: str = "drugmr",
    dashboard_script: str = "dashboard/mr_app.py",
    db_script: str = "bin/load_db_into_postgres.py",
    port_number: int = 5432,
):
    project_root = Path(__file__).resolve().parents[1]
    cfg = Config(project_root / config)
    pqtl_dataset = cfg.pqtl_dataset
    pheno_id = cfg.pheno_id
    mr_res = project_root / "results" / "cis-MR" / f"{pqtl_dataset}_{pheno_id}_all_MR.tsv"
    coloc_res = project_root / "results" / "coloc" / pqtl_dataset / f"{pqtl_dataset}_{pheno_id}_all_coloc.tsv"
    db_script = project_root / db_script
    dashboard_script = project_root / dashboard_script

    require_output(mr_res, "cis-MR", "PostgreSQL loading")
    require_output(coloc_res, "COLOC", "PostgreSQL loading")

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

def local(config: str = "assets/config.yaml"):
    project_root = Path(__file__).resolve().parents[1]
    cfg = Config(project_root / config)
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
    eqtl_dataset = cfg.eqtl_dataset ###### Only SingleBrain ATM
    genome_build = cfg.genome_build
    target_build = cfg.target_build
    out_dir = getattr(cfg, "out_dir", "results")
    maf = getattr(cfg, "maf", 0.01)
    info_threshold = getattr(cfg, "info_threshold", None)
    info_col = getattr(cfg, "info_col", None)
    mediators = getattr(cfg, "mediators", False)
    mediator_manifest = getattr(cfg, "mediator_manifest", "")
    remove_mhc = getattr(cfg, "remove_mhc", True)
    remove_apoe = getattr(cfg, "remove_apoe", False)
    overwrite = getattr(cfg, "overwrite", False)
    image_uri = getattr(cfg, "image_uri", "ghcr.io/guillermocomesanacimadevila/drugmr:latest")
    image_name = getattr(cfg, "image_name", "ghcr.io/guillermocomesanacimadevila/drugmr:latest")

    # set projectDir()
    project_root = Path(__file__).resolve().parents[1] # i.e. "Users/.../drugMR"

    # define all outputs first so pipeline knows what has already been ran
    qc_out = project_root / out_dir / "QC" / pheno_id / f"{pheno_id}.tsv"
    cis_dir = project_root / "dat" / "cis_regions" / pqtl_dataset
    mr_out = project_root / "results" / "cis-MR" / f"{pqtl_dataset}_{pheno_id}_all_MR.tsv"
    coloc_out = project_root / "results" / "coloc" / pqtl_dataset / f"{pqtl_dataset}_{pheno_id}_all_coloc.tsv"
    promising_smr_out = (
        project_root
        / "results"
        / "SMR"
        / eqtl_dataset
        / pheno_id
        / f"{pqtl_dataset}_{pheno_id}_promising_targets_SMR.tsv"
    )
    final_smr_out = (
        project_root
        / "results"
        / "SMR"
        / eqtl_dataset
        / pheno_id
        / f"{pqtl_dataset}_{pheno_id}_final_multi_omics_targets.tsv"
    )

    # coloc and moloc with multi-omics
    prepared_multi_omics_out = project_root / "results" / "SMR" / eqtl_dataset / pheno_id / f"{pqtl_dataset}_{pheno_id}_prepared_multi_omics_targets.tsv"
    eqtl_coloc_out = project_root / "results" / "eQTL_coloc" / pqtl_dataset / eqtl_dataset / pheno_id / f"{pqtl_dataset}_{pheno_id}_{eqtl_dataset}_all_eqtl_coloc.tsv"
    moloc_out = project_root / "results" / "QTL_moloc" / pqtl_dataset / eqtl_dataset / pheno_id / f"{pqtl_dataset}_{pheno_id}_{eqtl_dataset}_moloc_summary.tsv"

    # FINAL TABLE (NO MEDIATORS)
    summary_out = (
        project_root
        / "results"
        / "SMR"
        / eqtl_dataset
        / pheno_id
        / f"{pqtl_dataset}_{pheno_id}_multi_omics_overview.tsv"
    )

    snp_evidence_out = (
        project_root
        / "results"
        / "SMR"
        / eqtl_dataset
        / pheno_id
        / f"{pqtl_dataset}_{pheno_id}_multi_omics_snp_evidence.tsv"
    )
    
    # change this where NetworkMR saves its final compiled output
    network_mr_out = (
        project_root
        / "results"
        / "network-MR"
        / pqtl_dataset
        / f"{pqtl_dataset}_{pheno_id}_network_MR.tsv"
    )

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
        info_args += f" --info-col {info_col}"
    if info_threshold is not None:
        info_args += f" --info-threshold {info_threshold}"

    flag_args = ""
    if remove_mhc:
        flag_args += " --remove_mhc"
    if remove_apoe:
        flag_args += " --remove_apoe"

    # running individual modules
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

    if not check_output(qc_out, "GWAS QC", overwrite):
        print("[TRACKING] Running GWAS QC locally via Docker...")
        cmd_base(cmd_qc)

    require_output(qc_out, "GWAS QC", "cis-region preparation")

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
        # mediator preparation does not currently have one definite output file
        # so this reruns where mediators == True
        cmd_base(cmd_m_qc)

    else:
        print("[TRACKING] No mediators specificed, running drugMR without them then!")

    # cis-region module
    cmd_cis = f"""
set -euo pipefail
docker run --rm \\
  --platform linux/amd64 \\
  -v "{project_root}:/work" \\
  -w /work \\
  "{image_name}" \\
  python bin/prep_cis_regions.py \\
    --pqtl_dataset {pqtl_dataset} \\
    --pheno_id {pheno_id} \\
    --pqtl_dir {pqtl_dir}
"""

    if not check_cis_regions(cis_dir, overwrite):
        print("[TRACKING] Preparing cis-regions locally...")
        cmd_base(cmd_cis)

    print(f"[TRACKING] Checking cis-region output: {cis_dir}")

    if not cis_dir.exists():
        raise FileNotFoundError(f"cis-region directory not created: {cis_dir}")

    n_cis = len(list(cis_dir.glob("*/pqtl.parquet")))
    print(f"[TRACKING] cis-region loci generated: {n_cis}")

    if n_cis == 0:
        raise RuntimeError("No cis-region files generated. Check pqtl_dir path.")

    # cis-MR module 
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

    if not check_output(mr_out, "cis-MR", overwrite):
        print("[TRACKING] Running cis-MR locally via Docker...")
        cmd_base(cmd_mr)

    require_output(mr_out, "cis-MR", "COLOC")

    # CMD COLOC TARGETS
    # CMD RUN COLOC (Pairwise)
    # Need to test


    # networkMR (HERE)
        # networkMR
    if mediators:
        require_output(mr_out, "cis-MR", "NetworkMR")

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
    --pqtl_dir dat/cis_regions/{pqtl_dataset} \\
    --run_genomewide_mr \\
    --run_cis_mr_X_M \\
    --run_network_mr
"""

        if not check_output(network_mr_out, "NetworkMR", overwrite):
            print("[TRACKING] Running NetworkMR with mediators!")
            cmd_base(cmd_network_mr)
    else:
        print("[TRACKING] No mediators specified, skipping NetworkMR.")


    # coloc target module
    # ******** RE-DO -> AD if mediators:
    # cmd_coloc with mediators
    cmd_coloc_with_mediators = f"""
set -euo pipefail
docker run --rm \\
  -v "{project_root}:/work" \\
  -w /work \\
  -e PYTHONPATH=. \\
  "{image_name}" \\
  python bin/coloc_targets.py \\
    --pqtl_dataset {pqtl_dataset} \\
    --local_results_dir results/cis-MR \\
    --pqtl_dir dat/cis_regions/{pqtl_dataset} \\
    --pheno_id {pheno_id} \\
    --n_cases {n_cases} \\
    --n_controls {n_controls} \\
    --mediators \\
    --mediator_manifest {mediator_manifest}
"""

    # without mediators
    cmd_coloc_without_mediators = f"""
set -euo pipefail
docker run --rm \\
  -v "{project_root}:/work" \\
  -w /work \\
  -e PYTHONPATH=. \\
  "{image_name}" \\
  python bin/coloc_targets.py \\
    --pqtl_dataset {pqtl_dataset} \\
    --local_results_dir results/cis-MR \\
    --pqtl_dir dat/cis_regions/{pqtl_dataset} \\
    --pheno_id {pheno_id} \\
    --n_cases {n_cases} \\
    --n_controls {n_controls}
"""

    if not check_output(coloc_out, "COLOC", overwrite):
        print("[TRACKING] Running COLOC locally...")

        if mediators:
            cmd_base(cmd_coloc_with_mediators)
        else:
            cmd_base(cmd_coloc_without_mediators)

    require_output(coloc_out, "COLOC", "single-cell SMR")

    # Integration with other omics layers 
    # ------ sc-eQTL ------
    # single-ceLL SMR
    # ------ ------- ------

    cmd_smr = f"""
set -euo pipefail 
docker run --rm \\
  -v "{project_root}:/work" \\
  -w /work \\
  -e PYTHONPATH=. \\
  "{image_name}" \\
  python bin/sort_single_cell_smr.py \\
    --pqtl_dataset {pqtl_dataset} \\
    --pheno_id {pheno_id} \\
    --eqtl_dataset {eqtl_dataset} \\
    --sumstats {out_dir}/QC/{pheno_id}/{pheno_id}.tsv \\
    --ref_bfile {ref_bfile} \\
    --maf {maf}
"""

    # SMR depends on both cis-MR and COLOC because it extracts promising targets
    require_output(mr_out, "cis-MR", "single-cell SMR")
    require_output(coloc_out, "COLOC", "single-cell SMR")

    # check final compiled output first
    # where this exists -> no need to rerun genome-wide SMR and compilation
    if not check_output(final_smr_out, "single-cell SMR", overwrite):
        print(f"[TRACKING] Runnig single-cell SMR for {eqtl_dataset}!")
        cmd_base(cmd_smr)

    # check both intermediate and final SMR outputs
    if promising_smr_out.exists() and promising_smr_out.stat().st_size > 0:
        print(f"[TRACKING] Promising target SMR results found: {promising_smr_out}")
    else:
        print(f"[CONCERN] Promising target SMR file not found or empty: {promising_smr_out}")

    if final_smr_out.exists() and final_smr_out.stat().st_size > 0:
        print(f"[TRACKING] Final multi-omics targets found: {final_smr_out}")
    else:
        print(f"[CONCERN] Final multi-omics target file not found or empty: {final_smr_out}")

        # ------ -------------------- ------
    # Multi-omics QTL colocalisation
    # GWAS - sc-eQTL pairwise coloc
    # GWAS - pQTL - sc-eQTL MOLOC
    # ------ -------------------- ------

    require_output(final_smr_out, "single-cell SMR", "multi-omics QTL colocalisation")

    cmd_multi_omics_coloc = f"""
set -euo pipefail
docker run --rm \\
  --platform linux/amd64 \\
  -v "{project_root}:/work" \\
  -w /work \\
  -e PYTHONPATH=. \\
  "{image_name}" \\
  python bin/assort_moloc_for_sc_hits.py \\
    --pheno_id {pheno_id} \\
    --pqtl_dataset {pqtl_dataset} \\
    --eqtl_dataset {eqtl_dataset} \\
    --n_cases {n_cases} \\
    --n_controls {n_controls}
"""

    multi_omics_complete = (
        prepared_multi_omics_out.exists() and prepared_multi_omics_out.stat().st_size > 0
        and eqtl_coloc_out.exists() and eqtl_coloc_out.stat().st_size > 0
        and moloc_out.exists() and moloc_out.stat().st_size > 0
    )

    if overwrite or not multi_omics_complete:
        print("[TRACKING] Running multi-omics QTL colocalisation locally...")
        cmd_base(cmd_multi_omics_coloc)
    else:
        print("[TRACKING] Multi-omics QTL colocalisation already completed.")
        print("[TRACKING] Skipping multi-omics QTL colocalisation...")

    require_output(prepared_multi_omics_out, "multi-omics target preparation", "pipeline completion")
    require_output(eqtl_coloc_out, "GWAS - sc-eQTL COLOC", "pipeline completion")
    require_output(moloc_out, "GWAS - pQTL - sc-eQTL MOLOC", "pipeline completion")
    print(f"[TRACKING] Prepared multi-omics target manifest found: {prepared_multi_omics_out}")
    print(f"[TRACKING] GWAS - sc-eQTL COLOC results found: {eqtl_coloc_out}")
    print(f"[TRACKING] GWAS - pQTL - sc-eQTL MOLOC results found: {moloc_out}")

    # compile (for top SMR SNP and top pQTL SNP)
    cmd_summary = f"""
set -euo pipefail
docker run --rm \\
  -v "{project_root}:/work" \\
  -w /work \\
  -e PYTHONPATH=. \\
  "{image_name}" \\
  python bin/summarise_multi_omics.py \\
    --pheno_id {pheno_id} \\
    --pqtl_dataset {pqtl_dataset} \\
    --eqtl_dataset {eqtl_dataset}
"""
    
    require_output(prepared_multi_omics_out, "multi-omics target preparation", "final summary")
    require_output(eqtl_coloc_out, "GWAS - sc-eQTL COLOC", "final summary")
    require_output(moloc_out, "GWAS - pQTL - sc-eQTL MOLOC", "final summary")
    summary_complete = (check_output(summary_out, "multi-omics overview", overwrite) and check_output(snp_evidence_out, "multi-omics SNP evidence", overwrite))
    if not summary_complete:
        print("[TRACKING] Building final dashboard-ready multi-omics tables...")
        cmd_base(cmd_summary)

    require_output(summary_out, "multi-omics overview", "pipeline completion")
    require_output(snp_evidence_out, "multi-omics SNP evidence", "pipeline completion")
    print(f"[TRACKING] Multi-omics overview found: {summary_out}")
    print(f"[TRACKING] SNP evidence table found: {snp_evidence_out}")
    print("[DONE] Local Docker run completed.")