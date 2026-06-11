#!/usr/bin/env python3
import sys
import subprocess
from pathlib import Path

# LOCAL RESULTS / DASHBOARD STUFF
# To do's (After Greece)
# Need to create local running functions including a pulling docker from container function
# So then still the QC+MR runs in micromamba Docker env
# PostgreSQL db pulling and dashboard == jupyter (with .toml in ./)
# Goal == have a flagging variable (local/hpc)

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
    db_script = project_root / db_script
    dashboard_script = project_root / dashboard_script

    print("[TRACKING] Loading MR results into PostgreSQL...")

    subprocess.run(
        [
            sys.executable,
            str(db_script),
            "--mr_res",
            str(mr_res),
            "--db_id",
            db_id,
            "--pqtl_dataset",
            pqtl_dataset,
            "--pheno_id",
            pheno_id,
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