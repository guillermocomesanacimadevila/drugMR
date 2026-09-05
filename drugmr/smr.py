import os
import subprocess
from pathlib import Path

from drugmr import paths



class SMRUtils:


    """
    Take single parquet or a sum of parquets within a dir (/*.parquet)
    and transform to ESD/BESD/EPI format for later intake during eQTL-informed workflows
    * Stuff to bare in mind *
    - bulk vs single-cell
    - multiple regions within 1 dataset (e.g. GTEx)
    """

    def __init__(self):
        pass

    @staticmethod
    def harmonise_qtl():
        return

    @staticmethod
    def transform_to_esd():
        return

    @staticmethod
    def transform_to_flist():
        return

    @staticmethod
    def transform_to_besd():

        cmd = [
            "smr"
        ]
        return

    @staticmethod
    def run_smr():
        return





def SMR(
    pheno_id: str,
    sumstats: str,
    ref_bfile: str,
    beqtl_summary: str,
    eqtl_dataset: str,
    peqtl_smr: float,
    peqtl_heidi: float,
    thread_num: int,
    maf: float,
    out_dir: str = "synthesis"
):
    ref_bfile = Path(ref_bfile)
    sumstats = Path(sumstats)
    beqtl_summary = Path(beqtl_summary)

    # eqtl_dataset can be stuff like SingleBrain/Ast
    # use full path for directory but only cell name for output prefix
    eqtl_dataset = Path(eqtl_dataset)
    # SMR(GWAS x eQTL) result for a given (pheno_id, eqtl_dataset) never depends on
    # pqtl_dataset, so this defaults to the shared synthesis/ tree rather than a
    # per-run out_dir - every pqtl_dataset run reuses the same computation instead
    # of re-running the smr binary from scratch
    raw_out_dir = paths.smr_raw_dir(eqtl_dataset, pheno_id, out_dir)
    os.makedirs(raw_out_dir, exist_ok=True)
    out_file = paths.smr_raw_prefix(eqtl_dataset, pheno_id, out_dir)
    print(f"[TRACKING] Running SMR on {pheno_id} using {eqtl_dataset}")

    cmd_smr = [
        "smr",
        "--bfile", str(ref_bfile),
        "--gwas-summary", str(sumstats),
        "--beqtl-summary", str(beqtl_summary),
        "--maf", str(maf),
        "--peqtl-smr", str(peqtl_smr),
        "--peqtl-heidi", str(peqtl_heidi),
        "--thread-num", str(thread_num),
        "--out", str(out_file),
    ]

    subprocess.run(cmd_smr, check=True)