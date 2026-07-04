#!/usr/bin/env python3
import argparse
import polars as pl 
import subprocess
from pathlib import Path
import os

# TO DO'S
# ** INTEGRATE assort_networkMR.py onto drugMR/hpc.py and drugMR/local.py and test it


def qc_cmd(row, maf: float, remove_mhc: bool, remove_apoe: bool, overwrite: bool):
    # out_dir == already defined
    # already == polars dataframe
    # predefined colnames
    # pheno_id_col
    # snp_col
    # a1_col
    # a2_col
    # chr_col
    # pos_col
    # beta_col
    # pval_col
    # se_col
    # info_col
    # maf_col
    # genome_build
    # target_build
    # n_total - although the qc_gwas is cc focused so we need to do it like a way that when it calcs neff - neff == n_total

    out_dir = "./results/QC/mediators"
    os.makedirs(out_dir, exist_ok=True)
    mediator_id = row["pheno_id"]
    out_path = Path(out_dir) / f"{mediator_id}.tsv"

    if out_path.exists() and not overwrite:
        print(f"[TRACKING] {mediator_id} already QCed, skipping!")
        return

    cmd = [
        "python", "./bin/qc_gwas.py",
        "--pheno-id", str(row["pheno_id"]),
        "--sumstats", str(row["sumstats"]),
        "--out-dir", out_dir,
        "--maf", str(maf),
        "--snp-col", str(row["snp_col"]),
        "--a1-col", str(row["a1_col"]),
        "--a2-col", str(row["a2_col"]),
        "--beta-col", str(row["beta_col"]),
        "--se-col", str(row["se_col"]),
        "--p-col", str(row["p_col"]),
        "--pos-col", str(row["pos_col"]),
        "--chr-col", str(row["chr_col"]),
        "--af_col", str(row["af_col"]),
        "--genome_build", str(row["genome_build"]),
        "--target_build", str(row["target_build"]),
        "--falcon-user", str(row.get("falcon_user", "")),
        "--n_cases", str(row["n_total"]),
        "--n_controls", str(row["n_total"]),
    ]

    if row.get("info_col") is not None:
        cmd.extend(["--info-col", str(row["info_col"])])

    if row.get("info_threshold") is not None:
        cmd.extend(["--info-threshold", str(row["info_threshold"])])

    if remove_mhc:
        cmd.append("--remove_mhc")

    if remove_apoe:
        cmd.append("--remove_apoe")

    subprocess.run(cmd, check=True)
    print(f"[TRACKING] Saved mediator GWAS to: {out_path}")

def organise_mediators():
    qc_dir = Path("./results/QC/mediators")
    nested_qc_dir = qc_dir / "QC"

    if not nested_qc_dir.exists():
        return

    for pheno_dir in nested_qc_dir.iterdir():
        if not pheno_dir.is_dir():
            continue

        mediator_id = pheno_dir.name
        src = pheno_dir / f"{mediator_id}.tsv"
        dst = qc_dir / f"{mediator_id}.tsv"

        if src.exists():
            src.rename(dst)
            print(f"[TRACKING] Moved mediator GWAS to: {dst}")

        pheno_dir.rmdir()
    nested_qc_dir.rmdir()

def preprocess_mediators(mediators: bool, mediator_manifest: str, maf: float, remove_mhc: bool, remove_apoe: bool, overwrite: bool):
    if mediators:
        df = pl.read_csv(mediator_manifest)
        print("[TRACKING] Mediator manifest CSV file loaded!")
        for row in df.iter_rows(named=True):
            qc_cmd(
                row=row,
                maf=maf,
                remove_mhc=remove_mhc,
                remove_apoe=remove_apoe,
                overwrite=overwrite
            )
            organise_mediators()
    else:
        print("[TRACKING] No mediators specificed, running drugMR without them then!")

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--mediators", action="store_true")
    parser.add_argument("--mediator-manifest", default="")
    parser.add_argument("--maf", type=float, default=0.01)
    parser.add_argument("--remove_mhc", action="store_true")
    parser.add_argument("--remove_apoe", action="store_true")
    parser.add_argument("--overwrite", action="store_true")
    args = parser.parse_args()
    preprocess_mediators(
        mediators=args.mediators,
        mediator_manifest=args.mediator_manifest,
        maf=args.maf,
        remove_mhc=args.remove_mhc,
        remove_apoe=args.remove_apoe,
        overwrite=args.overwrite
    )

if __name__ == "__main__":
    main()