import subprocess
from pathlib import Path

import numpy as np
import polars as pl


# generic 3+ trait SNP matcher for multi-trait coloc-style analyses (HyPrColoc,
# MOLOC, ...) - each df in `datasets` needs SNP / A1 / A2 / BETA / SE columns and
# 1 row per SNP already (dedup by whatever the caller's significance criterion is,
# e.g. lowest P). Every non-reference dataset gets its effect allele aligned onto
# `reference`'s A1/A2 (BETA flipped where A1/A2 are swapped, SNP dropped where
# neither allele pairing resolves), then all datasets are reduced to the SNPs
# shared across the lot. Returns {name: df[SNP, BETA, SE]}, row-aligned by SNP.
def extract_common_snps(datasets: dict, reference: str):

    if reference not in datasets:
        raise ValueError(f"Reference dataset '{reference}' not found in datasets: {list(datasets.keys())}")

    ref_alleles = datasets[reference].select([
        "SNP",
        pl.col("A1").alias("ref_a1"),
        pl.col("A2").alias("ref_a2")
    ])

    aligned = {}
    for name, df in datasets.items():
        if name == reference:
            aligned[name] = df.select(["SNP", "BETA", "SE"])
            continue

        df = df.join(ref_alleles, on="SNP", how="inner")
        matched = (pl.col("A1") == pl.col("ref_a1")) & (pl.col("A2") == pl.col("ref_a2"))
        swapped = (pl.col("A1") == pl.col("ref_a2")) & (pl.col("A2") == pl.col("ref_a1"))

        aligned[name] = (
            df
            .filter(matched | swapped)
            .with_columns(
                pl.when(swapped).then(-pl.col("BETA")).otherwise(pl.col("BETA")).alias("BETA")
            )
            .select(["SNP", "BETA", "SE"])
        )

    shared_snps = None

    for df in aligned.values():
        snps = df.select("SNP")
        shared_snps = snps if shared_snps is None else shared_snps.join(snps, on="SNP", how="inner")

    return {name: shared_snps.join(df, on="SNP", how="inner").sort("SNP") for name, df in aligned.items()}


def impute_ld(ref_bfile, snp_1, snp_2):
    cmd = ["plink", "--bfile", ref_bfile, "--ld", snp_1, snp_2]
    return subprocess.run(cmd, check=False, capture_output=True, text=True)


def impute_ld_matrix(snps, out_prefix, ref_bfile):
    # snps -> list
    snp_list_path = f"{out_prefix}_snps.txt"
    with open(snp_list_path, "w") as f:
        f.write("\n".join(snps))

    cmd = [
        "plink",
        "--bfile", ref_bfile,
        "--extract", snp_list_path,
        "--r", "square",
        "--write-snplist",
        "--out", out_prefix,
    ]
    subprocess.run(cmd, check=True)
    ld = pl.read_csv(f"{out_prefix}.ld", separator="\t", has_header=False)
    with open(f"{out_prefix}.snplist") as f:
        snp_order = [line.strip() for line in f]
    return ld, snp_order


def grab_cis_mr_hits(csv_file, cochran_q_thresh: float, causal_thresh: float):
    targets = []
    df = pl.read_csv(csv_file, separator="\t")
    for row in df.iter_rows(named=True):
        protein = row["protein"]
        n_instruments = row["n_instruments"]
        ivw_fdr_q = row["IVW_FDR_q"]
        q_pval = row["Q_pval"]
        wald_fdr_q = row["Wald_FDR_q"]
        if int(n_instruments) == 1:
            if wald_fdr_q < causal_thresh:
                targets.append(protein)
        elif int(n_instruments) == 2:
            if ivw_fdr_q < causal_thresh:
                targets.append(protein)
        else:
            if int(n_instruments) >= 3:
                if ivw_fdr_q < causal_thresh and q_pval > cochran_q_thresh:
                    targets.append(protein)
    return targets


def extract_coloc_or_pwcoco_targets(coloc_csv_file, pwcoco_csv_file, pp4_thresh: float, method: tuple[str, ...] = ("pwcoco", "coloc")):
    targets = set()

    if "coloc" in method:
        df = pl.read_csv(coloc_csv_file, separator="\t")
        for row in df.iter_rows(named=True):
            protein = row["protein_id"]
            pp4 = row["PP.H4.abf"]
            if pp4 >= pp4_thresh:
                targets.add(protein)

    if "pwcoco" in method:
        df = pl.read_csv(pwcoco_csv_file, separator="\t")
        for row in df.iter_rows(named=True):
            protein = row["protein"]
            h4 = row["H4"]
            if h4 >= pp4_thresh:
                targets.add(protein)

    return list(targets)


def extract_smr_hits(bulk_smr_file, sc_smr_file,  p_heidi_thresh, p_smr_thresh, method: tuple[str, ...] = ("bulk", "sc")):
    targets_and_dataset = {} # target = [dataset, ...]

    if "bulk" in method:
        df = pl.read_csv(bulk_smr_file, separator="\t")
        for row in df.iter_rows(named=True):
            protein = row["protein"]
            p_smr = row["q_SMR"]
            p_heidi = row["p_HEIDI"]
            dataset = row["qtl_name"]
            if p_smr <= p_smr_thresh and p_heidi >= p_heidi_thresh:
                targets_and_dataset.setdefault(protein, []).append(dataset)

    if "sc" in method:
        df = pl.read_csv(sc_smr_file, separator="\t")
        for row in df.iter_rows(named=True):
            protein = row["protein"]
            p_smr = row["q_SMR"]
            p_heidi = row["p_HEIDI"]
            cell_type = row["cell_type"]
            if p_smr <= p_smr_thresh and p_heidi >= p_heidi_thresh:
                targets_and_dataset.setdefault(protein, []).append(cell_type)

    return targets_and_dataset


def find_bulk_eqtl(bulk_eqtl_dir, bulk_eqtl_dataset, eqtl_region=None,):
    bulk_eqtl_dir = Path(bulk_eqtl_dir)
    for dataset_dir in bulk_eqtl_dir.iterdir():
        if dataset_dir.is_dir() and bulk_eqtl_dataset in dataset_dir.name:
            for file_path in dataset_dir.glob("*.parquet"):
                if eqtl_region is None:
                    return file_path
                elif eqtl_region in file_path.name:
                    return file_path
    return None


def quick_f_statistic(beta_exposure, se_exposure):
    return (beta_exposure / se_exposure)**2


def lambda_sample_overlap(
        n_overlap,
        n_exposure_total,
        n_outcome_total):
    lambda_res = (
        n_overlap / np.sqrt(n_exposure_total * n_outcome_total)
    )
    return lambda_res


def sample_overlap_relative_bias(lambda_funct, f_statistic):
    raw = lambda_funct / f_statistic
    percent = raw * 100
    return percent