import subprocess
from pathlib import Path

import liftover
import numpy as np
import polars as pl

# generic 3+ trait SNP matcher for multi-trait coloc-style analyses (HyPrColoc,
# MOLOC, ...) - each df in `datasets` needs SNP / A1 / A2 / BETA / SE columns and
# 1 row per SNP already (dedup by whatever the caller's significance criterion is,
# e.g. lowest P). Every non-reference dataset gets its effect allele aligned onto
# `reference`'s A1/A2 (BETA flipped where A1/A2 are swapped, SNP dropped where
# neither allele pairing resolves), then all datasets are reduced to the SNPs
# shared across the lot. Returns {name: df[SNP, BETA, SE]}, row-aligned by SNP.

def quick_qc(sumstats: pl.DataFrame, a1_col: str, a2_col: str):
    bases = ["A", "C", "T", "G"]
    a1 = pl.col(a1_col).str.to_uppercase()
    a2 = pl.col(a2_col).str.to_uppercase()
    ok_len = (a1.str.len_chars() == 1) & (a2.str.len_chars() == 1)
    ok_bases = a1.is_in(bases) & a2.is_in(bases)
    df = (sumstats.drop_nulls().filter(ok_len & ok_bases))
    return df


def check_n_in_qtl(df: pl.DataFrame, n_total) -> pl.DataFrame:
    if "N" not in df.columns:
        df = df.with_columns(pl.lit(n_total).alias("N"))
    return df


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


# canonical cis-MR pass/fail rule - shared by standard COLOC's own protein
# selection (bin/coloc_targets.py's pairwise_coloc()) and PWCoCo's (below), so
# a target that fails cis-MR for one method can't still silently reach the
# other via 2 independently-drifted copies of this filter.
def select_cis_mr_passing_proteins(
    df: pl.DataFrame,
    wald_fdr_q: float = 0.05,
    ivw_fdr_q: float = 0.05,
    cochran_q_pval: float = 0.05,
    egger_intercept_pval_min: float = 0,
    min_instruments_for_ivw: int = 3,
) -> list[str]:
    return (
        df
        .filter(
            (
                (pl.col("n_instruments") >= min_instruments_for_ivw) &
                (pl.col("IVW_FDR_q") < ivw_fdr_q) &
                (pl.col("egger_intercept_pval") > egger_intercept_pval_min) &
                (pl.col("Q_pval") > cochran_q_pval)
            )
            |
            (
                (pl.col("n_instruments") == 1) &
                (pl.col("Wald_FDR_q") < wald_fdr_q)
            )
        )
        .select("protein")
        .unique()
        .get_column("protein")
        .to_list()
    )


def grab_cis_mr_hits(
    csv_file,
    wald_fdr_q: float = 0.05,
    ivw_fdr_q: float = 0.05,
    cochran_q_pval: float = 0.05,
    egger_intercept_pval_min: float = 0,
    min_instruments_for_ivw: int = 3,
):
    df = pl.read_csv(
        csv_file,
        separator="\t",
        infer_schema_length=None,
        schema_overrides={
            "n_instruments": pl.Int64,
            "IVW_FDR_q": pl.Float64,
            "egger_intercept_pval": pl.Float64,
            "Q_pval": pl.Float64,
            "Wald_FDR_q": pl.Float64,
        },
    )
    return select_cis_mr_passing_proteins(
        df,
        wald_fdr_q=wald_fdr_q,
        ivw_fdr_q=ivw_fdr_q,
        cochran_q_pval=cochran_q_pval,
        egger_intercept_pval_min=egger_intercept_pval_min,
        min_instruments_for_ivw=min_instruments_for_ivw,
    )


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


def extract_gene_coordinates(
        gene_id: str,
        ref: pl.DataFrame,
        genome_build: str = "hg38",
        converter=None,
    ):

    # current gencode == v50 so...

    start = 0
    end = 0
    chr = 0
    orientation = ""

    ensembl_id = None

    if genome_build == "hg38":
        for row in ref.iter_rows(named=True):
            if gene_id == row["Symbol"]:
                start = int(row["Begin"])
                end = int(row["End"])
                chr = int(row["Chromosome"])
                orientation = str(row["Orientation"])
                ensembl_id = row["Ensembl_ID"]
                break
        else:
            for row in ref.iter_rows(named=True):
                if gene_id in (row["Synonyms"] or "").split(","):
                    start = int(row["Begin"])
                    end = int(row["End"])
                    chr = int(row["Chromosome"])
                    orientation = str(row["Orientation"])
                    ensembl_id = row["Ensembl_ID"]
                    break
            else:
                raise ValueError(f"Gene '{gene_id}' not found in reference")

        accum_dict = {
            "ORIENTATION": orientation,
            "CHR": chr,
            "START": start,
            "END": end,
            "ENSEMBL_ID": ensembl_id,
        }
        df = pl.DataFrame(accum_dict)
        return df

    if genome_build == "hg19":
        if converter is None:
            converter = liftover.get_lifter("hg38", "hg19", one_based=True)
        for row in ref.iter_rows(named=True):
            if gene_id == row["Symbol"]:
                chr = str(row["Chromosome"])
                start = converter[chr][int(row["Begin"])][0][1]
                end = converter[chr][int(row["End"])][0][1]
                orientation = str(row["Orientation"])
                ensembl_id = row["Ensembl_ID"]
                break
        else:
            for row in ref.iter_rows(named=True):
                if gene_id in (row["Synonyms"] or "").split(","):
                    chr = str(row["Chromosome"])
                    start = converter[chr][int(row["Begin"])][0][1]
                    end = converter[chr][int(row["End"])][0][1]
                    orientation = str(row["Orientation"])
                    ensembl_id = row["Ensembl_ID"]
                    break
            else:
                raise ValueError(f"Gene '{gene_id}' not found in reference")

        accum_dict = {
            "ORIENTATION": orientation,
            "CHR": chr,
            "START": start,
            "END": end,
            "ENSEMBL_ID": ensembl_id,
        }
        df = pl.DataFrame(accum_dict)
        return df



#################################
#################################
#################################
# THIS SHIT NEEDS FIXING TO ENSURE 
# THAT SMR READY FILES FOR QTL DATASET X 
# ARE ACCOUNTED FOR!!!!



def needs_chr_split_etl(df: pl.DataFrame, chr_col: str) -> bool:
    return df[chr_col].n_unique() > 1


def detect_qtl_split(manifest_path: Path) -> dict:
    manifest_path = Path(manifest_path)
    literal_parts = []
    for part in manifest_path.parts:
        if any(ch in part for ch in "*?["):
            break
        literal_parts.append(part)
    search_dir = Path(*literal_parts) if literal_parts else Path(".")

    # a fully-literal path (no wildcard) is either an existing raw-data file, or
    # a target parquet that hasn't been built yet - either way, if it isn't
    # itself a directory, its parent is what we actually want to search
    if not search_dir.is_dir():
        search_dir = search_dir.parent

    if not search_dir.is_dir():
        return {}

    prefixes = {}
    for besd_file in search_dir.rglob("*.besd"):
        prefix = Path(str(besd_file)[: -len(".besd")])
        if Path(f"{prefix}.esi").exists() and Path(f"{prefix}.epi").exists():
            label = str(prefix.relative_to(search_dir))
            prefixes[label] = prefix

    return prefixes
