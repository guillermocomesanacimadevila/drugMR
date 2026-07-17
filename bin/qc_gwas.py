#!/usr/bin/env python3

import argparse
import polars as pl 
import os 
import sys
from pathlib import Path
from liftover import ChainFile


# include a liftover to GRCh38 - or keep GRCh38 if == same

# stuff to run before in sys
# module load Apptainer
# apptainer --version
# apptainer build genSEM.sif genSEM.def
# subprocess.run([
#             args.smr_bin,
#             "--eqtl-flist", str(flist_path),
#             "--make-besd",
#             "--out", str(Path(args.outdir) / output_label)
#         ], check=True)

def liftover_df_to_hg38(
    df: pl.DataFrame,
    chr_col: str,
    pos_col: str,
    chain_file: str,
) -> pl.DataFrame:
    
    converter = ChainFile(chain_file, one_based=True)
    chrs = df[chr_col].cast(pl.Utf8).to_list()
    poss = df[pos_col].cast(pl.Int64).to_list()
    lifted_chrs = []
    lifted_poss = []

    for chrom, pos in zip(chrs, poss):
        chrom = str(chrom)
        if chrom.startswith("chr"):
            query_chrom = chrom
        else:
            query_chrom = f"chr{chrom}"

        hits = converter[query_chrom][int(pos)]

        if len(hits) == 0:
            lifted_chrs.append(None)
            lifted_poss.append(None)
        else:
            lifted_chrs.append(hits[0][0].replace("chr", ""))
            lifted_poss.append(hits[0][1])

    df = df.with_columns([
        pl.Series("CHR_GRCh38", lifted_chrs),
        pl.Series("POS_GRCh38", lifted_poss),
    ])

    n_before = df.height
    df = df.drop_nulls(["CHR_GRCh38", "POS_GRCh38"])
    n_after = df.height
    print(f"Removed failed liftover variants: {n_before - n_after}")
    df = df.with_columns([
        pl.col("CHR_GRCh38").alias(chr_col),
        pl.col("POS_GRCh38").cast(pl.Int64).alias(pos_col),
    ])
    df = df.drop(["CHR_GRCh38", "POS_GRCh38"])
    return df

def liftover_df_to_hg19(
    df: pl.DataFrame,
    chr_col: str,
    pos_col: str,
    chain_file: str,
) -> pl.DataFrame:

    converter = ChainFile(chain_file, one_based=True)
    chrs = df[chr_col].cast(pl.Utf8).to_list()
    poss = df[pos_col].cast(pl.Int64).to_list()
    lifted_chrs = []
    lifted_poss = []

    for chrom, pos in zip(chrs, poss):
        chrom = str(chrom)
        if chrom.startswith("chr"):
            query_chrom = chrom
        else:
            query_chrom = f"chr{chrom}"

        hits = converter[query_chrom][int(pos)]

        if len(hits) == 0:
            lifted_chrs.append(None)
            lifted_poss.append(None)
        else:
            lifted_chrs.append(hits[0][0].replace("chr", ""))
            lifted_poss.append(hits[0][1])

    df = df.with_columns([
        pl.Series("CHR_GRCh37", lifted_chrs),
        pl.Series("POS_GRCh37", lifted_poss),
    ])

    n_before = df.height
    df = df.drop_nulls(["CHR_GRCh37", "POS_GRCh37"])
    n_after = df.height
    print(f"Removed failed liftover variants: {n_before - n_after}")
    df = df.with_columns([
        pl.col("CHR_GRCh37").alias(chr_col),
        pl.col("POS_GRCh37").cast(pl.Int64).alias(pos_col),
    ])
    df = df.drop(["CHR_GRCh37", "POS_GRCh37"])
    return df

def perform_qc(
    pheno_id: str,
    sumstats: str,
    out_dir: str,
    maf: float,
    info_threshold: float | None,
    info_col: str | None,
    snp_col: str,
    a1_col: str,
    a2_col: str,
    beta_col: str,
    se_col: str,
    af_col: str,
    p_col: str,
    pos_col: str,
    chr_col: str,
    remove_mhc: bool,
    remove_apoe: bool,
    genome_build: str,
    target_build: str,
    falcon_user: str,
    n_cases: int,
    n_controls: int):

    path = f"/work/{sumstats}"
    bases = ["A", "T", "C", "G"]
    out_dir = Path(out_dir)
    qc_dir = out_dir / "QC" / pheno_id
    qc_dir.mkdir(parents=True, exist_ok=True)
    df = pl.read_csv(path, separator="\t", comment_prefix="#", schema_overrides={chr_col: pl.Utf8})
    print(df.shape)

    # Remove empty stuff
    n_before = df.height
    df = df.filter(~pl.all_horizontal(pl.all().is_null()))
    n_after = df.height

    print(f"Removed empty rows: {n_before - n_after}")
    print(f"Rows after empty-row removal: {n_after}")

    # remove INDELs
    a1 = pl.col(a1_col).str.to_uppercase()
    a2 = pl.col(a2_col).str.to_uppercase()
    ok_len = (a1.str.len_chars() == 1) & (a2.str.len_chars() == 1)
    ok_bases = a1.is_in(bases) & a2.is_in(bases)
    no_gap = ~a1.str.contains("-") & ~a2.str.contains("-")
    df = df.filter(ok_len & ok_bases & no_gap)

    # remove palindromes
    pal = (
        ((a1 == "A") & (a2 == "T")) | ((a1 == "T") & (a2 == "A")) |
        ((a1 == "C") & (a2 == "G")) | ((a1 == "G") & (a2 == "C"))
    )
    df = df.filter(~pal)

    # remove rare variants anf info SCORE 
    if af_col not in [None, "", "None", "NA", "nan"]:
        df = df.with_columns(pl.col(af_col).cast(pl.Float64, strict=False))
        df = df.filter(
            pl.col(af_col).is_not_null()
            & (pl.col(af_col) >= maf)
            & (pl.col(af_col) <= 1 - maf)
        )
        print(f"[TRACKING] Applied MAF filter using {af_col}: {maf}")
    else:
        print("[TRACKING] No allele-frequency column provided, skipping MAF filtering...")
    
    if info_col is not None and info_threshold is not None:
        df = df.filter(pl.col(info_col) >= info_threshold)
    
    # remove duplicate SNPs 
    df = df.unique(subset=[snp_col], keep="first")

    # add N col 
    n_cases = int(n_cases)
    n_controls = int(n_controls)
    neff = (4 * n_cases * n_controls) / (n_cases + n_controls)
    df = df.with_columns(pl.lit(neff).alias("N"))

    # file is in /shared/home1/{falcon_user}/genSEM/ref/hg19_38
    if genome_build == target_build:
        print(f"Input GWAS already in {target_build}")
    elif genome_build == "GRCh38" and target_build == "GRCh37":
        print("Lifting GRCh38 -> GRCh37/hg19")
        chain_file = "./dat/ref/liftover/hg38ToHg19.over.chain"

        if not Path(chain_file).exists():
            raise FileNotFoundError(f"Missing chain file: {chain_file}")

        n_before = df.height

        df = liftover_df_to_hg19(
            df=df,
            chr_col=chr_col,
            pos_col=pos_col,
            chain_file=chain_file,
        )

        n_after = df.height

        print(f"Rows before liftover: {n_before}")
        print(f"Rows after liftover: {n_after}")

    elif genome_build == "GRCh37" and target_build == "GRCh38":
        print("Lifting GRCh37/hg19 -> GRCh38/hg38")
        chain_file = "./dat/ref/liftover/hg19ToHg38.over.chain"

        if not Path(chain_file).exists():
            raise FileNotFoundError(f"Missing chain file: {chain_file}")

        n_before = df.height

        df = liftover_df_to_hg38(
            df=df,
            chr_col=chr_col,
            pos_col=pos_col,
            chain_file=chain_file,
        )

        n_after = df.height

        print(f"Rows before liftover: {n_before}")
        print(f"Rows after liftover: {n_after}")

    else:
        raise ValueError(f"Unsupported genome_build/target_build: {genome_build} -> {target_build}")
    
    # fix chromosome ambiguity
    df = df.with_columns(pl.col(chr_col).cast(pl.Utf8).str.replace("^chr", "").alias(chr_col))
    df = df.filter(pl.col(chr_col).is_in([str(i) for i in range(1, 23)]))

    if target_build == "GRCh37":
        mhc_start = 25000000
        mhc_end = 34000000
        apoe_start = 44000000
        apoe_end = 46500000
    elif target_build == "GRCh38":
        mhc_start = 24999772
        mhc_end = 34032223
        apoe_start = 43495848
        apoe_end = 45996742
    else:
        raise ValueError(f"Unsupported target_build: {target_build}")

    # remove MHC
    if remove_mhc:
        n_before = df.height
        df = df.filter(
            ~(
                (pl.col(chr_col) == "6") &
                (pl.col(pos_col) >= mhc_start) &
                (pl.col(pos_col) <= mhc_end)
            )
        )
        n_after = df.height
        print(f"Removed MHC variants: {n_before - n_after}")

    # remove APOE
    if remove_apoe:
        n_before = df.height
        df = df.filter(
            ~(
                (pl.col(chr_col) == "19") &
                (pl.col(pos_col) >= apoe_start) &
                (pl.col(pos_col) <= apoe_end)
            )
        )
        n_after = df.height
        print(f"Removed APOE SNPs: {n_before - n_after}")

    # rename cols to LDSC / genSEM format
    if info_col is not None:
        ldsc = df.select(
            pl.col(snp_col).alias("SNP"),
            pl.col(a1_col).alias("A1"),
            pl.col(a2_col).alias("A2"),
            pl.col(af_col).alias("FRQ"),
            pl.col("N"),
            pl.col(beta_col).alias("BETA"),
            pl.col(se_col).alias("SE"),
            pl.col(p_col).alias("P"),
            pl.col(chr_col).alias("CHR"),
            pl.col(pos_col).alias("BP"),
            pl.col(info_col).alias("INFO"),
        )
    else:
        ldsc = df.select(
            pl.col(snp_col).alias("SNP"),
            pl.col(a1_col).alias("A1"),
            pl.col(a2_col).alias("A2"),
            pl.col(af_col).alias("FRQ"),
            pl.col("N"),
            pl.col(beta_col).alias("BETA"),
            pl.col(se_col).alias("SE"),
            pl.col(p_col).alias("P"),
            pl.col(chr_col).alias("CHR"),
            pl.col(pos_col).alias("BP"),
        )

    out_path = qc_dir / f"{pheno_id}.tsv"
    ldsc.write_csv(
        out_path,
        separator="\t"
    )

    print(f"Saved QCed GWAS to: {out_path}")

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--pheno-id", required=True)
    parser.add_argument("--sumstats", required=True)
    parser.add_argument("--out-dir", required=True)
    parser.add_argument("--maf", type=float, default=0.01)
    parser.add_argument("--info-threshold", type=float, default=None)
    parser.add_argument("--info-col", default=None)
    parser.add_argument("--snp-col", required=True)
    parser.add_argument("--a1-col", required=True)
    parser.add_argument("--a2-col", required=True)
    parser.add_argument("--beta-col", required=True)
    parser.add_argument("--se-col", required=True)
    parser.add_argument("--p-col", required=True)
    parser.add_argument("--pos-col", required=True)
    parser.add_argument("--chr-col", required=True)
    parser.add_argument("--falcon-user", required=True)
    parser.add_argument("--remove_mhc", action="store_true")
    parser.add_argument("--genome_build", required=True)
    parser.add_argument("--target_build", required=True)
    parser.add_argument("--n_cases", required=True, type=int)
    parser.add_argument("--n_controls", required=True, type=int)
    parser.add_argument("--af_col", required=True)
    parser.add_argument("--remove_apoe", action="store_true")
    args = parser.parse_args()
    perform_qc(
        pheno_id=args.pheno_id,
        sumstats=args.sumstats,
        out_dir=args.out_dir,
        maf=args.maf,
        info_threshold=args.info_threshold,
        info_col=args.info_col,
        snp_col=args.snp_col,
        a1_col=args.a1_col,
        a2_col=args.a2_col,
        beta_col=args.beta_col,
        se_col=args.se_col,
        p_col=args.p_col,
        af_col=args.af_col,
        pos_col=args.pos_col,
        chr_col=args.chr_col,
        falcon_user=args.falcon_user,
        remove_mhc=args.remove_mhc,
        genome_build=args.genome_build,
        target_build=args.target_build,
        n_cases=args.n_cases,
        n_controls=args.n_controls,
        remove_apoe=args.remove_apoe
    )

if __name__ == "__main__":
    main()