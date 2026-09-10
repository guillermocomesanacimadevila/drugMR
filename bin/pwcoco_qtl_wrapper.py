import argparse
from pathlib import Path

import polars as pl

from drugmr.paths import (
    DEFAULT_QTL_MANIFEST_PATH,
    pwcoco_eqtl_gwas_out,
    pwcoco_eqtl_pqtl_out,
    pwcoco_out,
    pwcoco_qtl_raw_prefix,
    pwcoco_qtl_shared_out,
    smr_final_targets_out,
)
from drugmr.pwcoco import PWCoCo
from drugmr.smr import SMRUtils

_smr = SMRUtils(manifest_path=DEFAULT_QTL_MANIFEST_PATH)


def resolve_maf_col(df):
    return "FRQ" if "FRQ" in df.columns else "MAF"


# 1 PWCoCo output df -> {protein: {snp: h4}}, keeping only rows that clear pp4_thresh -
# SNP1/SNP2 is literally "unconditioned" for PWCoCo's own unconditioned row, and a
# conditioned row's SNP carries a trailing "*" (PWCoCo's conditioning-SNP marker) -
# stripped here so the same variant matches across combos regardless of which row
# it came from
def snp_h4_map(df, pp4_thresh):
    m = {}
    for row in df.iter_rows(named=True):
        if row["H4"] < pp4_thresh:
            continue
        for snp in (row["SNP1"], row["SNP2"]):
            if snp and snp != "unconditioned":
                snp = snp.rstrip("*")
                protein_map = m.setdefault(row["protein"], {})
                protein_map[snp] = max(protein_map.get(snp, 0), row["H4"])
    return m


def pwcoco_qtl_wrapper(
        pqtl_dataset: str,
        pheno_id: str,
        ref_bfile: str,
        n_cases: int,
        n_controls: int,
        out_dir: str = "results",
        pp4_thresh: float = 0.7,
        cis_regions_dir: str | None = None
):

    """
    PWCoCo for each pair within eQTL-informed hits (pQTL-eQTL / GWAS-eQTL)
    - For any target situated in SMR (survies past SMR)
    - Grab pertaining GWAS and pQTL cis-reg - check dataset in which == hit in SMR ->
    - Extract cis-region for eQTL and overlap with the other two - then os.unlink() **
    - Run PWCoCo between eQTL-pQTL and GWAS-eQTL
    - For out_cols -> concordant - check original pQTL-GWAS to check and add that bool col to final df
    # ruff check --select I --fix .
    """

    pwcoco = PWCoCo()
    smr_targets = pl.read_csv(smr_final_targets_out(pqtl_dataset=pqtl_dataset, pheno_id=pheno_id, out_dir=out_dir), separator="\t") # data_type
    bulk = smr_targets.filter(pl.col("data_type") == "bulk")
    sc = smr_targets.filter(pl.col("data_type") == "single_cell")

    qtl_pqtl_rows = []
    qtl_gwas_rows = []

    # bulk hits first
    for row in bulk.iter_rows(named=True):
        p = row["protein"]
        dataset = row["qtl_dataset"]
        qtl_type = row.get("qtl_type", "eqtl")
        cell_type = row["cell_type"]
        probe = row["probe_id"]

        # 1 bad target must not lose every other already-computed target -
        # same reasoning as cis_mr.R's own per-protein tryCatch
        try:
            # cis regions
            dir = Path(cis_regions_dir) / p if cis_regions_dir else Path(f"./dat/cis_regions/{pqtl_dataset}/{p}")
            pqtl = dir / "pqtl.parquet"
            gwas = dir / "gwas.parquet"

            pqtl_df = pl.read_parquet(pqtl)

            base_gene = probe.split(".")[0]
            qtl_df = _smr.load_qtl_rows("bulk", dataset, cell_type, base_gene)

            if qtl_df is None or qtl_df.height == 0:
                continue

            gwas_df = pl.read_parquet(gwas)

            pqtl_h = pwcoco.harmonise_sumstats(pqtl_df, "SNP", "A1", "A2", resolve_maf_col(pqtl_df), "BETA", "SE", "P", "N")
            gwas_h = pwcoco.harmonise_sumstats(gwas_df, "SNP", "A1", "A2", resolve_maf_col(gwas_df), "BETA", "SE", "P", "N")
            qtl_h = pwcoco.harmonise_sumstats(qtl_df, "SNP", "A1", "A2", "FRQ", "BETA", "SE", "P", "N")

            qtl_n = int(qtl_h["n"][0])
            pqtl_n = int(pqtl_h["n"][0])
            qtl_source = f"{dataset}_{cell_type}"

            # eQTL - pQTL: both quantitative traits, n2_case=0
            out_ep = pwcoco_qtl_raw_prefix("eqtl_pqtl", pqtl_dataset, p, qtl_source, out_dir)
            out_ep.parent.mkdir(parents=True, exist_ok=True)
            pwcoco.pwcoco(
                ref_bfile=ref_bfile, sumstats_1=qtl_h, sumstats_2=pqtl_h,
                n_1=qtl_n, n_2=pqtl_n, n2_case=0, out_dir=str(out_ep), threads=8
            )
            if Path(f"{out_ep}.coloc").exists():
                qtl_pqtl_rows.append(
                    pl.read_csv(f"{out_ep}.coloc", separator="\t").with_columns(pl.lit(p).alias("protein"), pl.lit(dataset).alias("qtl_dataset"), pl.lit(cell_type).alias("cell_type"), pl.lit(qtl_type).alias("qtl_type"))
                )

            # eQTL - GWAS: GWAS is case-control
            out_eg = pwcoco_qtl_raw_prefix("eqtl_gwas", pqtl_dataset, p, qtl_source, out_dir)
            out_eg.parent.mkdir(parents=True, exist_ok=True)
            pwcoco.pwcoco(
                ref_bfile=ref_bfile, sumstats_1=qtl_h, sumstats_2=gwas_h,
                n_1=qtl_n, n_2=n_cases + n_controls, n2_case=n_cases, out_dir=str(out_eg), threads=8
            )
            if Path(f"{out_eg}.coloc").exists():
                qtl_gwas_rows.append(
                    pl.read_csv(f"{out_eg}.coloc", separator="\t").with_columns(pl.lit(p).alias("protein"), pl.lit(dataset).alias("qtl_dataset"), pl.lit(cell_type).alias("cell_type"), pl.lit(pheno_id).alias("outcome_trait"), pl.lit(qtl_type).alias("qtl_type"))
                )
        except Exception as error:
            print(f"[CONCERN] PWCoCo-QTL (bulk) failed for {p} x {cell_type} - continuing without it: {error}")

    # single-cell hits
    for row in sc.iter_rows(named=True):
        p = row["protein"]
        dataset = row["qtl_dataset"]
        qtl_type = row.get("qtl_type", "eqtl")
        cell_type = row["cell_type"]
        probe = row["probe_id"]

        # 1 bad target must not lose every other already-computed target -
        # same reasoning as cis_mr.R's own per-protein tryCatch
        try:
            dir = Path(cis_regions_dir) / p if cis_regions_dir else Path(f"./dat/cis_regions/{pqtl_dataset}/{p}")
            pqtl = dir / "pqtl.parquet"
            gwas = dir / "gwas.parquet"

            pqtl_df = pl.read_parquet(pqtl)
            gwas_df = pl.read_parquet(gwas)

            base_gene = probe.split(".")[0]
            qtl_df = _smr.load_qtl_rows("single_cell", dataset, cell_type, base_gene)

            if qtl_df is None or qtl_df.height == 0:
                continue

            pqtl_h = pwcoco.harmonise_sumstats(pqtl_df, "SNP", "A1", "A2", resolve_maf_col(pqtl_df), "BETA", "SE", "P", "N")
            gwas_h = pwcoco.harmonise_sumstats(gwas_df, "SNP", "A1", "A2", resolve_maf_col(gwas_df), "BETA", "SE", "P", "N")
            qtl_h = pwcoco.harmonise_sumstats(qtl_df, "SNP", "A1", "A2", "FRQ", "BETA", "SE", "P", "N")

            qtl_n = int(qtl_h["n"][0])
            pqtl_n = int(pqtl_h["n"][0])
            qtl_source = f"{dataset}_{cell_type}"

            out_ep = pwcoco_qtl_raw_prefix("eqtl_pqtl", pqtl_dataset, p, qtl_source, out_dir)
            out_ep.parent.mkdir(parents=True, exist_ok=True)
            pwcoco.pwcoco(
                ref_bfile=ref_bfile, sumstats_1=qtl_h, sumstats_2=pqtl_h,
                n_1=qtl_n, n_2=pqtl_n, n2_case=0, out_dir=str(out_ep), threads=8
            )
            if Path(f"{out_ep}.coloc").exists():
                qtl_pqtl_rows.append(
                    pl.read_csv(f"{out_ep}.coloc", separator="\t").with_columns(pl.lit(p).alias("protein"), pl.lit(dataset).alias("qtl_dataset"), pl.lit(cell_type).alias("cell_type"), pl.lit(qtl_type).alias("qtl_type"))
                )

            out_eg = pwcoco_qtl_raw_prefix("eqtl_gwas", pqtl_dataset, p, qtl_source, out_dir)
            out_eg.parent.mkdir(parents=True, exist_ok=True)
            pwcoco.pwcoco(
                ref_bfile=ref_bfile, sumstats_1=qtl_h, sumstats_2=gwas_h,
                n_1=qtl_n, n_2=n_cases + n_controls, n2_case=n_cases, out_dir=str(out_eg), threads=8
            )
            if Path(f"{out_eg}.coloc").exists():
                qtl_gwas_rows.append(
                    pl.read_csv(f"{out_eg}.coloc", separator="\t").with_columns(pl.lit(p).alias("protein"), pl.lit(dataset).alias("qtl_dataset"), pl.lit(cell_type).alias("cell_type"), pl.lit(pheno_id).alias("outcome_trait"), pl.lit(qtl_type).alias("qtl_type"))
                )
        except Exception as error:
            print(f"[CONCERN] PWCoCo-QTL (single-cell) failed for {p} x {cell_type} - continuing without it: {error}")

    qtl_pqtl_df = pl.concat(qtl_pqtl_rows, how="diagonal_relaxed") if qtl_pqtl_rows else pl.DataFrame()
    qtl_gwas_df = pl.concat(qtl_gwas_rows, how="diagonal_relaxed") if qtl_gwas_rows else pl.DataFrame()

    ep_out = pwcoco_eqtl_pqtl_out(pqtl_dataset, pheno_id, out_dir)
    ep_out.parent.mkdir(parents=True, exist_ok=True)
    if qtl_pqtl_df.height > 0:
        qtl_pqtl_df = qtl_pqtl_df.drop([c for c in ("Dataset1", "Dataset2") if c in qtl_pqtl_df.columns])
        qtl_pqtl_df.write_csv(ep_out, separator="\t")

    eg_out = pwcoco_eqtl_gwas_out(pqtl_dataset, pheno_id, out_dir)
    eg_out.parent.mkdir(parents=True, exist_ok=True)
    if qtl_gwas_df.height > 0:
        qtl_gwas_df = qtl_gwas_df.drop([c for c in ("Dataset1", "Dataset2") if c in qtl_gwas_df.columns])
        qtl_gwas_df.write_csv(eg_out, separator="\t")

    # triangulation across all 3 combos (pQTL-GWAS from pwcoco_out(), eQTL-pQTL,
    # eQTL-GWAS) - a target only counts when the SAME SNP clears pp4_thresh in ALL
    # 3, not just some of them. This is the eQTL-level equivalent of HyPrColoc: it
    # answers the same "do pQTL/eQTL/GWAS share 1 causal variant" question, but via
    # conditioning rather than HyPrColoc's single-causal-variant cluster assumption -
    # co-equal to HyPrColoc, not a downstream refinement of it.
    pg_file = pwcoco_out(pqtl_dataset, pheno_id, out_dir)
    pqtl_gwas_df = pl.read_csv(pg_file, separator="\t") if Path(pg_file).exists() else pl.DataFrame()

    pg_map = snp_h4_map(pqtl_gwas_df, pp4_thresh) if pqtl_gwas_df.height > 0 else {}
    ep_map = snp_h4_map(qtl_pqtl_df, pp4_thresh) if qtl_pqtl_df.height > 0 else {}
    eg_map = snp_h4_map(qtl_gwas_df, pp4_thresh) if qtl_gwas_df.height > 0 else {}

    shared_rows = []
    for protein in set(pg_map) & set(ep_map) & set(eg_map):
        shared_snps = set(pg_map[protein]) & set(ep_map[protein]) & set(eg_map[protein])
        for snp in shared_snps:
            shared_rows.append({
                "protein": protein,
                "snp": snp,
                "pqtl_gwas_h4": pg_map[protein][snp],
                "qtl_pqtl_h4": ep_map[protein][snp],
                "qtl_gwas_h4": eg_map[protein][snp],
            })

    if shared_rows:
        shared_file = pwcoco_qtl_shared_out(pqtl_dataset, pheno_id, out_dir)
        shared_file.parent.mkdir(parents=True, exist_ok=True)
        pl.DataFrame(shared_rows).write_csv(shared_file, separator="\t")


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--pqtl_dataset", required=True)
    p.add_argument("--pheno_id", required=True)
    p.add_argument("--ref_bfile", required=True)
    p.add_argument("--n_cases", required=True, type=int)
    p.add_argument("--n_controls", required=True, type=int)
    p.add_argument("--local_results_dir", default="results")
    p.add_argument("--pp4_threshold", type=float, default=0.7)
    p.add_argument("--manifest_path", default=DEFAULT_QTL_MANIFEST_PATH)
    p.add_argument("--repo_root", default=None)
    p.add_argument("--cis_regions_dir", default=None)
    args = p.parse_args()

    _smr.manifest_path = args.manifest_path
    if args.repo_root:
        _smr.base_dir = args.repo_root

    pwcoco_qtl_wrapper(
        pqtl_dataset=args.pqtl_dataset,
        pheno_id=args.pheno_id,
        ref_bfile=args.ref_bfile,
        n_cases=args.n_cases,
        n_controls=args.n_controls,
        out_dir=args.local_results_dir,
        pp4_thresh=args.pp4_threshold,
        cis_regions_dir=args.cis_regions_dir,
    )


if __name__ == "__main__":
    main()