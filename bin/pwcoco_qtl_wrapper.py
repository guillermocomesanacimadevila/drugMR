import argparse
from pathlib import Path

import polars as pl

from drugmr.paths import (
    pwcoco_eqtl_gwas_out,
    pwcoco_eqtl_pqtl_out,
    pwcoco_out,
    pwcoco_qtl_raw_prefix,
    pwcoco_qtl_shared_out,
    smr_final_targets_out,
)
from drugmr.pwcoco import PWCoCo
from drugmr.utils import find_bulk_eqtl


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
        pp4_thresh: float = 0.7
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

    eqtl_pqtl_rows = []
    eqtl_gwas_rows = []

    # bulk hits first
    for row in bulk.iter_rows(named=True):
        p = row["protein"]
        dataset = row["eqtl_dataset"]
        qtl_name = row["qtl_name"]
        cell_type = row["cell_type"]
        probe = row["probeID"]

        # cis regions
        dir = f"./dat/cis_regions/{pqtl_dataset}/{p}"
        dir = Path(dir)
        pqtl = dir / "pqtl.parquet"
        gwas = dir / "gwas.parquet"

        pqtl_df = pl.read_parquet(pqtl)
        chr = pqtl_df["CHR"].to_list()[0]
        pos = min(pqtl_df["BP"].to_list())
        start = min(pqtl_df["BP"].to_list())
        end = max(pqtl_df["BP"].to_list())

        # grab bulk eqtl dataset - full parquet - auto-detect flat (MetaBrain-style,
        # parquet sits directly in the dataset dir) vs nested-by-region (GTEx-style,
        # parquet sits one level down inside a tissue dir) instead of hardcoding either
        # dataset name
        eqtl_file = find_bulk_eqtl("./dat/bulk-eQTL", dataset)
        if eqtl_file is None:
            dataset_dir = Path(f"./dat/bulk-eQTL/{dataset}")
            for sub in dataset_dir.iterdir():
                if sub.is_dir() and sub.name in cell_type:
                    hits = list(sub.glob("*.parquet"))
                    if hits:
                        eqtl_file = hits[0]
                    break

        base_gene = probe.split(".")[0] # eqtl_file -> path to .parquet for that row

        if eqtl_file is None:
            continue

        eqtl_df = (
            pl.scan_parquet(eqtl_file)
            .filter(pl.col("Probe").str.split(".").list.first() == base_gene)
            .collect()
            .sort("p")
            .unique(subset="SNP", keep="first")
        )

        if eqtl_df.height == 0:
            continue

        gwas_df = pl.read_parquet(gwas)

        pqtl_h = pwcoco.harmonise_sumstats(pqtl_df, "SNP", "A1", "A2", resolve_maf_col(pqtl_df), "BETA", "SE", "P", "N")
        gwas_h = pwcoco.harmonise_sumstats(gwas_df, "SNP", "A1", "A2", resolve_maf_col(gwas_df), "BETA", "SE", "P", "N")
        eqtl_h = pwcoco.harmonise_sumstats(eqtl_df, "SNP", "A1", "A2", "Freq", "b", "SE", "p", "N")

        eqtl_n = int(eqtl_h["n"][0])
        pqtl_n = int(pqtl_h["n"][0])
        eqtl_source = f"{dataset}_{cell_type}"

        # eQTL - pQTL: both quantitative traits, n2_case=0
        out_ep = pwcoco_qtl_raw_prefix("eqtl_pqtl", pqtl_dataset, p, eqtl_source, out_dir)
        out_ep.parent.mkdir(parents=True, exist_ok=True)
        pwcoco.pwcoco(
            ref_bfile=ref_bfile, sumstats_1=eqtl_h, sumstats_2=pqtl_h,
            n_1=eqtl_n, n_2=pqtl_n, n2_case=0, out_dir=str(out_ep), threads=8
        )
        if Path(f"{out_ep}.coloc").exists():
            eqtl_pqtl_rows.append(
                pl.read_csv(f"{out_ep}.coloc", separator="\t").with_columns(pl.lit(p).alias("protein"), pl.lit(dataset).alias("eqtl_dataset"), pl.lit(cell_type).alias("cell_type"))
            )

        # eQTL - GWAS: GWAS is case-control
        out_eg = pwcoco_qtl_raw_prefix("eqtl_gwas", pqtl_dataset, p, eqtl_source, out_dir)
        out_eg.parent.mkdir(parents=True, exist_ok=True)
        pwcoco.pwcoco(
            ref_bfile=ref_bfile, sumstats_1=eqtl_h, sumstats_2=gwas_h,
            n_1=eqtl_n, n_2=n_cases + n_controls, n2_case=n_cases, out_dir=str(out_eg), threads=8
        )
        if Path(f"{out_eg}.coloc").exists():
            eqtl_gwas_rows.append(
                pl.read_csv(f"{out_eg}.coloc", separator="\t").with_columns(pl.lit(p).alias("protein"), pl.lit(dataset).alias("eqtl_dataset"), pl.lit(cell_type).alias("cell_type"))
            )

    # single-cell hits
    for row in sc.iter_rows(named=True):
        p = row["protein"]
        dataset = row["eqtl_dataset"]
        cell_type = row["cell_type"]
        probe = row["probeID"]

        dir = f"./dat/cis_regions/{pqtl_dataset}/{p}"
        dir = Path(dir)
        pqtl = dir / "pqtl.parquet"
        gwas = dir / "gwas.parquet"

        pqtl_df = pl.read_parquet(pqtl)
        gwas_df = pl.read_parquet(gwas)

        eqtl_file = Path(f"./dat/sc-eQTL/{dataset}/{cell_type}.parquet")
        if not eqtl_file.exists():
            continue

        base_gene = probe.split(".")[0]
        eqtl_df = (
            pl.scan_parquet(eqtl_file)
            .filter(pl.col("GENE").str.split(".").list.first() == base_gene)
            .collect()
            .sort("P")
            .unique(subset="SNP", keep="first")
        )

        if eqtl_df.height == 0:
            continue

        pqtl_h = pwcoco.harmonise_sumstats(pqtl_df, "SNP", "A1", "A2", resolve_maf_col(pqtl_df), "BETA", "SE", "P", "N")
        gwas_h = pwcoco.harmonise_sumstats(gwas_df, "SNP", "A1", "A2", resolve_maf_col(gwas_df), "BETA", "SE", "P", "N")
        eqtl_h = pwcoco.harmonise_sumstats(eqtl_df, "SNP", "A1", "A2", "FRQ", "BETA", "SE", "P", "N")

        eqtl_n = int(eqtl_h["n"][0])
        pqtl_n = int(pqtl_h["n"][0])
        eqtl_source = f"{dataset}_{cell_type}"

        out_ep = pwcoco_qtl_raw_prefix("eqtl_pqtl", pqtl_dataset, p, eqtl_source, out_dir)
        out_ep.parent.mkdir(parents=True, exist_ok=True)
        pwcoco.pwcoco(
            ref_bfile=ref_bfile, sumstats_1=eqtl_h, sumstats_2=pqtl_h,
            n_1=eqtl_n, n_2=pqtl_n, n2_case=0, out_dir=str(out_ep), threads=8
        )
        if Path(f"{out_ep}.coloc").exists():
            eqtl_pqtl_rows.append(
                pl.read_csv(f"{out_ep}.coloc", separator="\t").with_columns(pl.lit(p).alias("protein"), pl.lit(dataset).alias("eqtl_dataset"), pl.lit(cell_type).alias("cell_type"))
            )

        out_eg = pwcoco_qtl_raw_prefix("eqtl_gwas", pqtl_dataset, p, eqtl_source, out_dir)
        out_eg.parent.mkdir(parents=True, exist_ok=True)
        pwcoco.pwcoco(
            ref_bfile=ref_bfile, sumstats_1=eqtl_h, sumstats_2=gwas_h,
            n_1=eqtl_n, n_2=n_cases + n_controls, n2_case=n_cases, out_dir=str(out_eg), threads=8
        )
        if Path(f"{out_eg}.coloc").exists():
            eqtl_gwas_rows.append(
                pl.read_csv(f"{out_eg}.coloc", separator="\t").with_columns(pl.lit(p).alias("protein"), pl.lit(dataset).alias("eqtl_dataset"), pl.lit(cell_type).alias("cell_type"))
            )

    eqtl_pqtl_df = pl.concat(eqtl_pqtl_rows, how="diagonal_relaxed") if eqtl_pqtl_rows else pl.DataFrame()
    eqtl_gwas_df = pl.concat(eqtl_gwas_rows, how="diagonal_relaxed") if eqtl_gwas_rows else pl.DataFrame()

    ep_out = pwcoco_eqtl_pqtl_out(pqtl_dataset, pheno_id, out_dir)
    ep_out.parent.mkdir(parents=True, exist_ok=True)
    if eqtl_pqtl_df.height > 0:
        eqtl_pqtl_df.write_csv(ep_out, separator="\t")

    eg_out = pwcoco_eqtl_gwas_out(pqtl_dataset, pheno_id, out_dir)
    eg_out.parent.mkdir(parents=True, exist_ok=True)
    if eqtl_gwas_df.height > 0:
        eqtl_gwas_df.write_csv(eg_out, separator="\t")

    # triangulation across all 3 combos (pQTL-GWAS from pwcoco_out(), eQTL-pQTL,
    # eQTL-GWAS) - a target only counts when the SAME SNP clears pp4_thresh in ALL
    # 3, not just some of them. This is the eQTL-level equivalent of HyPrColoc: it
    # answers the same "do pQTL/eQTL/GWAS share 1 causal variant" question, but via
    # conditioning rather than HyPrColoc's single-causal-variant cluster assumption -
    # co-equal to HyPrColoc, not a downstream refinement of it.
    pg_file = pwcoco_out(pqtl_dataset, pheno_id, out_dir)
    pqtl_gwas_df = pl.read_csv(pg_file, separator="\t") if Path(pg_file).exists() else pl.DataFrame()

    pg_map = snp_h4_map(pqtl_gwas_df, pp4_thresh) if pqtl_gwas_df.height > 0 else {}
    ep_map = snp_h4_map(eqtl_pqtl_df, pp4_thresh) if eqtl_pqtl_df.height > 0 else {}
    eg_map = snp_h4_map(eqtl_gwas_df, pp4_thresh) if eqtl_gwas_df.height > 0 else {}

    shared_rows = []
    for protein in set(pg_map) & set(ep_map) & set(eg_map):
        shared_snps = set(pg_map[protein]) & set(ep_map[protein]) & set(eg_map[protein])
        for snp in shared_snps:
            shared_rows.append({
                "protein": protein,
                "snp": snp,
                "pqtl_gwas_h4": pg_map[protein][snp],
                "eqtl_pqtl_h4": ep_map[protein][snp],
                "eqtl_gwas_h4": eg_map[protein][snp],
            })

    if shared_rows:
        shared_file = pwcoco_qtl_shared_out(pqtl_dataset, pheno_id, out_dir)
        shared_file.parent.mkdir(parents=True, exist_ok=True)
        pl.DataFrame(shared_rows).write_csv(shared_file, separator="\t")


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--pqtl_dataset", required=True, choices=["ukb_ppp", "decode", "wu_csf", "wingo_brain"])
    p.add_argument("--pheno_id", required=True)
    p.add_argument("--ref_bfile", required=True)
    p.add_argument("--n_cases", required=True, type=int)
    p.add_argument("--n_controls", required=True, type=int)
    p.add_argument("--local_results_dir", default="results")
    p.add_argument("--pp4_threshold", type=float, default=0.7)
    args = p.parse_args()

    pwcoco_qtl_wrapper(
        pqtl_dataset=args.pqtl_dataset,
        pheno_id=args.pheno_id,
        ref_bfile=args.ref_bfile,
        n_cases=args.n_cases,
        n_controls=args.n_controls,
        out_dir=args.local_results_dir,
        pp4_thresh=args.pp4_threshold,
    )


if __name__ == "__main__":
    main()