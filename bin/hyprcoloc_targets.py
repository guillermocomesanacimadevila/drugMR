#!/usr/bin/env python3
import argparse
import subprocess
import polars as pl
from pathlib import Path
from drugmr import extract_common_snps
from drugmr import paths

# HyPrColoc for the final multi-omics targets
# -> for every target x cell-type/tissue hit that passed cis-MR + COLOC + SMR + HEIDI
#    in a bulk or single-cell eQTL dataset, run a 3-trait (pQTL / GWAS / eQTL) HyPrColoc
#    in that target's cis-region, restricted to the SNPs shared across all three
#    and aligned onto a common effect allele (drugmr.extract_common_snps aligns
#    onto the GWAS A1, same "align to GWAS" convention as compile_cis_hit_info.py
#    / sort_smr.py's align_to_risk_allele)
# both bulk (dat/bulk-eQTL/{GTEx_v10,MetaBrain}/...) and single-cell (dat/sc-eQTL/...)
# eQTL sources are supported - which one is used is driven entirely by the "data_type"
# tag already carried on each row of the combined SMR target table


# candidate_snp's own alleles/betas, aligned to the AD risk allele (b_GWAS > 0) -
# same "make A1 the GWAS risk allele" convention as bin/compile_cis_hit_info.py,
# extended from 2 traits (GWAS/pQTL) to 3 (GWAS/pQTL/eQTL). Looks the SNP up in the
# pre-alignment gwas/pqtl/eqtl tables (which still carry each trait's own A1/A2/P -
# extract_common_snps' output drops those) rather than the post-alignment `matched`
# tables. Returns None if the SNP is missing from any trait or its pQTL/eQTL allele
# pair doesn't match the GWAS one either way round (shouldn't happen for a SNP that
# HyPrColoc was run on, but guards against a silent mismatch if it ever does).
def resolve_candidate_snp_stats(snp: str, gwas: pl.DataFrame, pqtl: pl.DataFrame, eqtl: pl.DataFrame):
    gwas_row = gwas.filter(pl.col("SNP") == snp)
    pqtl_row = pqtl.filter(pl.col("SNP") == snp)
    eqtl_row = eqtl.filter(pl.col("SNP") == snp)

    if gwas_row.height == 0 or pqtl_row.height == 0 or eqtl_row.height == 0:
        return None

    gwas_row = gwas_row.row(0, named=True)
    pqtl_row = pqtl_row.row(0, named=True)
    eqtl_row = eqtl_row.row(0, named=True)

    a1, a2 = str(gwas_row["A1"]).upper(), str(gwas_row["A2"]).upper()
    gwas_beta, gwas_p = float(gwas_row["BETA"]), float(gwas_row["P"])

    if gwas_beta < 0:
        a1, a2 = a2, a1
        gwas_beta = -gwas_beta

    def realign(row):
        row_a1, row_a2 = str(row["A1"]).upper(), str(row["A2"]).upper()
        beta, p = float(row["BETA"]), float(row["P"])
        if row_a1 == a1 and row_a2 == a2:
            return beta, p
        if row_a1 == a2 and row_a2 == a1:
            return -beta, p
        return None, None

    pqtl_beta, pqtl_p = realign(pqtl_row)
    eqtl_beta, eqtl_p = realign(eqtl_row)

    if pqtl_beta is None or eqtl_beta is None:
        print(f"[CONCERN] Allele mismatch at candidate SNP {snp} against GWAS {a1}/{a2} - skipping SNP-level stats")
        return None

    return {
        "a1": a1,
        "a2": a2,
        "gwas_beta": gwas_beta,
        "gwas_p": 1e-300 if gwas_p == 0 else gwas_p,
        "pqtl_beta": pqtl_beta,
        "pqtl_p": 1e-300 if pqtl_p == 0 else pqtl_p,
        "eqtl_beta": eqtl_beta,
        "eqtl_p": 1e-300 if eqtl_p == 0 else eqtl_p,
    }


# attaches, for every result row, the candidate SNP's own aligned alleles/betas (see
# resolve_candidate_snp_stats) - null-filled where the SNP is missing or a row has no
# candidate_snp at all (e.g. a cluster HyPrColoc couldn't resolve to a single SNP)
CANDIDATE_SNP_STAT_COLS = ["a1", "a2", "gwas_beta", "gwas_p", "pqtl_beta", "pqtl_p", "eqtl_beta", "eqtl_p"]


def attach_candidate_snp_stats(result_df: pl.DataFrame, gwas: pl.DataFrame, pqtl: pl.DataFrame, eqtl: pl.DataFrame):
    if result_df.height == 0 or "candidate_snp" not in result_df.columns:
        return result_df

    stat_rows = [
        resolve_candidate_snp_stats(snp, gwas, pqtl, eqtl) if snp is not None else None
        for snp in result_df.get_column("candidate_snp").to_list()
    ]
    stats_df = pl.DataFrame(
        [row if row is not None else {col: None for col in CANDIDATE_SNP_STAT_COLS} for row in stat_rows],
        schema={
            "a1": pl.Utf8, "a2": pl.Utf8,
            "gwas_beta": pl.Float64, "gwas_p": pl.Float64,
            "pqtl_beta": pl.Float64, "pqtl_p": pl.Float64,
            "eqtl_beta": pl.Float64, "eqtl_p": pl.Float64
        }
    )
    return pl.concat([result_df, stats_df], how="horizontal")


# bulk eQTL parquet files are laid out differently per dataset: GTEx_v10 is
# tissue-resolved (1 file per tissue, cell_type carries a "GTEx_<Tissue>_v10" label
# that maps onto the dat/bulk-eQTL/GTEx_v10/<Tissue>/ directory name), MetaBrain is
# flat (1 genome-wide file regardless of the cell_type label). Returns None for any
# other/unrecognised bulk eqtl_dataset.
def resolve_bulk_eqtl_file(eqtl_dataset: str, cell_type: str):
    if eqtl_dataset == "GTEx_v10":
        tissue = cell_type.removeprefix("GTEx_").removesuffix("_v10")
        return Path(f"./dat/bulk-eQTL/GTEx_v10/{tissue}/{tissue}.parquet")
    if eqtl_dataset == "MetaBrain":
        return Path("./dat/bulk-eQTL/MetaBrain/BrainMeta_cis_eQTL.parquet")
    return None


# loads the 1 gene's eQTL rows (bulk or single-cell), aligned onto the pipeline's own
# SNP/A1/A2/BETA/SE/P convention - shared by the main HyPrColoc loop below and
# bin/backfill_hyprcoloc_snp_stats.py, which needs the exact same table to
# re-resolve a candidate SNP's stats for older result files without re-running
# HyPrColoc itself. Returns None (with a printed [CONCERN]) on any missing file or
# unrecognised data_type, same as the inline version this replaced.
def load_eqtl_table(data_type: str, eqtl_dataset: str, cell_type: str, base_gene_id: str):
    if data_type == "single_cell":
        eqtl_file = Path(f"./dat/sc-eQTL/{eqtl_dataset}/{cell_type}.parquet")

        if not eqtl_file.exists():
            print(f"[CONCERN] Missing {eqtl_dataset} eQTL file for {cell_type}: {eqtl_file}")
            return None

        # sc-eQTL files carry ref/alt as A1/A2 and the actual effect allele as EA -
        # re-point A1/A2 so A1 is always the effect allele the BETA belongs to,
        # same idea as sort_smr.py's pull_original_sc_eqtl_beta
        return (
            pl.scan_parquet(eqtl_file)
            .filter(pl.col("GENE").str.split(".").list.first() == base_gene_id)
            .select(["SNP", "A1", "A2", "EA", "BETA", "SE", "P"])
            .with_columns(
                pl.col("EA").alias("eqtl_a1"),
                pl.when(pl.col("EA") == pl.col("A2")).then(pl.col("A1")).otherwise(pl.col("A2")).alias("eqtl_a2")
            )
            .select(["SNP", pl.col("eqtl_a1").alias("A1"), pl.col("eqtl_a2").alias("A2"), "BETA", "SE", "P"])
            .sort("P")
            .unique(subset="SNP", keep="first")
            .collect()
        )

    if data_type == "bulk":
        eqtl_file = resolve_bulk_eqtl_file(eqtl_dataset, cell_type)

        if eqtl_file is None or not eqtl_file.exists():
            print(f"[CONCERN] Missing {eqtl_dataset} bulk eQTL file for {cell_type}: {eqtl_file}")
            return None

        # bulk eQTL parquets come straight from an SMR besd/esi/epi query, so A1 is
        # already the effect allele b belongs to (SMR's own convention) - no
        # re-pointing needed, just rename onto the pipeline's BETA/P convention
        return (
            pl.scan_parquet(eqtl_file)
            .filter(pl.col("Probe").str.split(".").list.first() == base_gene_id)
            .select(["SNP", "A1", "A2", pl.col("b").alias("BETA"), "SE", pl.col("p").alias("P")])
            .sort("P")
            .unique(subset="SNP", keep="first")
            .collect()
        )

    print(f"[CONCERN] Unrecognised data_type '{data_type}' for {cell_type} - skipping")
    return None


def hyprcoloc_targets(pqtl_dataset: str, pheno_id: str, eqtl_dataset: str, local_results_dir: str = "results"):
    hyprcoloc_script = "./bin/hyprcoloc.R"
    work_dir = Path(f"./work/hyprcoloc/{pqtl_dataset}/{pheno_id}")
    work_dir.mkdir(parents=True, exist_ok=True)
    out_dir = paths.hyprcoloc_dataset_out(pqtl_dataset, eqtl_dataset, pheno_id, local_results_dir).parent.parent
    out_dir.mkdir(parents=True, exist_ok=True)
    targets_file = paths.smr_final_targets_out(pqtl_dataset, pheno_id, local_results_dir)
    targets = pl.read_csv(targets_file, separator="\t", null_values=["NA"])
    targets = (
        targets
        .filter(pl.col("eqtl_dataset") == eqtl_dataset)
        .select(["protein", "cell_type", "probeID", "data_type"])
        .unique()
        .sort(["protein", "cell_type"])
    )

    print(f"[TRACKING] {targets.height} target x cell-type/tissue hit(s) found for HyPrColoc in {eqtl_dataset}")

    results = []

    for row in targets.iter_rows(named=True):
        protein = row["protein"]
        cell_type = row["cell_type"]
        probe_id = row["probeID"]
        data_type = row["data_type"]
        cis_region = Path(f"./dat/cis_regions/{pqtl_dataset}/{protein}")
        gwas_file = cis_region / "gwas.parquet"
        pqtl_file = cis_region / "pqtl.parquet"

        if not gwas_file.exists() or not pqtl_file.exists():
            print(f"[CONCERN] Missing cis-region parquet(s) for {protein}")
            continue

        # probeID (from the SMR .epi annotation) and the eQTL parquet's own gene ID
        # column can carry different Ensembl release versions for the same gene
        # (e.g. ENSG00000095585.20 vs ENSG00000095585.17) - match on the
        # version-stripped base ID rather than the raw string
        base_gene_id = probe_id.split(".")[0]

        eqtl = load_eqtl_table(data_type, eqtl_dataset, cell_type, base_gene_id)

        if eqtl is None:
            continue

        # 1 row per SNP, most significant kept - same pattern as coloc.R
        gwas = (
            pl.read_parquet(gwas_file)
            .select(["SNP", "A1", "A2", "BETA", "SE", "P"])
            .sort("P")
            .unique(subset="SNP", keep="first")
        )

        pqtl = (
            pl.read_parquet(pqtl_file)
            .select(["SNP", "A1", "A2", "BETA", "SE", "P"])
            .sort("P")
            .unique(subset="SNP", keep="first")
        )

        if eqtl.height == 0:
            print(f"[CONCERN] No {cell_type} eQTL rows found for probe {probe_id}")
            continue

        matched = extract_common_snps(
            {"pqtl": pqtl, "gwas": gwas, "eqtl": eqtl},
            reference="gwas"
        )

        n_shared = matched["gwas"].height
        print(f"[TRACKING] {protein} x {cell_type}: {n_shared} SNPs shared across pQTL/GWAS/eQTL")

        if n_shared < 2:
            print(f"[CONCERN] Fewer than 2 shared SNPs for {protein} x {cell_type} - skipping HyPrColoc")
            continue

        trio_dir = work_dir / f"{protein}__{cell_type}"
        trio_dir.mkdir(parents=True, exist_ok=True)

        for name, trait_df in matched.items():
            trait_df.write_parquet(trio_dir / f"{name}.parquet")

        cmd_hyprcoloc = [
            "Rscript", hyprcoloc_script,
            pqtl_dataset, protein, cell_type, pheno_id, str(trio_dir), local_results_dir
        ]
        print(f"[TRACKING] Running HyPrColoc for {protein} x {cell_type}")
        subprocess.run(cmd_hyprcoloc, check=True)

        result_file = out_dir / f"{pheno_id}_{protein}_{cell_type}_hyprcoloc.tsv"

        if result_file.exists():
            result_df = pl.read_csv(result_file, separator="\t").with_columns(
                pl.lit(eqtl_dataset).alias("eqtl_dataset"),
                pl.lit(data_type).alias("data_type")
            )
            result_df = attach_candidate_snp_stats(result_df, gwas, pqtl, eqtl)
            results.append(result_df)
            result_file.unlink()
        else:
            print(f"[CONCERN] Expected HyPrColoc output not found: {result_file}")

    if len(results) == 0:
        print(f"[CONCERN] No HyPrColoc results generated for any {eqtl_dataset} target")
        return

    dataset_results = pl.concat(results, how="diagonal_relaxed")

    # per-dataset output, used as the idempotency marker by local.py / hpc.py (mirrors
    # sort_smr.py's per-bulk_dataset promising_targets_SMR.tsv pattern) so re-running
    # for one eqtl_dataset doesn't require re-running every other one
    per_dataset_file = paths.hyprcoloc_dataset_out(pqtl_dataset, eqtl_dataset, pheno_id, local_results_dir)
    per_dataset_file.parent.mkdir(parents=True, exist_ok=True)
    dataset_results.write_csv(per_dataset_file, separator="\t")

    # canonical combined output (bulk + single-cell hits together) - upsert: drop any
    # stale rows for this eqtl_dataset, then append the fresh ones, so bulk and
    # single-cell runs (in either order) compose instead of overwriting each other
    master_file = paths.hyprcoloc_out(pqtl_dataset, pheno_id, local_results_dir)
    if master_file.exists() and master_file.stat().st_size > 0:
        existing = pl.read_csv(master_file, separator="\t", null_values=["NA"])
        if "eqtl_dataset" in existing.columns:
            existing = existing.filter(pl.col("eqtl_dataset") != eqtl_dataset)
        master = pl.concat([existing, dataset_results], how="diagonal_relaxed")
    else:
        master = dataset_results

    master.write_csv(master_file, separator="\t")
    print(f"[DONE] Saved master HyPrColoc table: {master_file}")


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--pqtl_dataset", required=True, choices=["ukb_ppp", "decode", "wu_csf", "wingo_brain"])
    p.add_argument("--pheno_id", required=True)
    p.add_argument("--eqtl_dataset", default="SingleBrain")
    p.add_argument("--local_results_dir", default="results")
    args = p.parse_args()
    hyprcoloc_targets(
        pqtl_dataset=args.pqtl_dataset,
        pheno_id=args.pheno_id,
        eqtl_dataset=args.eqtl_dataset,
        local_results_dir=args.local_results_dir
    )


if __name__ == "__main__":
    main()
