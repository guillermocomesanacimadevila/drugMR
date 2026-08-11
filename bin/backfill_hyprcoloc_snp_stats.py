#!/usr/bin/env python3
import argparse
import sys
from pathlib import Path

import polars as pl

sys.path.insert(0, str(Path(__file__).resolve().parent))
from hyprcoloc_targets import CANDIDATE_SNP_STAT_COLS, load_eqtl_table, resolve_candidate_snp_stats

# One-off backfill for HyPrColoc master tables written before attach_candidate_snp_stats
# existed (see hyprcoloc_targets.py) - those rows carry a candidate_snp but none of its
# aligned alleles/betas, which is what the dashboard's Final Targets tab needs to report
# the HyPrColoc SNP the same way it reports the top cis-pQTL SNP. Re-resolves each row's
# candidate SNP against the exact cis-region GWAS/pQTL parquets and eQTL file it was
# originally run against - HyPrColoc/Rscript itself never re-runs, since the candidate
# SNP it already found is taken as-is.
#
# Older master files don't carry every join key HyPrColoc itself now stamps on each row
# (probeID never made it in; some pre-multi-eqtl-dataset runs are also missing data_type/
# eqtl_dataset) - whichever of [protein, cell_type, data_type, eqtl_dataset] a file DOES
# have is used to look up the rest from that pQTL dataset's SMR final_multi_omics_targets
# table, which still carries probeID for every target x cell-type/tissue hit.


def load_smr_lookup(pqtl_dataset: str, pheno_id: str):
    smr_file = Path(f"./results/SMR/{pqtl_dataset}_{pheno_id}_final_multi_omics_targets.tsv")

    if not smr_file.exists():
        return None

    lookup_cols = ["protein", "cell_type", "data_type", "eqtl_dataset", "probeID"]
    return pl.read_csv(smr_file, separator="\t", null_values=["NA"]).select(lookup_cols).unique()


def load_trio(pqtl_dataset: str, protein: str, cell_type: str, data_type: str, eqtl_dataset: str, probe_id: str):
    cis_region = Path(f"./dat/cis_regions/{pqtl_dataset}/{protein}")
    gwas_file = cis_region / "gwas.parquet"
    pqtl_file = cis_region / "pqtl.parquet"

    if not gwas_file.exists() or not pqtl_file.exists():
        print(f"[CONCERN] Missing cis-region parquet(s) for {protein}")
        return None

    base_gene_id = str(probe_id).split(".")[0]
    eqtl = load_eqtl_table(data_type, eqtl_dataset, cell_type, base_gene_id)

    if eqtl is None or eqtl.height == 0:
        print(f"[CONCERN] No {cell_type} eQTL rows found for probe {probe_id}")
        return None

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

    return gwas, pqtl, eqtl


def backfill_dataset(pqtl_dataset: str, pheno_id: str):
    master_file = Path(f"./results/hyprcoloc/{pqtl_dataset}_{pheno_id}_all_hyprcoloc.tsv")

    if not master_file.exists():
        print(f"[SKIP] No HyPrColoc master file for {pqtl_dataset}: {master_file}")
        return

    hypr = pl.read_csv(master_file, separator="\t", null_values=["NA"])

    if hypr.height == 0:
        print(f"[SKIP] {pqtl_dataset} HyPrColoc master file is empty")
        return

    if set(CANDIDATE_SNP_STAT_COLS).issubset(hypr.columns):
        print(f"[SKIP] {pqtl_dataset} already carries candidate SNP stats")
        return

    smr_lookup = load_smr_lookup(pqtl_dataset, pheno_id)

    if smr_lookup is None:
        print(f"[CONCERN] No SMR final targets file for {pqtl_dataset} - cannot recover probeID, skipping")
        return

    join_keys = [key for key in ["protein", "cell_type", "data_type", "eqtl_dataset"] if key in hypr.columns]
    hypr = hypr.join(smr_lookup, on=join_keys, how="left")

    if "probeID" not in hypr.columns or hypr.get_column("probeID").null_count() == hypr.height:
        print(f"[CONCERN] Could not resolve probeID for any {pqtl_dataset} row via {join_keys} - skipping")
        return

    trio_cache = {}
    stat_rows = []

    for row in hypr.iter_rows(named=True):
        candidate_snp = row.get("candidate_snp")
        probe_id = row.get("probeID")

        if candidate_snp is None or probe_id is None:
            stat_rows.append(None)
            continue

        trio_key = (row["protein"], row["cell_type"], row.get("data_type"), row.get("eqtl_dataset"), probe_id)

        if trio_key not in trio_cache:
            trio_cache[trio_key] = load_trio(pqtl_dataset, *trio_key)

        trio = trio_cache[trio_key]
        stat_rows.append(resolve_candidate_snp_stats(candidate_snp, *trio) if trio is not None else None)

    stats_df = pl.DataFrame(
        [row if row is not None else {col: None for col in CANDIDATE_SNP_STAT_COLS} for row in stat_rows],
        schema={
            "a1": pl.Utf8, "a2": pl.Utf8,
            "gwas_beta": pl.Float64, "gwas_p": pl.Float64,
            "pqtl_beta": pl.Float64, "pqtl_p": pl.Float64,
            "eqtl_beta": pl.Float64, "eqtl_p": pl.Float64
        }
    )

    result = pl.concat([hypr.drop("probeID"), stats_df], how="horizontal")
    n_resolved = hypr.height - stat_rows.count(None)
    print(f"[DONE] {pqtl_dataset}: resolved candidate SNP stats for {n_resolved}/{hypr.height} rows")
    result.write_csv(master_file, separator="\t")


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--pqtl_dataset", nargs="+", default=["ukb_ppp", "decode", "wu_csf", "wingo_brain"])
    p.add_argument("--pheno_id", default="AD")
    args = p.parse_args()

    for pqtl_dataset in args.pqtl_dataset:
        backfill_dataset(pqtl_dataset, args.pheno_id)


if __name__ == "__main__":
    main()
