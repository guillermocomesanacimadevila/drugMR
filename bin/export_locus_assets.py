#!/usr/bin/env python3
import argparse
import json
import shutil
from pathlib import Path

import polars as pl

from drugmr.utils import compute_ld_to_lead


INDEX_SCHEMA = {
    "protein": pl.Utf8,
    "candidate_snp": pl.Utf8,
    "candidate_source": pl.Utf8,
    "pqtl_dataset": pl.Utf8,
    "status": pl.Utf8,
    "ld_available": pl.Boolean,
}


def _read_target_stats(path: Path) -> pl.DataFrame:
    if not path.is_file():
        raise FileNotFoundError(f"Target statistics file not found: {path}")
    if path.stat().st_size == 0:
        raise ValueError(f"Target statistics file is empty: {path}")
    targets = pl.read_csv(path, separator="\t")
    missing = {"protein", "SNP"}.difference(targets.columns)
    if missing:
        raise ValueError(
            f"Target statistics file is missing column(s): {', '.join(sorted(missing))}"
        )
    return targets


def export_locus_assets(
    pqtl_dataset: str,
    cis_regions_dir: Path,
    target_stats: Path,
    ref_bfile: Path,
    out_dir: Path,
    window_kb: int = 5000,
) -> pl.DataFrame:
    """Package regional GWAS/pQTL data and LD for portable dashboard use."""
    cis_regions_dir = Path(cis_regions_dir)
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    targets = _read_target_stats(Path(target_stats))
    index_rows = []

    for row in targets.unique(subset=["protein"], keep="first").iter_rows(named=True):
        protein = str(row["protein"])
        candidate_snp = str(row["SNP"])
        source_dir = cis_regions_dir / protein
        pqtl_file = source_dir / "pqtl.parquet"
        gwas_file = source_dir / "gwas.parquet"

        base_index = {
            "protein": protein,
            "candidate_snp": candidate_snp,
            "candidate_source": "top_cis_pqtl",
            "pqtl_dataset": pqtl_dataset,
        }
        if not pqtl_file.is_file() or not gwas_file.is_file():
            print(f"[CONCERN] Missing cis-region files for {protein}; skipping locus bundle")
            index_rows.append({**base_index, "status": "missing_cis_data", "ld_available": False})
            continue

        target_dir = out_dir / protein
        target_dir.mkdir(parents=True, exist_ok=True)
        shutil.copy2(pqtl_file, target_dir / "pqtl.parquet")
        shutil.copy2(gwas_file, target_dir / "gwas.parquet")
        pqtl = pl.read_parquet(pqtl_file)
        gwas = pl.read_parquet(gwas_file)
        candidate_rows = pl.concat(
            [pqtl.select(["SNP", "CHR", "BP"]), gwas.select(["SNP", "CHR", "BP"])],
            how="vertical",
        ).filter(pl.col("SNP").cast(pl.Utf8) == candidate_snp)

        ld_available = False
        ld_error = None
        chromosome = None
        candidate_bp = None
        if candidate_rows.height:
            candidate = candidate_rows.row(0, named=True)
            chromosome = candidate["CHR"]
            candidate_bp = candidate["BP"]
            try:
                compute_ld_to_lead(
                    ref_bfile=ref_bfile,
                    lead_snp=candidate_snp,
                    chromosome=chromosome,
                    out_file=target_dir / "ld.parquet",
                    window_kb=window_kb,
                )
                ld_available = True
            except (FileNotFoundError, RuntimeError, ValueError) as error:
                ld_error = str(error)
                print(f"[CONCERN] LD unavailable for {protein}: {error}")
        else:
            ld_error = f"Candidate SNP {candidate_snp} not found in the cis-region files"
            print(f"[CONCERN] {protein}: {ld_error}")

        metadata = {
            **base_index,
            "chromosome": chromosome,
            "candidate_bp": candidate_bp,
            "pqtl_rows": pqtl.height,
            "gwas_rows": gwas.height,
            "ld_available": ld_available,
            "ld_error": ld_error,
            "ld_window_kb": window_kb,
        }
        (target_dir / "metadata.json").write_text(
            json.dumps(metadata, indent=2, default=str) + "\n"
        )
        status = "complete" if ld_available else "complete_without_ld"
        index_rows.append({**base_index, "status": status, "ld_available": ld_available})

    index = pl.DataFrame(index_rows, schema=INDEX_SCHEMA)
    index.write_csv(out_dir / "index.tsv", separator="\t")
    print(f"[DONE] Exported {index.height} portable locus bundle(s) to {out_dir}")
    return index


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--pqtl-dataset", required=True)
    parser.add_argument("--cis-regions-dir", required=True, type=Path)
    parser.add_argument("--target-stats", required=True, type=Path)
    parser.add_argument("--ref-bfile", required=True, type=Path)
    parser.add_argument("--out-dir", required=True, type=Path)
    parser.add_argument("--window-kb", type=int, default=5000)
    args = parser.parse_args()
    export_locus_assets(
        pqtl_dataset=args.pqtl_dataset,
        cis_regions_dir=args.cis_regions_dir,
        target_stats=args.target_stats,
        ref_bfile=args.ref_bfile,
        out_dir=args.out_dir,
        window_kb=args.window_kb,
    )


if __name__ == "__main__":
    main()
