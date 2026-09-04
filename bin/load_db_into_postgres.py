import argparse
import json
from pathlib import Path

import polars as pl
from sqlalchemy import text

from drugmr import paths
from drugmr.db import get_engine

""" Simple script to load data from TSVs onto respective SQL table """

PROJECT_ROOT = Path(__file__).resolve().parents[1]

VALID_TABLES = [
    "cis_mr_results",
    "coloc_results",
    "pwcoco_results",
    "pwcoco_eqtl_pqtl_results",
    "pwcoco_eqtl_gwas_results",
    "smr_results",
    "hyprcoloc_results",
    "finngen_phewas_safety",
    "ukb_phewas_safety",
    "target_stats"
]


class PostgresLoader:

    """ Function accumulator-like object for connecting/updating to psql databse for run X """

    def __init__(
            self,
            run_id: str,
            db_id: str = None):

        self.run_id = run_id
        self.engine = get_engine(db_id=db_id)
        self.manifest = self._load_manifest(run_id)
        self._ensure_run_row()

    @staticmethod
    def _load_manifest(run_id: str) -> dict:
        manifest_path = paths.run_manifest_path(run_id, root=str(PROJECT_ROOT / "runs"))
        if not manifest_path.exists():
            raise FileNotFoundError(f"No manifest found for run '{run_id}': {manifest_path}")
        with open(manifest_path) as f:
            return json.load(f)

    def _ensure_run_row(self):
        with self.engine.begin() as conn:
            conn.execute(
                text("""
                    INSERT INTO runs (run_id, pheno_id, pqtl_dataset, git_sha7, mode, image_name, overwrite)
                    VALUES (:run_id, :pheno_id, :pqtl_dataset, :git_sha7, :mode, :image_name, :overwrite)
                    ON CONFLICT (run_id) DO NOTHING
                """),
                {
                    "run_id": self.run_id,
                    "pheno_id": self.manifest["pheno_id"],
                    "pqtl_dataset": self.manifest["pqtl_dataset"],
                    "git_sha7": self.manifest["git_sha7"],
                    "mode": self.manifest["mode"],
                    "image_name": self.manifest["image_name"],
                    "overwrite": self.manifest.get("overwrite"),
                },
            )

    def load_table(self, results_file: str, pqtl_dataset: str, table: str):
        if table not in VALID_TABLES:
            raise ValueError(f"Unknown table '{table}' - expected one of {VALID_TABLES}")

        df = pl.read_csv(Path(results_file), separator="\t")
        df.columns = [c.lower().replace(".", "_") for c in df.columns]
        df = df.with_columns(
            pl.lit(self.run_id).alias("run_id"),
            pl.lit(pqtl_dataset).alias("pqtl_dataset"),
        )

        # re-runs of the same run_id replace just that run's rows, not the whole table
        with self.engine.begin() as conn:
            conn.execute(text(f"DELETE FROM {table} WHERE run_id = :run_id"), {"run_id": self.run_id})

        df.to_pandas().to_sql(table, self.engine, if_exists="append", index=False)
        print(f"[DONE] Loaded {df.height:,} rows into {table} for run_id={self.run_id}")


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--results_file", required=True, type=str)
    p.add_argument("--run_id", required=True, type=str)
    p.add_argument("--pqtl_dataset", required=True, type=str)
    p.add_argument("--table", required=True, type=str)
    p.add_argument("--db_id", required=False, type=str, default=None)
    args = p.parse_args()

    loader = PostgresLoader(run_id=args.run_id, db_id=args.db_id)
    loader.load_table(
        results_file=args.results_file,
        pqtl_dataset=args.pqtl_dataset,
        table=args.table,
    )


if __name__ == "__main__":
    main()
