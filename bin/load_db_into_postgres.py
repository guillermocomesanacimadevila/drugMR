#!/usr/bin/env python3
from sqlalchemy import create_engine, text
import polars as pl 
import subprocess
import argparse
from pathlib import Path

# create db drugmr
# psql drugmr -c ""
# SELECT * FROM cis_mr_results LIMIT 5; 

def master_postgres(mr_res: str, db_id: str, pqtl_dataset: str, pheno_id: str, table: str):
    mr_res = Path(mr_res)
    analysis_id = f"{pheno_id}_{pqtl_dataset}"
    engine = create_engine("postgresql+psycopg2:///postgres")

    # baseline cmd command
    cmd = f"""
set -euo pipefail
createdb {db_id}
    """
    with engine.connect() as conn:
        result = conn.execute(
            text(
                """
                SELECT 1
                FROM pg_database
                WHERE datname = :db_name
                """
            ),
            {"db_name": db_id}
        )

        db_exists = result.scalar() is not None

    if db_exists:
        print(f"[TRACKING] PostgreSQL database '{db_id}' already exists.")
    else:
        print(f"[TRACKING] PostgreSQL database '{db_id}' does not exist.")
        print(f"[POSTGRESQL] Creating SQL database: '{db_id}'...")
        subprocess.run(cmd, shell=True, check=True, executable="/bin/bash")
        print(f"[POSTGRESQL] SQL database '{db_id}' created!")

    df = pl.read_csv(mr_res, separator="\t")
    df.columns = [c.lower() for c in df.columns]
    df = df.with_columns(
        pl.lit(analysis_id).alias("analysis_id"),
        pl.lit(pqtl_dataset).alias("pqtl_dataset")
    )

    engine = create_engine(f"postgresql+psycopg2:///{db_id}")
    df.to_pandas().to_sql(table, engine, if_exists="replace", index=False)
    print(f"[DONE] Loaded {df.height:,} rows into {db_id}.{table}")


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--mr_res", required=True, type=str)
    p.add_argument("--db_id", required=True, type=str)
    p.add_argument("--pqtl_dataset", required=True, type=str)
    p.add_argument("--pheno_id", required=True, type=str)
    p.add_argument("--table", default="cis_mr_results", type=str)
    args = p.parse_args()
    master_postgres(
        mr_res=args.mr_res,
        db_id=args.db_id,
        pqtl_dataset=args.pqtl_dataset,
        pheno_id=args.pheno_id,
        table=args.table
    )

if __name__ == "__main__":
    main()