import glob
from pathlib import Path

import polars as pl


class QTLManifest:

    def __init__(self, manifest_path: str):
        self.manifest_path = manifest_path
        self._manifest = pl.read_csv(manifest_path)

    def get_row(self, dataset: str) -> dict:
        rows = self._manifest.filter(pl.col("dataset") == dataset)
        if rows.height == 0:
            raise ValueError(f"'{dataset}' not found in manifest: {self.manifest_path}")
        return rows.row(0, named=True)

    @staticmethod
    def normalise_columns(df: pl.DataFrame, manifest_row: dict) -> pl.DataFrame:
        rename_map = {
            manifest_row["snp_col"]: "SNP",
            manifest_row["a1_col"]: "A1",
            manifest_row["a2_col"]: "A2",
            manifest_row["beta_col"]: "BETA",
            manifest_row["se_col"]: "SE",
            manifest_row["p_col"]: "P",
            manifest_row["chr_col"]: "CHR",
            manifest_row["pos_col"]: "BP",
        }
        return df.rename(rename_map)

    def resolve(self, dataset: str, gene: str):

        """
        One matched file  -> single pl.DataFrame, filtered to `gene`.
        Many matched files -> dict {file_label: pl.DataFrame}, one entry per
        matched file that actually has rows for `gene`
        """
        
        manifest_row = self.get_row(dataset)
        key_col = manifest_row.get("key_col") or None
        matched_files = sorted(glob.glob(manifest_row["path"]))
        if not matched_files:
            raise FileNotFoundError(f"No files matched path: {manifest_row['path']}")

        results = {}
        for f in matched_files:
            label = Path(f).stem
            df = self.normalise_columns(pl.read_parquet(f), manifest_row)

            if key_col:
                df = df.filter(pl.col(key_col) == gene)
                if df.height == 0:
                    continue  # gene not measured in this file/tissue - skip, not an error
            else:
                gene_from_filename = label.split("_")[0]
                if gene_from_filename != gene:
                    continue

            results[label] = df

        if not results:
            raise ValueError(f"No data for gene='{gene}' in dataset '{dataset}'")

        if len(matched_files) == 1:
            return next(iter(results.values()))  # single-file dataset: plain DataFrame
        return results  # multi-file dataset: {label: DataFrame}
