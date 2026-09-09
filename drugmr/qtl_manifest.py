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

    def get_rows_by_parent(self, parent_dataset: str) -> list[dict]:
        # case-insensitive on both sides - parent_dataset is stored mixed-case
        # (e.g. "GTEx_v10", matching the eqtl_dataset convention used elsewhere
        # in results/params), but every other manifest lookup (get_row) is
        # effectively case-insensitive too (callers lowercase before calling,
        # and "dataset" values are always stored lowercase) - this matches that
        rows = self._manifest.filter(pl.col("parent_dataset").str.to_lowercase() == parent_dataset.lower())
        return rows.to_dicts()

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
        eaf_col = manifest_row.get("eaf_col")
        if eaf_col:
            rename_map[eaf_col] = "FRQ"
        return df.rename(rename_map)

    def resolve(self, dataset: str, gene: str = None):

        """
        gene given -> One matched file  -> single pl.DataFrame, filtered to `gene`.
                      Many matched files -> dict {file_label: pl.DataFrame}, one
                      entry per matched file that actually has rows for `gene`.
        gene=None  -> whole-file mode, no filtering (e.g. for BESD conversion,
                      which needs every gene in the file at once). Same
                      single-file-vs-many-files return shape as above.
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

            if gene is not None:
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

        if len(results) == 1:
            return next(iter(results.values()))  # exactly one file matched (post gene-filter): plain DataFrame
        return results  # multiple files matched (post gene-filter): {label: DataFrame}
