import glob
from pathlib import Path

import polars as pl


class QTLManifest:

    def __init__(self, manifest_path: str, base_dir: str | Path | None = None):
        self.base_dir = Path(base_dir) if base_dir else None
        self.manifest_path = str(self.base_dir / manifest_path) if self.base_dir else manifest_path
        self._manifest = pl.read_csv(self.manifest_path)

    def _resolve_row_path(self, row: dict) -> dict:
        # every caller (get_row/get_rows_by_parent/get_rows_by_qtl_type) that hands back a manifest
        # row must return an already-resolved `path`, not the raw CSV value -
        # several call sites in drugmr/smr.py read row["path"] directly (glob.glob
        # or Path()) without ever calling resolve(), so applying base_dir only
        # inside resolve() left those call sites broken under a foreign CWD
        # (e.g. a Nextflow task sandbox). Resolving it once, here, fixes every
        # consumer at once instead of patching each call site individually.
        if self.base_dir and row.get("path"):
            row = dict(row)
            row["path"] = str(self.base_dir / row["path"])
        return row

    def get_row(self, dataset: str) -> dict:
        rows = self._manifest.filter(pl.col("dataset").str.to_lowercase() == dataset.lower())
        if rows.height == 0:
            raise ValueError(f"'{dataset}' not found in manifest: {self.manifest_path}")
        return self._resolve_row_path(rows.row(0, named=True))

    def get_rows_by_parent(self, parent_dataset: str) -> list[dict]:
        # case-insensitive on both sides, same as get_row(), since parent_dataset
        # is stored mixed-case (e.g. "GTEx_v10", matching the qtl_dataset
        # convention used elsewhere in results/params)
        rows = self._manifest.filter(pl.col("parent_dataset").str.to_lowercase() == parent_dataset.lower())
        return [self._resolve_row_path(row) for row in rows.to_dicts()]

    def get_rows_by_qtl_type(self, qtl_type: str) -> list[dict]:
        rows = self._manifest.filter(pl.col("qtl_type").str.to_lowercase() == qtl_type.lower())
        return [self._resolve_row_path(row) for row in rows.to_dicts()]

    @staticmethod
    def _read_qtl_file(path: str) -> pl.DataFrame:
        suffix = Path(path).suffix.lower()
        if suffix == ".parquet":
            return pl.read_parquet(path)
        if suffix == ".csv":
            return pl.read_csv(path, separator=",")
        if suffix in (".tsv", ".txt"):
            return pl.read_csv(path, separator="\t")
        raise ValueError(f"Unsupported QTL file extension '{suffix}' for {path}")

    @staticmethod
    def normalise_columns(df: pl.DataFrame, manifest_row: dict) -> pl.DataFrame:
        required_cols = ["snp_col", "a1_col", "a2_col", "beta_col", "se_col", "p_col", "chr_col", "pos_col"]
        blank = [col for col in required_cols if not manifest_row.get(col)]
        if blank:
            raise ValueError(
                f"Manifest row for dataset '{manifest_row.get('dataset')}' has blank "
                f"{', '.join(blank)} - fill these in assets/qtl_manifest.csv"
            )

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
            df = self.normalise_columns(self._read_qtl_file(f), manifest_row)

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
