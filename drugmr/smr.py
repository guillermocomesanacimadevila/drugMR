import glob
import os
import subprocess
from pathlib import Path

import liftover
import polars as pl

from drugmr import paths
from drugmr.qtl_manifest import QTLManifest
from drugmr.utils import (
    check_n_in_qtl,
    detect_qtl_split,
    extract_gene_coordinates,
    quick_qc,
)


class SMRUtils:


    """
    Take single parquet or a sum of parquets within a dir (/*.parquet)
    and transform to ESD/BESD/EPI format for later intake during QTL-informed workflows
    * Stuff to bare in mind *
    - bulk vs single-cell
    - multiple regions within 1 dataset (e.g. GTEx)
    """

    def __init__(self, manifest_path: str = None, ncbi_ref_path: str = None, base_dir: str = None, liftover_dir: str = None):
        self.manifest_path = manifest_path
        self.ncbi_ref_path = ncbi_ref_path
        self.base_dir = base_dir
        self.liftover_dir = liftover_dir
        self._qtl_manifest = None
        self._gene_positions = None
        self._hg38_to_hg19 = None

    @property
    def qtl_manifest(self):
        # lazy -> run_smr() doesn't need the manifest at all, only transform_to_esd() does
        if self._qtl_manifest is None:
            if self.manifest_path is None:
                raise ValueError("manifest_path not set - required for transform_to_esd()")
            self._qtl_manifest = QTLManifest(manifest_path=self.manifest_path, base_dir=self.base_dir)
        return self._qtl_manifest

    @property
    def gene_positions(self):
        # lazy - same reasoning, only transform_to_esd() needs the NCBI ref
        if self._gene_positions is None:
            if self.ncbi_ref_path is None:
                raise ValueError("ncbi_ref_path not set - required for transform_to_esd()")
            self._gene_positions = pl.read_csv(
                self.ncbi_ref_path,
                separator="\t",
                schema_overrides={"Chromosome": pl.Utf8},  # "X"/"Y"/"MT" break i64 inference
            )
        return self._gene_positions

    @property
    def hg38_to_hg19(self):
        if self._hg38_to_hg19 is None:
            if self.liftover_dir is None:
                raise ValueError("liftover_dir not set - required for GRCh37 QTL conversion")
            chain_dir = Path(self.liftover_dir)
            candidates = [
                chain_dir / "hg38ToHg19.over.chain",
                chain_dir / "hg38ToHg19.over.chain.gz",
            ]
            chain = next((path for path in candidates if path.exists()), None)
            if chain is None:
                raise FileNotFoundError(
                    f"No hg38-to-hg19 chain file found in {chain_dir}; expected one of "
                    f"{', '.join(path.name for path in candidates)}"
                )
            self._hg38_to_hg19 = liftover.ChainFile(str(chain), one_based=True)
        return self._hg38_to_hg19

    def besd_to_sumstats(self, file: Path, out_prefix: Path) -> pl.DataFrame: # != in wrappers -> convert to parquet and save as parquet

        """ Quick note for myself: we can hardcode A1 and A2 here given that file == from smr_out """

        file = Path(file)
        out_prefix = Path(out_prefix)

        cmd = [
            "smr",
            "--beqtl-summary", str(file),
            "--query", "1",
            "--out", str(out_prefix)
        ]
        try:
            subprocess.run(cmd, check=True, capture_output=True, text=True)
        except subprocess.CalledProcessError as error:
            raise RuntimeError(
                f"smr --query failed for {file} (exit {error.returncode}):\n"
                f"stdout: {error.stdout}\nstderr: {error.stderr}"
            ) from error

        # smr's --query output never has an N column at all (real production
        # BESD and self-built BESD alike) - N is filled in downstream from the
        # manifest's sample_size via check_n_in_qtl(), not detected here
        query_file = Path(f"{out_prefix}.txt")
        df = (
            pl.read_csv(query_file, separator="\t")
            .with_columns(
                pl.col("Chr").cast(pl.Int64),
                pl.col("BP").cast(pl.Int64),
                pl.col("Probe_Chr").cast(pl.Int64),
                pl.col("Probe_bp").cast(pl.Int64),
                pl.col("Freq").cast(pl.Float64),
                pl.col("b").cast(pl.Float64),
                pl.col("SE").cast(pl.Float64),
                pl.col("p").cast(pl.Float64),
            )
        )
        df = quick_qc(sumstats=df, a1_col="A1", a2_col="A2")
        query_file.unlink()
        return df

    def ensure_sumstats_parquet(self, besd_prefixes: dict, out_parquet: Path, sample_size: int = None) -> Path:
        out_parquet = Path(out_parquet)
        if out_parquet.exists():
            print(f"[TRACKING] '{out_parquet}' already exists, skipping BESD query")
            return out_parquet

        dfs = [self.besd_to_sumstats(prefix, prefix) for prefix in besd_prefixes.values()]
        combined = check_n_in_qtl(pl.concat(dfs), sample_size).sort(["Chr", "BP"])
        out_parquet.parent.mkdir(parents=True, exist_ok=True)
        combined.write_parquet(out_parquet)
        print(f"[TRACKING] Wrote {combined.height} rows to {out_parquet}")
        return out_parquet

    @staticmethod
    def _group_besd_prefixes(besd_prefixes: dict) -> dict:
        # groups by each BESD prefix's own parent directory - NOT by "/" in the
        # label string, which is ambiguous (a multi-tissue dataset covering only
        # 1 chromosome, e.g. a locus-slice toy fixture, never chr-splits, so its
        # per-tissue labels come back flat, indistinguishable by shape alone from
        # a genuinely single-file dataset's chr-split labels). The parent
        # directory is robust across both layouts actually produced by this
        # class: detect_qtl_split's real per-tissue/per-cell subdirectories
        # (GTEx_v10, SingleBrain) or single flat directory (MetaBrain), and
        # transform_to_esd's freshly-built one-directory-per-label layout.
        groups = {}
        for label, prefix in besd_prefixes.items():
            group = str(prefix.parent)
            groups.setdefault(group, {})[label] = prefix
        return groups

    def resolve_bulk_qtl_file(self, qtl_dataset: str, cell_type: str, esd_dir: Path = None) -> Path | None:
        # single-file datasets (MetaBrain) are registered under qtl_dataset alone;
        # multi-region datasets (GTEx_v10) are registered 1 manifest row per
        # region/tissue, named after the exact cell_type value results already
        # carry (e.g. "gtex_brain_amygdala_v10") - try the plain dataset name
        # first, fall back to cell_type. No hardcoded dataset names either way.
        try:
            manifest_row = self.qtl_manifest.get_row(qtl_dataset.lower())
        except ValueError:
            manifest_row = self.qtl_manifest.get_row(cell_type.lower())

        matched_files = sorted(glob.glob(manifest_row["path"]))
        if matched_files:
            if len(matched_files) == 1:
                return Path(matched_files[0])  # 1 file registered for this manifest row
            for f in matched_files:
                if Path(f).stem in cell_type:  # manifest row still covers >1 file
                    return Path(f)
            return None

        # nothing pre-built yet - fall back to building it from BESD
        if esd_dir is None:
            return None

        besd_prefixes = self.ensure_besd(manifest_row["dataset"], esd_dir)
        if not besd_prefixes:
            return None

        groups = self._group_besd_prefixes(besd_prefixes)
        if len(groups) == 1:
            _, group_prefixes = next(iter(groups.items()))
            out_parquet = Path(manifest_row["path"])  # literal path, no wildcard to fill in
        else:
            group_dir = next((g for g in groups if Path(g).name in cell_type or cell_type in Path(g).name), None)
            if group_dir is None:
                return None
            group_prefixes = groups[group_dir]
            group_name = Path(group_dir).name
            # substitute the matched tissue/cell name into the manifest's glob
            # pattern (e.g. "dat/bulk-eQTL/GTEx_v10/*/*.parquet" -> ".../Brain_Amygdala/Brain_Amygdala.parquet")
            parts = [p.replace("*", group_name) if "*" in p else p for p in Path(manifest_row["path"]).parts]
            out_parquet = Path(*parts)

        return self.ensure_sumstats_parquet(group_prefixes, out_parquet, sample_size=manifest_row.get("sample_size"))

    def resolve_sc_qtl_file(self, qtl_dataset: str, cell_type: str) -> Path | None:
        # manifest-driven (case-insensitive on qtl_dataset, matching the
        # convention everywhere else) instead of hardcoding "dat/sc-eQTL/{qtl_dataset}/"
        # - that hardcoded form only worked when the caller's casing happened to
        # match the real directory name (e.g. "SingleBrain") exactly
        manifest_row = self.qtl_manifest.get_row(qtl_dataset.lower())
        matched_files = sorted(glob.glob(manifest_row["path"]))
        for f in matched_files:
            if Path(f).stem.lower() == cell_type.lower():
                return Path(f)
        return None

    def load_qtl_rows(self, data_type: str, qtl_dataset: str, cell_type: str, base_gene_id: str, esd_dir: Path = None) -> pl.DataFrame | None:
        """
        1 dataset/cell-type/tissue's QTL rows for 1 gene, standardised onto
        SNP/A1(effect allele)/A2/BETA/SE/P/FRQ/N - single shared implementation for
        bin/hyprcoloc_targets.py and bin/pwcoco_qtl_wrapper.py, which each
        previously re-implemented this separately. Pass esd_dir to allow building
        the parquet from BESD on demand if it doesn't exist yet (rare - all 3
        registered QTL datasets already have it pre-built as of 2026-09).
        """
        if data_type == "single_cell":
            qtl_file = self.resolve_sc_qtl_file(qtl_dataset, cell_type)
            if qtl_file is None:
                if esd_dir is None:
                    print(f"[CONCERN] Missing {qtl_dataset} QTL file for {cell_type}")
                    return None

                manifest_row = self.qtl_manifest.get_row(qtl_dataset.lower())
                besd_prefixes = self.ensure_besd(qtl_dataset.lower(), esd_dir)
                group_prefixes = {label: prefix for label, prefix in besd_prefixes.items() if cell_type in label}
                if not group_prefixes:
                    print(f"[CONCERN] No BESD prefix found for {qtl_dataset}/{cell_type}")
                    return None
                # target path for the freshly-built parquet, derived from the
                # manifest's own registered directory rather than a hardcoded
                # dat/sc-eQTL/{qtl_dataset}/ guess
                qtl_file = Path(manifest_row["path"]).parent / f"{cell_type}.parquet"
                self.ensure_sumstats_parquet(group_prefixes, qtl_file, sample_size=manifest_row.get("sample_size"))

            # single-cell QTL files carry ref/alt as A1/A2 and the actual effect allele as EA -
            # re-point A1/A2 so A1 is always the effect allele BETA belongs to
            return (
                pl.scan_parquet(qtl_file)
                .filter(pl.col("GENE").str.split(".").list.first() == base_gene_id)
                .select(["SNP", "A1", "A2", "EA", "BETA", "SE", "P", "FRQ", "N"])
                .with_columns(
                    pl.col("EA").alias("qtl_a1"),
                    pl.when(pl.col("EA") == pl.col("A2")).then(pl.col("A1")).otherwise(pl.col("A2")).alias("qtl_a2")
                )
                .select(["SNP", pl.col("qtl_a1").alias("A1"), pl.col("qtl_a2").alias("A2"), "BETA", "SE", "P", "FRQ", "N"])
                .sort("P")
                .unique(subset="SNP", keep="first")
                .collect()
            )

        if data_type == "bulk":
            qtl_file = self.resolve_bulk_qtl_file(qtl_dataset, cell_type, esd_dir=esd_dir)
            if qtl_file is None or not qtl_file.exists():
                print(f"[CONCERN] Missing {qtl_dataset} bulk QTL file for {cell_type}: {qtl_file}")
                return None

            # bulk QTL parquets come straight from an SMR besd/esi/epi query, so A1 is
            # already the effect allele b belongs to (SMR's own convention) - no
            # re-pointing needed, just rename onto the pipeline's BETA/P/FRQ convention
            return (
                pl.scan_parquet(qtl_file)
                .filter(pl.col("Probe").str.split(".").list.first() == base_gene_id)
                .select(
                    "SNP", "A1", "A2",
                    pl.col("b").alias("BETA"), "SE", pl.col("p").alias("P"), pl.col("Freq").alias("FRQ"), "N",
                )
                .sort("P")
                .unique(subset="SNP", keep="first")
                .collect()
            )

        print(f"[CONCERN] Unrecognised data_type '{data_type}' for {cell_type} - skipping")
        return None

    def transform_to_esd(self, dataset: str, esd_dir: Path) -> dict:
        manifest_row = self.qtl_manifest.get_row(dataset)
        key_col = manifest_row.get("key_col") or None
        build = "hg38" if manifest_row["build"] == "GRCh38" else "hg19"

        all_flist_rows = {}
        for label, df in self._iter_qtl_partitions(manifest_row):

            gene_col = key_col
            if gene_col is None:
                gene_col = "_gene"
                df = df.with_columns(pl.lit(label.split("_")[0]).alias(gene_col))

            label_dir = esd_dir / label
            label_dir.mkdir(parents=True, exist_ok=True)
            flist_rows = []
            for gene_df in df.partition_by(gene_col, as_dict=False):
                gene = gene_df[gene_col][0]

                try:
                    coords = extract_gene_coordinates(
                        gene,
                        self.gene_positions,
                        genome_build=build,
                        converter=self.hg38_to_hg19 if build == "hg19" else None,
                    )
                except ValueError:
                    continue  # gene not in the NCBI ref - skip it, not a hard fail

                esd_df = (
                    gene_df
                    .select(["CHR", "SNP", "BP", "A1", "A2", "FRQ", "BETA", "SE", "P"])
                    .rename({"CHR": "Chr", "BP": "Bp", "FRQ": "Freq", "BETA": "Beta", "SE": "se", "P": "p"})
                    .drop_nulls()
                    .filter(
                        (pl.col("SNP") != "") &
                        (pl.col("A1") != pl.col("A2")) &
                        pl.col("Freq").is_finite() & pl.col("Freq").is_between(0, 1, closed="none") &
                        pl.col("Beta").is_finite() &
                        pl.col("se").is_finite() & (pl.col("se") > 0) &
                        pl.col("p").is_finite() & pl.col("p").is_between(0, 1, closed="right")
                    )
                    .unique(subset=["SNP", "Bp", "A1", "A2"])
                    .sort(["Chr", "Bp", "SNP"])
                )

                if esd_df.height == 0:
                    continue

                safe_gene = str(gene).replace("/", "_").replace(":", "_").replace(" ", "_")
                esd_path = label_dir / f"{safe_gene}.esd"
                esd_df.write_csv(esd_path, separator="\t")

                # ProbeID is the Ensembl ID (matching real production BESD's
                # native convention, which load_qtl_rows()'s Probe-based
                # filtering assumes) when the NCBI ref has one for this gene;
                # falls back to the symbol otherwise rather than hard-failing,
                # since ~30% of NCBI ref genes have no GENCODE symbol match
                probe_id = coords["ENSEMBL_ID"][0] or gene
                flist_rows.append({
                    "Chr": coords["CHR"][0],
                    "ProbeID": probe_id,
                    "GeneticDistance": 0,
                    "ProbeBp": coords["START"][0],
                    "Gene": gene,
                    "Orientation": coords["ORIENTATION"][0],
                    "PathOfEsd": str(esd_path.resolve()),
                })

            all_flist_rows[label] = flist_rows
        return all_flist_rows

    def _iter_qtl_partitions(self, manifest_row: dict):
        """Yield normalised QTL frames without expanding a whole dataset in RAM.

        Large genome-wide Parquet inputs such as MetaBrain are compressed to a
        few GB on disk but expand far beyond a normal SLURM allocation when read
        eagerly.  Discover the chromosome values cheaply, then rely on Parquet
        predicate pushdown to collect one chromosome at a time.
        """
        matched_files = sorted(glob.glob(manifest_row["path"]))
        if not matched_files:
            raise FileNotFoundError(f"No files matched path: {manifest_row['path']}")

        for file_name in matched_files:
            path = Path(file_name)
            suffix = path.suffix.lower()
            if suffix == ".parquet":
                frame = pl.scan_parquet(path)
            elif suffix == ".csv":
                frame = pl.scan_csv(path, separator=",")
            elif suffix in (".tsv", ".txt"):
                frame = pl.scan_csv(path, separator="\t")
            else:
                raise ValueError(f"Unsupported QTL file extension '{suffix}' for {path}")

            frame = self.qtl_manifest.normalise_columns(frame, manifest_row)
            label = path.stem
            key_col = manifest_row.get("key_col") or None
            if key_col is None:
                key_col = "_gene"
                frame = frame.with_columns(pl.lit(label.split("_")[0]).alias(key_col))

            chromosomes = (
                frame.select("CHR")
                .drop_nulls()
                .unique()
                .sort("CHR")
                .collect(engine="streaming")
                .get_column("CHR")
                .to_list()
            )
            split = len(chromosomes) > 1
            for chromosome in chromosomes:
                partition_label = f"{label}/chr{chromosome}" if split else label
                print(f"[TRACKING] Loading {manifest_row['dataset']} partition {partition_label}", flush=True)
                yield partition_label, frame.filter(pl.col("CHR") == chromosome).collect(engine="streaming")

    def transform_to_flist(self, all_flist_rows: dict, esd_dir: Path) -> dict:

        flist_paths = {}
        for label, rows in all_flist_rows.items():
            if not rows:
                raise ValueError(f"No ESD files/flist rows for '{label}'")

            flist = pl.DataFrame(rows).sort(["Chr", "ProbeBp"])
            label_dir = esd_dir / label
            flist_path = label_dir / f"{Path(label).name}.flist"
            flist.write_csv(flist_path, separator="\t")
            flist_paths[label] = flist_path

        return flist_paths

    def transform_to_besd(self, flist_paths: dict) -> dict:

        besd_prefixes = {}
        for label, flist_path in flist_paths.items():
            out_prefix = flist_path.parent / Path(label).name

            cmd = [
                "smr",
                "--eqtl-flist", str(flist_path),
                "--make-besd",
                "--out", str(out_prefix),
            ]
            subprocess.run(cmd, check=True)
            besd_prefixes[label] = out_prefix

        return besd_prefixes

    def ensure_besd(self, dataset: str, esd_dir: Path) -> dict:
        manifest_row = self.qtl_manifest.get_row(dataset)

        # detect_qtl_split handles both cases: a single besd/esi/epi prefix,
        # or several already split per-chromosome in the same directory (e.g.
        # MetaBrain's real downloaded bundle) - empty dict means neither exists
        already_ready = detect_qtl_split(Path(manifest_row["path"]))
        if already_ready:
            print(f"[TRACKING] '{dataset}' is already BESD-ready ({len(already_ready)} file(s)), skipping ETL")
            return already_ready

        all_flist_rows = self.transform_to_esd(dataset, esd_dir)
        flist_paths = self.transform_to_flist(all_flist_rows, esd_dir)
        return self.transform_to_besd(flist_paths)

    def run_smr(
            self,
            pheno_id: str,
            sumstats: str,
            ref_bfile: str,
            beqtl_summary: str,
            qtl_dataset: str,
            p_qtl_smr: float,
            p_qtl_heidi: float,
            thread_num: int,
            maf: float,
            out_dir: str = "synthesis"
    ):
        ref_bfile = Path(ref_bfile)
        sumstats = Path(sumstats)
        beqtl_summary = Path(beqtl_summary)

        # qtl_dataset can be stuff like SingleBrain/Ast
        # use full path for directory but only cell name for output prefix
        qtl_dataset = Path(qtl_dataset)
        # SMR(GWAS x QTL) result for a given (pheno_id, qtl_dataset) never depends on
        # pqtl_dataset, so this defaults to the shared synthesis/ tree rather than a
        # per-run out_dir - every pqtl_dataset run reuses the same computation instead
        # of re-running the smr binary from scratch
        raw_out_dir = paths.smr_raw_dir(qtl_dataset, pheno_id, out_dir)
        os.makedirs(raw_out_dir, exist_ok=True)
        out_file = paths.smr_raw_prefix(qtl_dataset, pheno_id, out_dir)
        print(f"[TRACKING] Running SMR on {pheno_id} using {qtl_dataset}")

        cmd_smr = [
            "smr",
            "--bfile", str(ref_bfile),
            "--gwas-summary", str(sumstats),
            "--beqtl-summary", str(beqtl_summary),
            "--maf", str(maf),
            "--peqtl-smr", str(p_qtl_smr),
            "--peqtl-heidi", str(p_qtl_heidi),
            "--thread-num", str(thread_num),
            "--out", str(out_file),
        ]

        subprocess.run(cmd_smr, check=True)
