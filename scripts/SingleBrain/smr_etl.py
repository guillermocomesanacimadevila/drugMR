#!/usr/bin/env python3
import argparse
import polars as pl
from pathlib import Path
import subprocess
import os
import shutil

# load parquet
# change format to SMR - i.e. appropiate colnames
# for each cell_type
# convert into smr format
# save onto .dat/sc-eQTL/SingleBrain/SMR_ready/...
# ALREADY IN GRH38 - no need for liftover
# write into work/temp/SMR/cell_type.ma (all of this within for loop)
# convert into smr format and save onto out_dir -> remove .ma from work/temp/...
# carry on with the next cell type


def load_gene_positions_from_gtf(gtf_path: str) -> pl.DataFrame:
    gtf = pl.read_csv(
        gtf_path,
        separator="\t",
        has_header=False,
        comment_prefix="#",
        new_columns=[
            "chrom", "source", "feature", "start", "end",
            "score", "strand", "frame", "attributes"
        ],
        truncate_ragged_lines=True,
    )

    gene_positions = (
        gtf
        .filter(pl.col("feature") == "gene")
        .with_columns([
            pl.col("attributes").str.extract(r'gene_id "([^"]+)"', 1).alias("ProbeID"),
            pl.col("attributes").str.extract(r'gene_name "([^"]+)"', 1).alias("Gene"),
        ])
        .with_columns(
            pl.col("ProbeID").str.replace(r"\.\d+$", "").alias("gene_id_base")
        )
        .select([
            "gene_id_base",
            "ProbeID",
            "Gene",
            pl.col("chrom").str.replace("^chr", "").alias("ProbeChr"),
            pl.col("start").cast(pl.Int64).alias("ProbeBp"),
            pl.col("strand").alias("Orientation"),
        ])
        .unique(subset=["gene_id_base"])
    )

    return gene_positions


def write_smr_esd_and_flist(smr_df: pl.DataFrame, out_dir: Path, cell_type: str) -> Path:
    cell_out_dir = out_dir / cell_type
    esd_dir = cell_out_dir / "esd"

    cell_out_dir.mkdir(parents=True, exist_ok=True)
    esd_dir.mkdir(parents=True, exist_ok=True)

    flist_rows = []

    for probe_df in smr_df.partition_by("Probe", as_dict=False):
        probe = probe_df["Probe"][0]
        probe_id = probe_df["ProbeID"][0]
        probe_chr = probe_df["ProbeChr"][0]
        probe_bp = probe_df["ProbeBp"][0]
        gene = probe_df["Gene"][0]
        orientation = probe_df["Orientation"][0]

        safe_probe = (
            str(probe)
            .replace("/", "_")
            .replace(":", "_")
            .replace(" ", "_")
        )

        esd_path = esd_dir / f"{safe_probe}.esd"

        esd_df = (
            probe_df
            .select(["Chr", "SNP", "Bp", "A1", "A2", "Freq", "Beta", "se", "p"])
            .drop_nulls()
            .filter(
                (pl.col("SNP") != "") &
                (pl.col("A1") != "") &
                (pl.col("A2") != "") &
                (pl.col("A1") != pl.col("A2")) &
                pl.col("Freq").is_finite() &
                pl.col("Beta").is_finite() &
                pl.col("se").is_finite() &
                pl.col("p").is_finite() &
                (pl.col("Freq") > 0) &
                (pl.col("Freq") < 1) &
                (pl.col("se") > 0) &
                (pl.col("p") > 0) &
                (pl.col("p") <= 1)
            )
            .unique(subset=["SNP", "Bp", "A1", "A2"])
            .sort(["Chr", "Bp", "SNP"])
        )

        if esd_df.height == 0:
            continue

        esd_df.write_csv(esd_path, separator="\t")

        flist_rows.append({
            "Chr": probe_chr,
            "ProbeID": probe_id,
            "GeneticDistance": 0,
            "ProbeBp": probe_bp,
            "Gene": gene,
            "Orientation": orientation if orientation is not None else "NA",
            "PathOfEsd": str(esd_path.resolve()),
        })

    if not flist_rows:
        raise ValueError(f"No valid ESD files written for {cell_type}")

    flist = pl.DataFrame(flist_rows).sort(["Chr", "ProbeBp", "Gene"])
    flist_path = cell_out_dir / f"{cell_type}.flist"
    flist.write_csv(flist_path, separator="\t")

    print(f"[TRACKING] Written {len(flist_rows)} ESD files for {cell_type}", flush=True)
    print(f"[TRACKING] flist: {flist_path}", flush=True)

    return flist_path


def run_smr_make_besd(flist_path: Path, out_dir: Path, cell_type: str, smr_bin: str):
    out_prefix = out_dir / cell_type / cell_type

    cmd = [
        smr_bin,
        "--eqtl-flist", str(flist_path),
        "--make-besd",
        "--out", str(out_prefix),
    ]

    print(f"[TRACKING] Running SMR for {cell_type}", flush=True)
    print("[TRACKING] " + " ".join(cmd), flush=True)

    subprocess.run(cmd, check=True)

    print(f"[TRACKING] SMR BESD done for {cell_type}", flush=True)


def preprocess_sc_eqtls_for_smr(
    input_dir: str,
    dataset: str,
    gtf_path: str,
    out_dir: str,
    smr_bin: str,
    run_smr: bool,
    clean_temp: bool,
):
    input_dir = Path(input_dir)
    out_dir = Path(out_dir)

    gene_positions = load_gene_positions_from_gtf(gtf_path)

    # input dir is ./dat/sc-eQTL/SingleBrain/...
    temp_out_dir = f"./work/temp/SMR/{dataset}" # temp out_dir
    temp_out_dir = Path(temp_out_dir) # save temp .ma file here then .unlin()
    os.makedirs(temp_out_dir, exist_ok=True)

    out_dir.mkdir(parents=True, exist_ok=True)

    for file in sorted(input_dir.glob("*.parquet")):
        cell_type = file.stem
        print(f"[TRACKING] Processing {cell_type}", flush=True)

        df = pl.read_parquet(file)

        print(df.columns, flush=True)
        print(df.head(), flush=True)

        # SNP    A1  A2  freq    b   se  p   n
        # in SingleBrain EA is almost always A2
        # so we have to rename EA - A1 because A1 == ref allele
        # SNP    A1  A2  Freq    Beta   se  p
        # SMR A1 must be the effect/coded allele
        smr_df = df.select([
            pl.col("CHR").cast(pl.Utf8).alias("Chr"),
            pl.col("SNP").cast(pl.Utf8).alias("SNP"),
            pl.col("BP").cast(pl.Int64).alias("Bp"),

            # EA is the effect allele, so this becomes SMR A1
            pl.col("EA").cast(pl.Utf8).alias("A1"),

            # A2 is whichever allele is not EA
            pl.when(pl.col("EA") == pl.col("A1"))
              .then(pl.col("A2"))
              .when(pl.col("EA") == pl.col("A2"))
              .then(pl.col("A1"))
              .otherwise(None)
              .cast(pl.Utf8)
              .alias("A2"),

            pl.col("FRQ").cast(pl.Float64).alias("Freq"),
            pl.col("BETA").cast(pl.Float64).alias("Beta"),
            pl.col("SE").cast(pl.Float64).alias("se"),
            pl.col("P").cast(pl.Float64).alias("p"),
            pl.col("GENE").cast(pl.Utf8).alias("Probe"),
        ])

        # if any weird FRQ was accidentally percentage-scaled, force it into 0-1
        smr_df = (
            smr_df
            .with_columns(
                pl.when(pl.col("Freq") > 1)
                .then(pl.col("Freq") / 100)
                .otherwise(pl.col("Freq"))
                .alias("Freq")
            )
        )

        print(smr_df.select([
            pl.col("A2").null_count().alias("missing_A2"),
            (pl.col("A1") == pl.col("A2")).sum().alias("same_A1_A2"),
            pl.col("Freq").min().alias("min_Freq"),
            pl.col("Freq").max().alias("max_Freq"),
            pl.len().alias("total")
        ]), flush=True)

        smr_df = (
            smr_df
            .with_columns(
                pl.col("Probe").str.replace(r"\.\d+$", "").alias("gene_id_base")
            )
            .join(gene_positions, on="gene_id_base", how="left")
        )

        print(smr_df.select([
            pl.col("ProbeID").null_count().alias("missing_ProbeID"),
            pl.col("ProbeChr").null_count().alias("missing_ProbeChr"),
            pl.col("ProbeBp").null_count().alias("missing_ProbeBp"),
            pl.len().alias("total")
        ]), flush=True)

        # SMR file processing
        smr_df = (
            smr_df
            .drop_nulls(["A1", "A2", "Freq", "Beta", "se", "p", "ProbeID", "ProbeChr", "ProbeBp"])
            .filter(
                (pl.col("A1") != pl.col("A2")) &
                pl.col("Freq").is_finite() &
                pl.col("Beta").is_finite() &
                pl.col("se").is_finite() &
                pl.col("p").is_finite() &
                (pl.col("Freq") > 0) &
                (pl.col("Freq") < 1) &
                (pl.col("se") > 0) &
                (pl.col("p") > 0) &
                (pl.col("p") <= 1)
            )
            .drop("gene_id_base")
        )

        print(smr_df.select([
            pl.col("Freq").min().alias("final_min_Freq"),
            pl.col("Freq").max().alias("final_max_Freq"),
            pl.len().alias("final_rows")
        ]), flush=True)

        flist_path = write_smr_esd_and_flist(
            smr_df=smr_df,
            out_dir=temp_out_dir,
            cell_type=cell_type,
        )

        if run_smr:
            run_smr_make_besd(
                flist_path=flist_path,
                out_dir=temp_out_dir,
                cell_type=cell_type,
                smr_bin=smr_bin,
            )

            final_cell_dir = out_dir / cell_type
            final_cell_dir.mkdir(parents=True, exist_ok=True)

            for ext in [".besd", ".esi", ".epi", ".summary", ".log"]:
                src = temp_out_dir / cell_type / f"{cell_type}{ext}"
                if src.exists():
                    shutil.copy2(src, final_cell_dir / src.name)

            shutil.copy2(flist_path, final_cell_dir / flist_path.name)

            print(f"[TRACKING] Saved SMR-ready files to {final_cell_dir}", flush=True)

    if clean_temp:
        shutil.rmtree(temp_out_dir, ignore_errors=True)
        print(f"[TRACKING] Removed temp dir: {temp_out_dir}", flush=True)

    return


# python scripts/SingleBrain/smr_etl.py --run-smr


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Convert processed SingleBrain parquet eQTLs into SMR-ready BESD files.")
    parser.add_argument("--input-dir", default="./dat/sc-eQTL/SingleBrain", help="Input directory containing processed SingleBrain parquet files.")
    parser.add_argument("--dataset", default="SingleBrain", help="Dataset name for temp SMR output.")
    parser.add_argument("--gtf", default="./dat/ref/GENCODE/gencode.v50.annotation.gtf", help="GENCODE GRCh38 GTF file.")
    parser.add_argument("--out-dir", default="./dat/sc-eQTL/SingleBrain/SMR_ready", help="Final output directory for SMR-ready files.")
    parser.add_argument("--smr-bin", default="smr", help="Path to SMR binary.")
    parser.add_argument("--run-smr", action="store_true", help="Actually run smr --make-besd after writing ESD/flist files.")
    parser.add_argument("--clean-temp", action="store_true", help="Remove temp SMR directory after finishing.")

    args = parser.parse_args()

    preprocess_sc_eqtls_for_smr(
        input_dir=args.input_dir,
        dataset=args.dataset,
        gtf_path=args.gtf,
        out_dir=args.out_dir,
        smr_bin=args.smr_bin,
        run_smr=args.run_smr,
        clean_temp=args.clean_temp,
    )