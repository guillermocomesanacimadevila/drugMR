import argparse
import os
from pathlib import Path

import polars as pl

from drugmr import paths
from drugmr.qtl_manifest import QTLManifest

# grab .parquet files from pQTLs
# add fixed N based on sample size, looked up from assets/qtl_manifest.csv rather than
# hardcoded here - this used to be a decode_n/ukb_ppp_n/wu_csf_n/wingo_brain_n if/elif
# chain duplicated in dashboard/mr_app.py too; single source of truth now
# for each parquet file - create a directory specific to it
# map grab exactly those same SNPs on the .parquet file and map the same SNPs on outcome GWAS
# save locus from GWAS specific to protein X onto the same dir as the protein
# make sure its harmonised to LDSC format
# these will be the ones used for MR and COLOC

# pre-established args from notebook
# * pheno_id
# * ref_bfile
# * pqtl_dataset
# * colnames (pQTL and GWAS)

def define_loci_from_cis_regions(pqtl_dataset: str, pheno_id: str, manifest_path: str = paths.DEFAULT_QTL_MANIFEST_PATH, qc_tsv: str | None = None, out_dir: str | None = None):
    gwas = pl.read_csv(qc_tsv or paths.qc_out(pheno_id), separator="\t")
    pqtl_dataset = pqtl_dataset.lower()

    # qtl_manifest = QTLManifest(manifest_path)
    repo_root = Path(__file__).resolve().parents[1]
    qtl_manifest = QTLManifest(str(Path(manifest_path).resolve()), base_dir=repo_root)

    manifest_row = qtl_manifest.get_row(pqtl_dataset)
    sample_size = int(manifest_row["sample_size"])

    resolved = qtl_manifest.resolve(pqtl_dataset)
    if isinstance(resolved, pl.DataFrame):
        # resolve() collapses a single matched file to a bare DataFrame and drops
        # its filename - re-derive the real protein label here (not pqtl_dataset),
        # since downstream stages need one dir per protein, not one dir total
        import glob
        matched_files = glob.glob(manifest_row["path"])
        protein_label = Path(matched_files[0]).stem
        resolved = {protein_label: resolved}

    for protein, df in resolved.items():
        protein_out_dir = Path(out_dir) / protein if out_dir else Path(f"./dat/cis_regions/{pqtl_dataset}/{protein}")
        os.makedirs(protein_out_dir, exist_ok=True)

        if df.height == 0:
            print(f"[SKIP] {protein}: empty pQTL parquet")
            continue

        df = df.with_columns(pl.lit(sample_size).alias("N"))

        # pos
        chr = df.select(pl.col("CHR").cast(pl.Int64).unique()).item()
        start = df.select(pl.col("BP").min()).item()
        end = df.select(pl.col("BP").max()).item()
        df2 = gwas.filter((pl.col("CHR").cast(pl.Int64) == chr) & (pl.col("BP").is_between(start, end)))

        # remove duplicates
        df = (df.sort("P").unique(subset=["SNP"], keep="first"))
        df2 = (df2.sort("P").unique(subset=["SNP"], keep="first"))

        # match SNPs with pQTL
        pqtl_matched = df.join(df2.select("SNP"), on="SNP", how="inner")
        gwas_matched = df2.join(df.select("SNP"), on="SNP", how="inner")

        # save onto protein_out_dir
        pqtl_matched.write_parquet(protein_out_dir / "pqtl.parquet")
        gwas_matched.write_parquet(protein_out_dir / "gwas.parquet")
        print(
            f"{protein}: "
            f"pQTL={df.height}, GWAS_region={df2.height}, matched={pqtl_matched.height}"
        )

def main():
    p = argparse.ArgumentParser()
    p.add_argument("--pqtl_dataset", required=True)
    p.add_argument("--pheno_id", required=True)
    p.add_argument("--manifest_path", default=paths.DEFAULT_QTL_MANIFEST_PATH)
    p.add_argument("--qc_tsv", default=None)
    p.add_argument("--out_dir", default=None)
    args = p.parse_args()
    define_loci_from_cis_regions(
        pqtl_dataset=args.pqtl_dataset,
        pheno_id=args.pheno_id,
        manifest_path=args.manifest_path,
        qc_tsv=args.qc_tsv,
        out_dir=args.out_dir,
    )

if __name__ == "__main__":
    main()
