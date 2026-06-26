#!/usr/bin/env python3
import polars as pl 
import argparse
from pathlib import Path
import os 
import subprocess

# ****** NEXT TO-DO'S
# * USE ukbppp_dl pakckage to grab full pQTLs for all proteins
# * Replace chr:start-end to rsIDs in the cloud

# ****** NEXT TO-DO'S
# python ukb_ppp.py --ukb_ppp_dir ../../dat/ukb_ppp --ukb_ppp_snps ../../dat/ukb_ppp/ukbb_ppp_snps
# re-tar dir: tar -cvf ...

# for each protein in dir X
# tar -xvf
# go into that new dir - 
# merge all chr into one single file per protein
# harmonise cols to LDSC/coloc standard
# overwrite previous files - and print whether any SNPs removed
# * tar -xvf A1BG_P04217_OID30771_v1_Inflammation_II.tar
# * rm -rf A1BG_P04217_OID30771_v1_Inflammation_II.tar
# * cd A1BG_P04217_OID30771_v1_Inflammation_II/
# * gunzip *
# data structure
# CHROM GENPOS ID ALLELE0 ALLELE1 A1FREQ INFO N TEST BETA SE CHISQ LOG10P EXTRA
# 1 17641 1:17641:G:A:imp:v1 G A 0.000845899 0.865845 33995 ADD 0.319378 0.12954 6.07854 1.86381 NA
# 1 55057 1:55057:A:G:imp:v1 A G 0.000948149 0.736728 33995 ADD -0.0596563 0.132627 0.202325 0.185186 NA
# 1 101551 1:101551:T:C:imp:v1 T C 0.000762861 0.771693 33995 ADD 0.0159175 0.144558 0.0121246 0.0398525 NA
# 1 108382 1:108382:C:A:imp:v1 C A 0.00103132 0.729595 33995 ADD -0.0525764 0.12788 0.169035 0.166871 NA
# discovery_chr1_A1BG:P04217:OID30771:v1:Inflammation_II - need to make sure we save teh protein ID

def preprocess_pqtls(ukb_ppp_dir: str, ukb_ppp_snps: str):
    ukb_ppp_dir = Path(ukb_ppp_dir)
    ukb_ppp_snps = Path(ukb_ppp_snps)
    proteins = [f.split("_")[0] for f in os.listdir(ukb_ppp_dir) if f.endswith(".tar")]
    if len(proteins) == 2923:
        print("[TRACKING] All good! We found all relevant proteins pertaining to the UKBB-PPP!")
    else:
        print(f"[CONCERN] Yowza! We're missing some stuff here as we have {len(proteins):,} / 2,923")

    # preprocessing cmd
    for tar_file in ukb_ppp_dir.glob("*.tar"):
        protein_dir = ukb_ppp_dir / tar_file.stem
        if protein_dir.exists():
            print(f"[TRACKING] {protein_dir.name} already extracted. Skipping.")
            continue
        print(f"[TRACKING] Extracting {tar_file.name}")

        cmd = f"""
set -euo pipefail
cd "{ukb_ppp_dir}"
tar -xvf "{tar_file.name}"
rm -f "{tar_file.name}"
cd "{protein_dir.name}"
for gz_file in *.gz; do
    [ -e "$gz_file" ] || continue
    gunzip "$gz_file"
done
"""
        subprocess.run(cmd, shell=True, check=True, executable="/bin/bash")

    # now within each dir for each protein
    for protein_dir in ukb_ppp_dir.iterdir():
        if not protein_dir.is_dir():
            continue
        if protein_dir == ukb_ppp_snps:
            continue
        if not any(protein_dir.glob("discovery_chr*")):
            continue

        protein = protein_dir.name.split("_")[0]
        out_file = protein_dir / f"{protein}.parquet"

        if out_file.exists():
            print(f"[TRACKING] {protein} already preprocessed. Skipping.")
            continue

        print(f"[TRACKING] Reading {protein_dir.name}")

        chr_dfs = []
        total_before_map = 0
        total_mapped = 0
        chr_files = sorted(protein_dir.glob("discovery_chr*"))
        for f in chr_files:
            chr_name = f.name.split("_")[1]
            snp_map_file = list(ukb_ppp_snps.glob(f"olink_rsid_map_*_{chr_name}_patched_v2.tsv*"))
            if len(snp_map_file) == 0:
                print(f"[CONCERN] No SNP map found for {chr_name}. Skipping {f.name}")
                continue

            snp_map_file = snp_map_file[0]
            print(f"[TRACKING] Mapping {protein_dir.name} {chr_name} using {snp_map_file.name}")
            df_chr = pl.read_csv(f, separator=" ", has_header=True)
            rsid_map = (pl.read_csv(snp_map_file, separator="\t", has_header=True).select(["ID", "rsid"]))
            n_before = df_chr.height
            df_chr = df_chr.join(rsid_map, on="ID", how="left")
            n_mapped = df_chr.filter(pl.col("rsid").is_not_null()).height
            total_before_map += n_before
            total_mapped += n_mapped
            df_chr = df_chr.with_columns(
                pl.when(pl.col("rsid").is_not_null())
                .then(pl.col("rsid"))
                .otherwise(pl.col("ID"))
                .alias("SNP")
            )

            chr_dfs.append(df_chr)

        if len(chr_dfs) == 0:
            print(f"[CONCERN] No chromosome files mapped for {protein_dir.name}. Skipping protein.")
            continue

        df = pl.concat(chr_dfs)
        pct_mapped = (total_mapped / total_before_map) * 100 if total_before_map > 0 else 0

        print(f"[TRACKING] {protein_dir.name}: {total_mapped:,} / {total_before_map:,} SNPs mapped to rsID ({pct_mapped:.2f}%)")
        print(df.head())

        # rename cols to Bayesian COLOC format
        # KEEP AL EXCEPT "EXTRA" 
        df = (
            df
            .drop(["EXTRA", "rsid"])
            .rename({
                "CHROM": "CHR",
                "GENPOS": "BP",
                "ALLELE1": "A1",
                "ALLELE0": "A2",
                "A1FREQ": "FRQ",
            })
            .with_columns(
                P = 10 ** (-pl.col("LOG10P"))
            )
            .drop("LOG10P")
        )
        
        n_before = df.height
        print(f"TOTAL SNPs (before removing missing IDs) = {n_before:,}")
        df = df.drop_nulls(subset=["SNP"])
        n_after = df.height
        print(f"TOTAL SNPs (after removing missing IDs) = {n_after:,}")
        print(f"SNPs removed = {n_before - n_after:,}")
        df.write_parquet(out_file)
        # pl.read_parquet('A1BG.parquet').head()

def main():
    p = argparse.ArgumentParser()
    p.add_argument("--ukb_ppp_dir", required=True, type=str)
    p.add_argument("--ukb_ppp_snps", required=True, type=str)
    args = p.parse_args()
    preprocess_pqtls(args.ukb_ppp_dir, args.ukb_ppp_snps)

if __name__ == "__main__":
    main()