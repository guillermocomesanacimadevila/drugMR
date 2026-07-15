#!/usr/bin/env python3
import argparse
import polars as pl
from pathlib import Path
import subprocess
from drugmr import SMR
import pandas as pd
import os 
from statsmodels.stats.multitest import fdrcorrection

# ------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------
# MAIN TO DO'S
# -> SLAP FUNCTION 3 (MAYBE 1 ONTO A DIFFERNT SCRIPT -> MAY CRASH SMR IF != RESULTS)
# -> CONSEQUENTLY UPDATE drugmr/local.py and drugmr/hpc.py
# ------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------

# need probs another function with mediators 
# -> load cis-MR results for pQTL dataset X  -> same for coloc -> check whether gene which passes coloc thresh and MR estimate FDR
# -> check on all cells - save the SMR output for that in 1 or > 1 cells onto results/...

# ADD TO DOCKER IMAGE 
# - STATSMODELS 
# - SMR PACKAGE (as part of ref/)

# COLOC and cis-MR filter for promising targets
def extract_promising_targets(pqtl_dataset: str, pheno_id: str):
    # extract stuff from here
    base_dir = "./results" # ukb_ppp_AD_all_MR.tsv
    base_dir = Path(base_dir) # coloc/ukb_ppp/ukb_ppp_AD_all_coloc.tsv
    cis_mr_res = base_dir / f"cis-MR/{pqtl_dataset}_{pheno_id}_all_MR.tsv"
    cis_mr_df = pl.read_csv(cis_mr_res, separator="\t")
    coloc_res = base_dir / f"coloc/{pqtl_dataset}/{pqtl_dataset}_{pheno_id}_all_coloc.tsv"
    coloc_df = pl.read_csv(coloc_res, separator="\t")
    
    # base parameters
    # wald_fdr = 0.05 
    # ivw_fdr = 0.05
    # cochran_q = 0.05
    # coloc_thresh = 0.75 #### subject to change  

    wald_hits = []
    ivw_hits = []
    coloc_hits = []

    for row in cis_mr_df.iter_rows(named=True):
        # separate where n_instruments == 1 or > 1
        n_instruments = row["n_instruments"]
        wald_fdr = row["Wald_FDR_q"]
        ivw_fdr = row["IVW_FDR_q"]
        cochran_q = row["Q_pval"]
        if n_instruments == 1 and wald_fdr is not None and wald_fdr < 0.05:
            wald_hits.append(row["protein"])
        elif (n_instruments > 1 and ivw_fdr is not None and ivw_fdr < 0.05 and cochran_q is not None and cochran_q > 0.05):
            ivw_hits.append(row["protein"])

    for row in coloc_df.iter_rows(named=True):
        pp4 = row["PP.H4.abf"]
        if pp4 is not None and pp4 > 0.75:
            coloc_hits.append(row["protein_id"])

    # compile final hits
    mr_hits = wald_hits + ivw_hits
    final_hits = [i for i in mr_hits if i in coloc_hits]
    print(f"[TRACKING] {len(wald_hits)} Wald ratio hits found")
    print(f"[TRACKING] {len(ivw_hits)} IVW hits found")
    print(f"[TRACKING] {len(coloc_hits)} coloc hits found")
    print(f"[TRACKING] {len(final_hits)} final promising targets found")
    return final_hits

    
def run_single_cell_smr(pqtl_dataset: str, eqtl_dataset: str, pheno_id: str, sumstats: str, ref_bfile: str, maf: float):
    ref_bfile = Path(ref_bfile)
    # temp dir to store .ma file per pheno
    temp_dir = "./work/SMR/"
    temp_dir = Path(temp_dir)
    os.makedirs(temp_dir, exist_ok=True)
    eqtl_temp = eqtl_dataset.lower()

    if eqtl_temp == "singlebrain":
        # Store sumstats (temp) within work/ as a .ma file 
        # which then delete after all cell types are ran -> just a temp file 
        # Store sumstats (temp) within work/ as a .ma file
        # which then delete after all cell types are ran -> just a temp file
        df = pl.read_csv(sumstats, separator="\t")
        n_before = df.height
        # rename and save as .ma
        # SMR GWAS format == SNP A1 A2 freq b se p N
        df = (
            df
            .select([
                pl.col("SNP"),
                pl.col("A1"),
                pl.col("A2"),
                pl.col("FRQ").alias("freq"),
                pl.col("BETA").alias("b"),
                pl.col("SE").alias("se"),
                pl.col("P").alias("p"),
                pl.col("N")
            ])
            .filter(
                pl.col("SNP").is_not_null(),
                pl.col("A1").is_not_null(),
                pl.col("A2").is_not_null(),
                pl.col("freq").is_not_null(),
                pl.col("b").is_not_null(),
                pl.col("se").is_not_null(),
                pl.col("p").is_not_null(),
                pl.col("N").is_not_null(),
                ~pl.col("SNP").str.contains(","),
                ~pl.col("SNP").str.contains(";"),
                ~pl.col("SNP").str.contains(" ")
            )
            .with_columns(
                pl.col("N").round(0).cast(pl.Int64)
            )
        )

        print(f"[TRACKING] Removed {n_before - df.height} invalid / incomplete GWAS rows for SMR")
        print(f"[TRACKING] {df.height} GWAS variants retained for SMR")
        temp_sumstats = temp_dir / f"{pheno_id}.ma"
        df.write_csv(temp_sumstats, separator="\t")
        # cell_types = ["Ast", "Ext", "MG", "OD", "OPC", "End", "IN"]
        cell_types = ["MG"]
        eqtls = "./dat/sc-eQTL/SingleBrain/SMR_ready"
        eqtls = Path(eqtls)
        for cell in cell_types:
            cell_dir = eqtls / cell
            besd_file = cell_dir / f"{cell}.besd"
            esi_file = cell_dir / f"{cell}.esi"
            epi_file = cell_dir / f"{cell}.epi"

            if not cell_dir.exists():
                print(f"[CONCERN] Cell type directory {cell_dir} not found")
                continue

            if not besd_file.exists():
                print(f"[CONCERN] {besd_file} not found")
                continue

            if not esi_file.exists():
                print(f"[CONCERN] {esi_file} not found")
                continue

            if not epi_file.exists():
                print(f"[CONCERN] {epi_file} not found")
                continue

            print(f"[TRACKING] Cell type {cell} found!")

            # continue
            # use prefix without .besd / .esi / .epi for SMR

            #####
            #####
            #####
            beqtl_summary = cell_dir / cell
            #####
            #####
            #####

            SMR(
                pheno_id=pheno_id,
                sumstats=temp_sumstats,
                ref_bfile=ref_bfile,
                beqtl_summary=beqtl_summary,
                eqtl_dataset=f"{eqtl_dataset}/{cell}",
                peqtl_smr=5.0e-8, #### change to default one 
                peqtl_heidi=1.57e-3, ###### change to real default 
                thread_num=8,
                maf=maf
            )

            # load SMR results
            # saving into out_dir 1 results file per cell type for trait X
            # results/SMR/SingleBrain/{cell}/{pheno_id}/...
            smr_res = f"./results/SMR/{eqtl_dataset}/{cell}/{pheno_id}"
            smr_res = Path(smr_res)
            for f in smr_res.glob("*.smr"):
                if pheno_id in f.name:
                    smr_df = pl.read_csv(f, separator="\t")
                    # now FDR correct p_SMR and add q_SMR col
                    if "p_SMR" in smr_df.columns:
                        valid_p = (
                            smr_df
                            .filter(pl.col("p_SMR").is_not_null())
                            .get_column("p_SMR")
                            .to_numpy()
                        )
                        if len(valid_p) > 0:
                            _, q_values = fdrcorrection(
                                valid_p,
                                alpha=0.05,
                                method="indep"
                            )
                            q_df = (
                                smr_df
                                .filter(pl.col("p_SMR").is_not_null())
                                .select(
                                    pl.int_range(pl.len()).alias("_row_id")
                                )
                                .with_columns(
                                    pl.Series("q_SMR", q_values)
                                )
                            )
                            smr_df = (
                                smr_df
                                .with_row_index("_row_id")
                                .join(q_df, on="_row_id", how="left")
                                .drop("_row_id")
                            )

                            # overwrite same SMR result with q_SMR added
                            smr_df.write_csv(f, separator="\t")
                            print(f"[TRACKING] FDR corrected SMR results saved for {pheno_id} in {cell}")
                        else:
                            print(f"[CONCERN] No valid p_SMR values found in {f.name}")

                    else:
                        print(f"[CONCERN] p_SMR not found in {f.name}")

        # delete temp GWAS .ma only after all cell types are ran
        if temp_sumstats.exists():
            temp_sumstats.unlink() 
        
        hits = extract_promising_targets(pheno_id=pheno_id, pqtl_dataset=pqtl_dataset)
        # now extract all of the SMR data from the results for each cell type pertaining to those targets and store as a dataframe in results/SMR/dataset
        # rows == 1 SMR result for target X on cell-type Y 
        # so 7 cell types x X targets in terms of rows
        all_target_smr = []
        for cell in cell_types:
            smr_res = Path(f"./results/SMR/{eqtl_dataset}/{cell}/{pheno_id}")
            for f in smr_res.glob("*.smr"):
                if pheno_id not in f.name:
                    continue
                smr_df = pl.read_csv(f, separator="\t")

                # SMR usually calls the gene / probe column Probe
                # match the gene part of GENE_UNIPROT targets to the SMR Probe column
                if "Gene" not in smr_df.columns:
                    print(f"[CONCERN] Gene column not found in {f.name}")
                    continue

                target_map = {
                    target.split("_")[0]: target for target in hits
                }

                target_genes = list(target_map.keys())

                target_smr = (
                    smr_df
                    .filter(pl.col("Gene").is_in(target_genes))
                    .with_columns(
                        pl.col("Gene").replace(target_map).alias("protein"),
                        pl.lit(cell).alias("cell_type"),
                        pl.lit(pheno_id).alias("phenotype"),
                        pl.lit(eqtl_dataset).alias("eqtl_dataset"),
                        pl.lit(pqtl_dataset).alias("pqtl_dataset")
                    )
                )

                if target_smr.height > 0:
                    all_target_smr.append(target_smr)

        out_dir = Path(f"./results/SMR/{eqtl_dataset}/{pheno_id}")
        os.makedirs(out_dir, exist_ok=True)
        out_file = out_dir / f"{pqtl_dataset}_{pheno_id}_promising_targets_SMR.tsv"

        if len(all_target_smr) > 0:
            final_smr_df = pl.concat(all_target_smr, how="diagonal_relaxed")
            final_smr_df.write_csv(out_file, separator="\t")
            print(f"[TRACKING] Compiled promising target SMR results saved to {out_file}")
        else:
            print(f"[CONCERN] No SMR results found for the promising {pqtl_dataset} targets")




def compile_multi_omics_targets(pheno_id: str, pqtl_dataset: str, eqtl_dataset: str):
    # out_dir = Path(f"./results/SMR/{eqtl_dataset}/{pheno_id}")
    # out_file = out_dir / f"{pqtl_dataset}_{pheno_id}_promising_targets_SMR.tsv"
    # necessary file -> "./results/SMR/{eqtl_dataset}/{pheno_id}/{pqtl_dataset}_{pheno_id}_promising_targets_SMR.tsv"
    out_dir = f"./results/SMR/{eqtl_dataset}/{pheno_id}"
    out_dir = Path(out_dir)
    os.makedirs(out_dir, exist_ok=True)
    targets_path = out_dir / f"{pqtl_dataset}_{pheno_id}_promising_targets_SMR.tsv"
    df = pl.read_csv(targets_path, separator="\t")

    # HEIDI CUT OFF = 0.01
    # SMR CUT OFF = 0.05
    # q_SMR - P_HEIDI

    heidi_col = "P_HEIDI"
    if "p_HEIDI" in df.columns:
        heidi_col = "p_HEIDI"

    final_targets_df = (
        df
        .filter(
            pl.col("q_SMR").is_not_null(),
            pl.col(heidi_col).is_not_null(),
            pl.col("q_SMR") < 0.05,
            pl.col(heidi_col) > 0.01
        )
        .sort(["protein", "cell_type", "q_SMR"])
    )

    if final_targets_df.height == 0:
        print(f"[CONCERN] No drug targets passed cis-MR (pQTLs) + COLOC + single-cell eQTL SMR")
        return []

    # save all rows - 1 row per target x cell type
    final_targets_file = out_dir / f"{pqtl_dataset}_{pheno_id}_final_multi_omics_targets.tsv"
    final_targets_df.write_csv(final_targets_file, separator="\t")
    # unique targets only for the next steps
    targets = (
        final_targets_df
        .select("protein")
        .unique()
        .sort("protein")
        .get_column("protein")
        .to_list()
    )

    print(f"[TRACKING] {final_targets_df.height} target x cell-type SMR hits found")
    print(f"[TRACKING] {len(targets)} unique drug targets passed cis-MR (pQTLs) + COLOC + single-cell eQTL SMR")
    print(f"[TRACKING] Final target x cell-type results saved to {final_targets_file}")
    print(f"[TRACKING] Drug targets: {targets}")
    return targets


# THEN
# -> For each prioritised target
# -> Check original cis-region (matched with GWAS)
# -> Match cis-region with sc-eQTL for cell type X
# -> RUN MOLOC / pairwise coloc

# sumstats: str, ref_bfile: str, maf: float
def main():
    p = argparse.ArgumentParser()
    p.add_argument("--pheno_id", required=True)
    p.add_argument("--sumstats", required=True)
    p.add_argument("--pqtl_dataset", required=True)
    p.add_argument("--eqtl_dataset", required=True)
    p.add_argument("--ref_bfile", required=True)
    p.add_argument("--maf", type=float, default=0.01)
    args = p.parse_args()

    # running SMR (genome-wide)
    run_single_cell_smr(
        pqtl_dataset=args.pqtl_dataset,
        eqtl_dataset=args.eqtl_dataset,
        pheno_id=args.pheno_id,
        sumstats=args.sumstats,
        maf=args.maf,
        ref_bfile=args.ref_bfile
    )

    # final hits
    compile_multi_omics_targets(
        pheno_id=args.pheno_id,
        pqtl_dataset=args.pqtl_dataset,
        eqtl_dataset=args.eqtl_dataset
    )

if __name__ == "__main__":
    main()