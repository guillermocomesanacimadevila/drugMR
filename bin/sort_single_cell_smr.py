#!/usr/bin/env python3
import argparse
import polars as pl
from pathlib import Path
import subprocess
from drugmr import SMR
import os 
from statsmodels.stats.multitest import fdrcorrection

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
        df = pl.read_csv(sumstats, separator="\t")

        # rename and save as .ma 
        df = (
            df.rename({
                "FRQ": "Freq",
                "SE": "se", 
                "P": "p",
                "BETA": "Beta",
                "CHR": "Chr",
                "BP": "Bp"
            })
        )

        temp_sumstats = temp_dir / f"{pheno_id}.ma"
        df.write_csv(temp_sumstats, separator="\t")
        cell_types = ["Ast", "Ext", "MG", "OD", "OPC", "End", "IN"]
        eqtls = "./dat/sc-eQTL/SingleBrain/SMR_ready"
        eqtls = Path(eqtls)
        for file in eqtls.glob("*.besd"):
            cell = file.stem
            if str(cell) in cell_types:
                print(f"[TRACKING] Cell type {cell} found!")
            else:
                print(f"Yowza! Cell type {cell} not found")
                continue
            
            # continue 
            # use prefix without .besd / .esi / .epi for SMR

            #####
            #####
            #####
            beqtl_summary = eqtls / cell
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
                if "Probe" not in smr_df.columns:
                    print(f"[CONCERN] Probe column not found in {f.name}")
                    continue

                target_genes = [target.split("_")[0] for target in hits]
                target_smr = (
                    smr_df
                    .filter(pl.col("Probe").is_in(target_genes))
                    .with_columns(
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

# sumstats: str, ref_bfile: str, maf: float
def main():
    p = argparse.ArgumentParser()
    p.add_argument("--pheno_id", required=True)
    p.add_argument("--sumstats", required=True)
    p.add_argument("--pqtl_dataset", required=True)
    p.add_argument("--eqtl_dataset", required=True)
    p.add_argument("--ref_bfile", required=True)
    p.add_argument("--maf", required=True, default=0.01)
    args = p.parse_args()
    run_single_cell_smr(
        pqtl_dataset=args.pqtl_dataset,
        eqtl_dataset=args.eqtl_dataset,
        pheno_id=args.pheno_id,
        sumstats=args.sumstats,
        maf=args.maf,
        ref_bfile=args.ref_bfile
    )

if __name__ == "__main__":
    main()