#!/usr/bin/env python3
import argparse
from pathlib import Path
import polars as pl

# FINAL MULTI-OMICS SUMMARY
# -> load all downstream results
# -> merge into one overview table
# -> extract SNP-level evidence
# -> save dashboard-ready tables

# load all downstream result files
def load_results(pheno_id: str, pqtl_dataset: str, eqtl_dataset: str):
    base_dir = Path("./results")
    cis_mr = pl.read_csv(base_dir / f"cis-MR/{pqtl_dataset}_{pheno_id}_all_MR.tsv", separator="\t")
    pqtl_coloc = pl.read_csv(base_dir / f"coloc/{pqtl_dataset}/{pqtl_dataset}_{pheno_id}_all_coloc.tsv", separator="\t")
    smr = pl.read_csv(base_dir / f"SMR/{eqtl_dataset}/{pheno_id}/{pqtl_dataset}_{pheno_id}_final_multi_omics_targets.tsv", separator="\t")
    # eqtl_coloc = pl.read_csv(base_dir / f"eqtl_coloc/{eqtl_dataset}/{pheno_id}/{pqtl_dataset}_{pheno_id}_eqtl_coloc.tsv", separator="\t")
    eqtl_coloc = pl.read_csv(base_dir / f"eQTL_coloc/{pqtl_dataset}/{eqtl_dataset}/{pheno_id}/{pqtl_dataset}_{pheno_id}_{eqtl_dataset}_all_eqtl_coloc.tsv", separator="\t")
    moloc = pl.read_csv(base_dir / f"QTL_moloc/{pqtl_dataset}/{eqtl_dataset}/{pheno_id}/{pqtl_dataset}_{pheno_id}_{eqtl_dataset}_moloc_summary.tsv",separator="\t")    
    # print(cis_mr.columns)
    # print(pqtl_coloc.columns)
    # print(smr.columns)
    # print(eqtl_coloc.columns)
    # print(moloc.columns)
    # mitmatch in colnames
    pqtl_coloc = pqtl_coloc.rename({"protein_id": "protein"})
    eqtl_coloc = eqtl_coloc.rename({"protein_id": "protein"})
    pqtl_coloc = pqtl_coloc.select(["protein", pl.col("top_snp").alias("pqtl_coloc_top_snp"), pl.col("PP.H4.abf").alias("pqtl_pp_h4"), "coloc_pass"])
    eqtl_coloc = eqtl_coloc.select(["protein", "cell_type", pl.col("top_snp").alias("eqtl_coloc_top_snp"), pl.col("PP.H4.abf").alias("eqtl_pp_h4"), "coloc_pass"])
    moloc = moloc.select(["protein", "cell_type", pl.col("model").alias("moloc_model"), pl.col("PPA").alias("moloc_ppa")])
    print(f"[TRACKING] cis-MR results loaded: {cis_mr.height}")
    print(f"[TRACKING] pQTL coloc results loaded: {pqtl_coloc.height}")
    print(f"[TRACKING] SMR results loaded: {smr.height}")
    print(f"[TRACKING] eQTL coloc results loaded: {eqtl_coloc.height}")
    print(f"[TRACKING] moloc results loaded: {moloc.height}")
    return cis_mr, pqtl_coloc, smr, eqtl_coloc, moloc


# merge everything into one final overview
def build_overview(cis_mr, pqtl_coloc, smr, eqtl_coloc, moloc):
    print("[TRACKING] Building multi-omics overview")
    overview = (
        smr
        .join(cis_mr, on="protein", how="left")
        .join(pqtl_coloc, on="protein", how="left")
        .join(eqtl_coloc, on=["protein","cell_type"], how="left")
        .join(moloc, on=["protein","cell_type"], how="left")
    )
    print(overview)
    print(f"[TRACKING] {overview.height} overview rows created")
    return overview

def build_snp_evidence(overview: pl.DataFrame, pheno_id: str, pqtl_dataset: str, eqtl_dataset: str):
    print("[TRACKING] Building SNP evidence table")
    pqtl_dir = Path(f"./dat/cis_regions/{pqtl_dataset}")
    eqtl_dir = Path(f"./dat/sc-eQTL/{eqtl_dataset}")
    all_snps = []
    for row in overview.iter_rows(named=True):
        protein = row["protein"]
        gene = row["probeID"].split(".")[0]
        cell_type = row["cell_type"]
        print(f"[TRACKING] Extracting SNP evidence for {protein} ({cell_type})")
        gwas = pl.read_parquet(pqtl_dir / protein / "gwas.parquet")
        pqtl = pl.read_parquet(pqtl_dir / protein / "pqtl.parquet")
        eqtl = (pl.read_parquet(eqtl_dir / f"{cell_type}.parquet").with_columns(pl.col("GENE").str.split(".").list.first().alias("_gene")).filter(pl.col("_gene") == gene).drop("_gene"))

        for snp, label in [(row["pqtl_coloc_top_snp"], "Lead pQTL SNP"), (row["topSNP"], "Lead sc-eQTL / SMR SNP")]:
            g = gwas.filter(pl.col("SNP") == snp)
            p = pqtl.filter(pl.col("SNP") == snp)
            e = eqtl.filter(pl.col("SNP") == snp)

            if g.height == 0 or p.height == 0:
                print(f"[CONCERN] {snp} not found for {protein}")
                continue

            g = g.row(0, named=True)
            p = p.row(0, named=True)
            e = e.row(0, named=True) if e.height else {}
            # orient everything to phenotype risk allele
            risk_a1 = g["A1"]
            risk_a2 = g["A2"]
            gwas_beta = g["BETA"]

            if gwas_beta < 0:
                risk_a1, risk_a2 = risk_a2, risk_a1
                gwas_beta *= -1

            pqtl_beta = p["BETA"]
            eqtl_beta = e.get("BETA")
            pqtl_flip = False
            eqtl_flip = False
            pqtl_a1 = p["A1"]
            pqtl_a2 = p["A2"]
            eqtl_effect_allele = e.get("EA")
            eqtl_a1 = risk_a1
            eqtl_a2 = risk_a2

            # align pQTL to GWAS risk allele
            if pqtl_a1 != risk_a1:
                pqtl_beta *= -1
                pqtl_a1, pqtl_a2 = pqtl_a2, pqtl_a1
                pqtl_flip = True

            # align sc-eQTL to GWAS risk allele (BETA is relative to EA, NOT A1)
            if eqtl_beta is not None:
                if eqtl_effect_allele == risk_a2:
                    eqtl_beta *= -1
                    eqtl_flip = True
                elif eqtl_effect_allele != risk_a1:
                    print(f"[CONCERN] eQTL EA does not match GWAS alleles for {snp}")
                    eqtl_beta = None
                    eqtl_a1 = None
                    eqtl_a2 = None

            all_snps.append({
                "protein": protein,
                "gene": gene,
                "cell_type": cell_type,
                "phenotype": pheno_id,
                "pqtl_dataset": pqtl_dataset,
                "eqtl_dataset": eqtl_dataset,
                "snp_type": label,
                "SNP": snp,
                "CHR": g["CHR"],
                "POS": g["BP"],
                # alleles (all oriented to AD risk allele)
                "A1": risk_a1,
                "A2": risk_a2,
                "pQTL_A1": pqtl_a1,
                "pQTL_A2": pqtl_a2,
                "sc_eQTL_A1": eqtl_a1,
                "sc_eQTL_A2": eqtl_a2,
                # GWAS
                "GWAS_beta": gwas_beta,
                "GWAS_SE": g["SE"],
                "GWAS_P": g["P"],
                # pQTL
                "pQTL_beta_raw": p["BETA"],
                "pQTL_beta": pqtl_beta,
                "pQTL_SE": p["SE"],
                "pQTL_P": p["P"],
                "pQTL_flipped": pqtl_flip,
                # sc-eQTL
                "sc_eQTL_beta_raw": e.get("BETA"),
                "sc_eQTL_beta": eqtl_beta,
                "sc_eQTL_SE": e.get("SE"),
                "sc_eQTL_P": e.get("P"),
                "sc_eQTL_flipped": eqtl_flip
            })

    if len(all_snps) == 0:
        raise RuntimeError("No SNP evidence could be extracted.")

    print(f"[TRACKING] {len(all_snps)} SNP evidence rows created")
    return pl.DataFrame(all_snps)

# save final tables
def save_results(overview: pl.DataFrame, snp_evidence: pl.DataFrame, pheno_id: str, pqtl_dataset: str, eqtl_dataset: str):
    out_dir = Path(f"./results/SMR/{eqtl_dataset}/{pheno_id}")
    overview.write_csv(out_dir / f"{pqtl_dataset}_{pheno_id}_multi_omics_overview.tsv", separator="\t")
    snp_evidence.write_csv(out_dir / f"{pqtl_dataset}_{pheno_id}_multi_omics_snp_evidence.tsv", separator="\t")
    print("[TRACKING] Multi-omics overview saved")
    print("[TRACKING] SNP evidence table saved")


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--pheno_id", required=True)
    p.add_argument("--pqtl_dataset", required=True)
    p.add_argument("--eqtl_dataset", required=True)
    args = p.parse_args()

    cis_mr, pqtl_coloc, smr, eqtl_coloc, moloc = load_results(
        args.pheno_id,
        args.pqtl_dataset,
        args.eqtl_dataset
    )

    overview = build_overview(
        cis_mr,
        pqtl_coloc,
        smr,
        eqtl_coloc,
        moloc
    )

    snp_evidence = build_snp_evidence(
        overview,
        args.pheno_id,
        args.pqtl_dataset,
        args.eqtl_dataset
    )

    save_results(
        overview,
        snp_evidence,
        args.pheno_id,
        args.pqtl_dataset,
        args.eqtl_dataset
    )

if __name__ == "__main__":
    main()