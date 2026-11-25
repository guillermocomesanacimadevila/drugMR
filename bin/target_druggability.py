import argparse
from datetime import date

import polars as pl
from drugmr.paths import (
    coloc_out,
    druggability_drugs_out,
    druggability_out,
    druggability_raw_dir,
    druggability_sources_out,
    mr_out,
    pwcoco_out,
)
from drugmr.utils import (
    extract_coloc_or_pwcoco_targets,
    grab_cis_mr_hits,
    map_symbol_to_ens_id,
)


def write_empty_outputs(hits, out_file, drugs_file, sources_file, status):
    """
    Write the 3 output tables when no druggability lookup was done, so the
    run still has every file it promises. status says why (no targets, or
    the error when a database could not be reached)
    """
    hits.write_csv(out_file, separator="\t")

    drugs = pl.DataFrame(schema={
        "ensembl_gene_id": pl.Utf8,
        "symbol": pl.Utf8,
        "uniprot": pl.Utf8,
        "required_action": pl.Utf8,
        "direction_vs_mr": pl.Utf8,
    })
    drugs.write_csv(drugs_file, separator="\t")

    sources = pl.DataFrame(
        {"source": ["none"], "release": [None], "accessed": [str(date.today())], "status": [status]},
        schema={
            "source": pl.Utf8,
            "release": pl.Utf8,
            "accessed": pl.Utf8,
            "status": pl.Utf8,
        },
    )
    sources.write_csv(sources_file, separator="\t")


def druggability_assessment(
        pheno_id: str,
        pqtl_dataset: str,
        gene_annotation,
        local_results_dir = "results",
        pchembl_thresh=6,
        wald_fdr_q=0.05,
        ivw_fdr_q=0.05,
        cochran_q_pval=0.05,
        egger_intercept_pval_min=0,
        min_instruments_for_ivw=3,
        pp4_threshold=0.7,
):

    """
    Steps for druggability assessment
    - Extract proteogenomic targets (function in drugmr/utils.py)
    - Compile them into a list 
    - Use map_symbol_to_ens_id to extract ENS_id
    - Run TargetDruggability()
    - Save main 3 tabs as TSVs and the rest as raw/*.parquet
    - Create three new tables in PosgreSQL db
    """

    cis_mr_res = mr_out(pheno_id=pheno_id, pqtl_dataset=pqtl_dataset, out_dir=local_results_dir)
    coloc_res = coloc_out(pheno_id=pheno_id, pqtl_dataset=pqtl_dataset, out_dir=local_results_dir)
    pwcoco_res = pwcoco_out(pheno_id=pheno_id, pqtl_dataset=pqtl_dataset, out_dir=local_results_dir)
    
    cis_mr_hits = grab_cis_mr_hits(
        wald_fdr_q=wald_fdr_q,
        csv_file=cis_mr_res,
        ivw_fdr_q=ivw_fdr_q,
        cochran_q_pval=cochran_q_pval,
        egger_intercept_pval_min=egger_intercept_pval_min,
        min_instruments_for_ivw=min_instruments_for_ivw,
    )

    coloc_pwcoco_hits = extract_coloc_or_pwcoco_targets(
        coloc_csv_file=coloc_res,
        pwcoco_csv_file=pwcoco_res,
        pp4_thresh=pp4_threshold,
    )

    final_targets = [i for i in cis_mr_hits if i in coloc_pwcoco_hits]
    symbols = []
    for target in final_targets:
        symbol = target.split("_")[0]
        if symbol not in symbols:
            symbols.append(symbol)

    symbol_to_ens_id_dict = map_symbol_to_ens_id(
        target_list=symbols,
        ncbi_ref=pl.read_csv(gene_annotation, separator="\t", infer_schema_length=None)
    )

    # MR direction: protein raises risk -> Inhibit, lowers risk -> Activate
    mr = pl.read_csv(cis_mr_res, separator="\t", infer_schema_length=None, schema_overrides={"primary_FDR_q": pl.Float64})
    mr = mr.filter(pl.col("protein").is_in(final_targets))
    
    # now look for direction of effect 
    rows = []
    for row in mr.iter_rows(named=True):
        symbol = row["protein"].split("_")[0]
        if float(row["primary_beta"]) > 0:
            required_action = "Inhibit"
        else:
            required_action = "Activate"

        rows.append({
            "protein": row["protein"],
            "symbol": str(symbol),
            "ensembl_gene_id": symbol_to_ens_id_dict.get(symbol),
            "mr_method": row["primary_method"],
            "mr_beta": row["primary_beta"],
            "mr_fdr_q": row["primary_FDR_q"],
            "required_action": required_action,
        })

    hits = pl.DataFrame(rows, schema={
        "protein": pl.Utf8,
        "symbol": pl.Utf8,
        "ensembl_gene_id": pl.Utf8,
        "mr_method": pl.Utf8,
        "mr_beta": pl.Float64,
        "mr_fdr_q": pl.Float64,
        "required_action": pl.Utf8,
    })

    hits = hits.with_columns(pl.col("ensembl_gene_id").is_null().alias("not_mapped"))

    # outputs in results/druggability/...
    out_file = druggability_out(pheno_id=pheno_id, pqtl_dataset=pqtl_dataset, out_dir=local_results_dir)
    drugs_file = druggability_drugs_out(pheno_id=pheno_id, pqtl_dataset=pqtl_dataset, out_dir=local_results_dir)
    sources_file = druggability_sources_out(pheno_id=pheno_id, pqtl_dataset=pqtl_dataset, out_dir=local_results_dir)
    raw_dir = druggability_raw_dir(pheno_id=pheno_id, pqtl_dataset=pqtl_dataset, out_dir=local_results_dir)
    raw_dir.mkdir(parents=True, exist_ok=True)

    # compilation of Ensembl IDs
    ens_ids = hits["ensembl_gene_id"].drop_nulls().to_list()
    if len(ens_ids) == 0:
        print("no targets to assess")
        write_empty_outputs(hits, out_file, drugs_file, sources_file, status="no targets")
        return

    # druggability, a database being down should not fail the whole run.
    # imported here because the ChEMBL client contacts ChEMBL on import.
    # a missing package (ImportError) still fails loudly
    try:
        from drugmr.druggability import TargetDruggability
        td = TargetDruggability()
        results = td.run(ens_ids, pchembl_thresh=pchembl_thresh)
    except ImportError:
        raise
    except Exception as error:
        print("druggability lookup failed:", error)
        write_empty_outputs(hits, out_file, drugs_file, sources_file, status=f"failed: {error}")
        return

    # 1 row per target, unmapped targets kept
    target_druggability = hits.join(
        results["target_druggability"].drop("symbol"),
        on="ensembl_gene_id",
        how="left",
    )

    # 1 row per target x drug, does the drug act the way the MR says we need?
    target_drugs = results["target_drugs"].join(
        hits.select(["ensembl_gene_id", "required_action"]).unique(),
        on="ensembl_gene_id",
        how="left",
    )
    target_drugs = target_drugs.with_columns(
        pl.when(pl.col("drug_direction") == pl.col("required_action")).then(pl.lit("Agrees"))
        .when(pl.col("drug_direction").is_in(["Inhibit", "Activate"])).then(pl.lit("Opposes"))
        .otherwise(pl.lit("Unclear"))
        .alias("direction_vs_mr")
    )

    # drugs that agree / oppose the MR direction, per target
    direction_counts = target_drugs.group_by("ensembl_gene_id").agg(
        (pl.col("direction_vs_mr") == "Agrees").sum().alias("n_drugs_agree"),
        (pl.col("direction_vs_mr") == "Opposes").sum().alias("n_drugs_oppose"),
    )
    target_druggability = target_druggability.join(direction_counts, on="ensembl_gene_id", how="left")
    target_druggability = target_druggability.with_columns(
        pl.col("n_drugs_agree").fill_null(0),
        pl.col("n_drugs_oppose").fill_null(0),
    )

    # the 3 main tables as TSVs, the untouched source tables as parquet
    target_druggability.write_csv(out_file, separator="\t")
    target_drugs.write_csv(drugs_file, separator="\t")
    sources = results["sources"].with_columns(pl.lit("ok").alias("status"))
    sources.write_csv(sources_file, separator="\t")
    for name, table in results["raw"].items():
        table.write_parquet(raw_dir / f"{name}.parquet")


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--pqtl_dataset", required=True)
    p.add_argument("--pheno_id", required=True)
    p.add_argument("--gene_annotation", required=True)
    p.add_argument("--pchembl_thresh", type=float, default=6)
    p.add_argument("--local_results_dir", default="results")
    p.add_argument("--wald_fdr_q", type=float, default=0.05)
    p.add_argument("--ivw_fdr_q", type=float, default=0.05)
    p.add_argument("--cochran_q_pval", type=float, default=0.05)
    p.add_argument("--egger_intercept_pval_min", type=float, default=0)
    p.add_argument("--min_instruments_for_ivw", type=int, default=3)
    p.add_argument("--pp4_threshold", type=float, default=0.7)
    args = p.parse_args()
    druggability_assessment(
        pheno_id=args.pheno_id,
        pqtl_dataset=args.pqtl_dataset,
        gene_annotation=args.gene_annotation,
        pchembl_thresh=args.pchembl_thresh,
        local_results_dir=args.local_results_dir,
        wald_fdr_q=args.wald_fdr_q,
        ivw_fdr_q=args.ivw_fdr_q,
        cochran_q_pval=args.cochran_q_pval,
        egger_intercept_pval_min=args.egger_intercept_pval_min,
        min_instruments_for_ivw=args.min_instruments_for_ivw,
        pp4_threshold=args.pp4_threshold,
    )

if __name__ == "__main__":
     main()