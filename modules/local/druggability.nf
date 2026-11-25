#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process DRUGGABILITY {

    tag "druggability_${meta.pheno_id}_${meta.pqtl_dataset}"
    // network bound (Pharos, ChEMBL, Open Targets, B3DB), little CPU or memory
    // but ChEMBL can be slow, so it gets the long time limit
    label "process_long"

    publishDir { "${params.runs_root}/${params.run_id}/results" }, mode: "copy"

    container {
        if (params.image_name) {
            return params.image_name
        } else {
            return 'ghcr.io/guillermocomesanacimadevila/drugmr:v1.1.0'
        }
    }

    input:
    tuple val(meta),
        path(mr_tsv, stageAs: "cis_mr/mr.tsv"),
        path(coloc_tsv, stageAs: "coloc/*"),
        path(pwcoco_tsv, stageAs: "pwcoco/summary/pwcoco.tsv")

    output:
    tuple val(meta),
        path("druggability/${meta.pqtl_dataset}_${meta.pheno_id}_target_druggability.tsv"),
        path("druggability/${meta.pqtl_dataset}_${meta.pheno_id}_target_drugs.tsv"),
        path("druggability/${meta.pqtl_dataset}_${meta.pheno_id}_sources.tsv"), emit: druggability_results
    path("druggability/raw/*.parquet"), emit: druggability_raw, optional: true

    script:
    def gene_annotation = meta.gene_annotation
    if (!gene_annotation.startsWith("/")) {
        gene_annotation = "${projectDir}/${meta.gene_annotation}"
    }
    """
    export PYTHONPATH=${projectDir}
    python ${projectDir}/bin/target_druggability.py \\
        --pqtl_dataset ${meta.pqtl_dataset} \\
        --pheno_id ${meta.pheno_id} \\
        --gene_annotation ${gene_annotation} \\
        --pchembl_thresh ${meta.gates.druggability.pchembl_thresh} \\
        --local_results_dir . \\
        --wald_fdr_q ${meta.gates.cis_mr.wald_fdr_q} \\
        --ivw_fdr_q ${meta.gates.cis_mr.ivw_fdr_q} \\
        --cochran_q_pval ${meta.gates.cis_mr.cochran_q_pval} \\
        --egger_intercept_pval_min ${meta.gates.cis_mr.egger_intercept_pval_min} \\
        --min_instruments_for_ivw ${meta.gates.cis_mr.min_instruments_for_ivw} \\
        --pp4_threshold ${meta.gates.coloc.pp4_threshold}
    """
}
