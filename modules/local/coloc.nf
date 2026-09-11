#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process PAIRWISE_COLOC {

    tag "pairwise_coloc_${meta.pheno_id}_${meta.pqtl_dataset}"
    label "process_low"

    // publishDir { "${params.runs_root}/${params.run_id}/results/coloc" }, mode: "copy"
    publishDir { "${params.runs_root}/${params.run_id}/results" }, mode: "copy"

    container {
        if (params.image_name) {
            return params.image_name
        } else {
            return 'ghcr.io/guillermocomesanacimadevila/drugmr:latest'
        }
    }

    input:
    tuple val(meta), path(mr_tsv, stageAs: "cis_mr/mr.tsv"), path(mr_instruments_tsv, stageAs: "cis_mr/instruments/mr_instruments.tsv"), path(protein_dirs, stageAs: "protein_dirs/*")

    output:
    tuple val(meta), path("coloc/${meta.pqtl_dataset}_${meta.pheno_id}_all_coloc.tsv"), path("coloc/coloc_sensitivity.tsv"), emit: coloc_results

    // no --repo_root/--manifest_path: coloc_targets.py never touches the QTL
    // manifest at all, it only reads already-extracted cis-region parquet
    // files passed in via protein_dirs.
    script:
    """
    export PYTHONPATH=${projectDir}
    python ${projectDir}/bin/coloc_targets.py \\
        --pqtl_dataset ${meta.pqtl_dataset} \\
        --local_results_dir . \\
        --pqtl_dir protein_dirs \\
        --pheno_id ${meta.pheno_id} \\
        --n_cases ${meta.n_cases} \\
        --n_controls ${meta.n_controls} \\
        --wald_fdr_q ${meta.gates.cis_mr.wald_fdr_q} \\
        --ivw_fdr_q ${meta.gates.cis_mr.ivw_fdr_q} \\
        --cochran_q_pval ${meta.gates.cis_mr.cochran_q_pval} \\
        --egger_intercept_pval_min ${meta.gates.cis_mr.egger_intercept_pval_min} \\
        --min_instruments_for_ivw ${meta.gates.cis_mr.min_instruments_for_ivw} \\
        --pp4_threshold ${meta.gates.coloc.pp4_threshold} \\
        --p1 ${meta.gates.coloc.p1} \\
        --p2 ${meta.gates.coloc.p2} \\
        --p12 ${meta.gates.coloc.p12}
    """

}
