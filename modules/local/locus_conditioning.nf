#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process LOCUS_CONDITIONING {

    tag "locus_conditioning_${meta.pheno_id}_${meta.pqtl_dataset}"
    label "process_massive"

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
        path(mr_instruments_tsv, stageAs: "cis_mr/instruments/mr_instruments.tsv"),
        path(coloc_tsv, stageAs: "coloc/*"),
        path(pwcoco_tsv, stageAs: "pwcoco/summary/pwcoco.tsv"),
        path(protein_dirs, stageAs: "protein_dirs/*")

    output:
    tuple val(meta), path("locus_conditioning/${meta.pqtl_dataset}_${meta.pheno_id}_locus_conditioning.tsv"), emit: locus_conditioning_results

    script:
    """
    export PYTHONPATH=${projectDir}
    python ${projectDir}/bin/condition_on_nearby_targets.py \\
        --pqtl_dataset ${meta.pqtl_dataset} \\
        --pheno_id ${meta.pheno_id} \\
        --ref_bfile ${projectDir}/${meta.ref_bfile} \\
        --n_cases ${meta.n_cases} \\
        --n_controls ${meta.n_controls} \\
        --local_results_dir . \\
        --cis_regions_dir protein_dirs \\
        --wald_fdr_q ${meta.gates.cis_mr.wald_fdr_q} \\
        --ivw_fdr_q ${meta.gates.cis_mr.ivw_fdr_q} \\
        --cochran_q_pval ${meta.gates.cis_mr.cochran_q_pval} \\
        --egger_intercept_pval_min ${meta.gates.cis_mr.egger_intercept_pval_min} \\
        --min_instruments_for_ivw ${meta.gates.cis_mr.min_instruments_for_ivw} \\
        --pp4_threshold ${meta.gates.coloc.pp4_threshold} \\
        --p1 ${meta.gates.coloc.p1} \\
        --p2 ${meta.gates.coloc.p2} \\
        --p12 ${meta.gates.coloc.p12} \\
        --gws_pval ${meta.gates.locus_conditioning.gws_pval} \\
        --clump_r2 ${meta.gates.locus_conditioning.clump_r2} \\
        --clump_kb ${meta.gates.locus_conditioning.clump_kb} \\
        --own_signal_r2 ${meta.gates.locus_conditioning.own_signal_r2} \\
        --flag_r2 ${meta.gates.locus_conditioning.flag_r2} \\
        --cojo_collinear ${meta.gates.locus_conditioning.cojo_collinear} \\
        --target_build ${meta.target_build}
    """
}
