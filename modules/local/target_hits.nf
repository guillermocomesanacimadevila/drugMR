#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

process EXTRACT_TWO_OMICS_HITS {

    tag "extract_targets_from_first_two_layers"

    publishDir { "${params.runs_root}/${params.run_id}/results/target_stats" }, mode: "copy"

    container {
        if (params.image_name) {
            return params.image_name
        } else {
            return 'ghcr.io/guillermocomesanacimadevila/drugmr:latest'
        }
    }

    input:
    tuple val(meta), path(coloc_tsv, stageAs: "coloc/coloc.tsv"), path(pwcoco_tsv, stageAs: "pwcoco/summary/pwcoco.tsv"), path(protein_dirs, stageAs: "protein_dirs/*")

    output:
    tuple val(meta), path("target_stats/top_cis_hits.tsv"), emit: target_stats

    script:
    """
    export PYTHONPATH=${projectDir}
    python ${projectDir}/bin/compile_cis_hit_info.py \\
        --pheno_id ${meta.pheno_id} \\
        --pqtl_dataset ${meta.pqtl_dataset} \\
        --local_results_dir . \\
        --cis_regions_dir protein_dirs
    """
}