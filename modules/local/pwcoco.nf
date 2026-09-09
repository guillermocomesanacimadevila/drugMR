#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process PWCOCO {

    tag "pwcoco_targets_in_phase1"

    publishDir { "${params.runs_root}/${params.run_id}/results/pwcoco" }, mode: "copy"

    container {
        if (params.image_name) {
            return params.image_name
        } else {
            return 'ghcr.io/guillermocomesanacimadevila/drugmr:latest'
        }
    }

    input:
    tuple val(meta), path(mr_tsv, stageAs: "cis_mr/mr.tsv"), path(protein_dirs, stageAs: "protein_dirs/*")

    output:
    tuple val(meta), path("pwcoco/summary/pwcoco.tsv"), emit: pwcoco_results

    script:
    """
    export PYTHONPATH=${projectDir}
    python ${projectDir}/bin/pwcoco_wrapper.py \\
        --pqtl_dataset ${meta.pqtl_dataset} \\
        --pheno_id ${meta.pheno_id} \\
        --ref_bfile ${projectDir}/${meta.ref_bfile} \\
        --n_cases ${meta.n_cases} \\
        --n_controls ${meta.n_controls} \\
        --local_results_dir . \\
        --cis_regions_dir protein_dirs
    """
}