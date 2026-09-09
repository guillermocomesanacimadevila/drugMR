#!/usr/bin/env nextflow
nextflow.enable.dsl=2


process CIS_MR {

    tag "cis_mr"

    publishDir { "${params.runs_root}/${params.run_id}/results/cis_mr" }, mode: "copy"

    container {
        if (params.image_name) {
            return params.image_name
        } else {
            return 'ghcr.io/guillermocomesanacimadevila/drugmr:latest'
        }
    }

    input:
    tuple val(meta), path(pheno_gwas), path(protein_dirs, stageAs: 'protein_dirs/*')

    output:
    tuple val(meta), path("cis_mr/mr.tsv"), path("cis_mr/instruments/mr_instruments.tsv"), emit: mr_results

    script:
    """
    Rscript ${projectDir}/bin/cis_mr.R \\
        ${meta.pqtl_dataset} \\
        protein_dirs \\
        ${meta.pheno_id} \\
        ${pheno_gwas} \\
        ${projectDir}/${meta.ref_bfile} \\
        . \\
    """
}