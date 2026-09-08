#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process EXTRACT_CIS_REGIONS {

    tag "extract_cis_regions"

    publishDir { "./dat/cis_regions/${params.pqtl_dataset}"}, mode: "copy"

    container {
        if (params.image_name) {
            return params.image_name
        } else {
            return 'ghcr.io/guillermocomesanacimadevila/drugmr:latest'
        }
    }

    input:
    val(meta)

    output:
    // we need the protein name
    // tuple val(meta), path("${meta.pheno_id}.tsv"), emit: qc_tsv

    script:
    """
    python bin/prep_cis_regions \\
        --pqtl_dataset ${meta.pqtl_dataset} \\
        --pheno_id ${meta.pheno_id}
    """
}