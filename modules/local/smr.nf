#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process SMR_BULK {
    
    tag "smr_bulk_${meta.pheno_id}_${eqtl_dataset}"

    publishDir { "${params.out_dir}/${meta.pheno_id}/${meta.pqtl_dataset}/smr/bulk/${eqtl_dataset}" }, mode: "copy"

    container {
        if (params.image_name) {
            return params.image_name
        } else {
            return 'ghcr.io/guillermocomesanacimadevila/drugmr:latest'
        }
    }

    input:
    tuple val(meta), 




}

process SMR_SC {
    
}