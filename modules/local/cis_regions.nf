#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process EXTRACT_CIS_REGIONS {

    tag "extract_cis_regions"

    publishDir { "./dat/cis_regions/${meta.pqtl_dataset}"}, mode: "copy"

    container {
        if (params.image_name) {
            return params.image_name
        } else {
            return 'ghcr.io/guillermocomesanacimadevila/drugmr:latest'
        }
    }

    input:
    tuple val(meta), path(qc_tsv), path(qtl_manifest)

    output:
    tuple val(meta), path("*", type: 'dir'), emit: protein_dirs

    script:
    """
    export PYTHONPATH=${projectDir}
    python ${projectDir}/bin/prep_cis_regions.py \\
        --pqtl_dataset ${meta.pqtl_dataset} \\
        --pheno_id ${meta.pheno_id} \\
        --manifest_path ${qtl_manifest} \\
        --qc_tsv ${qc_tsv}
    """
}