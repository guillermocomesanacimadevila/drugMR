#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

process EXPORT_LOCUS_ASSETS {

    tag "locus_assets_${meta.pheno_id}_${meta.pqtl_dataset}"
    label "process_low"

    publishDir { "${params.runs_root}/${params.run_id}/results" }, mode: "copy"

    container {
        if (params.image_name) {
            return params.image_name
        } else {
            return 'ghcr.io/guillermocomesanacimadevila/drugmr:latest'
        }
    }

    input:
    tuple val(meta), path(target_stats, stageAs: "target_stats/top_cis_hits.tsv"), path(protein_dirs, stageAs: "protein_dirs/*")

    output:
    tuple val(meta), path("locus_data"), emit: locus_assets

    script:
    def refBfile
    if (new File(meta.ref_bfile as String).isAbsolute()) {
        refBfile = meta.ref_bfile
    } else {
        refBfile = "${projectDir}/${meta.ref_bfile}"
    }

    """
    export PYTHONPATH=${projectDir}
    export PYTHONUNBUFFERED=1

    python ${projectDir}/bin/export_locus_assets.py \
        --pqtl-dataset ${meta.pqtl_dataset} \
        --cis-regions-dir protein_dirs \
        --target-stats target_stats/top_cis_hits.tsv \
        --ref-bfile ${refBfile} \
        --out-dir locus_data
    """
}
