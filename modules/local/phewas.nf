#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process PHEWAS_FINNGEN {

    tag "phewas_finngen_${meta.pheno_id}_${meta.pqtl_dataset}"
    label "process_low"
    label "process_long"

    publishDir { "${params.runs_root}/${params.run_id}/results/phewas" }, mode: "copy"

    container {
        if (params.image_name) {
            return params.image_name
        } else {
            return 'ghcr.io/guillermocomesanacimadevila/drugmr:latest'
        }
    }

    input:
    tuple val(meta), path(coloc_tsv, stageAs: "coloc/coloc.tsv"), path(pwcoco_tsv, stageAs: "pwcoco/summary/pwcoco.tsv"), path(mr_instruments_tsv, stageAs: "cis_mr/instruments/mr_instruments.tsv"), path(protein_dirs, stageAs: "protein_dirs/*")

    output:
    tuple val(meta), path("phewas/finngen/phewas_coverage.tsv"), emit: phewas_finngen_results
    tuple val(meta), path("phewas/finngen/phewas.tsv"), emit: phewas_finngen_hits, optional: true

    script:
    """
    export PYTHONPATH=${projectDir}
    python ${projectDir}/bin/phewas_cis_pqtls.py \\
        --pheno_id ${meta.pheno_id} \\
        --pqtl_dataset ${meta.pqtl_dataset} \\
        --local_results_dir . \\
        --cis_regions_dir protein_dirs \\
        --coloc_file coloc/coloc.tsv \\
        --coloc_threshold ${meta.gates.phewas.coloc_threshold} \\
        --bonferroni_alpha ${meta.gates.phewas.bonferroni_alpha}
    """
}

process PHEWAS_UKB {

    tag "phewas_ukb_${meta.pheno_id}_${meta.pqtl_dataset}"
    label "process_low"
    label "process_long"

    publishDir { "${params.runs_root}/${params.run_id}/results/phewas" }, mode: "copy"

    container {
        if (params.image_name) {
            return params.image_name
        } else {
            return 'ghcr.io/guillermocomesanacimadevila/drugmr:latest'
        }
    }

    input:
    tuple val(meta), path(coloc_tsv, stageAs: "coloc/coloc.tsv"), path(pwcoco_tsv, stageAs: "pwcoco/summary/pwcoco.tsv"), path(mr_instruments_tsv, stageAs: "cis_mr/instruments/mr_instruments.tsv"), path(protein_dirs, stageAs: "protein_dirs/*"), path(finngen_coverage_tsv, stageAs: "phewas/finngen/phewas_coverage.tsv")

    output:
    tuple val(meta), path("phewas/ukbb/phewas.tsv"), emit: phewas_ukb_results, optional: true

    script:
    """
    export PYTHONPATH=${projectDir}
    python ${projectDir}/bin/ukb_phewas.py \\
        --pheno_id ${meta.pheno_id} \\
        --pqtl_dataset ${meta.pqtl_dataset} \\
        --local_results_dir . \\
        --cis_regions_dir protein_dirs \\
        --coloc_file coloc/coloc.tsv \\
        --coloc_threshold ${meta.gates.phewas.coloc_threshold} \\
        --bonferroni_alpha ${meta.gates.phewas.bonferroni_alpha}
    """
}
