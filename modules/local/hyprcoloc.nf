#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process HYPRCOLOC {

    tag "multi_omics_hyprcoloc"

    publishDir { "${params.runs_root}/${params.run_id}/results/hyprcoloc" }, mode: "copy"

    container {
        if (params.image_name) {
            return params.image_name
        } else {
            return 'ghcr.io/guillermocomesanacimadevila/drugmr:latest'
        }
    } 

    input:
    tuple val(meta)

    output:
    tuple val(meta), path("hyprcoloc/hyprcoloc.tsv"), path("hyprcoloc/sensitivity_plots/*"), emit: hyprcoloc_results

    script:
    """
    export PYTHONPATH=${projectDir}
    python ${projectDir}/bin/hyprcoloc_targets.py \\
    
    """

}

