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
    tuple val(meta), val(eqtl_dataset), path(smr_final_targets, stageAs: "smr/final_multi_omics_targets.tsv"), path(protein_dirs, stageAs: "protein_dirs/*")

    output:
    tuple val(meta), val(eqtl_dataset), path("hyprcoloc/by_eqtl_source/${eqtl_dataset}/hyprcoloc.tsv"), path("hyprcoloc/by_eqtl_source/*.pdf"), emit: hyprcoloc_results

    script:
    """
    export PYTHONPATH=${projectDir}
    python ${projectDir}/bin/hyprcoloc_targets.py \\
        --pqtl_dataset ${meta.pqtl_dataset} \\
        --pheno_id ${meta.pheno_id} \\
        --eqtl_dataset ${eqtl_dataset} \\
        --local_results_dir . \\
        --repo_root ${projectDir} \\
        --manifest_path ${projectDir}/assets/qtl_manifest.csv \\
        --skip_merge
    """

}

process MERGE_HYPRCOLOC {

    tag "hyprcoloc_merge_${meta.pheno_id}_${meta.pqtl_dataset}"

    publishDir { "${params.runs_root}/${params.run_id}/results/hyprcoloc" }, mode: "copy"

    container {
        if (params.image_name) {
            return params.image_name
        } else {
            return 'ghcr.io/guillermocomesanacimadevila/drugmr:latest'
        }
    }

    input:
    tuple val(meta), val(labels), path(dataset_files, stageAs: "hc_target_?/hyprcoloc.tsv")

    output:
    tuple val(meta), path("hyprcoloc/hyprcoloc.tsv"), emit: hyprcoloc_final

    script:
    def input_flags = labels.withIndex(1).collect { label, i -> "--input '${label}:hc_target_${i}/hyprcoloc.tsv'" }.join(' \\\n        ')
    """
    export PYTHONPATH=${projectDir}
    python ${projectDir}/bin/hyprcoloc_targets.py \\
        --pheno_id ${meta.pheno_id} \\
        --pqtl_dataset ${meta.pqtl_dataset} \\
        --local_results_dir . \\
        --merge \\
        ${input_flags}
    """
}

