#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process HYPRCOLOC {

    tag "multi_omics_hyprcoloc"
    label "process_low"

    // publishDir { "${params.runs_root}/${params.run_id}/results/hyprcoloc" }, mode: "copy"
    publishDir { "${params.runs_root}/${params.run_id}/results" }, mode: "copy"

    container {
        if (params.image_name) {
            return params.image_name
        } else {
            return 'ghcr.io/guillermocomesanacimadevila/drugmr:latest'
        }
    } 

    input:
    tuple val(meta), val(qtl_dataset), path(smr_final_targets, stageAs: "smr/final_multi_omics_targets.tsv"), path(protein_dirs, stageAs: "protein_dirs/*")

    output:
    tuple val(meta), val(qtl_dataset), path("hyprcoloc/by_eqtl_source/${qtl_dataset}/hyprcoloc.tsv"), path("hyprcoloc/by_eqtl_source/*.pdf"), emit: hyprcoloc_results

    // --manifest_path is pinned explicitly here (belt and braces), since
    // hyprcoloc_targets.py also uses SMRUtils to load eQTL rows, same as
    // sort_smr.py - this happens to resolve to the same real file sort_smr.py
    // gets from its own default, just spelled out rather than relied upon.
    script:
    """
    export PYTHONPATH=${projectDir}
    python ${projectDir}/bin/hyprcoloc_targets.py \\
        --pqtl_dataset ${meta.pqtl_dataset} \\
        --pheno_id ${meta.pheno_id} \\
        --qtl_dataset ${qtl_dataset} \\
        --local_results_dir . \\
        --repo_root ${projectDir} \\
        --manifest_path ${projectDir}/assets/qtl_manifest.csv \\
        --prior_1 ${meta.gates.hyprcoloc.prior_1} \\
        --prior_c ${meta.gates.hyprcoloc.prior_c.join(',')} \\
        --reg_thresh ${meta.gates.hyprcoloc.reg_thresh.join(',')} \\
        --align_thresh ${meta.gates.hyprcoloc.align_thresh.join(',')} \\
        --equal_thresholds ${meta.gates.hyprcoloc.equal_thresholds} \\
        --cis_regions_dir protein_dirs \\
        --skip_merge
    """

}

process MERGE_HYPRCOLOC {

    tag "hyprcoloc_merge_${meta.pheno_id}_${meta.pqtl_dataset}"
    label "process_single"

    // publishDir { "${params.runs_root}/${params.run_id}/results/hyprcoloc" }, mode: "copy"
    publishDir { "${params.runs_root}/${params.run_id}/results" }, mode: "copy"

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
