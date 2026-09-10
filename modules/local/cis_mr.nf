#!/usr/bin/env nextflow
nextflow.enable.dsl=2


process CIS_MR {

    tag "cis_mr_${meta.pheno_id}_${meta.pqtl_dataset}"
    label "process_low"

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

    // no --repo_root/--manifest_path: cis_mr.R never touches the QTL
    // manifest at all, it only reads already-extracted cis-region parquet
    // files passed in via protein_dirs.
    script:
    """
    Rscript ${projectDir}/bin/cis_mr.R \\
        ${meta.pqtl_dataset} \\
        protein_dirs \\
        ${meta.pheno_id} \\
        ${pheno_gwas} \\
        ${projectDir}/${meta.ref_bfile} \\
        . \\
        ${meta.gates.cis_mr.clump_kb} \\
        ${meta.gates.cis_mr.clump_r2} \\
        ${meta.gates.cis_mr.instrument_pval_threshold} \\
        ${meta.gates.cis_mr.min_f_stat} \\
        ${meta.gates.cis_mr.apply_steiger_filter} \\
        ${meta.maf}
    """
}