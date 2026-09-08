#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

/*
Process mimicking local.py within drugmr/
for nextflow-informed module
*/

process GWAS_QC {

    tag "gwas_qc"
    
    publishDir { "${params.out_dir}/${meta.pheno_id}/qc" }, mode: 'copy'
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
    tuple val(meta), path("${meta.pheno_id}.tsv"), emit: qc_tsv

    script:

    // flag optional params for QC script
    def info_args = ""
    def mhc_flag = ""
    def apoe_flag = ""

    if (meta.info_col && meta.info_threshold != null) {
        info_args = "--info-col ${meta.info_col} --info_threshold ${meta.info_threshold}"
    }

    if (meta.remove_mhc != null) {
        mhc_flag = "--remove_mhc"
    }

    if (meta.remove_apoe != null) {
        apoe_flag = "--remove_apoe"
    }

    """
    python ${projectDir}/bin/qc_gwas.py \\
        --pheno-id ${meta.pheno_id} \\
        --sumstats ${meta.sumstats} \\
        --out-dir . \\
        --maf ${meta.maf} \\
        --snp-col ${meta.snp_col} \\
        --a1-col ${meta.a1_col} \\
        --a2-col ${meta.a2_col} \\
        --beta-col ${meta.beta_col} \\
        --se-col ${meta.se_col} \\
        --p-col ${meta.p_col} \\
        --pos-col ${meta.pos_col} \\
        --chr-col ${meta.chr_col} \\
        --af_col ${meta.af_col} \\
        --genome_build ${meta.genome_build} \\
        --target_build ${meta.target_build} \\
        --n_cases ${meta.n_cases} \\
        --n_controls ${meta.n_controls} \\
        --falcon-user nf \\
        ${info_args} \\
        ${mhc_flag} \\
        ${apoe_flag}
    """
}