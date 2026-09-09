#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process PWCOCO_QTL {

    tag "pwcoco_qtl_${meta.pheno_id}_${meta.pqtl_dataset}"

    publishDir { "${params.runs_root}/${params.run_id}/results/pwcoco" }, mode: "copy"

    container {
        if (params.image_name) {
            return params.image_name
        } else {
            return 'ghcr.io/guillermocomesanacimadevila/drugmr:latest'
        }
    }

    input:
    tuple val(meta), path(smr_final_targets, stageAs: "smr/final_multi_omics_targets.tsv"), path(pwcoco_tsv, stageAs: "pwcoco/summary/pwcoco.tsv"), path(protein_dirs, stageAs: "protein_dirs/*")

    output:
    tuple val(meta), path("pwcoco/summary/pwcoco_eqtl_pqtl.tsv"), emit: pwcoco_eqtl_pqtl_results, optional: true
    tuple val(meta), path("pwcoco/summary/pwcoco_eqtl_gwas.tsv"), emit: pwcoco_eqtl_gwas_results, optional: true
    tuple val(meta), path("pwcoco/summary/pwcoco_shared_snps.tsv"), emit: pwcoco_qtl_shared_results, optional: true

    script:
    """
    export PYTHONPATH=${projectDir}
    python ${projectDir}/bin/pwcoco_qtl_wrapper.py \\
        --pqtl_dataset ${meta.pqtl_dataset} \\
        --pheno_id ${meta.pheno_id} \\
        --ref_bfile ${projectDir}/${meta.ref_bfile} \\
        --n_cases ${meta.n_cases} \\
        --n_controls ${meta.n_controls} \\
        --local_results_dir . \\
        --repo_root ${projectDir} \\
        --manifest_path ${projectDir}/assets/qtl_manifest.csv \\
        --cis_regions_dir protein_dirs
    """
}
