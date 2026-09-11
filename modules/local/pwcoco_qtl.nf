#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process PWCOCO_QTL {

    tag "pwcoco_qtl_${meta.pheno_id}_${meta.pqtl_dataset}"
    label "process_massive"

    // publishDir { "${params.runs_root}/${params.run_id}/results/pwcoco" }, mode: "copy"
    publishDir { "${params.runs_root}/${params.run_id}/results" }, mode: "copy"

    container {
        if (params.image_name) {
            return params.image_name
        } else {
            return 'ghcr.io/guillermocomesanacimadevila/drugmr:latest'
        }
    }

    input:
    tuple val(meta), path(smr_final_targets, stageAs: "smr/final_multi_omics_targets.tsv"), path(pwcoco_tsv, stageAs: "pwcoco/summary/pwcoco.tsv"), path(protein_dirs, stageAs: "protein_dirs/*")

    // None of these 3 outputs are consumed by any further Nextflow stage,
    // deliberately, not a missing wiring step:
    // - pwcoco_eqtl_pqtl_results / pwcoco_eqtl_gwas_results are read directly
    //   from disk by the dashboard (dashboard/mr_app.py's own "PWCoCo-QTL
    //   (eQTL triangulation)" page).
    // - pwcoco_qtl_shared_results (pwcoco_shared_snps.tsv) is written but
    //   never read by anything: the dashboard recomputes the same shared-SNP
    //   triangulation live against its own interactive PP.H4 slider instead
    //   (see mr_app.py's compute_shared_snp_table), so it always agrees with
    //   whatever threshold the user has selected rather than being stuck at
    //   whatever value the pipeline used when this file was written.
    output:
    tuple val(meta), path("pwcoco/summary/pwcoco_eqtl_pqtl.tsv"), emit: pwcoco_eqtl_pqtl_results, optional: true
    tuple val(meta), path("pwcoco/summary/pwcoco_eqtl_gwas.tsv"), emit: pwcoco_eqtl_gwas_results, optional: true
    tuple val(meta), path("pwcoco/summary/pwcoco_shared_snps.tsv"), emit: pwcoco_qtl_shared_results, optional: true

    // --manifest_path is pinned explicitly here (belt and braces), since
    // pwcoco_qtl_wrapper.py also uses SMRUtils to load eQTL rows for
    // triangulation, same as sort_smr.py/hyprcoloc_targets.py.
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
        --cis_regions_dir protein_dirs \\
        --pp4_threshold ${meta.gates.pwcoco.pp4_threshold}
    """
}
