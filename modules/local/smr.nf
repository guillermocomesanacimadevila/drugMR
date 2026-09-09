#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process SMR_BULK {

    tag "smr_bulk_${meta.pheno_id}_${eqtl_dataset}"

    publishDir { "${params.runs_root}/${params.run_id}/results/smr/bulk/${eqtl_dataset}" }, mode: "copy"

    container {
        if (params.image_name) {
            return params.image_name
        } else {
            return 'ghcr.io/guillermocomesanacimadevila/drugmr:latest'
        }
    }

    input:
    tuple val(meta), val(eqtl_dataset), path(qc_tsv), path(mr_tsv, stageAs: "cis_mr/mr.tsv"), path(coloc_tsv, stageAs: "coloc/coloc.tsv"), path(pwcoco_tsv, stageAs: "pwcoco/summary/pwcoco.tsv")

    output:
    tuple val(meta), val("${eqtl_dataset}:bulk"), path("smr/bulk/${eqtl_dataset}/promising_targets.tsv"), emit: smr_bulk_results

    script:
    """
    export PYTHONPATH=${projectDir}
    python ${projectDir}/bin/sort_smr.py \\
        --pheno_id ${meta.pheno_id} \\
        --sumstats ${qc_tsv} \\
        --pqtl_dataset ${meta.pqtl_dataset} \\
        --eqtl_dataset ${eqtl_dataset} \\
        --eqtl_mode bulk \\
        --ref_bfile ${projectDir}/${meta.ref_bfile} \\
        --maf ${meta.maf} \\
        --local_results_dir . \\
        --repo_root ${projectDir} \\
        --synthesis_dir ${projectDir}/synthesis \\
        --skip_merge
    """
}

process SMR_SC {

    tag "smr_sc_${meta.pheno_id}_${meta.sc_eqtl_dataset}"

    publishDir { "${params.runs_root}/${params.run_id}/results/smr/sc/${meta.sc_eqtl_dataset}" }, mode: "copy"

    container {
        if (params.image_name) {
            return params.image_name
        } else {
            return 'ghcr.io/guillermocomesanacimadevila/drugmr:latest'
        }
    }

    input:
    tuple val(meta), path(qc_tsv), path(mr_tsv, stageAs: "cis_mr/mr.tsv"), path(coloc_tsv, stageAs: "coloc/coloc.tsv"), path(pwcoco_tsv, stageAs: "pwcoco/summary/pwcoco.tsv")

    output:
    tuple val(meta), val("${meta.sc_eqtl_dataset}:single_cell"), path("smr/sc/${meta.sc_eqtl_dataset}/promising_targets.tsv"), emit: smr_sc_results

    script:
    """
    export PYTHONPATH=${projectDir}
    python ${projectDir}/bin/sort_smr.py \\
        --pheno_id ${meta.pheno_id} \\
        --sumstats ${qc_tsv} \\
        --pqtl_dataset ${meta.pqtl_dataset} \\
        --eqtl_dataset ${meta.sc_eqtl_dataset} \\
        --eqtl_mode single_cell \\
        --ref_bfile ${projectDir}/${meta.ref_bfile} \\
        --maf ${meta.maf} \\
        --local_results_dir . \\
        --repo_root ${projectDir} \\
        --synthesis_dir ${projectDir}/synthesis \\
        --skip_merge
    """
}

// Fan-in: takes every SMR_BULK/SMR_SC task's own promising_targets.tsv (one
// per eqtl_dataset) for this (pheno_id, pqtl_dataset) run, plus the matching
// "eqtl_dataset:eqtl_mode" label for each (same order - built by the
// subworkflow's .mix().collect(), never reordered in between), and runs
// bin/sort_smr.py's --eqtl_mode merge (merge_multi_omics_targets_batch())
// once over the complete set. See project_nextflow_migration memory for why
// this can't be folded into SMR_BULK/SMR_SC themselves (each only ever sees
// its own one dataset - the union can only be built once every dataset's
// result already exists).
//
// UNVERIFIED: the stageAs "smr_target_?/promising_targets.tsv" pattern below
// disambiguates same-basename files from a collected list by auto-numbering
// each into its own subdirectory - this is the documented Nextflow idiom for
// exactly this problem, but there's no toy eQTL fixture yet to actually run
// this process against, unlike every other module this session. Verify for
// real (a tiny synthetic 2-dataset case, matching the merge_multi_omics_targets_batch
// unit test already run in chat) before trusting this in production.
process MERGE_MULTI_OMICS_TARGETS {

    tag "smr_merge_${meta.pheno_id}_${meta.pqtl_dataset}"

    publishDir { "${params.runs_root}/${params.run_id}/results/smr" }, mode: "copy"

    container {
        if (params.image_name) {
            return params.image_name
        } else {
            return 'ghcr.io/guillermocomesanacimadevila/drugmr:latest'
        }
    }

    input:
    tuple val(meta), val(labels), path(target_files, stageAs: "smr_target_?/promising_targets.tsv")

    output:
    tuple val(meta), path("smr/final_multi_omics_targets.tsv"), emit: smr_final_targets

    script:
    def input_flags = labels.withIndex(1).collect { label, i -> "--input '${label}:smr_target_${i}/promising_targets.tsv'" }.join(' \\\n        ')
    """
    export PYTHONPATH=${projectDir}
    python ${projectDir}/bin/sort_smr.py \\
        --pheno_id ${meta.pheno_id} \\
        --pqtl_dataset ${meta.pqtl_dataset} \\
        --eqtl_mode merge \\
        --local_results_dir . \\
        ${input_flags}
    """
}
