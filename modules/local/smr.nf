#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process SMR_BULK {

    tag "smr_bulk_${meta.pheno_id}_${qtl_dataset}"
    label "process_medium"

    // publishDir { "${params.runs_root}/${params.run_id}/results/smr/bulk/${qtl_dataset}" }, mode: "copy"
    publishDir { "${params.runs_root}/${params.run_id}/results" }, mode: "copy"

    container {
        if (params.image_name) {
            return params.image_name
        } else {
            return 'ghcr.io/guillermocomesanacimadevila/drugmr:latest'
        }
    }

    input:
    tuple val(meta), val(qtl_dataset), path(qc_tsv), path(mr_tsv, stageAs: "cis_mr/mr.tsv"), path(coloc_tsv, stageAs: "coloc/coloc.tsv"), path(pwcoco_tsv, stageAs: "pwcoco/summary/pwcoco.tsv")

    output:
    tuple val(meta), val("${qtl_dataset}:bulk"), path("smr/bulk/${qtl_dataset}/promising_targets.tsv"), emit: smr_bulk_results

    // no explicit --manifest_path here: sort_smr.py's own SMRUtils default
    // ("assets/qtl_manifest.csv", a relative path) combined with --repo_root
    // below already resolves to the same real file - only modules whose
    // underlying script has a different or no useful default need to pin it
    // explicitly (see hyprcoloc.nf/pwcoco_qtl.nf).
    script:
    """
    export PYTHONPATH=${projectDir}
    python ${projectDir}/bin/sort_smr.py \\
        --pheno_id ${meta.pheno_id} \\
        --sumstats ${qc_tsv} \\
        --pqtl_dataset ${meta.pqtl_dataset} \\
        --qtl_dataset ${qtl_dataset} \\
        --qtl_mode bulk \\
        --ref_bfile ${projectDir}/${meta.ref_bfile} \\
        --gene_annotation ${meta.gene_annotation ?: ''} \\
        --maf ${meta.maf} \\
        --local_results_dir . \\
        --repo_root ${projectDir} \\
        --synthesis_dir ${projectDir}/synthesis \\
        --coloc_file coloc/coloc.tsv \\
        --wald_fdr_q ${meta.gates.cis_mr.wald_fdr_q} \\
        --ivw_fdr_q ${meta.gates.cis_mr.ivw_fdr_q} \\
        --cochran_q_pval ${meta.gates.cis_mr.cochran_q_pval} \\
        --p_qtl_smr ${meta.gates.smr.p_qtl_smr} \\
        --p_qtl_heidi ${meta.gates.smr.p_qtl_heidi} \\
        --skip_merge
    """
}

process SMR_SC {

    tag "smr_sc_${meta.pheno_id}_${meta.sc_qtl_dataset}"
    label "process_medium"

    // publishDir { "${params.runs_root}/${params.run_id}/results/smr/sc/${meta.sc_qtl_dataset}" }, mode: "copy"
    publishDir { "${params.runs_root}/${params.run_id}/results" }, mode: "copy"

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
    tuple val(meta), val("${meta.sc_qtl_dataset}:single_cell"), path("smr/sc/${meta.sc_qtl_dataset}/promising_targets.tsv"), emit: smr_sc_results

    // no explicit --manifest_path here: sort_smr.py's own SMRUtils default
    // ("assets/qtl_manifest.csv", a relative path) combined with --repo_root
    // below already resolves to the same real file - only modules whose
    // underlying script has a different or no useful default need to pin it
    // explicitly (see hyprcoloc.nf/pwcoco_qtl.nf).
    script:
    """
    export PYTHONPATH=${projectDir}
    python ${projectDir}/bin/sort_smr.py \\
        --pheno_id ${meta.pheno_id} \\
        --sumstats ${qc_tsv} \\
        --pqtl_dataset ${meta.pqtl_dataset} \\
        --qtl_dataset ${meta.sc_qtl_dataset} \\
        --qtl_mode single_cell \\
        --ref_bfile ${projectDir}/${meta.ref_bfile} \\
        --gene_annotation ${meta.gene_annotation ?: ''} \\
        --maf ${meta.maf} \\
        --local_results_dir . \\
        --repo_root ${projectDir} \\
        --synthesis_dir ${projectDir}/synthesis \\
        --coloc_file coloc/coloc.tsv \\
        --wald_fdr_q ${meta.gates.cis_mr.wald_fdr_q} \\
        --ivw_fdr_q ${meta.gates.cis_mr.ivw_fdr_q} \\
        --cochran_q_pval ${meta.gates.cis_mr.cochran_q_pval} \\
        --p_qtl_smr ${meta.gates.smr.p_qtl_smr} \\
        --p_qtl_heidi ${meta.gates.smr.p_qtl_heidi} \\
        --skip_merge
    """
}

// Fan-in: takes every SMR_BULK/SMR_SC task's own promising_targets.tsv (one
// per qtl_dataset) for this (pheno_id, pqtl_dataset) run, plus the matching
// "qtl_dataset:qtl_mode" label for each (same order - built by the
// subworkflow's .mix().collect(), never reordered in between), and runs
// bin/sort_smr.py's --qtl_mode merge (merge_multi_omics_targets_batch())
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
    label "process_single"

    // publishDir { "${params.runs_root}/${params.run_id}/results/smr" }, mode: "copy"
    publishDir { "${params.runs_root}/${params.run_id}/results" }, mode: "copy"

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
        --qtl_mode merge \\
        --local_results_dir . \\
        --p_smr_threshold ${meta.gates.smr.p_smr_threshold} \\
        --p_heidi_threshold ${meta.gates.smr.p_heidi_threshold} \\
        ${input_flags}
    """
}
