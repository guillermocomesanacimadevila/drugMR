#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process PAIRWISE_COLOC {

    tag "pairwise_coloc"

    publishDir { "${params.runs_root}/${params.run_id}/results/coloc" }, mode: "copy"

    container {
        if (params.image_name) {
            return params.image_name
        } else {
            return 'ghcr.io/guillermocomesanacimadevila/drugmr:latest'
        }
    }

    input:
    tuple val(meta), path(mr_tsv, stageAs: "cis_mr/mr.tsv"), path(mr_instruments_tsv, stageAs: "cis_mr/instruments/mr_instruments.tsv"), path(protein_dirs, stageAs: "protein_dirs/*")

    output:
    tuple val(meta), path("coloc/${meta.pqtl_dataset}_${meta.pheno_id}_all_coloc.tsv"), path("coloc/coloc_sensitivity.tsv"), emit: coloc_results

    script:
    """
    export PYTHONPATH=${projectDir}
    python ${projectDir}/bin/coloc_targets.py \\
        --pqtl_dataset ${meta.pqtl_dataset} \\
        --local_results_dir . \\
        --pqtl_dir protein_dirs \\
        --pheno_id ${meta.pheno_id} \\
        --n_cases ${meta.n_cases} \\
        --n_controls ${meta.n_controls}
    """

}


/*
cmd_coloc = [
        "docker", "run", "--rm",
        "-v", f"{project_root}:/work",
        "-w", "/work",
        "-e", "PYTHONPATH=.",
        image_name,
        "python", "bin/coloc_targets.py",
        "--pqtl_dataset", pqtl_dataset, - meta
        "--local_results_dir", out_dir, - **
        "--pqtl_dir", f"dat/cis_regions/{pqtl_dataset}",
        "--pheno_id", pheno_id, - meta
        "--n_cases", str(n_cases), - meta
        "--n_controls", str(n_controls), - meta
        "--wald_fdr_q", str(wald_fdr_q), - mr_res
        "--ivw_fdr_q", str(ivw_fdr_q), - mr_res
        "--cochran_q_pval", str(cochran_q_pval),
        "--egger_intercept_pval_min", str(egger_intercept_pval_min), - mr_res
        "--min_instruments_for_ivw", str(min_instruments_for_ivw), - mr_res
    ]
*/
