#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

/*
ADD NF REPORT
*/

// nextflow run main.nf -profile docker -params-file tests/nf/params.qc_test.yaml --manifest_path tests/dat/qtl_manifest_toy.csv
// nextflow run main.nf -profile docker -params-file tests/nf/params.qc_test.yaml --manifest_path tests/dat/qtl_manifest_toy.csv -resume

include { QC_GWAS } from './subworkflows/qc_gwas/main.nf'
include { PREP_CIS_REGIONS } from './subworkflows/prep_cis_regions/main.nf'
include { MR_ON_CIS_REGIONS } from './subworkflows/cis_mr/main.nf'
include { COLOC } from './subworkflows/pairwise_coloc/main.nf'
include { PWCOCO_WF } from './subworkflows/pwcoco/main.nf'
include { TARGET_HITS } from './subworkflows/target_hits/main.nf'
include { SMR } from './subworkflows/smr/main.nf'
include { PWCOCO_QTL_WF } from './subworkflows/pwcoco_qtl/main.nf'
include { HYPRCOLOC_WF } from './subworkflows/hyprcoloc/main.nf'
include { PHEWAS_WF } from './subworkflows/phewas/main.nf'

workflow {

    main:

    QC_GWAS()
    PREP_CIS_REGIONS(QC_GWAS.out.qc_tsv)
    MR_ON_CIS_REGIONS(QC_GWAS.out.qc_tsv, PREP_CIS_REGIONS.out.protein_dirs)
    COLOC(MR_ON_CIS_REGIONS.out.mr_results, PREP_CIS_REGIONS.out.protein_dirs)
    PWCOCO_WF(MR_ON_CIS_REGIONS.out.mr_results, PREP_CIS_REGIONS.out.protein_dirs)
    TARGET_HITS(COLOC.out.coloc_results, PWCOCO_WF.out.pwcoco_results, PREP_CIS_REGIONS.out.protein_dirs)
    SMR(QC_GWAS.out.qc_tsv, MR_ON_CIS_REGIONS.out.mr_results, COLOC.out.coloc_results, PWCOCO_WF.out.pwcoco_results)
    PWCOCO_QTL_WF(SMR.out.smr_final_targets, PWCOCO_WF.out.pwcoco_results, PREP_CIS_REGIONS.out.protein_dirs)
    HYPRCOLOC_WF(SMR.out.smr_final_targets, PREP_CIS_REGIONS.out.protein_dirs)
    PHEWAS_WF(COLOC.out.coloc_results, PWCOCO_WF.out.pwcoco_results, MR_ON_CIS_REGIONS.out.mr_results, PREP_CIS_REGIONS.out.protein_dirs)
    // postgres loading + dashboard serving are deliberately NOT a nextflow
    // stage - see dm.results() in drugmr/local.py; run that manually
    // afterwards, pulling runs/<run_id> back first if this ran remotely

    onComplete:
    if (workflow.success) {
        def pheno_id = params.inputs[0].pheno_id
        def pqtl_dataset = params.inputs[0].pqtl_dataset
        def cmd = [
            "python3", "${projectDir}/bin/record_run.py",
            "--pheno_id", pheno_id,
            "--pqtl_dataset", pqtl_dataset,
            "--run_id", params.run_id,
            "--root", "${params.runs_root}"
        ]
        def proc = cmd.execute(["PYTHONPATH=${projectDir}"], new File("${projectDir}"))
        proc.waitFor()
        println proc.in.text
        if (proc.exitValue() != 0) {
            println "[CONCERN] Failed to record run in registry: ${proc.err.text}"
        }
    } else {
        println "[CONCERN] Pipeline did not complete successfully, registry NOT updated for run_id=${params.run_id}"
    }
}
