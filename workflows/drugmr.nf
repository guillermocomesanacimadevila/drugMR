#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include { QC_GWAS } from '../subworkflows/qc_gwas/main.nf'
include { PREP_CIS_REGIONS } from '../subworkflows/prep_cis_regions/main.nf'
include { MR_ON_CIS_REGIONS } from '../subworkflows/cis_mr/main.nf'
include { COLOC } from '../subworkflows/pairwise_coloc/main.nf'
include { PWCOCO_WF } from '../subworkflows/pwcoco/main.nf'
include { TARGET_HITS } from '../subworkflows/target_hits/main.nf'
include { LOCUS_ASSETS_WF } from '../subworkflows/locus_assets/main.nf'
include { SMR } from '../subworkflows/smr/main.nf'
include { PWCOCO_QTL_WF } from '../subworkflows/pwcoco_qtl/main.nf'
include { HYPRCOLOC_WF } from '../subworkflows/hyprcoloc/main.nf'
include { PHEWAS_WF } from '../subworkflows/phewas/main.nf'

workflow DRUGMR {

    main:
    // a run is always exactly 1 dataset x 1 trait (see subworkflows/target_hits)
    if (params.inputs.size() != 1) {
        error(
            "drugMR takes exactly 1 input per run, got ${params.inputs.size()}. " +
            "Each (pheno_id, pqtl_dataset) pair needs its own `nextflow run` " +
            "invocation -> see README's Configuration section."
        )
    }

    QC_GWAS()
    PREP_CIS_REGIONS(QC_GWAS.out.qc_tsv)
    MR_ON_CIS_REGIONS(QC_GWAS.out.qc_tsv, PREP_CIS_REGIONS.out.protein_dirs)
    COLOC(MR_ON_CIS_REGIONS.out.mr_results, PREP_CIS_REGIONS.out.protein_dirs)
    PWCOCO_WF(MR_ON_CIS_REGIONS.out.mr_results, PREP_CIS_REGIONS.out.protein_dirs)
    TARGET_HITS(COLOC.out.coloc_results, PWCOCO_WF.out.pwcoco_results, PREP_CIS_REGIONS.out.protein_dirs)
    LOCUS_ASSETS_WF(TARGET_HITS.out.target_stats, PREP_CIS_REGIONS.out.protein_dirs)
    SMR(QC_GWAS.out.qc_tsv, MR_ON_CIS_REGIONS.out.mr_results, COLOC.out.coloc_results, PWCOCO_WF.out.pwcoco_results)
    PWCOCO_QTL_WF(SMR.out.smr_final_targets, PWCOCO_WF.out.pwcoco_results, PREP_CIS_REGIONS.out.protein_dirs)
    HYPRCOLOC_WF(SMR.out.smr_final_targets, PREP_CIS_REGIONS.out.protein_dirs)
    PHEWAS_WF(COLOC.out.coloc_results, PWCOCO_WF.out.pwcoco_results, MR_ON_CIS_REGIONS.out.mr_results, PREP_CIS_REGIONS.out.protein_dirs)
    // postgres loading + dashboard serving are deliberately NOT a nextflow stage
    // see dm.results() in drugmr/local.py; run that manually
    // afterwards, pulling runs/<run_id> back first if this ran remotely
}
