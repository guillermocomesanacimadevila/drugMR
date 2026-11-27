#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include { QC_GWAS } from '../subworkflows/qc_gwas/main.nf'
include { PREP_CIS_REGIONS } from '../subworkflows/prep_cis_regions/main.nf'
include { MR_ON_CIS_REGIONS } from '../subworkflows/cis_mr/main.nf'
include { COLOC } from '../subworkflows/pairwise_coloc/main.nf'
include { PWCOCO_WF } from '../subworkflows/pwcoco/main.nf'
include { LOCUS_CONDITIONING_WF } from '../subworkflows/locus_conditioning/main.nf'
// include { PAV_CHECKS_WF } from '../subworkflows/pav_sensitivity/main.nf'
include { TARGET_HITS } from '../subworkflows/target_hits/main.nf'
include { LOCUS_ASSETS_WF } from '../subworkflows/locus_assets/main.nf'
include { SMR } from '../subworkflows/smr/main.nf'
include { PWCOCO_QTL_WF } from '../subworkflows/pwcoco_qtl/main.nf'
include { HYPRCOLOC_WF } from '../subworkflows/hyprcoloc/main.nf'
include { PHEWAS_WF } from '../subworkflows/phewas/main.nf'
include { DRUGGABILITY_WF } from '../subworkflows/druggability/main.nf'


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
    LOCUS_CONDITIONING_WF(MR_ON_CIS_REGIONS.out.mr_results, COLOC.out.coloc_results, PWCOCO_WF.out.pwcoco_results, PREP_CIS_REGIONS.out.protein_dirs)
    TARGET_HITS(COLOC.out.coloc_results, PWCOCO_WF.out.pwcoco_results, PREP_CIS_REGIONS.out.protein_dirs)
    LOCUS_ASSETS_WF(TARGET_HITS.out.target_stats, PREP_CIS_REGIONS.out.protein_dirs)
    SMR(QC_GWAS.out.qc_tsv, MR_ON_CIS_REGIONS.out.mr_results, COLOC.out.coloc_results, PWCOCO_WF.out.pwcoco_results)
    PWCOCO_QTL_WF(SMR.out.smr_final_targets, PWCOCO_WF.out.pwcoco_results, PREP_CIS_REGIONS.out.protein_dirs)
    HYPRCOLOC_WF(SMR.out.smr_final_targets, PREP_CIS_REGIONS.out.protein_dirs)
    PHEWAS_WF(COLOC.out.coloc_results, PWCOCO_WF.out.pwcoco_results, MR_ON_CIS_REGIONS.out.mr_results, PREP_CIS_REGIONS.out.protein_dirs)

    // druggability is the final step: it waits for every other stage to finish.
    // mix + collect rather than join, because some of these outputs are optional
    // and may never arrive. ifEmpty keeps the signal alive if all are empty
    ch_all_done = HYPRCOLOC_WF.out.hyprcoloc_final
        .mix(PHEWAS_WF.out.phewas_finngen_results)
        .mix(PHEWAS_WF.out.phewas_ukb_results)
        .mix(PWCOCO_QTL_WF.out.pwcoco_qtl_shared_results)
        .mix(LOCUS_ASSETS_WF.out.locus_assets)
        .mix(LOCUS_CONDITIONING_WF.out.locus_conditioning_results)
        .collect()
        .ifEmpty(["none"])

    // needs internet (Pharos, ChEMBL, Open Targets, GitHub), so it can be
    // switched off on offline compute nodes with --skip_druggability true
    // values given on the command line arrive as text ("false"), and any
    // non-empty text counts as true, so compare the text instead
    def skip_druggability = params.skip_druggability.toString().toLowerCase() == "true"
    if (!skip_druggability) {
        DRUGGABILITY_WF(MR_ON_CIS_REGIONS.out.mr_results, COLOC.out.coloc_results, PWCOCO_WF.out.pwcoco_results, ch_all_done)
    }
    // postgres loading + dashboard serving are deliberately NOT a nextflow stage
    // see dm.results() in drugmr/runners/local.py; run that manually
    // afterwards, pulling runs/<run_id> back first if this ran remotely
}
