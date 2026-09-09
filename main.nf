#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

/*
ADD NF REPORT
*/

// nextflow run main.nf -params-file tests/nf/params.qc_test.yaml --manifest_path tests/dat/qtl_manifest_toy.csv
// nextflow run main.nf -params-file tests/nf/params.qc_test.yaml --manifest_path tests/dat/qtl_manifest_toy.csv -resume

include { QC_GWAS } from './subworkflows/qc_gwas/main.nf'
include { PREP_CIS_REGIONS } from './subworkflows/prep_cis_regions/main.nf'
include { MR_ON_CIS_REGIONS } from './subworkflows/cis_mr/main.nf'
include { COLOC } from './subworkflows/pairwise_coloc/main.nf'
include { PWCOCO_WF } from './subworkflows/pwcoco/main.nf'

workflow {

    log.info("Welcome to the drugMR pipeline!")

    QC_GWAS()
    PREP_CIS_REGIONS(QC_GWAS.out.qc_tsv)
    MR_ON_CIS_REGIONS(QC_GWAS.out.qc_tsv, PREP_CIS_REGIONS.out.protein_dirs)
    COLOC(MR_ON_CIS_REGIONS.out.mr_results, PREP_CIS_REGIONS.out.protein_dirs)
    PWCOCO_WF(MR_ON_CIS_REGIONS.out.mr_results, PREP_CIS_REGIONS.out.protein_dirs)
}
