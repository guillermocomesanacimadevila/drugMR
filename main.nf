#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

// nextflow run main.nf -params-file tests/nf/params.qc_test.yaml
// nextflow run main.nf -params-file tests/nf/params.qc_test.yaml --manifest_path tests/dat/qtl_manifest_toy.csv

include { QC_GWAS } from './subworkflows/qc_gwas/main.nf'
include { PREP_CIS_REGIONS } from './subworkflows/prep_cis_regions/main.nf'

workflow {
    QC_GWAS()
    PREP_CIS_REGIONS(QC_GWAS.out.qc_tsv)
}
