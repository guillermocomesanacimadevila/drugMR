#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

// nextflow run main.nf -params-file tests/nf/params.qc_test.yaml

include { TEST_QC } from './tests/nf/qc_test'

workflow {
    TEST_QC()
}
