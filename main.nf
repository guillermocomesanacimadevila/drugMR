#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    drugMR/drugmr
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Github : https://github.com/guillermocomesanacimadevila/drugMR
----------------------------------------------------------------------------------------
*/

// nextflow run main.nf -profile docker -params-file tests/nf/params.qc_test.yaml --manifest_path tests/dat/qtl_manifest_toy.csv
// nextflow run main.nf -profile docker -params-file tests/nf/params.qc_test.yaml --manifest_path tests/dat/qtl_manifest_toy.csv -resume

include { DRUGMR } from './workflows/drugmr.nf'

workflow {

    main:
    DRUGMR()

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
