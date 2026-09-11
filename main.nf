#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    drugMR/drugmr
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Github : https://github.com/guillermocomesanacimadevila/drugMR
----------------------------------------------------------------------------------------
*/

include { DRUGMR } from './workflows/drugmr.nf'

workflow {

    main:
    DRUGMR()

    onComplete:
    if (workflow.preview) {
        println "[TRACKING] -preview run - registry NOT updated for run_id=${params.run_id}"
    } else if (workflow.success) {
        def pheno_id = params.inputs[0].pheno_id
        def pqtl_dataset = params.inputs[0].pqtl_dataset
        def host
        try {
            host = InetAddress.getLocalHost().getHostName()
        } catch (Exception e) {
            host = "unknown"
        }
        // Run bookkeeping with the same pinned Python environment validated by
        // env/activate.sh.  Groovy's execute(envp, dir) does not reliably retain
        // the activated shell PATH on HPC, so a bare `python3` can resolve to a
        // system interpreter without PyYAML or the rest of drugMR's dependencies.
        def projectPython = "${projectDir}/.venv/bin/python"
        def cmd = [
            projectPython, "${projectDir}/bin/record_run.py",
            "--pheno_id", pheno_id,
            "--pqtl_dataset", pqtl_dataset,
            "--run_id", params.run_id,
            "--root", "${params.runs_root}",
            "--git_sha7", "${workflow.commitId ?: 'unknown'}",
            "--image_name", "${params.image_name}",
            "--host", host,
            "--params_json", groovy.json.JsonOutput.toJson(params.inputs[0])
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
