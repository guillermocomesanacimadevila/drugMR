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
    } else {
        def pheno_id = params.inputs[0].pheno_id
        def pqtl_dataset = params.inputs[0].pqtl_dataset
        def host
        try {
            host = InetAddress.getLocalHost().getHostName()
        } catch (Exception e) {
            host = "unknown"
        }
        def status = workflow.success ? "success" : "failed"
        // git state read once at launch in nextflow.config (the same value run_id
        // uses); workflow.commitId is only set when run from a remote repository
        def git_sha7 = workflow.commitId ? workflow.commitId.take(7) : "${params.git_sha7}"
        def git_dirty = "${params.git_dirty}"
        // every process runs params.image_name (workflow.container is an empty map
        // when containers are set per process, so it cannot be used here)
        def container = "${params.image_name}"
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
            "--git_sha7", git_sha7,
            "--git_dirty", git_dirty,
            "--image_name", "${params.image_name}",
            "--container", container,
            "--host", host,
            "--status", status,
            "--params_json", groovy.json.JsonOutput.toJson(params.inputs[0] + [manifest_path: params.manifest_path])
        ]
        def proc = cmd.execute(["PYTHONPATH=${projectDir}"], new File("${projectDir}"))
        proc.waitFor()
        println proc.in.text
        if (proc.exitValue() != 0) {
            println "[CONCERN] Failed to record run in registry: ${proc.err.text}"
        } else if (!workflow.success) {
            println "[CONCERN] Pipeline did not complete successfully, run recorded as failed for run_id=${params.run_id}"
        }
    }
}
