#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include { PWCOCO } from '../../modules/local/pwcoco.nf'

workflow PWCOCO_WF {

    take:
    mr_out        // tuple(meta, mr.tsv, mr_instruments.tsv)
    protein_dirs  // tuple(meta, dirs)

    main:
    ch_in = mr_out
        .map { meta, mr_tsv, mr_instruments_tsv -> tuple(meta.pheno_id, meta, mr_tsv) }
        .join(protein_dirs.map { meta, dirs -> tuple(meta.pheno_id, dirs) })
        .map { pheno_id, meta, mr_tsv, dirs -> tuple(meta, mr_tsv, dirs) }

    PWCOCO(ch_in)

    emit:
    pwcoco_results = PWCOCO.out.pwcoco_results
}
