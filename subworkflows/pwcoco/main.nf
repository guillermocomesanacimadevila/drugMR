#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include { PWCOCO } from '../../modules/local/pwcoco.nf'

workflow PWCOCO_WF {

    take:
    mr_out        // tuple(meta, mr.tsv, mr_instruments.tsv)
    protein_dirs  // tuple(meta, dirs)

    main:
    // keyed on (pheno_id, pqtl_dataset), not pheno_id alone -> a run is always
    // exactly 1 dataset x 1 trait, and pheno_id-only would wrongly cross-join
    // two different pqtl_datasets sharing the same outcome trait
    ch_in = mr_out
        .map { meta, mr_tsv, mr_instruments_tsv -> tuple([meta.pheno_id, meta.pqtl_dataset], meta, mr_tsv) }
        .join(protein_dirs.map { meta, dirs -> tuple([meta.pheno_id, meta.pqtl_dataset], dirs) })
        .map { key, meta, mr_tsv, dirs -> tuple(meta, mr_tsv, dirs) }

    PWCOCO(ch_in)

    emit:
    pwcoco_results = PWCOCO.out.pwcoco_results
}
