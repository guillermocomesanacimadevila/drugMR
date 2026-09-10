#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include { CIS_MR } from '../../modules/local/cis_mr.nf'

workflow MR_ON_CIS_REGIONS {

    take:
    qc_out
    protein_dirs

    main:

    // keyed on (pheno_id, pqtl_dataset), not pheno_id alone -> a run is always
    // exactly 1 dataset x 1 trait, and pheno_id-only would wrongly cross-join
    // two different pqtl_datasets sharing the same outcome trait
    ch_in = qc_out
        .map { meta, qc_tsv -> tuple([meta.pheno_id, meta.pqtl_dataset], meta, qc_tsv) }
        .join(protein_dirs.map { meta, dirs -> tuple([meta.pheno_id, meta.pqtl_dataset], dirs) })
        .map { key, meta, qc_tsv, dirs -> tuple(meta, qc_tsv, dirs) }

    CIS_MR(ch_in)

    emit:
    mr_results = CIS_MR.out.mr_results
}