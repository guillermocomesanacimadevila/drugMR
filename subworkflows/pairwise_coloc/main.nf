#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include { PAIRWISE_COLOC } from '../../modules/local/coloc.nf'

workflow COLOC {

    take:
    mr_out
    protein_dirs

    main:
    // we need
    // meta -> pqtl_dataset, pheno_id, n_cases, n_controls,
    // protein_dirs -> pqtl_dir, local_results_dir
    //
    // keyed on (pheno_id, pqtl_dataset), not pheno_id alone -> a run is always
    // exactly 1 dataset x 1 trait, and pheno_id-only would wrongly cross-join
    // two different pqtl_datasets sharing the same outcome trait
    ch_in = mr_out
        .map { meta, mr_tsv, mr_instruments_tsv -> tuple([meta.pheno_id, meta.pqtl_dataset], meta, mr_tsv, mr_instruments_tsv) }
        .join(protein_dirs.map { meta, dirs -> tuple([meta.pheno_id, meta.pqtl_dataset], dirs) })
        .map { key, meta, mr_tsv, mr_instruments_tsv, dirs -> tuple(meta, mr_tsv, mr_instruments_tsv, dirs) }

    PAIRWISE_COLOC(ch_in)

    emit:
    coloc_results = PAIRWISE_COLOC.out.coloc_results
}