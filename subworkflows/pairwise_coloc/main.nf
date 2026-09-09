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
    ch_in = mr_out
        .map { meta, mr_tsv, mr_instruments_tsv -> tuple(meta.pheno_id, meta, mr_tsv, mr_instruments_tsv) }
        .join(protein_dirs.map { meta, dirs -> tuple(meta.pheno_id, dirs) })
        .map { pheno_id, meta, mr_tsv, mr_instruments_tsv, dirs -> tuple(meta, mr_tsv, mr_instruments_tsv, dirs) }

    PAIRWISE_COLOC(ch_in)

    emit:
    coloc_results = PAIRWISE_COLOC.out.coloc_results
}