#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include { DRUGGABILITY } from '../../modules/local/druggability.nf'

workflow DRUGGABILITY_WF {

    take:
    mr_results      // tuple(meta, mr.tsv, mr_instruments.tsv)
    coloc_results   // tuple(meta, coloc.tsv, coloc_sensitivity.tsv)
    pwcoco_results  // tuple(meta, pwcoco.tsv)
    all_done        // one list holding every final output of the run

    main:
    // keyed on (pheno_id, pqtl_dataset), same composite-key discipline as phewas/target_hits
    ch_mr = mr_results.map { meta, mr_tsv, mr_instruments_tsv -> tuple([meta.pheno_id, meta.pqtl_dataset], meta, mr_tsv) }
    ch_coloc = coloc_results.map { meta, coloc_tsv, coloc_sensitivity_tsv -> tuple([meta.pheno_id, meta.pqtl_dataset], coloc_tsv) }
    ch_pwcoco = pwcoco_results.map { meta, pwcoco_tsv -> tuple([meta.pheno_id, meta.pqtl_dataset], pwcoco_tsv) }

    // all_done only carries the order: druggability waits until every other
    // stage has finished. Mapped to a plain string so the collected list of
    // files does not spread into the tuple
    ch_done = all_done.map { files -> "done" }

    ch_in = ch_mr
        .join(ch_coloc)
        .join(ch_pwcoco)
        .combine(ch_done)
        .map { key, meta, mr_tsv, coloc_tsv, pwcoco_tsv, done -> tuple(meta, mr_tsv, coloc_tsv, pwcoco_tsv) }

    DRUGGABILITY(ch_in)

    emit:
    druggability_results = DRUGGABILITY.out.druggability_results
}
