#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include { EXTRACT_TWO_OMICS_HITS } from '../../modules/local/target_hits.nf'

workflow TARGET_HITS {

    take:
    coloc_results   // tuple(meta, coloc.tsv, coloc_sensitivity.tsv)
    pwcoco_results  // tuple(meta, pwcoco.tsv)
    protein_dirs    // tuple(meta, dirs)

    main:
    // keyed on (pheno_id, pqtl_dataset), not pheno_id alone -> a run is always
    // exactly 1 dataset x 1 trait, and pheno_id-only would wrongly cross-join
    // two different pqtl_datasets sharing the same outcome trait
    ch_coloc = coloc_results.map { meta, coloc_tsv, coloc_sensitivity_tsv -> tuple([meta.pheno_id, meta.pqtl_dataset], meta, coloc_tsv) }
    ch_pwcoco = pwcoco_results.map { meta, pwcoco_tsv -> tuple([meta.pheno_id, meta.pqtl_dataset], pwcoco_tsv) }
    ch_dirs = protein_dirs.map { meta, dirs -> tuple([meta.pheno_id, meta.pqtl_dataset], dirs) }

    ch_in = ch_coloc
        .join(ch_pwcoco)
        .join(ch_dirs)
        .map { key, meta, coloc_tsv, pwcoco_tsv, dirs -> tuple(meta, coloc_tsv, pwcoco_tsv, dirs) }

    EXTRACT_TWO_OMICS_HITS(ch_in)

    emit:
    target_stats = EXTRACT_TWO_OMICS_HITS.out.target_stats
}
