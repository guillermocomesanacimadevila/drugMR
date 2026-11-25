#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include { LOCUS_CONDITIONING } from '../../modules/local/locus_conditioning.nf'

workflow LOCUS_CONDITIONING_WF {

    take:
    mr_results      // tuple(meta, mr.tsv, mr_instruments.tsv)
    coloc_results   // tuple(meta, coloc.tsv, coloc_sensitivity.tsv)
    pwcoco_results  // tuple(meta, pwcoco.tsv)
    protein_dirs    // tuple(meta, dirs)

    main:
    // keyed on (pheno_id, pqtl_dataset), same composite-key discipline as phewas/target_hits
    ch_mr = mr_results.map { meta, mr_tsv, mr_instruments_tsv -> tuple([meta.pheno_id, meta.pqtl_dataset], meta, mr_tsv, mr_instruments_tsv) }
    ch_coloc = coloc_results.map { meta, coloc_tsv, coloc_sensitivity_tsv -> tuple([meta.pheno_id, meta.pqtl_dataset], coloc_tsv) }
    ch_pwcoco = pwcoco_results.map { meta, pwcoco_tsv -> tuple([meta.pheno_id, meta.pqtl_dataset], pwcoco_tsv) }
    ch_dirs = protein_dirs.map { meta, dirs -> tuple([meta.pheno_id, meta.pqtl_dataset], dirs) }

    ch_in = ch_mr
        .join(ch_coloc)
        .join(ch_pwcoco)
        .join(ch_dirs)
        .map { key, meta, mr_tsv, mr_instruments_tsv, coloc_tsv, pwcoco_tsv, dirs -> tuple(meta, mr_tsv, mr_instruments_tsv, coloc_tsv, pwcoco_tsv, dirs) }

    LOCUS_CONDITIONING(ch_in)

    emit:
    locus_conditioning_results = LOCUS_CONDITIONING.out.locus_conditioning_results
}
