#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include { SMR_BULK; SMR_SC; MERGE_MULTI_OMICS_TARGETS } from '../../modules/local/smr.nf'

workflow SMR {

    take:
    qc_out          // tuple(meta, qc_tsv)
    mr_out          // tuple(meta, mr.tsv, mr_instruments.tsv)
    coloc_results   // tuple(meta, coloc.tsv, coloc_sensitivity.tsv)
    pwcoco_results  // tuple(meta, pwcoco.tsv)

    main:
    // keyed on (pheno_id, pqtl_dataset), same composite-key discipline as
    ch_qc = qc_out.map { meta, qc_tsv -> tuple([meta.pheno_id, meta.pqtl_dataset], meta, qc_tsv) }
    ch_mr = mr_out.map { meta, mr_tsv, mr_instruments_tsv -> tuple([meta.pheno_id, meta.pqtl_dataset], mr_tsv) }
    ch_coloc = coloc_results.map { meta, coloc_tsv, coloc_sensitivity_tsv -> tuple([meta.pheno_id, meta.pqtl_dataset], coloc_tsv) }
    ch_pwcoco = pwcoco_results.map { meta, pwcoco_tsv -> tuple([meta.pheno_id, meta.pqtl_dataset], pwcoco_tsv) }

    ch_joined = ch_qc
        .join(ch_mr)
        .join(ch_coloc)
        .join(ch_pwcoco)
        .map { key, meta, qc_tsv, mr_tsv, coloc_tsv, pwcoco_tsv -> tuple(meta, qc_tsv, mr_tsv, coloc_tsv, pwcoco_tsv) }

    // bulk: one task per entry in meta.bulk_eqtl_datasets, matching local.py's
    ch_bulk_in = ch_joined
        .filter { meta, qc_tsv, mr_tsv, coloc_tsv, pwcoco_tsv -> meta.run_smr && meta.bulk_eqtl_datasets }
        .flatMap { meta, qc_tsv, mr_tsv, coloc_tsv, pwcoco_tsv ->
            meta.bulk_eqtl_datasets.collect { ds -> tuple(meta, ds, qc_tsv, mr_tsv, coloc_tsv, pwcoco_tsv) }
        }

    SMR_BULK(ch_bulk_in)

    // single-cell: at most 1 task, only when meta.sc_eqtl_dataset is actually set
    ch_sc_in = ch_joined
        .filter { meta, qc_tsv, mr_tsv, coloc_tsv, pwcoco_tsv -> meta.run_smr && meta.sc_eqtl_dataset }

    SMR_SC(ch_sc_in)

    ch_merge_in = SMR_BULK.out.smr_bulk_results
        .mix(SMR_SC.out.smr_sc_results)
        .map { meta, label, file -> tuple([meta.pheno_id, meta.pqtl_dataset], meta, label, file) }
        .groupTuple(by: 0)
        .map { key, metas, labels, files -> tuple(metas[0], labels, files) }

    MERGE_MULTI_OMICS_TARGETS(ch_merge_in)

    emit:
    smr_bulk_results  = SMR_BULK.out.smr_bulk_results
    smr_sc_results    = SMR_SC.out.smr_sc_results
    smr_final_targets = MERGE_MULTI_OMICS_TARGETS.out.smr_final_targets
}
