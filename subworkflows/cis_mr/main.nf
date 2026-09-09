#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include { CIS_MR } from '../../modules/local/cis_mr.nf'

workflow MR_ON_CIS_REGIONS {

    take:
    qc_out
    protein_dirs

    main:

    // phase 1 - tuple("AD", meta, qc_tsv, dirs) -> join on AD
    ch_in = qc_out
        .map { meta, qc_tsv -> tuple(meta.pheno_id, meta, qc_tsv) } // meta, qc_tsv -> AD, meta, qc_tsv
        .join(protein_dirs.map { meta, dirs -> tuple(meta.pheno_id, dirs) }) // so like meta, [dirs], becomes AD, [dirs] -> then joins above...
        .map { pheno_id, meta, qc_tsv, dirs -> tuple(meta, qc_tsv, dirs) } // drop AD -> just useful for join

    CIS_MR(ch_in)

    emit:
    mr_results = CIS_MR.out.mr_results
}