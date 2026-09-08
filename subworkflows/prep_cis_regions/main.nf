#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include { EXTRACT_CIS_REGIONS } from '../../modules/local/cis_regions.nf'

/*
We need the following...
- pheno_id (from meta)
- pqtl_dataset (from meta)
- qtl manifest
- the output from bloody qc_gwas.py
*/

workflow PREP_CIS_REGIONS {

    take:
    qc_out

    main:
    ch_in = qc_out.map { meta, qc_tsv -> 
        tuple(meta, qc_tsv, file(params.manifest_path)) 
    }

    EXTRACT_CIS_REGIONS(ch_in)

    emit:
    protein_dirs = EXTRACT_CIS_REGIONS.out.protein_dirs
}