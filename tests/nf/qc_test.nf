#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

// run via the repo-root launcher so ${projectDir} inside GWAS_QC resolves correctly:

include { GWAS_QC } from '../../modules/local/qc_gwas.nf'

workflow TEST_QC {
    ch_in = Channel
        .fromList(params.inputs)
        .map { row -> 
            
            // defs for optional params
            
            // maf
            def maf_def = 0.01
            def maf

            // mhc
            def remove_mhc_def = true
            def remove_mhc

            // apoe
            def remove_apoe_def = false
            def remove_apoe

            if (row.maf) {
                maf = row.maf
            } else {
                maf = maf_def
            }

            if (row.containsKey("remove_mhc")) {
                remove_mhc = row.remove_mhc
            } else {
                remove_mhc = remove_mhc_def
            }

            if (row.containsKey("remove_apoe")) {
                remove_apoe = row.remove_apoe
            } else {
                remove_apoe = remove_apoe_def
            }
            
            // define meta tuple as input channel for QC module
            def meta = [
                pheno_id      : row.pheno_id,
                sumstats      : file(row.sumstats),
                maf           : maf,
                info_col      : row.info_col,
                info_threshold: row.info_threshold,
                snp_col       : row.snp_col,
                a1_col        : row.a1_col,
                a2_col        : row.a2_col,
                beta_col      : row.beta_col,
                se_col        : row.se_col,
                p_col         : row.p_col,
                pos_col       : row.pos_col,
                chr_col       : row.chr_col,
                af_col        : row.af_col,
                genome_build  : row.genome_build,
                target_build  : row.target_build,
                n_cases       : row.n_cases,
                n_controls    : row.n_controls,
                remove_mhc    : remove_mhc,
                remove_apoe   : remove_apoe,
            ]
        }

    // now run QC module
    GWAS_QC(ch_in)
}