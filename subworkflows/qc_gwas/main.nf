#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

// run via the repo-root launcher so ${projectDir} inside GWAS_QC resolves correctly:

include { GWAS_QC } from '../../modules/local/qc_gwas.nf'

workflow QC_GWAS {

    main:
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

            // SMR (matches params/schema.json's own default: true)
            def run_smr_def = true
            def run_smr

            // bulk_qtl_datasets / sc_qtl_dataset are genuinely optional -
            // no bulk fan-out / no single-cell task at all when unset,
            // matching local.py's own "no bulk_qtl_datasets specified,
            // skipping bulk SMR" / "no sc_qtl_dataset specified" behaviour
            def bulk_qtl_datasets_def = []
            def bulk_qtl_datasets

            def sc_qtl_dataset_def = null
            def sc_qtl_dataset

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

            if (row.containsKey("run_smr")) {
                run_smr = row.run_smr
            } else {
                run_smr = run_smr_def
            }

            if (row.containsKey("bulk_qtl_datasets")) {
                bulk_qtl_datasets = row.bulk_qtl_datasets
            } else {
                bulk_qtl_datasets = bulk_qtl_datasets_def
            }

            if (row.containsKey("sc_qtl_dataset")) {
                sc_qtl_dataset = row.sc_qtl_dataset
            } else {
                sc_qtl_dataset = sc_qtl_dataset_def
            }

            def gates_cis_mr
            def gates_coloc
            def gates_smr
            def gates_hyprcoloc
            def gates_pwcoco
            def gates_phewas

            if (row.containsKey("gates") && row.gates.containsKey("cis_mr")) {
                gates_cis_mr = row.gates.cis_mr
            } else {
                gates_cis_mr = [:]
            }

            if (row.containsKey("gates") && row.gates.containsKey("coloc")) {
                gates_coloc = row.gates.coloc
            } else {
                gates_coloc = [:]
            }

            if (row.containsKey("gates") && row.gates.containsKey("smr")) {
                gates_smr = row.gates.smr
            } else {
                gates_smr = [:]
            }

            if (row.containsKey("gates") && row.gates.containsKey("hyprcoloc")) {
                gates_hyprcoloc = row.gates.hyprcoloc
            } else {
                gates_hyprcoloc = [:]
            }

            if (row.containsKey("gates") && row.gates.containsKey("pwcoco")) {
                gates_pwcoco = row.gates.pwcoco
            } else {
                gates_pwcoco = [:]
            }

            if (row.containsKey("gates") && row.gates.containsKey("phewas")) {
                gates_phewas = row.gates.phewas
            } else {
                gates_phewas = [:]
            }

            def wald_fdr_q_gate
            def ivw_fdr_q_gate
            def cochran_q_pval_gate
            def egger_intercept_pval_min_gate
            def min_instruments_for_ivw_gate
            def apply_steiger_filter_gate
            def clump_kb_gate
            def clump_r2_gate
            def instrument_pval_threshold_gate
            def min_f_stat_gate
            def pp4_threshold_gate
            def p1_gate
            def p2_gate
            def p12_gate
            def p_qtl_smr_gate
            def p_qtl_heidi_gate
            def p_smr_threshold_gate
            def p_heidi_threshold_gate
            def prior_1_gate
            def prior_c_gate
            def reg_thresh_gate
            def align_thresh_gate
            def equal_thresholds_gate
            def pwcoco_pp4_threshold_gate
            def bonferroni_alpha_gate
            def coloc_threshold_gate

            if (gates_cis_mr.containsKey("wald_fdr_q")) {
                wald_fdr_q_gate = gates_cis_mr.wald_fdr_q
            } else {
                wald_fdr_q_gate = 0.05
            }

            if (gates_cis_mr.containsKey("ivw_fdr_q")) {
                ivw_fdr_q_gate = gates_cis_mr.ivw_fdr_q
            } else {
                ivw_fdr_q_gate = 0.05
            }

            if (gates_cis_mr.containsKey("cochran_q_pval")) {
                cochran_q_pval_gate = gates_cis_mr.cochran_q_pval
            } else {
                cochran_q_pval_gate = 0.05
            }

            if (gates_cis_mr.containsKey("egger_intercept_pval_min")) {
                egger_intercept_pval_min_gate = gates_cis_mr.egger_intercept_pval_min
            } else {
                egger_intercept_pval_min_gate = 0
            }

            if (gates_cis_mr.containsKey("min_instruments_for_ivw")) {
                min_instruments_for_ivw_gate = gates_cis_mr.min_instruments_for_ivw
            } else {
                min_instruments_for_ivw_gate = 3
            }

            if (gates_cis_mr.containsKey("apply_steiger_filter")) {
                apply_steiger_filter_gate = gates_cis_mr.apply_steiger_filter
            } else {
                apply_steiger_filter_gate = false
            }

            if (gates_cis_mr.containsKey("clump_kb")) {
                clump_kb_gate = gates_cis_mr.clump_kb
            } else {
                clump_kb_gate = 10000
            }

            if (gates_cis_mr.containsKey("clump_r2")) {
                clump_r2_gate = gates_cis_mr.clump_r2
            } else {
                clump_r2_gate = 0.001
            }

            if (gates_cis_mr.containsKey("instrument_pval_threshold")) {
                instrument_pval_threshold_gate = gates_cis_mr.instrument_pval_threshold
            } else {
                instrument_pval_threshold_gate = 5.0e-8
            }

            if (gates_cis_mr.containsKey("min_f_stat")) {
                min_f_stat_gate = gates_cis_mr.min_f_stat
            } else {
                min_f_stat_gate = 10
            }

            if (gates_coloc.containsKey("pp4_threshold")) {
                pp4_threshold_gate = gates_coloc.pp4_threshold
            } else {
                pp4_threshold_gate = 0.7
            }

            if (gates_coloc.containsKey("p1")) {
                p1_gate = gates_coloc.p1
            } else {
                p1_gate = 1.0e-4
            }

            if (gates_coloc.containsKey("p2")) {
                p2_gate = gates_coloc.p2
            } else {
                p2_gate = 1.0e-4
            }

            if (gates_coloc.containsKey("p12")) {
                p12_gate = gates_coloc.p12
            } else {
                p12_gate = 1.0e-5
            }

            if (gates_smr.containsKey("p_qtl_smr")) {
                p_qtl_smr_gate = gates_smr.p_qtl_smr
            } else {
                p_qtl_smr_gate = 5.0e-8
            }

            if (gates_smr.containsKey("p_qtl_heidi")) {
                p_qtl_heidi_gate = gates_smr.p_qtl_heidi
            } else {
                p_qtl_heidi_gate = 1.57e-3
            }

            if (gates_smr.containsKey("p_smr_threshold")) {
                p_smr_threshold_gate = gates_smr.p_smr_threshold
            } else {
                p_smr_threshold_gate = 0.05
            }

            if (gates_smr.containsKey("p_heidi_threshold")) {
                p_heidi_threshold_gate = gates_smr.p_heidi_threshold
            } else {
                p_heidi_threshold_gate = 0.01
            }

            if (gates_hyprcoloc.containsKey("prior_1")) {
                prior_1_gate = gates_hyprcoloc.prior_1
            } else {
                prior_1_gate = 1.0e-4
            }

            if (gates_hyprcoloc.containsKey("prior_c")) {
                prior_c_gate = gates_hyprcoloc.prior_c
            } else {
                prior_c_gate = [0.05, 0.02, 0.01, 0.005]
            }

            if (gates_hyprcoloc.containsKey("reg_thresh")) {
                reg_thresh_gate = gates_hyprcoloc.reg_thresh
            } else {
                reg_thresh_gate = [0.5, 0.6, 0.7]
            }

            if (gates_hyprcoloc.containsKey("align_thresh")) {
                align_thresh_gate = gates_hyprcoloc.align_thresh
            } else {
                align_thresh_gate = [0.5, 0.6, 0.7]
            }

            if (gates_hyprcoloc.containsKey("equal_thresholds")) {
                equal_thresholds_gate = gates_hyprcoloc.equal_thresholds
            } else {
                equal_thresholds_gate = true
            }

            if (gates_pwcoco.containsKey("pp4_threshold")) {
                pwcoco_pp4_threshold_gate = gates_pwcoco.pp4_threshold
            } else {
                pwcoco_pp4_threshold_gate = 0.7
            }

            if (gates_phewas.containsKey("bonferroni_alpha")) {
                bonferroni_alpha_gate = gates_phewas.bonferroni_alpha
            } else {
                bonferroni_alpha_gate = 0.05
            }

            if (gates_phewas.containsKey("coloc_threshold")) {
                coloc_threshold_gate = gates_phewas.coloc_threshold
            } else {
                coloc_threshold_gate = 0
            }

            def gates = [
                cis_mr: [
                    wald_fdr_q                : wald_fdr_q_gate,
                    ivw_fdr_q                 : ivw_fdr_q_gate,
                    cochran_q_pval            : cochran_q_pval_gate,
                    egger_intercept_pval_min  : egger_intercept_pval_min_gate,
                    min_instruments_for_ivw   : min_instruments_for_ivw_gate,
                    apply_steiger_filter      : apply_steiger_filter_gate,
                    clump_kb                  : clump_kb_gate,
                    clump_r2                  : clump_r2_gate,
                    instrument_pval_threshold : instrument_pval_threshold_gate,
                    min_f_stat                : min_f_stat_gate
                ],
                coloc: [
                    pp4_threshold: pp4_threshold_gate,
                    p1           : p1_gate,
                    p2           : p2_gate,
                    p12          : p12_gate
                ],
                smr: [
                    p_qtl_smr        : p_qtl_smr_gate,
                    p_qtl_heidi      : p_qtl_heidi_gate,
                    p_smr_threshold  : p_smr_threshold_gate,
                    p_heidi_threshold: p_heidi_threshold_gate
                ],
                hyprcoloc: [
                    prior_1         : prior_1_gate,
                    prior_c         : prior_c_gate,
                    reg_thresh      : reg_thresh_gate,
                    align_thresh    : align_thresh_gate,
                    equal_thresholds: equal_thresholds_gate
                ],
                pwcoco: [
                    pp4_threshold: pwcoco_pp4_threshold_gate
                ],
                phewas: [
                    bonferroni_alpha: bonferroni_alpha_gate,
                    coloc_threshold : coloc_threshold_gate
                ]
            ]

            // define meta tuple as input channel for QC module
            def meta = [
                pheno_id      : row.pheno_id,
                pqtl_dataset  : row.pqtl_dataset,
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
                ref_bfile     : row.ref_bfile,
                liftover_dir  : row.liftover_dir ?: 'dat/ref/liftover',
                gene_annotation: row.gene_annotation,
                run_smr       : run_smr,
                bulk_qtl_datasets: bulk_qtl_datasets,
                sc_qtl_dataset   : sc_qtl_dataset,
                gates             : gates
            ]
        }

    // now run QC module
    GWAS_QC(ch_in)

    emit:
    qc_tsv = GWAS_QC.out.qc_tsv
}
