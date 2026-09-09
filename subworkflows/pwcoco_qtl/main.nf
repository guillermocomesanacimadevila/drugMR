#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include { PWCOCO_QTL } from '../../modules/local/pwcoco_qtl.nf'

workflow PWCOCO_QTL_WF {

    take:
    smr_final_targets  // tuple(meta, final_multi_omics_targets.tsv) - empty when SMR didn't run
    pwcoco_results      // tuple(meta, pwcoco.tsv)
    protein_dirs         // tuple(meta, dirs)

    main:
    // keyed on (pheno_id, pqtl_dataset), same composite-key discipline as target_hits/hyprcoloc/phewas.
    // Joining against smr_final_targets (empty unless run_smr && a bulk/sc eqtl
    // dataset was configured) is what makes this correctly no-op exactly like
    // HYPRCOLOC_WF, matching local.py's own "if run_smr:" gate.
    ch_targets = smr_final_targets.map { meta, targets_file -> tuple([meta.pheno_id, meta.pqtl_dataset], meta, targets_file) }
    ch_pwcoco = pwcoco_results.map { meta, pwcoco_tsv -> tuple([meta.pheno_id, meta.pqtl_dataset], pwcoco_tsv) }
    ch_dirs = protein_dirs.map { meta, dirs -> tuple([meta.pheno_id, meta.pqtl_dataset], dirs) }

    ch_in = ch_targets
        .join(ch_pwcoco)
        .join(ch_dirs)
        .map { key, meta, targets_file, pwcoco_tsv, dirs -> tuple(meta, targets_file, pwcoco_tsv, dirs) }

    PWCOCO_QTL(ch_in)

    emit:
    pwcoco_eqtl_pqtl_results  = PWCOCO_QTL.out.pwcoco_eqtl_pqtl_results
    pwcoco_eqtl_gwas_results  = PWCOCO_QTL.out.pwcoco_eqtl_gwas_results
    pwcoco_qtl_shared_results = PWCOCO_QTL.out.pwcoco_qtl_shared_results
}
