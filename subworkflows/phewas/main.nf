#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include { PHEWAS_FINNGEN; PHEWAS_UKB } from '../../modules/local/phewas.nf'

workflow PHEWAS_WF {

    take:
    coloc_results   // tuple(meta, coloc.tsv, coloc_sensitivity.tsv)
    pwcoco_results  // tuple(meta, pwcoco.tsv)
    mr_results      // tuple(meta, mr.tsv, mr_instruments.tsv)
    protein_dirs    // tuple(meta, dirs)

    main:
    // keyed on (pheno_id, pqtl_dataset), same composite-key discipline as target_hits/hyprcoloc
    ch_coloc = coloc_results.map { meta, coloc_tsv, coloc_sensitivity_tsv -> tuple([meta.pheno_id, meta.pqtl_dataset], meta, coloc_tsv) }
    ch_pwcoco = pwcoco_results.map { meta, pwcoco_tsv -> tuple([meta.pheno_id, meta.pqtl_dataset], pwcoco_tsv) }
    ch_instruments = mr_results.map { meta, mr_tsv, mr_instruments_tsv -> tuple([meta.pheno_id, meta.pqtl_dataset], mr_instruments_tsv) }
    ch_dirs = protein_dirs.map { meta, dirs -> tuple([meta.pheno_id, meta.pqtl_dataset], dirs) }

    ch_keyed = ch_coloc
        .join(ch_pwcoco)
        .join(ch_instruments)
        .join(ch_dirs)
        // key, meta, coloc_tsv, pwcoco_tsv, mr_instruments_tsv, dirs

    PHEWAS_FINNGEN(ch_keyed.map { key, meta, coloc_tsv, pwcoco_tsv, mr_instruments_tsv, dirs -> tuple(meta, coloc_tsv, pwcoco_tsv, mr_instruments_tsv, dirs) })

    // UKB fallback runs strictly after FinnGen: it reads FinnGen's coverage
    // manifest (always written, even when FinnGen finds zero hits) to restrict
    // itself to targets with zero FinnGen instrument coverage
    ch_finngen_coverage = PHEWAS_FINNGEN.out.phewas_finngen_results
        .map { meta, coverage_tsv -> tuple([meta.pheno_id, meta.pqtl_dataset], coverage_tsv) }

    ch_ukb_in = ch_keyed
        .join(ch_finngen_coverage)
        .map { key, meta, coloc_tsv, pwcoco_tsv, mr_instruments_tsv, dirs, coverage_tsv -> tuple(meta, coloc_tsv, pwcoco_tsv, mr_instruments_tsv, dirs, coverage_tsv) }

    PHEWAS_UKB(ch_ukb_in)

    emit:
    phewas_finngen_results = PHEWAS_FINNGEN.out.phewas_finngen_results
    phewas_finngen_hits    = PHEWAS_FINNGEN.out.phewas_finngen_hits
    phewas_ukb_results     = PHEWAS_UKB.out.phewas_ukb_results
}
