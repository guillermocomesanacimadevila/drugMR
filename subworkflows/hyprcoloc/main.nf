#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include { HYPRCOLOC; MERGE_HYPRCOLOC } from '../../modules/local/hyprcoloc.nf'

workflow HYPRCOLOC_WF {

    take:
    smr_final_targets   // tuple(meta, final_multi_omics_targets.tsv)
    protein_dirs        // tuple(meta, dirs)

    main:
    // keyed on (pheno_id, pqtl_dataset), same composite-key discipline as
    ch_targets = smr_final_targets.map { meta, targets_file -> tuple([meta.pheno_id, meta.pqtl_dataset], meta, targets_file) }
    ch_dirs = protein_dirs.map { meta, dirs -> tuple([meta.pheno_id, meta.pqtl_dataset], dirs) }

    ch_joined = ch_targets
        .join(ch_dirs)
        .map { key, meta, targets_file, dirs -> tuple(meta, targets_file, dirs) }

    ch_hyprcoloc_in = ch_joined
        .filter { meta, targets_file, dirs -> meta.run_smr && (meta.bulk_qtl_datasets || meta.sc_qtl_dataset) }
        .flatMap { meta, targets_file, dirs ->
            def qtl_datasets = (meta.bulk_qtl_datasets ?: []) + (meta.sc_qtl_dataset ? [meta.sc_qtl_dataset] : [])
            qtl_datasets.collect { ds -> tuple(meta, ds, targets_file, dirs) }
        }

    HYPRCOLOC(ch_hyprcoloc_in)


    ch_merge_in = HYPRCOLOC.out.hyprcoloc_results
        .map { meta, qtl_dataset, dataset_file, pdfs -> tuple([meta.pheno_id, meta.pqtl_dataset], meta, qtl_dataset, dataset_file) }
        .groupTuple(by: 0)
        .map { key, metas, labels, files -> tuple(metas[0], labels, files) }

    MERGE_HYPRCOLOC(ch_merge_in)

    emit:
    hyprcoloc_results = HYPRCOLOC.out.hyprcoloc_results
    hyprcoloc_final    = MERGE_HYPRCOLOC.out.hyprcoloc_final
}
