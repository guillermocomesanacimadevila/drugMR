#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { EXPORT_LOCUS_ASSETS } from '../../modules/local/locus_assets.nf'

workflow LOCUS_ASSETS_WF {

    take:
    target_stats
    protein_dirs

    main:
    ch_targets = target_stats.map { meta, target_stats_file ->
        tuple([meta.pheno_id, meta.pqtl_dataset], meta, target_stats_file)
    }

    ch_regions = protein_dirs.map { meta, dirs ->
        tuple([meta.pheno_id, meta.pqtl_dataset], dirs)
    }

    ch_export = ch_targets
        .join(ch_regions)
        .map { key, meta, target_stats_file, dirs ->
            tuple(meta, target_stats_file, dirs)
        }

    EXPORT_LOCUS_ASSETS(ch_export)

    emit:
    locus_assets = EXPORT_LOCUS_ASSETS.out.locus_assets
}
