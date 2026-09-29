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
    // hand the process only the header plus this dataset's row, so editing an
    // unrelated manifest row (an eQTL path, a new dataset) doesn't invalidate
    // the -resume cache for cis regions and everything downstream of it
    // same resolution as smr.nf / hyprcoloc.nf / pwcoco_qtl.nf: absolute stays,
    // relative is taken from the repo root, not wherever nextflow was launched
    def manifest_file = new File(params.manifest_path as String).isAbsolute()
        ? file(params.manifest_path)
        : file("${projectDir}/${params.manifest_path}")
    if (!manifest_file.exists()) {
        error "QTL manifest not found: ${manifest_file}"
    }
    def manifest_lines = manifest_file.readLines().findAll { line -> line.trim() }
    ch_in = qc_out.map { meta, qc_tsv ->
        def row = manifest_lines.drop(1).find { line ->
            line.split(',', -1)[0].trim().equalsIgnoreCase(meta.pqtl_dataset as String)
        }
        if (!row) {
            error "pqtl_dataset '${meta.pqtl_dataset}' not found in ${manifest_file}"
        }
        tuple(meta, qc_tsv, "${manifest_lines[0]}\n${row}\n")
    }

    EXTRACT_CIS_REGIONS(ch_in)

    emit:
    protein_dirs = EXTRACT_CIS_REGIONS.out.protein_dirs
}