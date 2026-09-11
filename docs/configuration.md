# Configure a run

Create one YAML file under `params/` for each outcome GWAS and pQTL dataset pair. The file is validated against `params/schema.json`.

```yaml
pheno_id: AD
sumstats: dat/gwas/AD.tsv
n_cases: 128681
n_controls: 849833
genome_build: GRCh38
target_build: GRCh38

snp_col: SNP
a1_col: A1
a2_col: A2
beta_col: BETA
se_col: SE
p_col: P
pos_col: BP
chr_col: CHR
af_col: FRQ

pqtl_dataset: test_ukb_ppp
reference_data:
  ref_bfile: ref/1000G_EUR_Phase3_plink/1000G.EUR.QC.ALL
  liftover_dir: ref/liftover
  gene_annotation: ref/NCBI/NCBI_genes_grch38_with_synonyms.tsv

run_smr: true
bulk_qtl_datasets: [metabrain]
sc_qtl_dataset: ""

maf: 0.01
overwrite: false
remove_mhc: true
remove_apoe: false
```

## Register QTL datasets

`pqtl_dataset`, each entry in `bulk_qtl_datasets`, and `sc_qtl_dataset` must match an ID in `assets/qtl_manifest.csv`. Each manifest row maps the dataset to its files, column names, genome build, sample size, and optional parent dataset.

Supported source formats are Parquet, CSV, TSV, and TXT. Large QTL inputs are processed chromosome by chromosome during conversion to SMR format. `synthesis/qtl_esd/` is the conversion workspace. Once a complete `.besd`, `.esi`, and `.epi` triple has been built, drugMR moves it into a dataset-specific subdirectory beside the source file declared in the manifest. Later phenotypes discover those files from the manifest location and reuse them without repeating conversion. The declared QTL location must therefore be writable when SMR-format files need to be generated.

The statistical thresholds live in the optional `gates` block. Start from `params/AD.test.yaml` when creating a new configuration.
