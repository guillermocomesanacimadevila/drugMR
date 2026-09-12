# Configure a run

One params file represents exactly one outcome trait and one pQTL dataset.

For example:

```text
params/AD.ukb_ppp.yaml    = AD outcome GWAS  + ukb_ppp pQTL panel
params/SCZ.wingo.yaml     = SCZ outcome GWAS + wingo_brain pQTL panel
```

Run each pair with its own `nextflow run` command. drugMR rejects a params file containing more than one input pair because the run ID, results, provenance, and dashboard record must all refer to one unambiguous analysis.

## Params file

Start from an existing file under `params/` and change the outcome, pQTL dataset, paths, column names, and optional downstream datasets.

```yaml
pheno_id: SCZ
sumstats: ../../data/GWAS/SCZ_Trubetskoy_2022_EUR_postQC.tsv
n_cases: 67390
n_controls: 94015
genome_build: GRCh37
target_build: GRCh38

snp_col: SNP
a1_col: A1
a2_col: A2
beta_col: BETA
se_col: SE
p_col: P
pos_col: POS
chr_col: CHR
af_col: FRQ

pqtl_dataset: wingo_brain

reference_data:
  ref_bfile: ref/1000G_EUR_Phase3_plink/1000G.EUR.QC.ALL
  liftover_dir: ref/liftover
  gene_annotation: ref/NCBI/NCBI_genes_grch38_with_synonyms.tsv

run_smr: true
bulk_qtl_datasets: [metabrain, gtex_v10]
sc_qtl_dataset: singlebrain

maf: 0.01
overwrite: false
remove_mhc: true
remove_apoe: false
```

The outcome GWAS is read as a tab separated summary statistics file. Its column names are declared in the params file. The file can use any supported column names because drugMR maps them using `snp_col`, `a1_col`, `a2_col`, `beta_col`, `se_col`, `p_col`, `pos_col`, `chr_col`, and `af_col`.

The params file is validated against `params/schema.json` when Nextflow starts.

## The QTL manifest

`assets/qtl_manifest.csv` is the dataset registry. It keeps dataset specific paths and column names out of the analysis code.

The params values `pqtl_dataset`, `bulk_qtl_datasets`, and `sc_qtl_dataset` refer to IDs in the `dataset` column. The corresponding manifest row tells drugMR:

1. Where the files live.
2. Whether the dataset is a pQTL or another QTL type.
3. Which genome build it uses.
4. Which column contains the gene or probe identifier.
5. How its SNP, allele, effect, standard error, p value, chromosome, position, and frequency columns are named.
6. Whether several rows belong to one parent collection, such as the tissues in `gtex_v10`.

A simplified manifest looks like this:

```text
dataset,path,key_col,qtl_type,build,n_col,sample_size,snp_col,a1_col,a2_col,beta_col,se_col,p_col,chr_col,pos_col,eaf_col,parent_dataset
wingo_brain,../../data/pQTL/mass-spec/wingo_brain/*.parquet,,pqtl,GRCh38,,1013,SNP,A1,A2,BETA,SE,P,CHR,BP,FRQ,
metabrain,../../data/eQTL/bulk-eQTL/MetaBrain/BrainMeta_cis_eQTL.parquet,Gene,eqtl,GRCh37,N,2865,SNP,A1,A2,b,SE,p,Chr,BP,Freq,
gtex_brain_amygdala_v10,../../data/eQTL/bulk-eQTL/GTEx_v10/Brain_Amygdala/Brain_Amygdala.parquet,Gene,eqtl,GRCh38,N,180,SNP,A1,A2,b,SE,p,Chr,BP,Freq,gtex_v10
singlebrain,../../data/eQTL/sc-eQTL/SingleBrain/*.parquet,GENE,eqtl,GRCh38,N,983,SNP,A1,A2,BETA,SE,P,CHR,BP,FRQ,
```

Paths may name one file or use a wildcard to register many files. They may be repository relative or absolute. On Linux, path case matters. `data/GWAS` and `data/gwas` are different directories.

## Supported QTL inputs

Raw QTL summary statistics can be supplied as:

1. Parquet files ending in `.parquet`.
2. Comma separated files ending in `.csv`.
3. Tab separated files ending in `.tsv`.
4. Tab separated text files ending in `.txt`.

The source format and original column names do not change the downstream analysis. The manifest normalises each dataset to the columns used by drugMR.

## Existing SMR files

drugMR can use a QTL dataset that is already in SMR format. A complete SMR dataset prefix contains:

```text
dataset.besd
dataset.esi
dataset.epi
```

Keep complete triples under the directory containing the source path declared in the manifest. Nested directories are supported. For example:

```text
MetaBrain/
├── BrainMeta_cis_eQTL.parquet
└── BrainMeta_cis_eQTL/
    ├── chr1/
    │   ├── chr1.besd
    │   ├── chr1.esi
    │   └── chr1.epi
    └── chr2/
        ├── chr2.besd
        ├── chr2.esi
        └── chr2.epi
```

drugMR searches this location recursively and reuses every complete triple. An incomplete prefix is ignored.

If no complete SMR triples are found, drugMR reads the registered Parquet, CSV, TSV, or TXT source, converts it chromosome by chromosome in `synthesis/qtl_esd/`, and moves the completed triples into a dataset specific directory beside the declared source. The source location must be writable when this conversion is required.

## Check paths before a run

Run path checks from the repository root because relative paths are resolved from there:

```bash
cd /path/to/drugMR
ls -lh ../../data/GWAS/SCZ_Trubetskoy_2022_EUR_postQC.tsv
ls ../../data/pQTL/mass-spec/wingo_brain/*.parquet | head
ls -lh ref/liftover/hg19ToHg38.over.chain
```
