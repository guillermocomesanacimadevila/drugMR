# Configure a run

One params file represents exactly one outcome trait and one pQTL dataset.

For example:

```text
params/AD.ukb_ppp.yaml    = AD outcome GWAS  + ukb_ppp pQTL panel
params/SCZ.wingo.yaml     = SCZ outcome GWAS + wingo_brain pQTL panel
```

Run each pair with its own `nextflow run` command. drugMR rejects a params file containing more than one input pair because the run ID, results, provenance, and dashboard record must all refer to one unambiguous analysis.

## Params file and analysis gates

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

The outcome GWAS is read as a tab separated (`separator=\t`) summary statistics file. Its column names are declared in the params file. The file can use any supported column names because drugMR maps them using `snp_col`, `a1_col`, `a2_col`, `beta_col`, `se_col`, `p_col`, `pos_col`, `chr_col`, and `af_col`. The params file is validated against `params/schema.json` when Nextflow starts.

The `gates` block records the statistical thresholds used by the run. Keep it in the params file, so the values are copied into `params.lock.yaml` and stay attached to the results.

- `cis_mr`: instrument selection, F statistic, Steiger filtering, Wald and IVW FDR thresholds, heterogeneity tests.
- `coloc`: coloc priors and the minimum PP4.
- `smr`: SNP selection for SMR and HEIDI, then the final SMR FDR and HEIDI thresholds.
- `hyprcoloc`: priors and sensitivity grid.
- `pwcoco`: conditional colocalisation PP4 threshold.
- `phewas`: Bonferroni alpha for the safety screen.

```yaml
gates:
  cis_mr:
    wald_fdr_q: 0.05
    ivw_fdr_q: 0.05
    cochran_q_pval: 0.05
    egger_intercept_pval_min: 0
    min_instruments_for_ivw: 3
    apply_steiger_filter: true
    clump_kb: 10000
    clump_r2: 0.001
    instrument_pval_threshold: 5.0e-8
    min_f_stat: 10
  coloc:
    pp4_threshold: 0.7
    p1: 1.0e-4
    p2: 1.0e-4
    p12: 1.0e-5
  smr:
    p_qtl_smr: 5.0e-8
    p_qtl_heidi: 1.57e-3
    p_smr_threshold: 0.05
    p_heidi_threshold: 0.01
  hyprcoloc:
    prior_1: 1.0e-4
    prior_c: [0.05, 0.02, 0.01, 0.005]
    reg_thresh: [0.5, 0.6, 0.7]
    align_thresh: [0.5, 0.6, 0.7]
    equal_thresholds: true
  pwcoco:
    pp4_threshold: 0.7
  phewas:
    bonferroni_alpha: 0.05
```

These values affect which targets proceed. Changing a gate and running with `-resume` invalidates tasks whose commands contain that value while allowing unrelated completed tasks to remain cached.

## The QTL manifest

`assets/qtl_manifest.csv` is the dataset registry. It keeps dataset specific paths and column names out of the analysis code.

The params values (within `params/*.yml`) `pqtl_dataset`, `bulk_qtl_datasets`, and `sc_qtl_dataset` refer to IDs in the `dataset` column. The corresponding manifest row tells drugMR:

1. Where the files live.
2. Whether the dataset is a pQTL or another QTL type.
3. Which genome build it uses.
4. Which column contains the gene or probe identifier.
5. How its SNP, allele, effect, standard error, p value, chromosome, position, and frequency columns are named.
6. Whether several rows belong to one parent collection, such as the tissues in `gtex_v10`.

Dataset matching is case insensitive. `metabrain`, `MetaBrain`, and `METABRAIN` resolve to the same `dataset` value. Parent matching is also case insensitive. Use one consistent lowercase spelling in new params files and manifest rows because that keeps run IDs and output directory names predictable.

`pqtl_dataset` must resolve to one pQTL dataset row. `sc_qtl_dataset` must resolve to one dataset row, which may use a wildcard to cover several cell type files. Each value in `bulk_qtl_datasets` may resolve in either of two ways. It can match one dataset row such as `metabrain`, or it can match a shared `parent_dataset` value such as `gtex_v10`.

For a grouped resource, give every tissue its own unique `dataset` value and give every row the same `parent_dataset`. Then place the parent name in `bulk_qtl_datasets`:

```yaml
bulk_qtl_datasets: [metabrain, gtex_v10]
```

```text
dataset,path,...,parent_dataset
gtex_brain_amygdala_v10,../../data/eQTL/bulk-eQTL/GTEx_v10/Brain_Amygdala/Brain_Amygdala.parquet,...,gtex_v10
gtex_brain_cortex_v10,../../data/eQTL/bulk-eQTL/GTEx_v10/Brain_Cortex/Brain_Cortex.parquet,...,gtex_v10
```

Selecting `gtex_v10` runs the registered tissues as one bulk QTL collection. Selecting `gtex_brain_amygdala_v10` addresses that individual manifest row instead.

A simplified manifest looks like this:

```text
dataset,path,key_col,qtl_type,build,n_col,sample_size,snp_col,a1_col,a2_col,beta_col,se_col,p_col,chr_col,pos_col,eaf_col,parent_dataset
wingo_brain,../../data/pQTL/mass-spec/wingo_brain/*.parquet,,pqtl,GRCh38,,1013,SNP,A1,A2,BETA,SE,P,CHR,BP,FRQ,
metabrain,../../data/eQTL/bulk-eQTL/MetaBrain/BrainMeta_cis_eQTL.parquet,Gene,eqtl,GRCh37,N,2865,SNP,A1,A2,b,SE,p,Chr,BP,Freq,
gtex_brain_amygdala_v10,../../data/eQTL/bulk-eQTL/GTEx_v10/Brain_Amygdala/Brain_Amygdala.parquet,Gene,eqtl,GRCh38,N,180,SNP,A1,A2,b,SE,p,Chr,BP,Freq,gtex_v10
singlebrain,../../data/eQTL/sc-eQTL/SingleBrain/*.parquet,GENE,eqtl,GRCh38,N,983,SNP,A1,A2,BETA,SE,P,CHR,BP,FRQ,
```

Paths may name one file or use a wildcard to register many files. They may be repository relative or absolute. On Linux, path case matters. `data/GWAS` and `data/gwas` are different directories.

Raw QTL summary statistics can be supplied as:

1. Parquet files ending in `.parquet`.
2. Comma separated files ending in `.csv`.
3. Tab separated files ending in `.tsv`.
4. Tab separated text files ending in `.txt`.

The source format and original column names do not change the downstream analysis. The manifest normalises each dataset to the columns used by drugMR.

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

## The synthesis workspace

`synthesis/` holds reusable intermediate data that is expensive to create but is not part of any one portable run.

- `synthesis/qtl_esd/`: temporary, restartable workspace where tabular QTL files are converted to ESD, FLIST, and BESD form. Completed `.besd`, `.esi`, and `.epi` triples are moved beside the manifest source so later runs can reuse them.
- `synthesis/SMR/`: reusable SMR calculations, stored by QTL dataset and outcome.
- Other subdirectories hold derived target summaries and manifests shared across stages.

For example, an SCZ run using the `wingo_brain` pQTL panel and MetaBrain as a bulk QTL dataset produces an SCZ and MetaBrain SMR calculation that does not depend on `wingo_brain`. A later SCZ run using `ukb_ppp` and MetaBrain again reuses that same calculation from `synthesis/SMR/`, instead of repeating the chromosome level SMR analysis. A different outcome GWAS or QTL dataset needs its own calculation.

```text
Run 1
params/SCZ.wingo_brain.yaml
outcome: SCZ
pQTL panel: wingo_brain
bulk QTL dataset: metabrain

Run 2
params/SCZ.ukb_ppp.yaml
outcome: SCZ
pQTL panel: ukb_ppp
bulk QTL dataset: metabrain

Shared reusable calculation
SCZ outcome GWAS + MetaBrain QTL
synthesis/SMR/bulk/MetaBrain/<SCZ MetaBrain SMR output>

Run 2 finds the completed SCZ + MetaBrain calculation and reuses it.
```

Do not treat `synthesis/` as the final results directory. Final run outputs are copied to `runs/<run_id>/results/`. Do not routinely delete `synthesis/` between runs because doing so can force expensive conversion or SMR work to run again. It can be rebuilt from the registered inputs, but only if those inputs remain available and their destination directories are writable.

## Validate the configuration

Run path checks from the repository root because relative paths are resolved from there:

```bash
cd /path/to/drugMR
ls -lh ../../data/GWAS/SCZ_Trubetskoy_2022_EUR_postQC.tsv
ls ../../data/pQTL/mass-spec/wingo_brain/*.parquet | head
ls -lh ref/liftover/hg19ToHg38.over.chain
```
