# Results schema

Column reference for the TSVs under `results/`. Column names are copied verbatim from real output files (`results/cis-MR/`, `results/coloc/`, `results/SMR/`, `results/PheWAS*/` for `wingo_brain`/`AD`); columns aren't renamed here for consistency with the dashboard.

## `results/cis-MR/<pqtl_dataset>_<pheno_id>_all_MR.tsv`

One row per protein. `primary_*` mirrors whichever method (`Wald` or `IVW`) was actually used for that protein, based on `n_instruments`.

| Column | Meaning |
| --- | --- |
| `protein` | Protein identifier (`GENE_UNIPROT` or `GENE_ENSEMBL`, dataset-dependent) |
| `pqtl_dataset`, `outcome_trait` | Which pQTL panel and outcome this row is from |
| `n_instruments` | Number of independent cis-instruments for this protein |
| `methods_run` | Which MR methods were run given `n_instruments` |
| `primary_method`, `primary_beta`, `primary_se`, `primary_pval`, `primary_FDR_q` | The method actually used to gate this protein: Wald ratio if `n_instruments == 1`, IVW otherwise |
| `MDE_OR` | Minimum detectable effect (odds ratio), for power context |
| `IVW_beta/se/pval/FDR_q` | Inverse-variance-weighted estimate (`n_instruments > 1`) |
| `Egger_beta/se/pval/FDR_q`, `egger_intercept`, `egger_intercept_pval` | MR-Egger estimate and pleiotropy intercept test |
| `WME_beta/se/pval/FDR_q` | Weighted median estimator |
| `Wald_beta/se/pval/FDR_q` | Wald ratio estimate (`n_instruments == 1`) |
| `Q`, `Q_df`, `Q_pval` | Cochran's Q heterogeneity statistic, degrees of freedom, p-value |

**Gate to next stage:** `Wald_FDR_q < 0.05` (1 instrument) OR (`IVW_FDR_q < 0.05` AND `Q_pval > 0.05`).

## `results/coloc/<pqtl_dataset>/<pqtl_dataset>_<pheno_id>_all_coloc.tsv`

One row per protein, from `bin/coloc.R`.

| Column | Meaning |
| --- | --- |
| `protein_id`, `outcome_trait` | Protein and outcome this coloc run covers |
| `top_snp`, `nsnps` | Lead SNP in the tested region and number of SNPs used |
| `PP.H0.abf` … `PP.H4.abf` | Posterior probabilities for coloc's five hypotheses (H4 = shared causal variant) |
| `coloc_pass` | Pre-computed boolean: `PP.H4.abf > pp4_threshold` |
| `pp4_threshold` | The PP.H4 cutoff actually applied (currently `0.7` — read from the data rather than assumed, since `bin/coloc_targets.py` has some threshold values still marked for post-testing cleanup) |
| `n_pqtl_snps`, `n_gwas_snps`, `n_gwas`, `s_gwas` | SNP counts and outcome GWAS sample size / case fraction used in the coloc run |
| `n_cases`, `n_controls` | Outcome GWAS case/control counts |

**Gate to next stage:** `coloc_pass == True` (i.e. `PP.H4.abf > 0.7`).

## `results/SMR/<pqtl_dataset>_<pheno_id>_final_multi_omics_targets.tsv`

The canonical combined SMR output the dashboard reads — one row per target × cell-type/tissue, bulk and single-cell rows appended together. Columns up to `Lead_SNP_BP` come straight from the `smr` binary (bulk files are pre-computed and harmonised onto this same shape); everything from `q_SMR` onward is added by `bin/sort_smr.py`.

| Column | Meaning |
| --- | --- |
| `Gene`, `probeID`, `gene_id` | Gene symbol / probe / Ensembl-style ID for the eQTL probe |
| `qtl_name` | Source eQTL sub-dataset label (e.g. `eQTL_GTEx_Brain_Cortex_v10`) |
| `ProbeChr`, `Probe_bp`, `chr`, `start`, `end`, `strand` | Probe genomic coordinates |
| `topSNP`, `topSNP_chr`, `topSNP_bp` | SMR's top associated SNP and its position |
| `A1`, `A2`, `Freq` | Alleles and effect-allele frequency, aligned so `A1` is always the AD risk-increasing allele (see `align_to_risk_allele` in `bin/sort_smr.py`) |
| `b_GWAS`, `se_GWAS`, `p_GWAS` | Outcome GWAS effect at the top SNP |
| `b_eQTL`, `se_eQTL`, `p_eQTL` | eQTL effect at the top SNP (for single-cell targets, re-sourced from the original per-cell-type eQTL file rather than SMR's own report — see `pull_original_sc_eqtl_beta`) |
| `b_SMR`, `se_SMR`, `p_SMR`, `p_SMR_multi` | SMR effect estimate (`b_GWAS / b_eQTL`) and test p-values |
| `p_HEIDI`, `nsnp_HEIDI` | HEIDI heterogeneity test p-value and number of SNPs used |
| `GWAS_LOCUS`, `Lead_SNP`, `Lead_SNP_BP` | Outcome GWAS locus this probe falls in |
| `q_SMR` | FDR-corrected `p_SMR`, computed per-file across every probe (not just promising targets) by `fdr_correct_smr_file` |
| `protein` | Target identifier, mapped back from `Gene` |
| `cell_type` | Cell type (single-cell) or tissue/sub-dataset label (bulk) |
| `data_type` | `bulk` or `single_cell` |
| `phenotype`, `eqtl_dataset`, `pqtl_dataset` | Run identity columns |

**Gate to "Final Targets":** `q_SMR < 0.05` AND `p_HEIDI > 0.01`.

## `results/PheWAS-FinnGen/<pqtl_dataset>/<pheno_id>/*_PheWAS-FinnGen.tsv` and `results/PheWAS_UKBB/...` (same shape)

One row per protein × PheWAS outcome phenotype. `results/PheWAS-FinnGen` is the FinnGen screen, `results/PheWAS_UKBB` is the UK Biobank screen — same columns in both.

| Column | Meaning |
| --- | --- |
| `protein`, `pqtl_dataset`, `method` | Target and which MR method produced `beta_mr` |
| `n_instruments_original/available`, `missing_instruments`, `instrument_completeness` | How many of the protein's cis-instruments were found in this PheWAS outcome's GWAS |
| `rsid`, `snp`, `A1`, `A2` | Instrument SNP and alleles, aligned to the AD risk allele |
| `ad_effect_allele_original`, `ad_other_allele_original`, `beta_ad_original`, `beta_ad`, `ad_A1_flipped` | Original vs. risk-allele-aligned AD effect |
| `exposure_effect_allele`, `exposure_other_allele`, `beta_exposure(_original)`, `se_exposure`, `p_exposure`, `exposure_A1_flipped` | pQTL (exposure) effect, original vs. aligned |
| `finngen_ref`, `finngen_alt`, `phewas_A1_flipped` | PheWAS outcome's own allele coding and whether it was flipped to match |
| `beta_mr`, `se_mr`, `p_mr` | MR estimate of the protein's effect on this PheWAS outcome |
| `beta_phewas`, `se_phewas`, `p_phewas` | Raw PheWAS outcome GWAS association at the instrument SNP |
| `PHENOCODE`, `PHENOSTRING`, `CATEGORY` | PheWAS outcome phenotype identity |
| `p_bonferroni`, `bonferroni_significant` | Bonferroni-corrected p-value and significance flag |

**Safety gate (used by the dashboard and `analysis/pipeline_dag.dot`'s "Safety screen" stage):** a target *fails* only if it has a `bonferroni_significant` association with `beta_mr >= 0` (the same allele that raises AD risk also raises the PheWAS outcome). No significant hit, a protective (`beta_mr < 0`) significant hit, or no coverage at all — each counts as passing.
