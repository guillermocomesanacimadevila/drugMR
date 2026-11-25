# Methods

This page states what each stage tests and the rule a target must meet to move to the next stage. Every threshold below is a default from the `gates` block of the params file. The values used by a run are saved in its `params.lock.yaml`.

## Stages and gates

| Stage | What it does | Pass rule (default) | Params key |
| --- | --- | --- | --- |
| GWAS QC | Lifts the outcome GWAS to the target build, removes variants below the MAF and INFO cut-offs, and optionally removes the MHC and APOE regions. | Variant level filter only. | `maf`, `info_threshold`, `remove_mhc`, `remove_apoe` |
| cis regions | Cuts the outcome GWAS to the positions covered by each protein's pQTL summary statistics. | Proteins with no overlapping variants are not tested. | |
| cis-MR | Selects instruments and estimates the effect of protein level on the outcome. | See [cis-MR](#cis-mr). | `gates.cis_mr` |
| Colocalisation | Runs `coloc.abf` and conditional PWCoCo between the pQTL and the outcome GWAS for every cis-MR target. | PP.H4 ≥ 0.7 in coloc **or** H4 ≥ 0.7 in PWCoCo. | `gates.coloc`, `gates.pwcoco` |
| Locus conditioning | Conditions the outcome GWAS on other independent signals in each colocalised target's cis window, then reruns cis-MR and coloc. | Reported, not a gate. See [Locus conditioning](#locus-conditioning). | `gates.locus_conditioning` |
| PheWAS | Tests colocalised targets against FinnGen R13 endpoints. Targets with no FinnGen instrument fall back to UK Biobank. | Reported as a flag, not a gate. See [PheWAS](#phewas). | `gates.phewas` |
| SMR and HEIDI | Tests colocalised targets against bulk and single cell QTL datasets. | q_SMR ≤ 0.05 and p_HEIDI ≥ 0.01. | `gates.smr` |
| QTL PWCoCo | Conditional colocalisation across the pQTL, QTL and outcome GWAS for SMR targets. | The same SNP reaches H4 ≥ 0.7 in all three pairs (pQTL and GWAS, QTL and pQTL, QTL and GWAS). | `gates.pwcoco` |
| HyPrColoc | Tests whether the pQTL, QTL and outcome GWAS share one causal variant. | A cluster containing all three traits with posterior probability ≥ 0.7. The threshold is set in the dashboard. | `gates.hyprcoloc` |
| Druggability | Looks up every target that passes cis-MR and colocalisation in Pharos, ChEMBL, Open Targets and B3DB. Runs last and needs internet access. | Reported, not a gate. See [Druggability](#druggability). | `gates.druggability`, `--skip_druggability` |

## cis-MR

**Instruments.** SNPs with pQTL P < 5 × 10⁻⁸ are clumped against the LD reference (r² < 0.001 within 10,000 kb). Instruments with F < 10 are removed. Steiger filtering is off by default; when on, instruments explaining more variance in the outcome than in the protein are removed.

**Estimators.**

| Instruments | Primary method | Sensitivity methods |
| --- | --- | --- |
| 1 | Wald ratio | None |
| 2 | IVW | Cochran's Q |
| 3 or more | IVW | Cochran's Q, weighted median, MR Egger |

**FDR.** `primary_FDR_q` is the Benjamini Hochberg q value of the primary method's P value, corrected across every protein tested in the panel.

**Gate.**

- 1 instrument: `primary_FDR_q` < 0.05.
- 2 or more instruments: `primary_FDR_q` < 0.05 and Cochran's Q P > 0.05.
- 3 or more instruments: also MR Egger intercept P > `egger_intercept_pval_min` (0 by default, so the intercept is reported but does not filter).

## Locus conditioning

A target can pass cis-MR and colocalisation because a neighbouring gene's signal sits in the same window. For every target that passes cis-MR and coloc or PWCoCo, genome-wide significant outcome GWAS SNPs in the cis window are clumped (r² 0.1, 1,000 kb) into independent signal leads. Each lead is compared with the target's instruments: r² ≥ 0.8 means it is the target's own signal; anything lower makes it a rival, flagged when r² is between 0.3 and 0.8. The outcome GWAS is then conditioned on the rivals with GCTA-COJO, and cis-MR and coloc are rerun with the pQTL unchanged.

A target holds when its MR estimate and PP.H4 survive conditioning. The APOE region is not conditioned, because its signal is too strong to model with the LD reference. Results are reported alongside each target and do not remove targets.

## PheWAS

Each target's cis instruments are aligned to the protein increasing allele and tested against 2,511 FinnGen R13 ICD endpoints. Significance is Bonferroni corrected across all 2,511 endpoints, because the FinnGen API only returns endpoints with P < 0.05. Targets with no FinnGen instrument are tested in UK Biobank instead.

A significant association is labelled by comparing its direction with the cis-MR estimate for the outcome:

- **Additional indication**: same direction. Changing the protein in the direction that lowers outcome risk would also lower risk of this trait.
- **Adverse effect**: opposite direction. The same change would raise risk of this trait.

PheWAS results do not remove targets. They are shown alongside each target in the dashboard.

## Druggability

Every target that passes cis-MR and coloc or PWCoCo is looked up in four public resources. Pharos/TCRD gives the target development level (Tclin, Tchem, Tbio or Tdark) and protein family. ChEMBL gives the drugs with a recorded mechanism against the target, their highest clinical phase and indications, and the number of compounds with pChEMBL ≥ `pchembl_thresh` (6 by default, 1 µM). Potency is counted only against the single protein, not against families or complexes that contain it. Open Targets gives tractability for small molecules, antibodies and oligonucleotides. B3DB labels whether a small molecule drug crosses the blood brain barrier.

Each drug's mechanism is classed as Inhibit or Activate and compared with the direction the MR estimate calls for: Inhibit when higher protein raises outcome risk, Activate when it lowers it. A drug acting the required way agrees with the MR, one acting the other way opposes it, and a drug recorded with both directions or neither is unclear.

The release and access date of each resource are written to the sources table. If a resource cannot be reached, the run still completes with empty druggability tables and the error in the sources table. Druggability does not remove targets.

## Effect direction

Most tables report effects for the **outcome GWAS risk allele** (the allele with a positive GWAS beta). On that allele:

- a positive pQTL beta means higher protein level raises outcome risk, so the target is a candidate for **inhibition**;
- a negative pQTL beta means higher protein level lowers outcome risk, so the target is a candidate for **activation**.

PheWAS tables are the exception: they use the protein increasing allele. sQTL betas measure a splicing ratio, so their sign cannot be compared with pQTL or eQTL betas.

## Assumptions and limitations

- **Genome build.** All positions are GRCh38 after GWAS QC. Set `genome_build` to the build the GWAS file is actually in. A GWAS lifted twice shows a constant position offset from the pQTL data (see [Troubleshooting](troubleshooting.md#gwas-and-pqtl-positions-do-not-match)).
- **Ancestry.** The default LD reference is 1000 Genomes Phase 3 European. Clumping, conditional analysis and HEIDI assume the GWAS and QTL samples match it.
- **One pQTL panel per run.** Run each panel separately and compare runs in the dashboard.
- **Variant matching.** Variants are matched across datasets on rsID. Before SMR, rsIDs that carry more than one allele pair in a QTL dataset are dropped.
- **Sample overlap.** No correction is applied for overlap between the pQTL and outcome GWAS samples.
- **Colocalisation priors.** coloc uses p1 = p2 = 10⁻⁴ and p12 = 10⁻⁵. `results/coloc/coloc_sensitivity.tsv` shows how PP.H4 changes across other priors.
