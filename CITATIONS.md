# drugMR: Citations

## Pipeline

If you use drugMR, please cite it. See [`CITATION.cff`](CITATION.cff).

## Nextflow

> Di Tommaso P, Chatzou M, Floden EW, Barja PP, Palumbo E, Notredame C. Nextflow enables reproducible computational workflows. Nat Biotechnol. 2017 Apr 11;35(4):316 to 319. doi: 10.1038/nbt.3820.

## Statistical methods

### Mendelian randomisation (TwoSampleMR / MR Base)

> Hemani G, Zheng J, Elsworth B, et al. The MR Base platform supports systematic causal inference across the human phenome. eLife. 2018;7:e34408. doi: 10.7554/eLife.34408.

### Colocalisation (coloc)

> Giambartolomei C, Vukcevic D, Schadt EE, et al. Bayesian test for colocalisation between pairs of genetic association studies using summary statistics. PLoS Genet. 2014;10(5):e1004383. doi: 10.1371/journal.pgen.1004383.

### HyPrColoc

> Foley CN, Staley JR, Breen PG, et al. A fast and efficient colocalization algorithm for identifying shared genetic risk factors across multiple traits. Nat Commun. 2021;12:764. doi: 10.1038/s41467-020-20885-8.

### SMR and HEIDI

> Zhu Z, Zhang F, Hu H, et al. Integration of summary data from GWAS and eQTL studies predicts complex trait gene targets. Nat Genet. 2016;48:481 to 487. doi: 10.1038/ng.3538.

### PWCoCo (pairwise conditional and colocalisation)

Bundled in this repository as a git submodule (`tools/pwcoco`). See the tool's own repository for methodological details: https://github.com/jwr-git/pwcoco

### GCTA (used here for PWCoCo's conditional analysis)

> Yang J, Lee SH, Goddard ME, Visscher PM. GCTA: a tool for genome wide complex trait analysis. Am J Hum Genet. 2011;88(1):76 to 82. doi: 10.1016/j.ajhg.2010.11.011.

### PLINK

> Purcell S, Neale B, Todd Brown K, et al. PLINK: a tool set for whole genome association and population based linkage analyses. Am J Hum Genet. 2007;81(3):559 to 575. doi: 10.1086/519795.

## Software and infrastructure

- [Docker](https://www.docker.com/): Merkel D. Docker: lightweight linux containers for consistent development and deployment. Linux Journal. 2014;2014(239):2.
- [Apptainer](https://apptainer.org/) / [Singularity](https://sylabs.io/singularity/): Kurtzer GM, Sochat V, Bauer MW. Singularity: Scientific containers for mobility of compute. PLoS One. 2017;12(5):e0177459. doi: 10.1371/journal.pone.0177459.
- [PostgreSQL](https://www.postgresql.org/)
- [Streamlit](https://streamlit.io/)
- [Synapse](https://www.synapse.org/), for distribution of some pQTL cohorts.

## Data sources

This pipeline is dataset agnostic. Each cohort's own citation requirement is set by its data provider, not listed here.
