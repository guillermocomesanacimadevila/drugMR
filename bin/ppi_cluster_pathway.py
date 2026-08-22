import polars as pl 
from pathlib import Path
from drugmr.string_ppi import StringPPI
from drugmr.unsupervised_algorithms import MarkovClustering
from drugmr.enrichr import EnrichR
from drugmr.utils import filter_mr_targets, filter_coloc_targets, filter_phewas
from drugmr.paths import mr_out, coloc_out, phewas_out


# REMINDER -> remove apptamer/ens ID

# inherit objects from drugmr
StringPpi        = StringPPI()
Markovclustering = MarkovClustering()
Enrichr          = EnrichR()


def extract_candidate_targets(pqtl_dataset: str, pheno_id: str):

    """
    Reqs:
    - cis-MR beta P_FDR (wald/IVW) < 0.05 
    - Coloc (default priors) PP.H4 > tresh
    - Pass FinnGen safety -> UKBB (not needed)
    """

    mr_df = pl.DataFrame(mr_out(pqtl_dataset=pqtl_dataset, pheno_id=pheno_id))
    coloc_df = pl.DataFrame(coloc_out(pqtl_dataset=pqtl_dataset, pheno_id=pheno_id))
    phewas_df = pl.DataFrame(phewas_out(pqtl_dataset=pqtl_dataset, pheno_id=pheno_id))
    mr_targets = {target.split("_", 1)[0] for target in filter_mr_targets(mr_df)}
    coloc_targets = {target.split("_", 1)[0] for target in filter_coloc_targets(coloc_df)}
    phewas_targets = {target.split("_", 1)[0] for target in filter_phewas(phewas_df)}
    final_targets = sorted(mr_targets & coloc_targets & phewas_targets)
    return final_targets


def pathway_enrichment_pipeline(pqtl_dataset: str, pheno_id: str):

    """
    Grab candidate targets from previous function ->
    - build STRING PPI
    - Filter by score > 0.40 or X
    - Markov clutering with r = 2.0 ** subject to change 
    - Submit to EnrichR API
    """

    # out_dirs -> need to slap them onto paths/
    targets = extract_candidate_targets(pqtl_dataset=pqtl_dataset, pheno_id=pheno_id)
    ppi_network = StringPpi.query_ppi(targets) # just a df for STRING results
    return ppi_network