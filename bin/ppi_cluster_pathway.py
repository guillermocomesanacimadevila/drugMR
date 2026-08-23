import argparse
import polars as pl 
from pathlib import Path
from drugmr.string_ppi import StringPPI
from drugmr.unsupervised_algorithms import MarkovClustering
from drugmr.enrichr import EnrichR
from drugmr.utils import filter_mr_targets, filter_coloc_targets, filter_phewas, strip_protein_id
from drugmr.paths import mr_out, coloc_out, phewas_out


# REMINDER -> remove apptamer/ens ID

class PathwayEnrichmentPipeline:

    """
    Pathway enrichment object comprising STRING PPI -> MCL clustering -> EnrichR
    """

    def __init__(self):

        self.StringPPI = StringPPI()
        self.MarkovClustering = MarkovClustering()
        self.EnrichR = EnrichR()


    def extract_candidate_targets(self, pqtl_dataset: str, pheno_id: str):

        """
        Reqs:
        - cis-MR beta P_FDR (wald/IVW) < 0.05 
        - Coloc (default priors) PP.H4 > tresh
        - Pass FinnGen safety -> UKBB (not needed)
        """

        mr_df = pl.DataFrame(mr_out(pqtl_dataset=pqtl_dataset, pheno_id=pheno_id))
        coloc_df = pl.DataFrame(coloc_out(pqtl_dataset=pqtl_dataset, pheno_id=pheno_id))
        phewas_df = pl.DataFrame(phewas_out(pqtl_dataset=pqtl_dataset, pheno_id=pheno_id))
        mr_targets = set(strip_protein_id(filter_mr_targets(mr_df)))
        coloc_targets = set(strip_protein_id(filter_coloc_targets(coloc_df)))
        phewas_targets = set(strip_protein_id(filter_phewas(phewas_df)))
        final_targets = sorted(mr_targets & coloc_targets & phewas_targets)
        return final_targets


    def query_and_filter_ppi(self, targets: list[str], string_score: float):

        """
        Query STRING PPI for the given targets and filter by score.
        """

        ppi_network = self.StringPPI.query_ppi(targets) # just a df for STRING results
        # add score gait into params/ for STRING ##################
        ppi_network = ppi_network.filter(pl.col("score") > string_score)
        for row in ppi_network.iter_rows(named=True):
            gene1 = row["preferredName_A"]
            gene2 = row["preferredName_B"]
            sc = row["score"]
            print(f"[TRACKING] STRING score between {gene1} and {gene2}: {sc}")
        return ppi_network


    def build_ppi_network(self, pqtl_dataset: str, pheno_id: str, string_score: float):

        """
        Grab candidate targets from previous function -> build + filter STRING PPI.
        """

        # out_dirs -> need to slap them onto paths/ ##################
        targets = self.extract_candidate_targets(pqtl_dataset=pqtl_dataset, pheno_id=pheno_id)
        return self.query_and_filter_ppi(targets=targets, string_score=string_score)


    def perform_markov_clustering(self):
        return 


if __name__ == "__main__":
    pipeline = PathwayEnrichmentPipeline()
    gene_list = ["ADAM10_Q1234", "APP_928D36", "SPI1_183937"]
    targets = strip_protein_id(gene_list)
    ppi_network = pipeline.query_and_filter_ppi(targets=targets, string_score=0.4)
    print(ppi_network)
