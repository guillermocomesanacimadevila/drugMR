import pyenrichr as pye
import polars as pl 


class EnrichR:

    """
    EnrichR module for drugmr package
    """

    def __init__(self, library1: str = "GO_Biological_Process_2023", library2: str = "KEGG_2021_Human"):
        self.lib_1 = pye.libraries.get_library(library1)
        self.lib_2 = pye.libraries.get_library(library2)


    def run_enrichr(self, gene_list: list[str]) -> tuple[pl.DataFrame, pl.DataFrame]:
        gene_set = set(gene_list)

        result = pye.enrichment.fisher(
            gene_set,
            self.lib_1
        )

        res = pye.enrichment.fisher(
            gene_set,
            self.lib_2
        )

        return (pl.from_pandas(result), pl.from_pandas(res))