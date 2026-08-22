import stringdb
import polars as pl


class StringPPI:

    """
    Building a PPI network from the StringDB API 
    """

    def __init__(self):
        pass


    def query_string_enrichment(self, gene_list: list[str]) -> pl.DataFrame:
        gene_list = [i.upper() for i in gene_list]
        string_ids = stringdb.get_string_ids(gene_list)
        enrichment_df = stringdb.get_enrichment(string_ids.queryItem)
        return pl.from_pandas(enrichment_df)


    def query_ppi(self, gene_list: list[str]) -> pl.DataFrame:
        gene_list = [i.upper() for i in gene_list]
        string_ids = stringdb.get_string_ids(gene_list)
        network_df = stringdb.get_network(string_ids.queryItem)       
        return pl.from_pandas(network_df)
