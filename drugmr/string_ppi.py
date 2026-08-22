import stringdb
import polars as pl


class StringPPI:

    """
    Building a PPI network from the StringDB API 
    """

    def __init__(self):
        pass


    def query_ppi(self, gene_list: list[str]) -> pl.DataFrame:
        string_ids = stringdb.get_string_ids(gene_list)
        enrichment_df = stringdb.get_enrichment(string_ids.queryItem)
        return pl.from_pandas(enrichment_df)


if __name__ == "__main__":
    ppi = StringPPI()
    aa = ppi.query_ppi(
        ["ADAM10", "APP"]
    )
    print(aa)