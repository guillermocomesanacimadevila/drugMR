from gprofiler import GProfiler
import polars as pl


class RunGProfiler:

    """Run g:Profiler through its API."""

    def __init__(self, organism: str = "hsapiens"):
        self.organism = organism
        self.gp = GProfiler(return_dataframe=True)


    def run_gprofiler(self, gene_list: list[str]) -> pl.DataFrame:
        result = self.gp.profile(
            organism=self.organism,
            query=gene_list,
        )
        return pl.from_pandas(result)
