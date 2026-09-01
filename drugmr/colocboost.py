import subprocess
import polars as pl 


class ColocBoost:


    """
    Python wrapper for colocboost package inside R/
    for usage within bin/ and snakemake wrappers
    """

    def __init__(self):
        pass


    @staticmethod    
    def harmonise_data(
        df: pl.DataFrame,
        snp_col: str,
        se_col: str,
        beta_col: str,
        n_col: str
    ) -> pl.DataFrame:

        """ Simple funct to harmonise sumstats to colocboost format """
        
        return df.select(
            pl.col(snp_col).alias("variant"),
            pl.col(beta_col).alias("beta"),
            pl.col(se_col).alias("sebeta"),
            pl.col(n_col).alias("n")
        )

    @staticmethod
    def run_colocboost():

        """
        Python wrapper function in-charge of running coloboost.R 
        args:
        - dict[str, df] -> whereby pheno_id: df
        - out_dir
        - ld_matrix for locus (overlapping SNPs across dfs???)
        """


        
        return