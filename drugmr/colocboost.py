import os 
import subprocess
import tempfile 
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
    def run_colocboost(
        sumstats: dict[str, pl.DataFrame],
        ld: pl.DataFrame,
        snp_order: list[str],
        out_file: str,
        r_script: str = "R/colocboost.R",
    ) -> subprocess.CompletedProcess:


        """
        Python wrapper function in-charge of running coloboost.R 
        args:
        - dict[str, df] -> whereby pheno_id: df
        - out_dir
        - ld_matrix for locus (overlapping SNPs across dfs???)
        """

        work_dir = tempfile.mkdtemp(prefix="colocboost_")
        sumstat_paths = {trait: os.path.join(work_dir, f"{trait}.txt") for trait in sumstats}
        ld_path = os.path.join(work_dir, "ld.txt")

        try:
            for trait, df in sumstats.items():
                if df["variant"] != snp_order:
                    raise ValueError(f"sumstat '{trait}' is not reindexed to snp_order -> reindex before calling run()")
                df.write_csv(sumstat_paths[trait], separator="\t")
            ld.write_csv(ld_path, separator="\t", include_header=False)
            cmd = ["Rscript", r_script, ld_path, out_file, *sumstat_paths.values()]
            return subprocess.run(cmd, check=True)
        finally:
            for path in sumstat_paths.values():
                if os.path.exists(path):
                    os.unlink(path)
                if os.path.exists(ld_path):
                    os.unlink(ld_path)
                os.rmdir(work_dir) 