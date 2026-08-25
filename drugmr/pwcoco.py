import os
import tempfile
import polars as pl
import subprocess


class PWCoCo:

    """
    Object to run PWCoCo in parallel to standard pairwise coloc
    """
    
    def __init__(self):
        pass

    @staticmethod
    def harmonise_sumstats(
            df: pl.DataFrame,
            snp_col: str,
            a1_col: str,
            a2_col: str,
            maf_col: str,
            beta_col: str,
            se_col: str,
            p_col: str,
            n_col: str
    ) -> pl.DataFrame:
        
        return df.select(
            pl.col(snp_col).alias("SNP"),
            pl.col(a1_col).alias("A1"),
            pl.col(a2_col).alias("A2"),
            pl.col(maf_col).alias("A1_freq"),
            pl.col(beta_col).alias("beta"),
            pl.col(se_col).alias("se"),
            pl.col(p_col).alias("p"),
            pl.col(n_col).alias("n"),
        )

    @staticmethod
    def pwcoco(
            ref_bfile: str,
            sumstats_1: pl.DataFrame,
            sumstats_2: pl.DataFrame,
            n_1: int,
            n_2: int,
            n2_case: int,
            out_dir: str,
            threads: int
            ):

        work_dir = tempfile.mkdtemp(prefix="pwcoco_")
        sumstats_1_path = os.path.join(work_dir, "sumstats1.txt")
        sumstats_2_path = os.path.join(work_dir, "sumstats2.txt")

        try:
            sumstats_1.write_csv(sumstats_1_path, separator="\t")
            sumstats_2.write_csv(sumstats_2_path, separator="\t")

            cmd = [
                "pwcoco",
                "--bfile", ref_bfile,
                "--sum_stats1", sumstats_1_path,
                "--sum_stats2", sumstats_2_path,
                "--n1", str(n_1),
                "--n2", str(n_2),
                "--n2_case", str(n2_case),
                "--out", out_dir,
                "--out_cond",
                "--threads", str(threads),
            ]

            return subprocess.run(cmd, check=True)
        finally:
            os.unlink(sumstats_1_path)
            os.unlink(sumstats_2_path)
            os.rmdir(work_dir)