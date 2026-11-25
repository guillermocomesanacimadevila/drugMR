import os
import shutil
import subprocess
from pathlib import Path

import polars as pl

GCTA_BINARY = shutil.which("gcta") or shutil.which("gcta64") or "gcta"


class COJO:

    """ 
    Need: collinear = 0.99, ref_bfile, work_dir, gcta
    - on workflow:
    -> exclude cis-MR instruments + the ones in LD with it (i.e. r2 >= 0.80)
    """


    def __init__(
            self,
            ref_bfile,
            work_dir,
            collinear=0.99,
            gcta=GCTA_BINARY
    ):

        self.ref_bfile = str(ref_bfile)
        self.work_dir = Path(work_dir)
        os.makedirs(str(self.work_dir), exist_ok=True)
        self.collinear = collinear
        self.gcta = gcta

    @staticmethod
    def to_ma(
        sumstats: pl.DataFrame,
        snp_col: str = "SNP",
        a1_col: str = "A1",
        a2_col: str = "A2",
        freq_col: str = "FRQ",
        beta_col: str = "BETA",
        se_col: str = "SE",
        p_col: str = "P",
        n_col: str = "N"
    ) -> pl.DataFrame:
    
        """ Convert current summary stats (1 input) -> .ma file within $work_dir"""

        # ensure data types == ok
        ma = sumstats.select(
            pl.col(snp_col).alias("SNP").cast(pl.String),
            pl.col(a1_col).alias("A1").cast(pl.String),
            pl.col(a2_col).alias("A2").cast(pl.String),
            pl.col(freq_col).alias("freq").cast(pl.Float64), # freq
            pl.col(beta_col).alias("b").cast(pl.Float64), # b
            pl.col(se_col).alias("se").cast(pl.Float64), # se
            pl.col(p_col).alias("p").cast(pl.Float64), # p
            pl.col(n_col).alias("N") # n_col can vary
        )

        # filter 
        ma = ma.filter(pl.col("SNP").is_not_null())
        ma = ma.unique(subset="SNP", keep="first", maintain_order=True)
        return ma # converted sumstats

    def conditional(
            self,
            ma: pl.DataFrame,
            condition_on,
            chromosome,
            tag

    ):

        """ Function to run cond GCTA-COJO on SNPs """
        ma_file = self.work_dir / f"{tag}.ma"
        cond_file = self.work_dir / f"{tag}.cond.txt" ######
        out_prefix = self.work_dir / tag

        ma.write_csv(ma_file, separator="\t")
        
        # extract from "extract_relevant_snps() within utils.py"
        cond_file.write_text("\n".join(condition_on) + "\n") #### SNPs which are being conditioned on 

        cmd = [
            self.gcta,
            "--bfile", self.ref_bfile,
            "--chr", str(chromosome).removeprefix("chr"),
            "--cojo-file", str(ma_file),
            "--cojo-cond", str(cond_file), # just a list of SNPs (i.e. a .txt file of rsIDs to condition on)
            "--cojo-collinear", str(self.collinear),
            "--out", str(out_prefix)
        ]

        result = subprocess.run(cmd, capture_output=True, text=True)
        
        cma_file = Path(f"{out_prefix}.cma.cojo")
        if result.returncode != 0 or not cma_file.exists():

            # GCTA writes its error messages to stdout, not stderr
            message = (result.stdout + result.stderr).strip()
            print(f"[CONCERN] COJO gave no result for {tag}: {message[-300:]}")
            return None
        
        # COJO writes NA for the conditioning SNPs themselves (and any it drops)
        return pl.read_csv(cma_file, separator="\t", null_values="NA", infer_schema_length=None)