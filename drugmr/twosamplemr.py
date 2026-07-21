#!/usr/bin/env python3
import polars as pl 
import numpy as np 
from scipy.stats import norm 


class PyTwoSampleMR:

    """
    Assortment of Mendelian Randomisation functions for generalised usability in drugMR pipeline
    # 1. IVW
    # -> need to add cochran Q
    """

    def __init__(self):
        pass

    def IVW(
            self,
            exposure_df,
            outcome_df,
            exposure_snp_col,
            exposure_beta_col,
            exposure_se_col,
            outcome_snp_col,
            outcome_beta_col,
            outcome_se_col):  # IVW + delta method
        
        # compile inputs and match SNPs between exposure and outcome
        exposure = exposure_df.select(
            pl.col(exposure_snp_col).cast(pl.Utf8).alias("SNP"),
            pl.col(exposure_beta_col).cast(pl.Float64).alias("BETA_EXPOSURE"),
            pl.col(exposure_se_col).cast(pl.Float64).alias("SE_EXPOSURE")
        )

        outcome = outcome_df.select(
            pl.col(outcome_snp_col).cast(pl.Utf8).alias("SNP"),
            pl.col(outcome_beta_col).cast(pl.Float64).alias("BETA_OUTCOME"),
            pl.col(outcome_se_col).cast(pl.Float64).alias("SE_OUTCOME")
        )

        df = (
            exposure
            .join(
                outcome,
                on="SNP",
                how="inner"
            )
            .drop_nulls()
            .filter(
                pl.col("BETA_EXPOSURE").is_finite() &
                pl.col("SE_EXPOSURE").is_finite() &
                pl.col("BETA_OUTCOME").is_finite() &
                pl.col("SE_OUTCOME").is_finite() &
                (pl.col("BETA_EXPOSURE") != 0) &
                (pl.col("SE_EXPOSURE") > 0) &
                (pl.col("SE_OUTCOME") > 0)
            )
        )

        if df.height == 0:
            raise ValueError(
                "No valid overlapping SNPs were found between exposure and outcome."
            )

        # compile matched inputs
        snps           = df["SNP"].to_numpy()
        betas_exposure = df["BETA_EXPOSURE"].to_numpy()
        ses_exposure   = df["SE_EXPOSURE"].to_numpy()
        betas_outcome  = df["BETA_OUTCOME"].to_numpy()
        ses_outcome    = df["SE_OUTCOME"].to_numpy()
        
        # WALD RATIO
        THETA_I = betas_outcome / betas_exposure
        
        # variance of theta_i using second-order delta method
        # assuming no sample overlap and therefore covariance = 0
        VAR_THETA_I = (
            ((ses_outcome ** 2) / (betas_exposure ** 2)) + 
            (
                (betas_outcome ** 2) *
                (ses_exposure ** 2) /
                (betas_exposure ** 4)
            )
        )

        SE_THETA_I = np.sqrt(VAR_THETA_I)

        # inverse-variance weights
        W_I = 1 / VAR_THETA_I

        # IVW estimate
        BETA_IVW = (np.sum(W_I * THETA_I) / np.sum(W_I) )

        # variance and standard error
        VAR_BETA_IVW = 1 / np.sum(W_I)
        SE_BETA_IVW = np.sqrt(VAR_BETA_IVW)

        # hypothesis test with respect to Z-score
        Z_IVW = BETA_IVW / SE_BETA_IVW
        P_IVW = 2 * norm.sf(abs(Z_IVW))

        # 95% confidence interval
        CI_LOW_IVW = BETA_IVW - 1.96 * SE_BETA_IVW
        CI_HIGH_IVW = BETA_IVW + 1.96 * SE_BETA_IVW

        # SNP-specific estimates
        snp_results = pl.DataFrame({
            "SNP": snps,
            "BETA_EXPOSURE": betas_exposure,
            "SE_EXPOSURE": ses_exposure,
            "BETA_OUTCOME": betas_outcome,
            "SE_OUTCOME": ses_outcome,
            "WALD_RATIO": THETA_I,
            "SE_WALD_RATIO": SE_THETA_I,
            "VAR_WALD_RATIO": VAR_THETA_I,
            "WEIGHT": W_I
        })

        # pooled IVW result
        ivw_result = pl.DataFrame({
            "METHOD": ["IVW_delta"],
            "N_SNPS": [df.height],
            "BETA": [BETA_IVW],
            "SE": [SE_BETA_IVW],
            "Z": [Z_IVW],
            "P": [P_IVW],
            "CI_LOW": [CI_LOW_IVW],
            "CI_HIGH": [CI_HIGH_IVW]
        })

        return ivw_result, snp_results