import numpy as np
import polars as pl
from scipy.stats import chi2, norm


class PyTwoSampleMR:

    """
    Assortment of Mendelian Randomisation functions for generalised usability in drugMR pipeline
    # 1. IVW + DELTA METHOD 
    """

    def __init__(self):
        pass

    def wald_ratio(
            self,
            exposure_df,
            outcome_df,
            exposure_snp_col,
            exposure_beta_col,
            exposure_se_col,
            outcome_snp_col,
            outcome_beta_col,
            outcome_se_col
    ):  # single-instrument Wald ratio + delta method

        # compile inputs and match SNPs between exposure and outcome
        exposure = exposure_df.select(pl.col(exposure_snp_col).cast(pl.Utf8).alias("SNP"), pl.col(exposure_beta_col).cast(pl.Float64).alias("BETA_EXPOSURE"), pl.col(exposure_se_col).cast(pl.Float64).alias("SE_EXPOSURE"))
        outcome = outcome_df.select(pl.col(outcome_snp_col).cast(pl.Utf8).alias("SNP"), pl.col(outcome_beta_col).cast(pl.Float64).alias("BETA_OUTCOME"), pl.col(outcome_se_col).cast(pl.Float64).alias("SE_OUTCOME"))
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

        # a Wald ratio uses exactly one instrument; use IVW() for more
        if df.height > 1:
            raise ValueError(
                f"Wald ratio needs exactly 1 SNP, got {df.height}; use IVW() instead."
            )

        snp            = df["SNP"][0]
        beta_exposure  = df["BETA_EXPOSURE"][0]
        se_exposure    = df["SE_EXPOSURE"][0]
        beta_outcome   = df["BETA_OUTCOME"][0]
        se_outcome     = df["SE_OUTCOME"][0]

        # WALD RATIO
        BETA_WR = beta_outcome / beta_exposure

        # variance using second-order delta method (same as IVW())
        # assuming no sample overlap and therefore covariance = 0
        VAR_WR = (
            ((se_outcome ** 2) / (beta_exposure ** 2)) +
            (
                (beta_outcome ** 2) *
                (se_exposure ** 2) /
                (beta_exposure ** 4)
            )
        )
        SE_WR = np.sqrt(VAR_WR)

        # hypothesis test with respect to Z-score
        Z_WR = BETA_WR / SE_WR
        P_WR = 2 * norm.sf(abs(Z_WR))

        # 95% confidence interval
        CI_LOW_WR = BETA_WR - 1.96 * SE_WR
        CI_HIGH_WR = BETA_WR + 1.96 * SE_WR

        # same layout as IVW()'s pooled result
        wald_result = pl.DataFrame({
            "METHOD": ["Wald_ratio_delta"],
            "SNP": [snp],
            "N_SNPS": [1],
            "BETA": [BETA_WR],
            "SE": [SE_WR],
            "Z": [Z_WR],
            "P": [P_WR],
            "CI_LOW": [CI_LOW_WR],
            "CI_HIGH": [CI_HIGH_WR]
        })
        return wald_result

    # standard IVW
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
        exposure = exposure_df.select(pl.col(exposure_snp_col).cast(pl.Utf8).alias("SNP"), pl.col(exposure_beta_col).cast(pl.Float64).alias("BETA_EXPOSURE"), pl.col(exposure_se_col).cast(pl.Float64).alias("SE_EXPOSURE"))
        outcome = outcome_df.select(pl.col(outcome_snp_col).cast(pl.Utf8).alias("SNP"), pl.col(outcome_beta_col).cast(pl.Float64).alias("BETA_OUTCOME"), pl.col(outcome_se_col).cast(pl.Float64).alias("SE_OUTCOME"))
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

        # first-order weights (1/se_outcome^2), same as TwoSampleMR mr_heterogeneity,
        N_SNPS = df.height
        if N_SNPS >= 2:
            W_FIRST = 1 / (ses_outcome ** 2)
            BETA_FIRST = np.sum(W_FIRST * betas_exposure * betas_outcome) / np.sum(W_FIRST * betas_exposure ** 2)
            Q = np.sum(W_FIRST * (betas_outcome - BETA_FIRST * betas_exposure) ** 2)
            Q_DF = N_SNPS - 1
            Q_P = chi2.sf(Q, Q_DF)
        else:
            # Q does not exist with one instrument
            Q = None
            Q_DF = None
            Q_P = None

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
            "N_SNPS": [N_SNPS],
            "BETA": [BETA_IVW],
            "SE": [SE_BETA_IVW],
            "Z": [Z_IVW],
            "P": [P_IVW],
            "CI_LOW": [CI_LOW_IVW],
            "CI_HIGH": [CI_HIGH_IVW],
            "Q": [Q],
            "Q_DF": [Q_DF],
            "Q_P": [Q_P]
        })
        return ivw_result, snp_results
