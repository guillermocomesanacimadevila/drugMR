import polars as pl

from drugmr.utils import (
    lambda_sample_overlap,
    quick_f_statistic,
    sample_overlap_relative_bias,
)


class SampleOverlapBias:

    """
    Object to flag cis-MR / PheWAS targets at risk of two-sample MR sample-overlap bias
    """

    def __init__(self):
        pass

    @staticmethod
    def per_protein_f_statistic(mr_instruments: pl.DataFrame) -> pl.DataFrame:
        return (
            mr_instruments
            .with_columns(
                quick_f_statistic(pl.col("beta.exposure"), pl.col("se.exposure")).alias("f_statistic")
            )
            .group_by("protein")
            .agg(pl.col("f_statistic").mean())
        )

    @staticmethod
    def relative_bias(mr_instruments: pl.DataFrame, n_overlap: int, n_exposure_total: int, n_outcome_total: int,) -> pl.DataFrame:
        f_stats = SampleOverlapBias.per_protein_f_statistic(mr_instruments)
        lambda_val = float(lambda_sample_overlap(n_overlap, n_exposure_total, n_outcome_total))
        return f_stats.with_columns(
            pl.lit(lambda_val).alias("lambda_overlap"),
            sample_overlap_relative_bias(lambda_val, pl.col("f_statistic")).alias("relative_bias_pct"),
        )
