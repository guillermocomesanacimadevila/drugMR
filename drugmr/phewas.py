#!/usr/bin/env python3
import argparse
import polars as pl
import numpy as np 
from scipy.stats import norm

# WE ARE NOT USING UKB PheWAS
# AS WALD RATIO H:0 TESTING -> WITH DELTA METHOD 
# SO Cov(B_JX, B_JY) != 0
# WE SHALL USE FinnGen then 

# Receive an input 
# Calc wald ratio
# Delta method 
# FDR correction here? - Maybe as a bool? 

# WE ASSUME NO SAMPLE OVERLAP 
def PheWAS(B_X, B_Y, SE_X, SE_Y):

    # Wald ratio
    B_XY = B_Y / B_X

    # Variance of the Wald ratio (delta method)
    VAR_XY = ((SE_Y ** 2) / (B_X ** 2) + ((B_Y ** 2) * (SE_X ** 2)) / (B_X ** 4))  # assumes Cov(B_X, B_Y) == 0

    # Standard error of the Wald ratio
    SE_XY = np.sqrt(VAR_XY)

    # Z-score it 
    Z_XY = B_XY / SE_XY

    # H0: Wald ratio == 0, assuming Z ~ N(0, 1)
    P = 2 * norm.sf(abs(Z_XY))

    # GENERATE SOME PRINTS HERE 
    print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
    print("PheWAS MENDELIAN RANDOMISATION")
    print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")

    print(f"Causal effect (BETA) between exposure X and outcome Y   = {B_XY:.3f}")
    print(f"Variance of B_XY (between exposure X and outcome Y) = {VAR_XY:.3f}")
    print(f"Standard error of B_XY (between exposure X and outcome Y) = {SE_XY:.3f}")
    print(f"Z-score of B_XY (between exposure X and outcome Y) = {Z_XY:.3f}")
    print(f"P-value of B_XY (between exposure X and outcome Y) = {P:.3f}")

    return {
        "wald_ratio": B_XY,
        "se_wald_ratio": SE_XY,
        "z_score": Z_XY,
        "P_nominal": P
    }