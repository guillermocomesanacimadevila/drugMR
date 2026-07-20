#!/usr/bin/env python3
import numpy as np
from scipy.stats import norm 

# pheno_id
# mediator_id
# ivw_mr_X_M_Y_results (fixed dir)
# pqtl_dataset
# B_XY -> Total effect of X on Y (including known and unknown Ms)
# B_D -> Pure effect of X on Y (excluding known and unknown Ms)

# network MR workflow
# X = protein (e.g. PU.1)
# M = mediator (e.g. pTau)
# Y = outcome (e.g. AD)
# 1. TwoSampleMR - X -> M && M -> Y
# 2. Calculate indirect effect (B_I = B_XM * B_MY) - hyp test under gaussian
# 3. Proportion of XY effect -> B_I / B_XY
# 4. If PROP == HIGH + B_I p < 0.05 -> MOLOC at pQTL locus (X - M - Y)

# DS NetworkMR pipeline
# dictionary in jupyter notebook {M_id: 'User/Path/...'}
# FROM NOTEBOOK -> MAKE A MEDIATOR MANIFEST
# For each protein part of dataset X
# Run cis-MR (twice) -> For each X -> M
# Also run M -> Y (whole genome) 
# results/networkMR/ 3 subdirectories
# results/networkMR/M_Y/....csv (Genome-wide - one CSV with MR outputs where 1 entry == univariable MR from a mediator M on Y)
# results/networkMR/X_M/mediator_1/....csv (1 entry == univariable cis-MR - 1 protein vs that mediator)
# results/networkMR/X_M/mediator_2/....csv (1 entry == univariable cis-MR - 1 protein vs that mediator)
# results/networkMR/X_M/mediator_N/....csv (1 entry == univariable cis-MR - 1 protein vs that mediator)
# results/networkMR/mediation_estimates/...csv (massive CSV with a given protein that FDR significant in X->M and X->Y and also if IVW_p < 0.05 in X->Y run NetworkMR package - here the output of NetworkMR package)

# MR object

# IVW for  
# class PyTwoSampleMR(self):


# CREATE A MEDIATOR MANIFEST FUNCTION

def NetworkMR(
    B_XM: float,   # protein -> mediator
    SE_XM: float,
    B_XY: float,   # protein -> outcome = total effect
    SE_XY: float,
    B_MY: float,   # mediator -> outcome
    SE_MY: float,
):
    # indirect effect
    B_I = B_XM * B_MY

    # delta-method SE for indirect effect
    SE_I = np.sqrt((B_MY ** 2) * (SE_XM ** 2) + (B_XM ** 2) * (SE_MY ** 2))

    Z_I = B_I / SE_I
    P_I = 2 * norm.sf(abs(Z_I))

    # direct effect
    B_D = B_XY - B_I

    # approximate SE for direct effect
    SE_D = np.sqrt(SE_XY ** 2 + SE_I ** 2)

    Z_D = B_D / SE_D
    P_D = 2 * norm.sf(abs(Z_D))

    # proportion mediated
    P_M = B_I / B_XY if B_XY != 0 else np.nan

    SE_PM = np.sqrt(
        (SE_I ** 2 / B_XY ** 2) +
        ((B_I ** 2 * SE_XY ** 2) / B_XY ** 4)
    ) if B_XY != 0 else np.nan

    Z_PM = P_M / SE_PM if B_XY != 0 else np.nan
    P_PM = 2 * norm.sf(abs(Z_PM)) if B_XY != 0 else np.nan

    # CIs for direct and indirect effects
    # indirect effect
    CI_LOW_I = B_I - 1.96 * SE_I
    CI_HIGH_I = B_I + 1.96 * SE_I

    # direct effect
    CI_LOW_D = B_D - 1.96 * SE_D
    CI_HIGH_D = B_D + 1.96 * SE_D

    # proportion mediated
    CI_LOW_PM = P_M - 1.96 * SE_PM if B_XY != 0 else np.nan
    CI_HIGH_PM = P_M + 1.96 * SE_PM if B_XY != 0 else np.nan

    # directionality checks
    consistent_direction = np.sign(B_I) == np.sign(B_XY)

    # prints here:
    print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
    print("NETWORK MENDELIAN RANDOMISATION")
    print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")

    print(f"Causal effect (BETA) between exposure X and outcome Y   = {B_XY:.3f}")
    print(f"Causal effect (BETA) between exposure X and mediator M  = {B_XM:.3f}")
    print(f"Causal effect (BETA) between mediator M and outcome Y   = {B_MY:.3f}")

    print("----------------------------------------------------------------------")

    print("Indirect causal effect (X → M → Y)")
    print(f"  Beta      : {B_I:.3f}")
    print(f"  SE        : {SE_I:.3f}")
    print(f"  95% CI    : ({CI_LOW_I:.3f}, {CI_HIGH_I:.3f})")
    print(f"  Z-score   : {Z_I:.3f}")
    print(f"  P-value   : {P_I:.3e}")

    print("----------------------------------------------------------------------")

    print("Direct causal effect (X → Y)")
    print(f"  Beta      : {B_D:.3f}")
    print(f"  SE        : {SE_D:.3f}")
    print(f"  95% CI    : ({CI_LOW_D:.3f}, {CI_HIGH_D:.3f})")
    print(f"  Z-score   : {Z_D:.3f}")
    print(f"  P-value   : {P_D:.3e}")

    print("----------------------------------------------------------------------")

    print("Proportion mediated")
    print(f"  Estimate  : {P_M:.3f}")
    print(f"  Percent   : {P_M * 100:.1f}%")
    print(f"  SE        : {SE_PM:.3f}")
    print(f"  95% CI    : ({CI_LOW_PM:.3f}, {CI_HIGH_PM:.3f})")
    print(f"  Z-score   : {Z_PM:.3f}")
    print(f"  P-value   : {P_PM:.3e}")

    print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")

    print(f"Consistent direction? => {consistent_direction} ")

    print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")


    return {
        "beta_total": B_XY,
        "se_total": SE_XY,
        "beta_indirect": B_I,
        "se_indirect": SE_I,
        "ci_low_indirect": CI_LOW_I,
        "ci_high_indirect": CI_HIGH_I,
        "z_indirect": Z_I,
        "p_indirect": P_I,
        "beta_direct": B_D,
        "se_direct": SE_D,
        "ci_low_direct": CI_LOW_D,
        "ci_high_direct": CI_HIGH_D,
        "z_direct": Z_D,
        "p_direct": P_D,
        "prop_mediated": P_M,
        "percent_mediated": P_M * 100,
        "se_prop_mediated": SE_PM,
        "ci_low_prop_mediated": CI_LOW_PM,
        "ci_high_prop_mediated": CI_HIGH_PM,
        "z_prop_mediated": Z_PM,
        "p_prop_mediated": P_PM,
        "consistent_direction": consistent_direction
    }