import numpy as np
from scipy.stats import norm

# pheno_id
# mediator_id
# ivw_mr_X_M_Y_results (fixed dir)
# pqtl_dataset
# B_XY -> Total effect of X on Y (including known and unknown Ms)
# B_D -> Pure effect of X on Y (excluding known and unknown Ms)

# network MR workflow
# X = protein (e.g. PU.1)
# M = mediator (e.g. pTau)
# Y = outcome (e.g. AD)
# 1. TwoSampleMR - X -> M && M -> Y
# 2. Calculate indirect effect (B_I = B_XM * B_MY) - hyp test under gaussian
# 3. Proportion of XY effect -> B_I / B_XY
# 4. If PROP == HIGH + B_I p < 0.05 -> MOLOC at pQTL locus (X - M - Y)

# DS NetworkMR pipeline
# dictionary in jupyter notebook {M_id: 'User/Path/...'}
# FROM NOTEBOOK -> MAKE A MEDIATOR MANIFEST
# For each protein part of dataset X
# Run cis-MR (twice) -> For each X -> M
# Also run M -> Y (whole genome)
# results/networkMR/ 3 subdirectories
# results/networkMR/M_Y/....csv (Genome-wide - one CSV with MR outputs where 1 entry == univariable MR from a mediator M on Y)
# results/networkMR/X_M/mediator_1/....csv (1 entry == univariable cis-MR - 1 protein vs that mediator)
# results/networkMR/X_M/mediator_2/....csv (1 entry == univariable cis-MR - 1 protein vs that mediator)
# results/networkMR/X_M/mediator_N/....csv (1 entry == univariable cis-MR - 1 protein vs that mediator)
# results/networkMR/mediation_estimates/...csv (massive CSV with a given protein that FDR significant in X->M and X->Y and also if IVW_p < 0.05 in X->Y run NetworkMR package - here the output of NetworkMR package)


class NetworkMR:

    """ Object to estimate indirect/direct/mediated effects for a
    protein (X) -> mediator (M) -> outcome (Y) triple """

    def __init__(
            self,
            B_XM: float,   # protein -> mediator
            SE_XM: float,
            B_XY: float,   # protein -> outcome = total effect
            SE_XY: float,
            B_MY: float,   # mediator -> outcome
            SE_MY: float):

        self.B_XM = B_XM
        self.SE_XM = SE_XM
        self.B_XY = B_XY
        self.SE_XY = SE_XY
        self.B_MY = B_MY
        self.SE_MY = SE_MY

        self.indirect = self.estimate_indirect_effect()
        self.direct = self.estimate_direct_effect()
        self.mediated = self.estimate_proportion_mediated()

    def estimate_indirect_effect(self) -> dict:
        B_I = self.B_XM * self.B_MY
        SE_I = np.sqrt((self.B_MY ** 2) * (self.SE_XM ** 2) + (self.B_XM ** 2) * (self.SE_MY ** 2))
        Z_I = B_I / SE_I
        P_I = 2 * norm.sf(abs(Z_I))

        return {
            "beta": B_I,
            "se": SE_I,
            "ci_low": B_I - 1.96 * SE_I,
            "ci_high": B_I + 1.96 * SE_I,
            "z": Z_I,
            "p": P_I,
        }

    def estimate_direct_effect(self) -> dict:
        B_I = self.indirect["beta"]
        SE_I = self.indirect["se"]

        B_D = self.B_XY - B_I
        SE_D = np.sqrt(self.SE_XY ** 2 + SE_I ** 2)
        Z_D = B_D / SE_D
        P_D = 2 * norm.sf(abs(Z_D))

        return {
            "beta": B_D,
            "se": SE_D,
            "ci_low": B_D - 1.96 * SE_D,
            "ci_high": B_D + 1.96 * SE_D,
            "z": Z_D,
            "p": P_D,
        }

    def estimate_proportion_mediated(self) -> dict:
        B_I = self.indirect["beta"]
        SE_I = self.indirect["se"]
        B_XY, SE_XY = self.B_XY, self.SE_XY

        if B_XY == 0:
            return {"beta": np.nan, "percent": np.nan, "se": np.nan, "ci_low": np.nan, "ci_high": np.nan, "z": np.nan, "p": np.nan}

        P_M = B_I / B_XY
        SE_PM = np.sqrt((SE_I ** 2 / B_XY ** 2) + ((B_I ** 2 * SE_XY ** 2) / B_XY ** 4))
        Z_PM = P_M / SE_PM
        P_PM = 2 * norm.sf(abs(Z_PM))

        return {
            "beta": P_M,
            "percent": P_M * 100,
            "se": SE_PM,
            "ci_low": P_M - 1.96 * SE_PM,
            "ci_high": P_M + 1.96 * SE_PM,
            "z": Z_PM,
            "p": P_PM,
        }

    def is_consistent_direction(self) -> bool:
        return np.sign(self.indirect["beta"]) == np.sign(self.B_XY)

    def report(self) -> dict:
        I, D, M = self.indirect, self.direct, self.mediated
        consistent_direction = self.is_consistent_direction()

        print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
        print("NETWORK MENDELIAN RANDOMISATION")
        print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")

        print(f"Causal effect (BETA) between exposure X and outcome Y   = {self.B_XY:.3f}")
        print(f"Causal effect (BETA) between exposure X and mediator M  = {self.B_XM:.3f}")
        print(f"Causal effect (BETA) between mediator M and outcome Y   = {self.B_MY:.3f}")

        print("----------------------------------------------------------------------")

        print("Indirect causal effect (X → M → Y)")
        print(f"  Beta      : {I['beta']:.3f}")
        print(f"  SE        : {I['se']:.3f}")
        print(f"  95% CI    : ({I['ci_low']:.3f}, {I['ci_high']:.3f})")
        print(f"  Z-score   : {I['z']:.3f}")
        print(f"  P-value   : {I['p']:.3e}")

        print("----------------------------------------------------------------------")

        print("Direct causal effect (X → Y)")
        print(f"  Beta      : {D['beta']:.3f}")
        print(f"  SE        : {D['se']:.3f}")
        print(f"  95% CI    : ({D['ci_low']:.3f}, {D['ci_high']:.3f})")
        print(f"  Z-score   : {D['z']:.3f}")
        print(f"  P-value   : {D['p']:.3e}")

        print("----------------------------------------------------------------------")

        print("Proportion mediated")
        print(f"  Estimate  : {M['beta']:.3f}")
        print(f"  Percent   : {M['percent']:.1f}%")
        print(f"  SE        : {M['se']:.3f}")
        print(f"  95% CI    : ({M['ci_low']:.3f}, {M['ci_high']:.3f})")
        print(f"  Z-score   : {M['z']:.3f}")
        print(f"  P-value   : {M['p']:.3e}")

        print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")

        print(f"Consistent direction? => {consistent_direction} ")

        print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")

        return {
            "beta_total": self.B_XY,
            "se_total": self.SE_XY,
            "beta_indirect": I["beta"],
            "se_indirect": I["se"],
            "ci_low_indirect": I["ci_low"],
            "ci_high_indirect": I["ci_high"],
            "z_indirect": I["z"],
            "p_indirect": I["p"],
            "beta_direct": D["beta"],
            "se_direct": D["se"],
            "ci_low_direct": D["ci_low"],
            "ci_high_direct": D["ci_high"],
            "z_direct": D["z"],
            "p_direct": D["p"],
            "prop_mediated": M["beta"],
            "percent_mediated": M["percent"],
            "se_prop_mediated": M["se"],
            "ci_low_prop_mediated": M["ci_low"],
            "ci_high_prop_mediated": M["ci_high"],
            "z_prop_mediated": M["z"],
            "p_prop_mediated": M["p"],
            "consistent_direction": consistent_direction,
        }
