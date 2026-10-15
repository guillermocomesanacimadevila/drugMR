import json
import math

import pandas as pd
import polars as pl
import yaml

from drugmr import paths
from drugmr.utils import extract_coloc_or_pwcoco_targets, grab_cis_mr_hits


# rows shown in the default preview of a long table (half at the top, half at the bottom)
PREVIEW_ROWS = 10
# tables up to this many rows are always shown in full
FULL_TABLE_ROWS = 20


def _format_number(value):

    """ 3 significant figures: 8.73e-26, 0.774, 1.22 """

    return f"{value:.3g}"


class ResultsTable(pd.DataFrame):

    """
    A pandas DataFrame that displays cleanly in a notebook: numbers to 3
    significant figures and every column. Long tables show a preview of their
    rows unless show_all is True. The data itself keeps full precision.
    """

    _metadata = ["show_all"]
    show_all = False

    @property
    def _constructor(self):
        return ResultsTable

    def _display_settings(self):
        if self.show_all:
            max_rows = None
        else:
            max_rows = FULL_TABLE_ROWS
        return [
            "display.float_format", _format_number,
            "display.max_columns", None,
            "display.max_rows", max_rows,
            "display.min_rows", PREVIEW_ROWS,
            "display.max_colwidth", 60,
        ]

    def _repr_html_(self):
        with pd.option_context(*self._display_settings()):
            return super()._repr_html_()

    def __repr__(self):
        with pd.option_context(*self._display_settings()):
            return super().__repr__()


class ResultsExplorer:

    """ Object compiling functions for sample notebook """

    def __init__(
            self,
            run_id: str,
            root="runs"
    ):

        self.run_id = run_id
        self.root = root
        self.results_dir = paths.run_results_dir(run_id, root)

        with open(paths.run_manifest_path(run_id, root)) as fh:
            self.manifest = json.load(fh)

        with open(paths.run_params_lock_path(run_id, root)) as fh:
            self.params = yaml.safe_load(fh)

        self.pheno_id = self.params["pheno_id"]
        self.pqtl_dataset = self.params["pqtl_dataset"]

    def run_info(self):

        """ One table describing what this run is and the thresholds it used """

        manifest = self.manifest
        params = self.params
        cis_mr = params["gates"]["cis_mr"]
        coloc = params["gates"]["coloc"]

        # (setting, value)
        rows = []

        # identity
        rows.append(("Run ID", self.run_id))
        rows.append(("Outcome", self.pheno_id))
        rows.append(("pQTL dataset", self.pqtl_dataset))
        rows.append(("Status", manifest.get("status")))
        rows.append(("Created", manifest.get("created_at")))

        # provenance: what code and container produced it
        git_sha = manifest.get("git_sha7")
        if manifest.get("git_dirty"):
            git_sha = f"{git_sha} (uncommitted changes)"
        rows.append(("Git commit", git_sha))
        # container: image name plus a short digest, enough to identify it
        container = manifest.get("image_name")
        if container and "@sha256:" in container:
            image, digest = container.split("@sha256:")
            image_name = image.split("/")[-1]
            container = f"{image_name} (sha256:{digest[:12]})"
        rows.append(("Container", container))
        rows.append(("Host", manifest.get("host")))

        # outcome GWAS inputs
        rows.append(("Cases", params.get("n_cases")))
        rows.append(("Controls", params.get("n_controls")))
        rows.append(("Genome build", params.get("target_build")))
        rows.append(("MAF filter", params.get("maf")))
        rows.append(("MHC removed", params.get("remove_mhc")))
        rows.append(("APOE removed", params.get("remove_apoe")))
        rows.append(("LD reference", params.get("ref_bfile")))

        # cis-MR gates
        rows.append(("Instrument P threshold", cis_mr.get("instrument_pval_threshold")))
        rows.append(("Clumping r2 / kb", f"{cis_mr.get('clump_r2')} / {cis_mr.get('clump_kb')}"))
        rows.append(("Minimum F statistic", cis_mr.get("min_f_stat")))
        rows.append(("Steiger filter", cis_mr.get("apply_steiger_filter")))
        rows.append(("Wald ratio FDR q", cis_mr.get("wald_fdr_q")))
        rows.append(("IVW FDR q", cis_mr.get("ivw_fdr_q")))
        rows.append(("Minimum instruments for IVW", cis_mr.get("min_instruments_for_ivw")))
        rows.append(("Cochran's Q P (heterogeneity)", cis_mr.get("cochran_q_pval")))

        # colocalisation gates
        rows.append(("PP4 threshold", coloc.get("pp4_threshold")))
        rows.append(("Priors p1 / p2 / p12", f"{coloc.get('p1')} / {coloc.get('p2')} / {coloc.get('p12')}"))

        # xQTL triangulation
        if params.get("run_smr"):
            bulk = ", ".join(params.get("bulk_qtl_datasets") or [])
            rows.append(("Bulk QTL datasets", bulk))
            rows.append(("Single cell QTL dataset", params.get("sc_qtl_dataset")))
        else:
            rows.append(("SMR / HyPrColoc", "not run"))

        # values are mixed types, so show them all as text
        settings = []
        values = []
        for setting, value in rows:
            settings.append(setting)
            values.append(str(value))

        table = pl.DataFrame(
            {
                "setting": settings,
                "value": values,
            }
        )
        # pandas renders as a clean table in Jupyter, one row per setting
        table = ResultsTable(table.to_pandas().set_index("setting"))
        table.index.name = None
        table.show_all = True
        return table

    def funnel(self):

        """ How many proteins survive each stage of this run """

        dataset = self.pqtl_dataset
        pheno = self.pheno_id
        results = self.results_dir
        pp4 = self.params["gates"]["coloc"]["pp4_threshold"]

        # (stage, n, base count the share is taken from, note)
        rows = []

        # cis-MR: everything tested
        mr_file = paths.mr_out(dataset, pheno, results)
        mr = pl.read_csv(mr_file, separator="\t", infer_schema_length=None)
        n_tested = mr.height
        n_wald = mr.filter(pl.col("primary_method") == "Wald ratio").height
        note = f"{n_wald} Wald ratio, {n_tested - n_wald} multi-instrument"
        rows.append(("Tested in cis-MR", n_tested, None, note))

        # instruments behind them (a SNP count, so no share)
        instruments = pl.read_csv(paths.mr_instruments_out(dataset, pheno, results), separator="\t", infer_schema_length=None)
        instruments = instruments.filter(pl.col("mr_keep"))
        per_protein = instruments.group_by("protein").len()
        note = f"median {per_protein['len'].median():.0f} per protein, min F {instruments['F'].min():.1f}"
        rows.append(("Instruments (SNPs)", instruments.height, None, note))

        # cis-MR gate, exactly as the pipeline applies it
        cis_hits = self._cis_mr_hits()
        rows.append(("Passed cis-MR", len(cis_hits), n_tested, "of tested, FDR and heterogeneity gates"))

        # colocalisation among cis-MR hits
        coloc_file = paths.coloc_out(dataset, pheno, results)
        pwcoco_file = paths.pwcoco_out(dataset, pheno, results)
        coloc_hits = set(extract_coloc_or_pwcoco_targets(coloc_file, pwcoco_file, pp4, ("coloc",)))
        pwcoco_hits = set(extract_coloc_or_pwcoco_targets(coloc_file, pwcoco_file, pp4, ("pwcoco",)))
        coloc_hits = coloc_hits & cis_hits
        pwcoco_hits = pwcoco_hits & cis_hits
        prioritised = coloc_hits | pwcoco_hits
        rows.append(("Colocalised (COLOC)", len(coloc_hits), len(cis_hits), f"of cis-MR hits, PP4 >= {pp4}"))
        rows.append(("Colocalised (PWCoCo)", len(pwcoco_hits), len(cis_hits), f"of cis-MR hits, PP4 >= {pp4}"))
        rows.append(("Prioritised", len(prioritised), len(cis_hits), "of cis-MR hits, COLOC or PWCoCo"))

        # xQTL support among prioritised, bulk tissue and single cell kept apart
        # (data_type column: "bulk" or "single_cell")
        smr_file = paths.smr_final_targets_out(dataset, pheno, results)
        if smr_file.exists():
            smr = pl.read_csv(smr_file, separator="\t", infer_schema_length=None)
            smr_bulk = set(smr.filter(pl.col("data_type") == "bulk")["protein"].to_list())
            smr_sc = set(smr.filter(pl.col("data_type") == "single_cell")["protein"].to_list())
            smr_any = smr_bulk | smr_sc
            rows.append(("SMR support (bulk)", len(smr_bulk & prioritised), len(prioritised), "of prioritised, bulk tissue eQTL"))
            rows.append(("SMR support (single cell)", len(smr_sc & prioritised), len(prioritised), "of prioritised, cell type eQTL"))
            rows.append(("SMR support (either)", len(smr_any & prioritised), len(prioritised), "of prioritised"))
        else:
            rows.append(("SMR support", None, None, "not run"))

        # HyPrColoc only runs on SMR targets, so a row alone is not support:
        # same rule as the dashboard, one cluster holding pQTL + GWAS + an xQTL
        # together with posterior >= 0.7
        hyprcoloc_file = paths.hyprcoloc_out(dataset, pheno, results)
        if hyprcoloc_file.exists():
            hyprcoloc = pl.read_csv(hyprcoloc_file, separator="\t", infer_schema_length=None)
            hyprcoloc_bulk = set()
            hyprcoloc_sc = set()
            for row in hyprcoloc.iter_rows(named=True):
                posterior = row["posterior_prob"]
                if posterior is None or posterior < 0.7:
                    continue
                has_pqtl = False
                has_gwas = False
                has_xqtl = False
                for trait in str(row["traits"]).lower().split(","):
                    trait = trait.strip()
                    if trait.startswith("pqtl_"):
                        has_pqtl = True
                    elif trait.startswith("gwas_"):
                        has_gwas = True
                    elif "qtl_" in trait:
                        has_xqtl = True
                if has_pqtl and has_gwas and has_xqtl:
                    if row["data_type"] == "single_cell":
                        hyprcoloc_sc.add(row["protein"])
                    else:
                        hyprcoloc_bulk.add(row["protein"])
            hyprcoloc_any = hyprcoloc_bulk | hyprcoloc_sc
            note = "of prioritised, pQTL + GWAS + xQTL, posterior >= 0.7"
            rows.append(("HyPrColoc support (bulk)", len(hyprcoloc_bulk & prioritised), len(prioritised), note))
            rows.append(("HyPrColoc support (single cell)", len(hyprcoloc_sc & prioritised), len(prioritised), note))
            rows.append(("HyPrColoc support (either)", len(hyprcoloc_any & prioritised), len(prioritised), note))
        else:
            rows.append(("HyPrColoc support", None, None, "not run"))

        # PheWAS: prioritised proteins with any Bonferroni significant hit (FinnGen or UKB)
        phewas_files = [paths.phewas_out(dataset, pheno, results), paths.phewas_ukbb_out(dataset, pheno, results)]
        phewas_ran = False
        flagged = set()
        for phewas_file in phewas_files:
            if phewas_file.exists():
                phewas_ran = True
                phewas = pl.read_csv(phewas_file, separator="\t", infer_schema_length=None)
                phewas = phewas.filter(pl.col("bonferroni_significant"))
                flagged = flagged | set(phewas["protein"].to_list())
        if phewas_ran:
            rows.append(("PheWAS hit", len(flagged & prioritised), len(prioritised), "of prioritised, any Bonferroni significant trait"))
        else:
            rows.append(("PheWAS hit", None, None, "not run"))

        stages = []
        counts = []
        shares = []
        notes = []
        for stage, n, base, note in rows:
            stages.append(stage)
            notes.append(note)
            if n is None:
                counts.append("")
            else:
                counts.append(str(n))
            if n is not None and base:
                shares.append(f"{100 * n / base:.0f}%")
            else:
                shares.append("")

        table = pl.DataFrame(
            {
                "stage": stages,
                "n": counts,
                "share": shares,
                "note": notes,
            }
        )
        # pandas renders as a clean table in Jupyter, one row per stage
        table = ResultsTable(table.to_pandas().set_index("stage"))
        table.index.name = None
        table.show_all = True
        return table

    def cis_mr(self, protein=None, hits_only=True, detailed=False):

        """
        One row per protein: the cis-MR estimate on the outcome.

        protein: optional symbol (CTSH) or full protein id to show just that protein
        hits_only: True shows only proteins passing the run's cis-MR gate
        """

        mr_file = paths.mr_out(self.pqtl_dataset, self.pheno_id, self.results_dir)
        mr = pl.read_csv(mr_file, separator="\t", infer_schema_length=None)
        cis_hits = self._cis_mr_hits()

        records = []
        for row in mr.iter_rows(named=True):
            protein_id = row["protein"]

            if hits_only and protein_id not in cis_hits:
                continue
            if protein is not None and not self._matches(protein_id, protein):
                continue

            beta = row["primary_beta"]
            se = row["primary_se"]

            # binary outcome: beta is a log odds ratio per unit protein
            odds_ratio = math.exp(beta)
            or_low = math.exp(beta - 1.96 * se)
            or_high = math.exp(beta + 1.96 * se)

            # higher protein raising risk means a drug should lower it (inhibit)
            if beta > 0:
                effect = "higher protein raises risk"
                drug_direction = "inhibit"
            else:
                effect = "higher protein lowers risk"
                drug_direction = "activate"

            # heterogeneity only exists from 2 instruments, Egger intercept from 3
            n_instruments = row["n_instruments"]
            q_pval = row["Q_pval"] if n_instruments >= 2 else None
            egger_pval = row["egger_intercept_pval"] if n_instruments >= 3 else None

            records.append({
                "symbol": self._symbol(protein_id),
                "protein": protein_id,
                "method": row["primary_method"],
                "n_instruments": n_instruments,
                "beta": round(beta, 4),
                "se": round(se, 4),
                "OR": round(odds_ratio, 3),
                "OR 95% CI": f"{or_low:.3f} to {or_high:.3f}",
                "P": row["primary_pval"],
                "FDR q": row["primary_FDR_q"],
                "Cochran's Q P": q_pval,
                "Egger intercept P": egger_pval,
                "effect": effect,
                "drug direction": drug_direction,
                "passed cis-MR": protein_id in cis_hits,
            })

        if len(records) == 0:
            print("No proteins match.")
            return None

        # numeric even when every value is empty, so blanks always show the same way
        table = pl.DataFrame(
            records,
            infer_schema_length=None,
            schema_overrides={"Cochran's Q P": pl.Float64, "Egger intercept P": pl.Float64},
        )

        # strongest evidence first
        table = table.sort("P")

        # hits_only already means every row passed, so the flag adds nothing
        if hits_only:
            table = table.drop("passed cis-MR")

        table = self._finish(table, ["method", "n_instruments", "OR", "OR 95% CI", "P", "FDR q", "drug direction", "passed cis-MR"], detailed)
        return table

    def colocalisation(self, protein=None, hits_only=True, detailed=False):

        """
        One row per cis-MR hit: does the protein share a causal variant with the outcome?

        COLOC assumes one causal variant per trait; PWCoCo conditions on each
        independent signal and colocalises every pair, so it can rescue windows
        with more than one signal. Best PWCoCo row = highest PP4 over all pairs.

        protein: optional symbol (CTSH) or full protein id
        hits_only: True shows only proteins passing COLOC or PWCoCo (the prioritised set)
        """

        pp4_threshold = self.params["gates"]["coloc"]["pp4_threshold"]
        cis_hits = self._cis_mr_hits()

        coloc_file = paths.coloc_out(self.pqtl_dataset, self.pheno_id, self.results_dir)
        coloc = pl.read_csv(coloc_file, separator="\t", infer_schema_length=None)

        # COLOC: one row per protein
        coloc_by_protein = {}
        for row in coloc.iter_rows(named=True):
            coloc_by_protein[row["protein"]] = row

        # PWCoCo: many rows per protein (unconditioned + each conditioned pair)
        pwcoco_file = paths.pwcoco_out(self.pqtl_dataset, self.pheno_id, self.results_dir)
        pwcoco_unconditioned = {}
        pwcoco_best = {}
        pwcoco_n_pairs = {}
        if pwcoco_file.exists():
            pwcoco = pl.read_csv(pwcoco_file, separator="\t", infer_schema_length=None)
            for row in pwcoco.iter_rows(named=True):
                protein_id = row["protein"]
                is_unconditioned = row["SNP1"] == "unconditioned" and row["SNP2"] == "unconditioned"

                if is_unconditioned:
                    pwcoco_unconditioned[protein_id] = row["H4"]
                else:
                    pwcoco_n_pairs[protein_id] = pwcoco_n_pairs.get(protein_id, 0) + 1

                best = pwcoco_best.get(protein_id)
                if best is None or row["H4"] > best["H4"]:
                    pwcoco_best[protein_id] = row

        records = []
        for protein_id in sorted(cis_hits):

            if protein is not None and not self._matches(protein_id, protein):
                continue

            coloc_row = coloc_by_protein.get(protein_id)
            best_row = pwcoco_best.get(protein_id)

            coloc_pp4 = None
            coloc_pp3 = None
            coloc_snp = None
            coloc_nsnps = None
            if coloc_row is not None:
                coloc_pp4 = coloc_row["PP.H4.abf"]
                coloc_pp3 = coloc_row["PP.H3.abf"]
                coloc_snp = coloc_row["top_snp"]
                coloc_nsnps = coloc_row["nsnps"]

            best_pp4 = None
            best_pp3 = None
            best_pair = None
            best_kind = None
            if best_row is not None:
                best_pp4 = best_row["H4"]
                best_pp3 = best_row["H3"]
                best_pair = f"{best_row['SNP1']} / {best_row['SNP2']}"

                # SNP1 is the pQTL side, SNP2 the outcome side; either can stay unconditioned
                pqtl_conditioned = best_row["SNP1"] != "unconditioned"
                outcome_conditioned = best_row["SNP2"] != "unconditioned"
                if pqtl_conditioned and outcome_conditioned:
                    best_kind = "conditioned (both)"
                elif pqtl_conditioned:
                    best_kind = "conditioned (pQTL only)"
                elif outcome_conditioned:
                    best_kind = "conditioned (outcome only)"
                else:
                    best_kind = "unconditioned"

            passed_coloc = coloc_pp4 is not None and coloc_pp4 >= pp4_threshold
            passed_pwcoco = best_pp4 is not None and best_pp4 >= pp4_threshold

            # which method carries the colocalisation
            if passed_coloc and passed_pwcoco:
                supported_by = "COLOC and PWCoCo"
            elif passed_coloc:
                supported_by = "COLOC only"
            elif passed_pwcoco:
                supported_by = "PWCoCo only"
            else:
                supported_by = "not colocalised"

            if hits_only and not (passed_coloc or passed_pwcoco):
                continue

            records.append({
                "symbol": self._symbol(protein_id),
                "protein": protein_id,
                "COLOC PP4": coloc_pp4,
                "COLOC PP3": coloc_pp3,
                "COLOC top SNP": coloc_snp,
                "COLOC n SNPs": coloc_nsnps,
                "PWCoCo unconditioned PP4": pwcoco_unconditioned.get(protein_id),
                "PWCoCo best PP4": best_pp4,
                "PWCoCo best PP3": best_pp3,
                "PWCoCo best row": best_kind,
                "PWCoCo best pair (pQTL / outcome SNP)": best_pair,
                "PWCoCo conditioned pairs": pwcoco_n_pairs.get(protein_id, 0),
                "supported by": supported_by,
            })

        if len(records) == 0:
            print("No proteins match.")
            return None

        table = pl.DataFrame(
            records,
            infer_schema_length=None,
            schema_overrides={
                "COLOC PP4": pl.Float64,
                "COLOC PP3": pl.Float64,
                "PWCoCo unconditioned PP4": pl.Float64,
                "PWCoCo best PP4": pl.Float64,
                "PWCoCo best PP3": pl.Float64,
                "PWCoCo best pair (pQTL / outcome SNP)": pl.String,
            },
        )

        # strongest colocalisation first (best of the two methods)
        table = table.with_columns(
            pl.max_horizontal(pl.col("COLOC PP4").fill_null(0), pl.col("PWCoCo best PP4").fill_null(0)).alias("_best")
        )
        table = table.sort("_best", descending=True).drop("_best")

        # round probabilities for reading
        for column in ["COLOC PP4", "COLOC PP3", "PWCoCo unconditioned PP4", "PWCoCo best PP4", "PWCoCo best PP3"]:
            table = table.with_columns(pl.col(column).round(3))

        table = self._finish(table, ["COLOC PP4", "COLOC PP3", "PWCoCo best PP4", "PWCoCo best row", "supported by"], detailed)
        return table

    def smr(self, protein=None, data_type=None, detailed=False):

        """
        One row per protein x tissue or cell type with significant SMR and a passing HEIDI test.

        SMR tests whether the xQTL (molecular trait of the same gene) and the outcome share a causal
        variant; HEIDI checks it is one variant rather than two in LD. Direction compares
        the SMR effect (xQTL trait on outcome) with the cis-MR effect (protein on outcome).

        protein: optional symbol (CTSH) or full protein id
        data_type: None for all, "bulk" for tissue eQTL, "single_cell" for cell type eQTL
        """

        smr_file = paths.smr_final_targets_out(self.pqtl_dataset, self.pheno_id, self.results_dir)
        if not smr_file.exists():
            print("SMR was not run for this run.")
            return None

        smr = pl.read_csv(smr_file, separator="\t", infer_schema_length=None)

        # cis-MR effect of each protein, to compare directions
        mr_file = paths.mr_out(self.pqtl_dataset, self.pheno_id, self.results_dir)
        mr = pl.read_csv(mr_file, separator="\t", infer_schema_length=None)
        mr_beta = {}
        for row in mr.iter_rows(named=True):
            mr_beta[row["protein"]] = row["primary_beta"]

        records = []
        for row in smr.iter_rows(named=True):
            protein_id = row["protein"]

            if protein is not None and not self._matches(protein_id, protein):
                continue
            if data_type is not None and row["data_type"] != data_type:
                continue

            # only eQTL effects share the protein's direction;
            # splicing QTL betas are not comparable (see the sQTL note in the docs)
            protein_beta = mr_beta.get(protein_id)
            if row["qtl_type"] != "eqtl" or protein_beta is None:
                direction = "not comparable"
            elif (row["b_SMR"] > 0) == (protein_beta > 0):
                direction = "agrees"
            else:
                direction = "opposes"

            if row["data_type"] == "single_cell":
                kind = "single cell"
            else:
                kind = "bulk"

            records.append({
                "symbol": self._symbol(protein_id),
                "protein": protein_id,
                "data type": kind,
                "QTL dataset": row["qtl_dataset"],
                "tissue / cell type": row["cell_type"],
                "top SNP": row["top_snp"],
                "b_SMR": round(row["b_SMR"], 4),
                "se_SMR": round(row["se_SMR"], 4),
                "P_SMR": row["p_SMR"],
                "q_SMR": row["q_SMR"],
                "P_HEIDI": row["p_HEIDI"],
                "HEIDI SNPs": row["nsnp_HEIDI"],
                "cis-MR beta": None if protein_beta is None else round(protein_beta, 4),
                "direction vs cis-MR": direction,
            })

        if len(records) == 0:
            print("No SMR results match.")
            return None

        table = pl.DataFrame(records, infer_schema_length=None)

        # group each protein's tissues together, strongest first
        table = table.sort(["symbol", "P_SMR"])
        table = self._finish(table, ["data type", "QTL dataset", "tissue / cell type", "top SNP", "b_SMR", "P_SMR", "P_HEIDI", "direction vs cis-MR"], detailed)
        return table

    def hyprcoloc(self, protein=None, data_type=None, supported_only=True, detailed=False):

        """
        One row per HyPrColoc cluster: protein x tissue or cell type x cluster.

        HyPrColoc tests whether the pQTL, the outcome GWAS and an xQTL share one
        causal variant. A cluster supports the target only if it holds all three
        traits with posterior >= 0.7 (same rule as the dashboard). Betas are all on
        the same allele (effect allele column), so directions compare directly.

        protein: optional symbol (CTSH) or full protein id
        data_type: None for all, "bulk" or "single_cell"
        supported_only: True shows only supporting clusters
        """

        hyprcoloc_file = paths.hyprcoloc_out(self.pqtl_dataset, self.pheno_id, self.results_dir)
        if not hyprcoloc_file.exists():
            print("HyPrColoc was not run for this run.")
            return None

        hyprcoloc = pl.read_csv(hyprcoloc_file, separator="\t", infer_schema_length=None)

        records = []
        for row in hyprcoloc.iter_rows(named=True):
            protein_id = row["protein"]

            if protein is not None and not self._matches(protein_id, protein):
                continue
            if data_type is not None and row["data_type"] != data_type:
                continue

            # which kinds of trait sit in this cluster
            has_pqtl = False
            has_gwas = False
            has_xqtl = False
            labels = []
            for trait in str(row["traits"]).split(","):
                trait = trait.strip()
                lower = trait.lower()
                if lower.startswith("pqtl_"):
                    has_pqtl = True
                    labels.append("pQTL")
                elif lower.startswith("gwas_"):
                    has_gwas = True
                    labels.append("GWAS")
                elif "qtl_" in lower:
                    has_xqtl = True
                    labels.append(trait.split("_")[0])

            all_three = has_pqtl and has_gwas and has_xqtl
            posterior = row["posterior_prob"]
            supported = all_three and posterior is not None and posterior >= 0.7

            if supported_only and not supported:
                continue

            # xQTL and pQTL on the same allele: same sign means the xQTL trait and protein move together
            pqtl_beta = row["pqtl_beta"]
            qtl_beta = row["qtl_beta"]
            if row["qtl_type"] != "eqtl" or pqtl_beta is None or qtl_beta is None:
                direction = "not comparable"
            elif (pqtl_beta > 0) == (qtl_beta > 0):
                direction = "agrees"
            else:
                direction = "opposes"

            if row["data_type"] == "single_cell":
                kind = "single cell"
            else:
                kind = "bulk"

            records.append({
                "symbol": self._symbol(protein_id),
                "protein": protein_id,
                "data type": kind,
                "QTL dataset": row["qtl_dataset"],
                "tissue / cell type": row["cell_type"],
                "cluster": row["iteration"],
                "traits in cluster": " + ".join(labels),
                "posterior": row["posterior_prob"],
                "regional probability": row["regional_prob"],
                "candidate SNP": row["candidate_snp"],
                "posterior explained by SNP": row["posterior_explained_by_snp"],
                "effect allele": row["a1"],
                "GWAS beta": row["gwas_beta"],
                "pQTL beta": pqtl_beta,
                "xQTL beta": qtl_beta,
                "xQTL vs pQTL direction": direction,
                "supported": supported,
            })

        if len(records) == 0:
            print("No HyPrColoc clusters match.")
            return None

        table = pl.DataFrame(
            records,
            infer_schema_length=None,
            schema_overrides={
                "posterior": pl.Float64,
                "regional probability": pl.Float64,
                "posterior explained by SNP": pl.Float64,
                "GWAS beta": pl.Float64,
                "pQTL beta": pl.Float64,
                "xQTL beta": pl.Float64,
            },
        )

        for column in ["posterior", "regional probability", "posterior explained by SNP"]:
            table = table.with_columns(pl.col(column).round(3))
        for column in ["GWAS beta", "pQTL beta", "xQTL beta"]:
            table = table.with_columns(pl.col(column).round(4))

        # supported_only already means every row passed, so the flag adds nothing
        if supported_only:
            table = table.drop("supported")

        # group each protein's tissues together, strongest first
        table = table.sort(["symbol", "posterior"], descending=[False, True])
        table = self._finish(table, ["data type", "tissue / cell type", "traits in cluster", "posterior", "candidate SNP", "xQTL vs pQTL direction", "supported"], detailed)
        return table

    def phewas(self, protein=None, significant_only=True, detailed=False):

        """
        One row per protein x trait from the phenome wide MR (FinnGen, then UK Biobank).

        A Bonferroni significant trait with the same sign as the protein's cis-MR effect
        on the outcome is a possible additional indication (the same drug direction helps
        both); the opposite sign is a possible adverse effect (same rule as the dashboard).

        protein: optional symbol (TREM2) or full protein id
        significant_only: True shows only Bonferroni significant traits
        """

        sources = [
            ("FinnGen", paths.phewas_out(self.pqtl_dataset, self.pheno_id, self.results_dir)),
            ("UK Biobank", paths.phewas_ukbb_out(self.pqtl_dataset, self.pheno_id, self.results_dir)),
        ]
        mr_beta = self._mr_betas()

        records = []
        any_source = False
        for source, phewas_file in sources:
            phewas = self._read(phewas_file)
            if phewas is None:
                continue
            any_source = True

            for row in phewas.iter_rows(named=True):
                protein_id = row["protein"]

                if protein is not None and not self._matches(protein_id, protein):
                    continue

                significant = bool(row["bonferroni_significant"])
                if significant_only and not significant:
                    continue

                beta = row["beta_mr"]
                se = row["se_mr"]
                protein_beta = mr_beta.get(protein_id)

                if not significant:
                    classification = "not significant"
                elif protein_beta is None or beta is None:
                    classification = "not classified"
                elif (beta > 0) == (protein_beta > 0):
                    classification = "additional indication"
                else:
                    classification = "possible adverse effect"

                # code column differs between sources
                code = row.get("PHENOCODE")

                records.append({
                    "symbol": self._symbol(protein_id),
                    "protein": protein_id,
                    "source": source,
                    "trait": row["PHENOSTRING"],
                    "code": code,
                    "category": row["CATEGORY"],
                    "method": row["method"],
                    "n_instruments": row["n_instruments"],
                    "MR beta": None if beta is None else round(beta, 4),
                    "OR": None if beta is None else round(math.exp(beta), 3),
                    "OR 95% CI": None if beta is None or se is None else f"{math.exp(beta - 1.96 * se):.3f} to {math.exp(beta + 1.96 * se):.3f}",
                    "P": row["p_mr"],
                    "Bonferroni P": row["p_bonferroni"],
                    "cis-MR beta": None if protein_beta is None else round(protein_beta, 4),
                    "class": classification,
                })

        if not any_source:
            print("PheWAS was not run for this run.")
            return None
        if len(records) == 0:
            if significant_only and protein is not None:
                print(f"No Bonferroni significant PheWAS traits for {protein}.")
            else:
                print("No PheWAS results match.")
            return None

        table = pl.DataFrame(
            records,
            infer_schema_length=None,
            schema_overrides={"code": pl.String, "MR beta": pl.Float64, "OR": pl.Float64, "cis-MR beta": pl.Float64},
        )
        table = table.sort(["symbol", "P"])
        table = self._finish(table, ["source", "trait", "category", "OR", "P", "Bonferroni P", "class"], detailed)
        return table

    def instruments(self, protein=None, hits_only=True, detailed=False):

        """
        One row per instrument SNP used in cis-MR.

        Effect allele is the allele that raises the protein's measured level direction
        in the pQTL (exposure), and outcome betas are aligned to it. Wald ratio is that
        single SNP's estimate (outcome beta / pQTL beta). Steiger checks the SNP explains
        more variance in the protein than in the outcome (correct causal direction).

        protein: optional symbol (TREM2) or full protein id
        hits_only: True restricts to proteins passing cis-MR
        """

        instruments = self._read(paths.mr_instruments_out(self.pqtl_dataset, self.pheno_id, self.results_dir))
        if instruments is None:
            print("No instruments file for this run.")
            return None

        cis_hits = self._cis_mr_hits()

        records = []
        for row in instruments.iter_rows(named=True):
            protein_id = row["protein"]

            if not row["mr_keep"]:
                continue
            if hits_only and protein_id not in cis_hits:
                continue
            if protein is not None and not self._matches(protein_id, protein):
                continue

            wald_ratio = None
            if row["beta.exposure"] not in (None, 0):
                wald_ratio = row["beta.outcome"] / row["beta.exposure"]

            records.append({
                "symbol": self._symbol(protein_id),
                "protein": protein_id,
                "SNP": row["SNP"],
                "effect allele": row["effect_allele.exposure"],
                "other allele": row["other_allele.exposure"],
                "EAF": row["eaf.exposure"],
                "pQTL beta": row["beta.exposure"],
                "pQTL P": row["pval.exposure"],
                "F": row["F"],
                "outcome beta": row["beta.outcome"],
                "outcome P": row["pval.outcome"],
                "Wald ratio": wald_ratio,
                "Steiger correct direction": row["steiger_dir"],
                "Steiger P": row["steiger_pval"],
                "palindromic": row["palindromic"],
            })

        if len(records) == 0:
            print("No instruments match.")
            return None

        table = pl.DataFrame(records, infer_schema_length=None, schema_overrides={"EAF": pl.Float64, "Wald ratio": pl.Float64})
        for column in ["EAF", "pQTL beta", "outcome beta", "Wald ratio"]:
            table = table.with_columns(pl.col(column).round(4))
        table = table.with_columns(pl.col("F").round(1))

        table = table.sort(["symbol", "pQTL P"])
        table = self._finish(table, ["SNP", "effect allele", "EAF", "pQTL beta", "pQTL P", "F", "outcome beta", "outcome P", "Steiger correct direction"], detailed)
        return table

    def sensitivity(self, protein=None, hits_only=True, verbose=True, detailed=False):

        """
        One row per multi instrument protein: does the estimate hold across MR methods?

        IVW assumes all instruments are valid; MR Egger allows directional pleiotropy
        (intercept tests it); weighted median holds if at least half the weight is valid.
        Egger and weighted median need 3 or more instruments. Single instrument proteins
        (Wald ratio) have no sensitivity methods and are left out.

        protein: optional symbol (TREM2) or full protein id
        hits_only: True restricts to proteins passing cis-MR
        verbose: False hides the note about single instrument proteins
        """

        mr = self._read(paths.mr_out(self.pqtl_dataset, self.pheno_id, self.results_dir))
        cis_hits = self._cis_mr_hits()

        records = []
        n_single = 0
        for row in mr.iter_rows(named=True):
            protein_id = row["protein"]

            if hits_only and protein_id not in cis_hits:
                continue
            if protein is not None and not self._matches(protein_id, protein):
                continue

            n_instruments = row["n_instruments"]
            if n_instruments < 2:
                n_single = n_single + 1
                continue

            # sign of every method that could be run
            # IVW from 2 instruments, Egger and weighted median from 3
            signs = []
            if row["IVW_beta"] is not None:
                signs.append(row["IVW_beta"] > 0)
            if n_instruments >= 3:
                for column in ["Egger_beta", "WME_beta"]:
                    if row[column] is not None:
                        signs.append(row[column] > 0)

            if len(signs) <= 1:
                agreement = "IVW only"
            elif all(signs) or not any(signs):
                agreement = "all agree"
            else:
                agreement = "disagree"

            egger_beta = row["Egger_beta"] if n_instruments >= 3 else None
            egger_p = row["Egger_pval"] if n_instruments >= 3 else None
            wme_beta = row["WME_beta"] if n_instruments >= 3 else None
            wme_p = row["WME_pval"] if n_instruments >= 3 else None
            intercept_p = row["egger_intercept_pval"] if n_instruments >= 3 else None

            records.append({
                "symbol": self._symbol(protein_id),
                "protein": protein_id,
                "n_instruments": n_instruments,
                "IVW beta": row["IVW_beta"],
                "IVW P": row["IVW_pval"],
                "Egger beta": egger_beta,
                "Egger P": egger_p,
                "Weighted median beta": wme_beta,
                "Weighted median P": wme_p,
                "Egger intercept P": intercept_p,
                "Cochran's Q P": row["Q_pval"],
                "method signs": agreement,
            })

        if n_single > 0 and protein is None and verbose:
            print(f"{n_single} single instrument proteins (Wald ratio) have no sensitivity methods and are not shown.")
        if len(records) == 0:
            if verbose and protein is not None and n_single > 0:
                print(f"{protein} has a single instrument (Wald ratio), so there are no sensitivity methods.")
            elif verbose:
                print("No multi instrument proteins match.")
            return None

        numeric = ["IVW beta", "IVW P", "Egger beta", "Egger P", "Weighted median beta", "Weighted median P", "Egger intercept P", "Cochran's Q P"]
        overrides = {}
        for column in numeric:
            overrides[column] = pl.Float64
        table = pl.DataFrame(records, infer_schema_length=None, schema_overrides=overrides)
        for column in ["IVW beta", "Egger beta", "Weighted median beta"]:
            table = table.with_columns(pl.col(column).round(4))

        table = table.sort("IVW P")
        table = self._finish(table, ["n_instruments", "IVW beta", "IVW P", "Egger beta", "Weighted median beta", "Egger intercept P", "method signs"], detailed)
        return table

    def directions(self, protein=None, detailed=False):

        """
        One row per prioritised target: do all lines of evidence point the same way?

        cis-MR gives the protein's direction on the outcome. Each other check counts how
        many results agree with it: SMR (xQTL trait on outcome) split into bulk and
        single cell, HyPrColoc (xQTL and pQTL on the same allele), and the MR
        sensitivity methods. Overall is "consistent" when nothing opposes, "mixed" when
        something does.

        protein: optional symbol (TREM2) or full protein id
        """

        prioritised = self._prioritised()
        mr_beta = self._mr_betas()

        smr = self.smr(detailed=True)
        hyprcoloc = self.hyprcoloc(detailed=True)
        sensitivity = self.sensitivity(verbose=False, detailed=True)

        # count agree / total per protein from each table
        smr_bulk = self._count_directions(smr, "direction vs cis-MR", "data type", "bulk")
        smr_sc = self._count_directions(smr, "direction vs cis-MR", "data type", "single cell")
        hyprcoloc_counts = self._count_directions(hyprcoloc, "xQTL vs pQTL direction", None, None)

        sensitivity_call = {}
        if sensitivity is not None:
            for row in sensitivity.reset_index(names="symbol").to_dict("records"):
                sensitivity_call[row["protein"]] = row["method signs"]

        records = []
        for protein_id in sorted(prioritised):

            if protein is not None and not self._matches(protein_id, protein):
                continue

            beta = mr_beta.get(protein_id)
            if beta is not None and beta > 0:
                effect = "raises risk"
                drug_direction = "inhibit"
            else:
                effect = "lowers risk"
                drug_direction = "activate"

            any_opposes = False
            any_check = False
            cells = {}
            for label, counts in [("SMR bulk", smr_bulk), ("SMR single cell", smr_sc), ("HyPrColoc", hyprcoloc_counts)]:
                agree, total = counts.get(protein_id, (0, 0))
                if total == 0:
                    cells[label] = ""
                else:
                    cells[label] = f"{agree}/{total} agree"
                    any_check = True
                    if agree < total:
                        any_opposes = True

            method_signs = sensitivity_call.get(protein_id, "")
            if method_signs == "disagree":
                any_opposes = True
            if method_signs in ("all agree", "disagree"):
                any_check = True

            if not any_check:
                overall = "cis-MR only"
            elif any_opposes:
                overall = "mixed"
            else:
                overall = "consistent"

            records.append({
                "symbol": self._symbol(protein_id),
                "protein": protein_id,
                "higher protein": effect,
                "drug direction": drug_direction,
                "MR methods": method_signs,
                "SMR bulk": cells["SMR bulk"],
                "SMR single cell": cells["SMR single cell"],
                "HyPrColoc": cells["HyPrColoc"],
                "overall": overall,
            })

        if len(records) == 0:
            print("No prioritised targets match.")
            return None

        table = pl.DataFrame(records, infer_schema_length=None)
        table = table.sort("symbol")
        table = self._finish(table, ["higher protein", "drug direction", "MR methods", "SMR bulk", "SMR single cell", "HyPrColoc", "overall"], detailed)
        return table

    def targets(self, detailed=False):

        """
        One row per prioritised target, every stage summarised: the evidence table.
        """

        prioritised = self._prioritised()
        if len(prioritised) == 0:
            print("No prioritised targets in this run.")
            return None

        cis_mr = self.cis_mr(detailed=True).reset_index(names="symbol")
        coloc = self.colocalisation(detailed=True).reset_index(names="symbol")
        smr = self.smr(detailed=True)
        hyprcoloc = self.hyprcoloc(detailed=True)
        phewas = self.phewas(detailed=True)
        directions = self.directions(detailed=True).reset_index(names="symbol")

        mr_by_protein = {}
        for row in cis_mr.to_dict("records"):
            mr_by_protein[row["protein"]] = row
        coloc_by_protein = {}
        for row in coloc.to_dict("records"):
            coloc_by_protein[row["protein"]] = row
        direction_by_protein = {}
        for row in directions.to_dict("records"):
            direction_by_protein[row["protein"]] = row["overall"]

        smr_bulk = self._count_rows(smr, "data type", "bulk")
        smr_sc = self._count_rows(smr, "data type", "single cell")
        hyprcoloc_n = self._count_rows(hyprcoloc, None, None)
        indications = self._count_rows(phewas, "class", "additional indication")
        adverse = self._count_rows(phewas, "class", "possible adverse effect")

        records = []
        for protein_id in prioritised:
            mr_row = mr_by_protein[protein_id]
            coloc_row = coloc_by_protein[protein_id]

            best_pp4 = coloc_row["COLOC PP4"]
            pwcoco_pp4 = coloc_row["PWCoCo best PP4"]
            if best_pp4 is None or (pwcoco_pp4 is not None and pwcoco_pp4 > best_pp4):
                best_pp4 = pwcoco_pp4

            records.append({
                "symbol": self._symbol(protein_id),
                "protein": protein_id,
                "OR": mr_row["OR"],
                "OR 95% CI": mr_row["OR 95% CI"],
                "P": mr_row["P"],
                "FDR q": mr_row["FDR q"],
                "drug direction": mr_row["drug direction"],
                "colocalisation": coloc_row["supported by"],
                "best PP4": best_pp4,
                "SMR bulk tissues": smr_bulk.get(protein_id, 0),
                "SMR cell types": smr_sc.get(protein_id, 0),
                "HyPrColoc clusters": hyprcoloc_n.get(protein_id, 0),
                "directions": direction_by_protein.get(protein_id, ""),
                "PheWAS indications": indications.get(protein_id, 0),
                "PheWAS adverse": adverse.get(protein_id, 0),
            })

        table = pl.DataFrame(records, infer_schema_length=None)
        table = table.sort("P")
        table = self._finish(table, ["OR", "P", "drug direction", "colocalisation", "best PP4", "SMR bulk tissues", "SMR cell types", "HyPrColoc clusters", "directions", "PheWAS adverse"], detailed)
        return table

    def target(self, protein):

        """
        Every table above for one protein, in reading order. Displays each section
        and returns them as a dictionary {section: table}.
        """

        titles = [
            "cis-MR",
            "Instruments",
            "MR sensitivity",
            "Lead variant",
            "Colocalisation",
            "COLOC prior sensitivity",
            "SMR",
            "HyPrColoc",
            "PWCoCo xQTL",
            "Directions",
            "PheWAS",
        ]

        results = {}
        for title in titles:
            # heading first, so any message from the method lands under it
            self._heading(title)

            if title == "cis-MR":
                table = self.cis_mr(protein, hits_only=False)
            elif title == "Instruments":
                table = self.instruments(protein, hits_only=False)
            elif title == "MR sensitivity":
                table = self.sensitivity(protein, hits_only=False)
            elif title == "Lead variant":
                table = self.lead_variants(protein)
            elif title == "Colocalisation":
                table = self.colocalisation(protein, hits_only=False)
            elif title == "COLOC prior sensitivity":
                table = self.coloc_priors(protein)
            elif title == "SMR":
                table = self.smr(protein)
            elif title == "HyPrColoc":
                table = self.hyprcoloc(protein)
            elif title == "PWCoCo xQTL":
                table = self.pwcoco_qtl(protein)
            elif title == "Directions":
                table = self.directions(protein)
            else:
                table = self.phewas(protein)

            results[title] = table
            if table is not None:
                self._show(table)
        return results

    def coloc_priors(self, protein=None, detailed=False):

        """
        One row per protein: COLOC PP4 under each prior set the pipeline tried.

        COLOC's answer depends on p12, the prior that one variant affects both traits.
        If PP4 stays above the threshold under every prior, the colocalisation is robust.

        protein: optional symbol (TREM2) or full protein id
        """

        sensitivity = self._read(paths.coloc_sensitivity_out(self.pqtl_dataset, self.pheno_id, self.results_dir))
        if sensitivity is None:
            print("No COLOC prior sensitivity file for this run.")
            return None

        pp4_threshold = self.params["gates"]["coloc"]["pp4_threshold"]

        # keep the order the pipeline wrote the scenarios in
        scenarios = []
        for scenario in sensitivity["scenario"].to_list():
            if scenario not in scenarios:
                scenarios.append(scenario)

        by_protein = {}
        for row in sensitivity.iter_rows(named=True):
            protein_id = row["protein_id"]
            if protein is not None and not self._matches(protein_id, protein):
                continue
            if protein_id not in by_protein:
                by_protein[protein_id] = {}
            by_protein[protein_id][row["scenario"]] = row["PP.H4.abf"]

        if len(by_protein) == 0:
            print("No COLOC prior results match.")
            return None

        records = []
        for protein_id in by_protein:
            record = {"symbol": self._symbol(protein_id), "protein": protein_id}
            robust = True
            for scenario in scenarios:
                pp4 = by_protein[protein_id].get(scenario)
                record[f"PP4 {scenario.replace('_', ' ')}"] = None if pp4 is None else round(pp4, 3)
                if pp4 is None or pp4 < pp4_threshold:
                    robust = False
            record["passes under all priors"] = robust
            records.append(record)

        table = pl.DataFrame(records, infer_schema_length=None)
        table = table.sort("symbol")
        # compact view: every PP4 scenario, without the full protein id
        table = self._finish(table, None, detailed)
        return table

    def lead_variants(self, protein=None, detailed=False):

        """
        One row per protein: the top cis variant with outcome and pQTL effects on the same allele.

        GWAS beta / pQTL beta at this SNP gives the implied direction of the protein on the
        outcome, which should match the cis-MR estimate.

        protein: optional symbol (TREM2) or full protein id
        """

        lead = self._read(paths.target_stats_out(self.pqtl_dataset, self.pheno_id, self.results_dir))
        if lead is None:
            print("No lead variant file for this run.")
            return None

        mr_beta = self._mr_betas()

        records = []
        for row in lead.iter_rows(named=True):
            protein_id = row["protein"]
            if protein is not None and not self._matches(protein_id, protein):
                continue

            gwas_beta = row["GWAS_BETA"]
            pqtl_beta = row["pQTL_BETA"]
            protein_beta = mr_beta.get(protein_id)

            if gwas_beta is None or pqtl_beta is None or pqtl_beta == 0 or protein_beta is None:
                direction = "not comparable"
            elif (gwas_beta / pqtl_beta > 0) == (protein_beta > 0):
                direction = "agrees"
            else:
                direction = "opposes"

            records.append({
                "symbol": self._symbol(protein_id),
                "protein": protein_id,
                "SNP": row["SNP"],
                "effect allele": row["A1"],
                "other allele": row["A2"],
                "EAF": row["FRQ"],
                "GWAS beta": gwas_beta,
                "GWAS P": row["GWAS_P"],
                "pQTL beta": pqtl_beta,
                "pQTL P": row["pQTL_P"],
                "implied direction vs cis-MR": direction,
            })

        if len(records) == 0:
            print("No lead variants match.")
            return None

        table = pl.DataFrame(records, infer_schema_length=None, schema_overrides={"GWAS beta": pl.Float64, "pQTL beta": pl.Float64})
        for column in ["GWAS beta", "pQTL beta"]:
            table = table.with_columns(pl.col(column).round(4))
        table = table.sort("symbol")
        table = self._finish(table, ["SNP", "effect allele", "GWAS beta", "GWAS P", "pQTL beta", "pQTL P", "implied direction vs cis-MR"], detailed)
        return table

    def pwcoco_qtl(self, protein=None, comparison=None, detailed=False):

        """
        One row per protein x xQTL dataset x tissue or cell type x comparison.

        PWCoCo between the xQTL and (a) the pQTL and (b) the outcome GWAS.
        Both colocalising means the xQTL, the protein and the outcome share a signal.
        Best = highest PP4 over unconditioned and conditioned pairs.

        protein: optional symbol (TREM2) or full protein id
        comparison: None for both, "eQTL vs pQTL" or "eQTL vs GWAS"
        """

        files = [
            ("eQTL vs pQTL", paths.pwcoco_eqtl_pqtl_out(self.pqtl_dataset, self.pheno_id, self.results_dir)),
            ("eQTL vs GWAS", paths.pwcoco_eqtl_gwas_out(self.pqtl_dataset, self.pheno_id, self.results_dir)),
        ]
        pp4_threshold = self.params["gates"]["coloc"]["pp4_threshold"]

        # best and unconditioned PP4 per (comparison, protein, dataset, cell type)
        best = {}
        unconditioned = {}
        any_file = False
        for label, pwcoco_file in files:
            if comparison is not None and label != comparison:
                continue
            pwcoco = self._read(pwcoco_file)
            if pwcoco is None:
                continue
            any_file = True

            for row in pwcoco.iter_rows(named=True):
                protein_id = row["protein"]
                if protein is not None and not self._matches(protein_id, protein):
                    continue

                key = (label, protein_id, row["qtl_dataset"], row["cell_type"])
                if row["SNP1"] == "unconditioned" and row["SNP2"] == "unconditioned":
                    unconditioned[key] = row["H4"]
                current = best.get(key)
                if current is None or row["H4"] > current["H4"]:
                    best[key] = row

        if not any_file:
            print("PWCoCo xQTL was not run for this run.")
            return None
        if len(best) == 0:
            print("No PWCoCo xQTL results match.")
            return None

        records = []
        for key in best:
            label, protein_id, qtl_dataset, cell_type = key
            row = best[key]
            records.append({
                "symbol": self._symbol(protein_id),
                "protein": protein_id,
                "comparison": label,
                "QTL dataset": qtl_dataset,
                "tissue / cell type": cell_type,
                "unconditioned PP4": unconditioned.get(key),
                "best PP4": row["H4"],
                "best PP3": row["H3"],
                "best pair": f"{row['SNP1']} / {row['SNP2']}",
                "colocalised": row["H4"] >= pp4_threshold,
            })

        table = pl.DataFrame(records, infer_schema_length=None, schema_overrides={"unconditioned PP4": pl.Float64})
        for column in ["unconditioned PP4", "best PP4", "best PP3"]:
            table = table.with_columns(pl.col(column).round(3))
        table = table.sort(["symbol", "comparison", "best PP4"], descending=[False, False, True])
        table = self._finish(table, ["comparison", "QTL dataset", "tissue / cell type", "unconditioned PP4", "best PP4", "colocalised"], detailed)
        return table

    def regional_plot(self, protein):

        """
        LocusZoom style plot of the cis window: outcome GWAS on top, pQTL below,
        each SNP coloured by LD (r2) with the lead variant. Returns the figure.
        """

        import matplotlib.pyplot as plt

        protein_id = self._find_protein(protein)
        if protein_id is None:
            print(f"No protein matches {protein}.")
            return None

        locus_dir = self.results_dir / "locus_data" / protein_id
        if not (locus_dir / "gwas.parquet").exists():
            print(f"No locus data for {protein_id}; locus plots are only exported for prioritised targets.")
            return None

        gwas = pl.read_parquet(locus_dir / "gwas.parquet")
        pqtl = pl.read_parquet(locus_dir / "pqtl.parquet")
        with open(locus_dir / "metadata.json") as fh:
            metadata = json.load(fh)
        lead_snp = metadata.get("candidate_snp")
        chromosome = metadata.get("chromosome")
        symbol = self._symbol(protein_id)

        ld = None
        if (locus_dir / "ld.parquet").exists():
            ld = pl.read_parquet(locus_dir / "ld.parquet")

        # LocusZoom r2 bins: (upper bound, colour, legend label); weakest drawn first
        ld_bins = [
            (0.2, "#2c4a8a", "< 0.2"),
            (0.4, "#5fb4e0", "0.2 to 0.4"),
            (0.6, "#4dac26", "0.4 to 0.6"),
            (0.8, "#f1a340", "0.6 to 0.8"),
            (1.01, "#d7301f", "0.8 to 1"),
        ]
        genome_wide = -math.log10(5e-8)

        fig, axes = plt.subplots(2, 1, figsize=(9, 6.2), dpi=150, sharex=True)
        panels = [
            (axes[0], gwas, f"{self.pheno_id} GWAS"),
            (axes[1], pqtl, f"{symbol} pQTL ({self.pqtl_dataset})"),
        ]

        for axis, data, title in panels:
            data = data.filter(pl.col("P") > 0).with_columns((-pl.col("P").log10()).alias("log10p"))

            if ld is None:
                axis.scatter((data["BP"] / 1e6).to_list(), data["log10p"].to_list(),
                             color="#7f7f7f", s=9, alpha=0.7, linewidths=0)
            else:
                data = data.join(ld, on="SNP", how="left").with_columns(pl.col("R2").fill_null(0))
                lower = -1
                for upper, colour, label in ld_bins:
                    in_bin = data.filter((pl.col("R2") > lower) & (pl.col("R2") <= upper))
                    # the many unlinked SNPs are drawn smaller and fainter
                    if upper <= 0.2:
                        size = 7
                        alpha = 0.45
                    else:
                        size = 14
                        alpha = 0.9
                    axis.scatter((in_bin["BP"] / 1e6).to_list(), in_bin["log10p"].to_list(),
                                 color=colour, s=size, alpha=alpha, linewidths=0, label=label)
                    lower = upper

            axis.axhline(genome_wide, color="#9e9e9e", linestyle="--", linewidth=0.8)

            lead = data.filter(pl.col("SNP") == lead_snp)
            if lead.height > 0:
                x = lead["BP"][0] / 1e6
                y = lead["log10p"][0]
                axis.scatter([x], [y], marker="D", s=55, color="#7b3294", edgecolor="black",
                             linewidths=0.6, zorder=3, label="Lead variant")
                axis.annotate(lead_snp, (x, y), xytext=(6, 4), textcoords="offset points", fontsize=8)
                axis.set_ylim(0, y * 1.15)

            axis.set_ylabel(r"$-\log_{10}$ P")
            axis.text(0.01, 0.95, title, transform=axis.transAxes, fontsize=10, fontweight="bold", va="top")
            self._style_axis(axis)

        if ld is not None:
            legend = axes[0].legend(loc="upper right", frameon=False, fontsize=8, markerscale=1.3,
                                    title=f"r$^2$ with {lead_snp}", title_fontsize=8, handletextpad=0.3)
            legend._legend_box.align = "left"

        axes[1].set_xlabel(f"Position on chromosome {chromosome} (Mb)")
        fig.suptitle(f"{symbol}: cis window", x=0.07, ha="left", fontsize=12, fontweight="bold")
        fig.text(0.07, 0.925, f"Dashed line: P = 5e-8. Lead variant from {metadata.get('candidate_source', 'the pipeline')}.",
                 fontsize=8, color="#555555")
        fig.subplots_adjust(left=0.07, right=0.98, top=0.9, bottom=0.09, hspace=0.12)
        plt.close(fig)
        return fig

    def forest_plot(self, proteins=None):

        """
        Odds ratio with 95% CI for prioritised targets (or a chosen list), strongest
        at the top, coloured by whether higher protein raises or lowers risk.
        Returns the figure.
        """

        import matplotlib.pyplot as plt
        from matplotlib.lines import Line2D
        from matplotlib.ticker import NullLocator

        table = self.cis_mr(detailed=True).reset_index(names="symbol")
        prioritised = self._prioritised()

        rows = []
        for row in table.to_dict("records"):
            if proteins is None:
                if row["protein"] not in prioritised:
                    continue
            else:
                keep = False
                for query in proteins:
                    if self._matches(row["protein"], query):
                        keep = True
                if not keep:
                    continue
            rows.append(row)

        if len(rows) == 0:
            print("No proteins to plot.")
            return None

        # strongest evidence at the top of the plot
        rows.reverse()

        # symbols shared by several aptamers or probes get the full protein id instead
        symbol_counts = {}
        for row in rows:
            symbol_counts[row["symbol"]] = symbol_counts.get(row["symbol"], 0) + 1

        raises_colour = "#b2182b"
        lowers_colour = "#2166ac"

        labels = []
        odds = []
        lows = []
        highs = []
        colours = []
        estimates = []
        for row in rows:
            if symbol_counts[row["symbol"]] > 1:
                # e.g. PILRA_Q9UKJ1_10816-150 -> PILRA (10816-150)
                labels.append(f"{row['symbol']} ({row['protein'].split('_')[-1]})")
            else:
                labels.append(row["symbol"])
            odds_ratio = math.exp(row["beta"])
            low = math.exp(row["beta"] - 1.96 * row["se"])
            high = math.exp(row["beta"] + 1.96 * row["se"])
            odds.append(odds_ratio)
            lows.append(odds_ratio - low)
            highs.append(high - odds_ratio)
            colours.append(raises_colour if odds_ratio > 1 else lowers_colour)
            estimates.append(f"{odds_ratio:.2f} ({low:.2f} to {high:.2f})")

        longest_label = 0
        for label in labels:
            longest_label = max(longest_label, len(label))
        left_inches = 0.3 + 0.075 * longest_label
        width = left_inches + 4.6 + 1.9
        height = 1.5 + 0.3 * len(rows)
        fig, axis = plt.subplots(figsize=(width, height), dpi=150)
        positions = list(range(len(rows)))

        axis.axvline(1, color="#9e9e9e", linestyle="--", linewidth=0.8, zorder=1)
        for index in positions:
            axis.errorbar(odds[index], index, xerr=[[lows[index]], [highs[index]]], fmt="s", markersize=5,
                          color=colours[index], ecolor=colours[index], elinewidth=1.2, capsize=0, zorder=2)

        axis.set_yticks(positions)
        axis.set_yticklabels(labels, fontsize=9)
        axis.set_ylim(-0.7, len(rows) - 0.3)

        # log scale with plain, readable tick labels
        axis.set_xscale("log")
        smallest = min(odd - low for odd, low in zip(odds, lows))
        largest = max(odd + high for odd, high in zip(odds, highs))
        candidates = [0.05, 0.1, 0.2, 0.3, 0.5, 0.67, 0.8, 0.9, 1.1, 1.25, 1.5, 2, 3, 5, 10, 20, 50]
        span = max(largest, 1) / min(smallest, 1)
        min_step = max(1.08, span ** (1 / 7))
        ticks = [1]
        last = 1
        for tick in candidates:
            if 1 < tick <= largest * 1.1 and tick / last >= min_step:
                ticks.append(tick)
                last = tick
        last = 1
        for tick in reversed(candidates):
            if smallest * 0.9 <= tick < 1 and last / tick >= min_step:
                ticks.append(tick)
                last = tick
        ticks.sort()
        axis.set_xticks(ticks)
        axis.set_xticklabels([f"{tick:g}" for tick in ticks], fontsize=8)
        axis.xaxis.set_minor_locator(NullLocator())
        axis.set_xlim(min(smallest, 1) * 0.92, max(largest, 1) * 1.08)
        axis.set_xlabel(f"Odds ratio for {self.pheno_id} per unit higher protein (95% CI)", fontsize=9)

        # estimates printed to the right of the plot
        for index in positions:
            axis.text(1.02, index, estimates[index], transform=axis.get_yaxis_transform(),
                      va="center", fontsize=8, family="monospace")
        axis.text(1.02, len(rows) - 0.3, "OR (95% CI)", transform=axis.get_yaxis_transform(),
                  va="bottom", fontsize=8, fontweight="bold")

        self._style_axis(axis)
        axis.spines["left"].set_visible(False)
        axis.tick_params(axis="y", length=0)

        # legend lists only the colours used
        legend_items = []
        if raises_colour in colours:
            legend_items.append(Line2D([0], [0], marker="s", color=raises_colour, linestyle="", label="Higher protein raises risk (inhibit)"))
        if lowers_colour in colours:
            legend_items.append(Line2D([0], [0], marker="s", color=lowers_colour, linestyle="", label="Higher protein lowers risk (activate)"))
        fig.legend(handles=legend_items, loc="lower left", bbox_to_anchor=(left_inches / width, 0.02),
                   ncol=2, frameon=False, fontsize=8, handletextpad=0.3, columnspacing=1.5)

        axis.set_title(f"cis-MR estimates, prioritised targets ({len(rows)})", loc="left", fontsize=11, fontweight="bold")
        fig.subplots_adjust(left=left_inches / width, right=(left_inches + 4.6) / width,
                            top=1 - 0.55 / height, bottom=1.0 / height)
        plt.close(fig)
        return fig

    @staticmethod
    def _style_axis(axis):

        """ Clean axis: no top or right frame, light grid behind the data """

        axis.spines["top"].set_visible(False)
        axis.spines["right"].set_visible(False)
        axis.grid(axis="x", color="#e6e6e6", linewidth=0.6)
        axis.set_axisbelow(True)

    def _read(self, path):

        """ Read a results table, or None when that stage did not run """

        if not path.exists():
            return None
        return pl.read_csv(path, separator="\t", infer_schema_length=None)

    def _mr_betas(self):

        """ protein id -> cis-MR beta (primary method) """

        mr = self._read(paths.mr_out(self.pqtl_dataset, self.pheno_id, self.results_dir))
        betas = {}
        for row in mr.iter_rows(named=True):
            betas[row["protein"]] = row["primary_beta"]
        return betas

    def _prioritised(self):

        """ cis-MR hits that also colocalise (COLOC or PWCoCo): the run's final target set """

        pp4 = self.params["gates"]["coloc"]["pp4_threshold"]
        coloc_file = paths.coloc_out(self.pqtl_dataset, self.pheno_id, self.results_dir)
        pwcoco_file = paths.pwcoco_out(self.pqtl_dataset, self.pheno_id, self.results_dir)
        coloc_hits = set(extract_coloc_or_pwcoco_targets(coloc_file, pwcoco_file, pp4, ("coloc", "pwcoco")))
        return self._cis_mr_hits() & coloc_hits

    def _find_protein(self, query):

        """ First full protein id among cis-MR results matching a symbol or id """

        for protein_id in self._mr_betas():
            if self._matches(protein_id, query):
                return protein_id
        return None

    @staticmethod
    def _count_rows(table, column, value):

        """ protein id -> number of rows (optionally where column == value) """

        counts = {}
        if table is None:
            return counts
        for row in table.reset_index(names="symbol").to_dict("records"):
            if column is not None and row[column] != value:
                continue
            counts[row["protein"]] = counts.get(row["protein"], 0) + 1
        return counts

    @staticmethod
    def _count_directions(table, direction_column, filter_column, filter_value):

        """ protein id -> (rows that agree, comparable rows) """

        counts = {}
        if table is None:
            return counts
        for row in table.reset_index(names="symbol").to_dict("records"):
            if filter_column is not None and row[filter_column] != filter_value:
                continue
            if row[direction_column] == "not comparable":
                continue
            agree, total = counts.get(row["protein"], (0, 0))
            if row[direction_column] == "agrees":
                agree = agree + 1
            counts[row["protein"]] = (agree, total + 1)
        return counts

    @staticmethod
    def _finish(table, compact_columns, detailed):

        """
        polars table -> pandas indexed by symbol, ready to display.

        detailed=False keeps only the key columns (compact_columns, or every column
        except the full protein id when None) and previews long tables;
        detailed=True keeps every column and shows every row. The full protein id
        is kept whenever one symbol has several aptamers or probes.
        """

        frame = table.to_pandas().set_index("symbol")

        if not detailed:
            keep = []

            # symbols measured by more than one aptamer or probe need the full id
            if "protein" in frame.columns:
                ids_per_symbol = frame.groupby(level=0)["protein"].nunique()
                if (ids_per_symbol > 1).any():
                    keep.append("protein")

            if compact_columns is None:
                for column in frame.columns:
                    if column != "protein":
                        keep.append(column)
            else:
                for column in compact_columns:
                    if column in frame.columns:
                        keep.append(column)

            frame = frame[keep]

        # no separate header row for the index name
        frame = ResultsTable(frame)
        frame.index.name = None
        frame.show_all = detailed
        return frame

    @staticmethod
    def _heading(title):

        """ A section heading: Markdown in Jupyter, plain text elsewhere """

        try:
            from IPython import get_ipython
            from IPython.display import Markdown, display
        except ImportError:
            print(f"\n=== {title} ===")
            return

        if get_ipython() is None:
            print(f"\n=== {title} ===")
        else:
            display(Markdown(f"#### {title}"))

    @staticmethod
    def _show(table):

        """ Display a table nicely in Jupyter, plain text elsewhere """

        try:
            from IPython.display import display
            display(table)
        except ImportError:
            print(table.to_string())

    def _cis_mr_hits(self):

        """ Proteins passing the run's cis-MR gate, exactly as the pipeline applies it """

        cis_gates = self.params["gates"]["cis_mr"]
        hits = grab_cis_mr_hits(
            csv_file=paths.mr_out(self.pqtl_dataset, self.pheno_id, self.results_dir),
            wald_fdr_q=cis_gates["wald_fdr_q"],
            ivw_fdr_q=cis_gates["ivw_fdr_q"],
            cochran_q_pval=cis_gates["cochran_q_pval"],
            egger_intercept_pval_min=cis_gates["egger_intercept_pval_min"],
            min_instruments_for_ivw=cis_gates["min_instruments_for_ivw"],
        )
        return set(hits)

    @staticmethod
    def _symbol(protein_id):

        """ CTSH_ENSG00000103811 or CTSH_P09668_8465-52 -> CTSH """

        return str(protein_id).split("_")[0]

    def _matches(self, protein_id, query):

        """ True if query is the full protein id or its symbol (case insensitive) """

        if str(protein_id) == str(query):
            return True
        return self._symbol(protein_id).upper() == str(query).upper()
