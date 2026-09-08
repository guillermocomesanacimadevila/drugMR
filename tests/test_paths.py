import os
from pathlib import Path

import pytest

from drugmr import paths

PQTL = "ukb_ppp"
PHENO = "AD"
BULK = "GTEx_v10"
SC = "SingleBrain"


def test_qc_out():
    assert paths.qc_out(PHENO) == Path("dat/derived/AD/qc_gwas/AD.tsv")


def test_mr_out():
    assert paths.mr_out(PQTL, PHENO) == Path("results/cis_mr/mr.tsv")


def test_mr_instruments_out():
    assert paths.mr_instruments_out(PQTL, PHENO) == Path(
        "results/cis_mr/instruments/mr_instruments.tsv"
    )


def test_coloc_out():
    assert paths.coloc_out(PQTL, PHENO) == Path("results/coloc/coloc.tsv")


def test_coloc_sensitivity_out():
    assert paths.coloc_sensitivity_out(PQTL, PHENO) == Path("results/coloc/coloc_sensitivity.tsv")


def test_coloc_susie_out():
    assert paths.coloc_susie_out(PQTL, PHENO) == Path("results/coloc/coloc_susie.tsv")


def test_target_stats_out():
    assert paths.target_stats_out(PQTL, PHENO) == Path("results/target_stats/top_cis_hits.tsv")


def test_smr_bulk_out():
    assert paths.smr_bulk_out(PQTL, PHENO, BULK) == Path(
        "results/smr/bulk/GTEx_v10/promising_targets.tsv"
    )


def test_smr_sc_out():
    assert paths.smr_sc_out(PQTL, PHENO, SC) == Path(
        "results/smr/sc/SingleBrain/promising_targets.tsv"
    )


def test_smr_bulk_dir():
    # UNCHANGED - synthesis/ tier, dataset-independent (see paths.py module docstring)
    assert paths.smr_bulk_dir(BULK) == Path("results/SMR/bulk/GTEx_v10")


def test_smr_final_targets_out():
    assert paths.smr_final_targets_out(PQTL, PHENO) == Path(
        "results/smr/final_multi_omics_targets.tsv"
    )


def test_hyprcoloc_out():
    assert paths.hyprcoloc_out(PQTL, PHENO) == Path("results/hyprcoloc/hyprcoloc.tsv")


def test_hyprcoloc_dataset_out():
    assert paths.hyprcoloc_dataset_out(PQTL, BULK, PHENO) == Path(
        "results/hyprcoloc/by_eqtl_source/GTEx_v10/hyprcoloc.tsv"
    )


def test_phewas_out():
    assert paths.phewas_out(PQTL, PHENO) == Path("results/phewas/finngen/phewas.tsv")


def test_phewas_ukbb_out():
    assert paths.phewas_ukbb_out(PQTL, PHENO) == Path("results/phewas/ukbb/phewas.tsv")


def test_phewas_finngen_coverage_out():
    assert paths.phewas_finngen_coverage_out(PQTL, PHENO) == Path(
        "results/phewas/finngen/phewas_coverage.tsv"
    )


def test_pwcoco_raw_prefix():
    assert paths.pwcoco_raw_prefix(PQTL, PHENO, "ANKRD54_Q6NXT1") == Path(
        "results/pwcoco/cis_pqtl/ANKRD54_Q6NXT1/ANKRD54_Q6NXT1"
    )


def test_pwcoco_out():
    assert paths.pwcoco_out(PQTL, PHENO) == Path("results/pwcoco/summary/pwcoco.tsv")


def test_pwcoco_qtl_raw_prefix():
    assert paths.pwcoco_qtl_raw_prefix("eqtl_gwas", PQTL, "ANKRD54_Q6NXT1", "GTEx_v10_Cerebellum") == Path(
        "results/pwcoco/eqtl_gwas/ANKRD54_Q6NXT1_GTEx_v10_Cerebellum/ANKRD54_Q6NXT1_GTEx_v10_Cerebellum"
    )


def test_pwcoco_eqtl_pqtl_out():
    assert paths.pwcoco_eqtl_pqtl_out(PQTL, PHENO) == Path("results/pwcoco/summary/pwcoco_eqtl_pqtl.tsv")


def test_pwcoco_eqtl_gwas_out():
    assert paths.pwcoco_eqtl_gwas_out(PQTL, PHENO) == Path("results/pwcoco/summary/pwcoco_eqtl_gwas.tsv")


def test_pwcoco_qtl_shared_out():
    assert paths.pwcoco_qtl_shared_out(PQTL, PHENO) == Path("results/pwcoco/summary/pwcoco_shared_snps.tsv")


def test_smr_raw_prefix_matches_smr_binary_out_prefix():
    # UNCHANGED - synthesis/ tier, dataset-independent (see paths.py module docstring).
    # drugmr.smr.SMR()'s --out prefix, and the .smr file bin/sort_smr.py reads back
    # (it appends ".smr" itself - that's the smr binary's own naming convention)
    eqtl_dataset = "bulk_raw/GTEx_v10/eQTL_GTEx_Brain_Cortex_v10/chr1"
    assert paths.smr_raw_prefix(eqtl_dataset, PHENO) == Path(
        "results/SMR/bulk_raw/GTEx_v10/eQTL_GTEx_Brain_Cortex_v10/chr1/AD/AD_chr1"
    )


def test_make_run_id():
    assert paths.make_run_id(PHENO, PQTL, "20260811", "a3a2aa1") == "AD_ukb_ppp_20260811_a3a2aa1"


def test_run_dirs():
    rid = "AD_ukb_ppp_20260811_a3a2aa1"
    assert paths.run_dir(rid) == Path("runs/AD_ukb_ppp_20260811_a3a2aa1")
    assert paths.run_results_dir(rid) == Path("runs/AD_ukb_ppp_20260811_a3a2aa1/results")
    assert paths.run_work_dir(rid) == Path("runs/AD_ukb_ppp_20260811_a3a2aa1/work")
    assert paths.run_logs_dir(rid) == Path("runs/AD_ukb_ppp_20260811_a3a2aa1/logs")
    assert paths.run_manifest_path(rid) == Path("runs/AD_ukb_ppp_20260811_a3a2aa1/manifest.json")
    assert paths.run_params_lock_path(rid) == Path("runs/AD_ukb_ppp_20260811_a3a2aa1/params.lock.yaml")


def test_registry_path():
    assert paths.registry_path() == Path("runs/registry.json")


def test_work_dir_for_results_dir():
    rid = "AD_ukb_ppp_20260811_a3a2aa1"
    local_results_dir = str(paths.run_results_dir(rid))
    assert paths.work_dir_for_results_dir(local_results_dir) == paths.run_work_dir(rid)
    assert paths.work_dir_for_results_dir(local_results_dir) == Path(
        "runs/AD_ukb_ppp_20260811_a3a2aa1/work"
    )


def test_synthesis_paths():
    assert paths.synthesis_dir(PHENO) == Path("synthesis/AD")
    assert paths.synthesis_target_stats_out(PHENO) == Path(
        "synthesis/AD/target_stats/all_datasets_mined_targets.tsv"
    )
    assert paths.synthesis_manifest_path(PHENO) == Path("synthesis/AD/manifest.json")


def test_tier2_out_dir_composes_with_run_results_dir():
    # this is the actual Phase 3 wiring: out_dir passed to every Tier-2 function
    # is now run_results_dir(run_id), not a bare "results" default
    rid = "AD_ukb_ppp_20260811_a3a2aa1"
    run_root = str(paths.run_results_dir(rid))
    assert paths.mr_out(PQTL, PHENO, run_root) == Path(
        "runs/AD_ukb_ppp_20260811_a3a2aa1/results/cis_mr/mr.tsv"
    )


# real run_id from a completed local pipeline run, migrated to the Phase 4 schema
# by analysis/migrate_results_schema.py - see that script if this run is ever
# cleaned up/archived
REAL_RUN_ID = "AD_ukb_ppp_20260811_149fc55"


@pytest.mark.skipif(
    os.environ.get("GITHUB_ACTIONS") == "true",
    reason="requires runs/ output from a real local pipeline run, not present on a fresh CI checkout",
)
def test_paths_match_real_files_on_disk():
    """Sanity check against actual files already on disk from a real prior run -
    not just the literals, but that paths.py resolves to files that exist."""
    repo_root = Path(__file__).resolve().parents[1]
    out_dir = str(paths.run_results_dir(REAL_RUN_ID))
    for p in [
        paths.mr_out(PQTL, PHENO, out_dir),
        paths.coloc_out(PQTL, PHENO, out_dir),
        paths.hyprcoloc_out(PQTL, PHENO, out_dir),
        paths.smr_final_targets_out(PQTL, PHENO, out_dir),
        paths.phewas_out(PQTL, PHENO, out_dir),
        paths.phewas_ukbb_out(PQTL, PHENO, out_dir),
        paths.target_stats_out(PQTL, PHENO, out_dir),
        paths.pwcoco_out(PQTL, PHENO, out_dir),
    ]:
        assert (repo_root / p).exists(), f"{p} should exist from a real prior run"
