from pathlib import Path

from drugmr import paths

PQTL = "ukb_ppp"
PHENO = "AD"
BULK = "GTEx_v10"
SC = "SingleBrain"


def test_qc_out():
    assert paths.qc_out(PHENO) == Path("dat/derived/AD/qc_gwas/AD.tsv")


def test_qc_mediator_dir():
    assert paths.qc_mediator_dir() == Path("dat/derived/mediators/qc_gwas")


def test_qc_mediator_out():
    assert paths.qc_mediator_out("AB42") == Path("dat/derived/mediators/qc_gwas/AB42.tsv")


def test_mr_out():
    assert paths.mr_out(PQTL, PHENO) == Path("results/cis-MR/ukb_ppp_AD_all_MR.tsv")


def test_mr_instruments_out():
    assert paths.mr_instruments_out(PQTL, PHENO) == Path(
        "results/cis-MR/instruments/ukb_ppp_AD_all_MR_instruments.tsv"
    )


def test_coloc_out():
    assert paths.coloc_out(PQTL, PHENO) == Path("results/coloc/ukb_ppp/ukb_ppp_AD_all_coloc.tsv")


def test_network_mr_m_y_out_is_tier1_pheno_scoped():
    # M_Y is mediator -> outcome genome-wide MR, shared across every pqtl_dataset
    assert paths.network_mr_m_y_out(PHENO) == Path(
        "dat/derived/AD/network_mr/M_Y/AD_mediator_genomewide_MR.tsv"
    )


def test_network_mr_x_m_out_is_tier2_per_dataset():
    assert paths.network_mr_x_m_out(PQTL, "AB42") == Path(
        "results/network_mr/X_M/ukb_ppp/ukb_ppp_AB42_all_MR.tsv"
    )


def test_network_mr_mediation_estimates_out_matches_actual_writer():
    # This is the literal bin/assort_network_mr.py's perform_network_mr() and
    # bin/coloc_targets.py's coloc_with_mediators() both read/write, and what
    # drugmr/local.py + drugmr/hpc.py now gate the NetworkMR step on directly.
    assert paths.network_mr_mediation_estimates_out(PQTL, PHENO) == Path(
        "results/network_mr/mediation_estimates/ukb_ppp/ukb_ppp_AD_networkMR.tsv"
    )


def test_target_stats_out():
    assert paths.target_stats_out(PQTL, PHENO) == Path(
        "results/target_stats/ukb_ppp/AD/ukb_ppp_AD_top_cis_hits.tsv"
    )


def test_smr_bulk_out():
    assert paths.smr_bulk_out(PQTL, PHENO, BULK) == Path(
        "results/SMR/bulk/GTEx_v10/AD/ukb_ppp_AD_promising_targets_SMR.tsv"
    )


def test_smr_sc_out():
    assert paths.smr_sc_out(PQTL, PHENO, SC) == Path(
        "results/SMR/sc/SingleBrain/AD/ukb_ppp_AD_promising_targets_SMR.tsv"
    )


def test_smr_bulk_dir():
    assert paths.smr_bulk_dir(BULK) == Path("results/SMR/bulk/GTEx_v10")


def test_smr_final_targets_out():
    assert paths.smr_final_targets_out(PQTL, PHENO) == Path(
        "results/SMR/ukb_ppp_AD_final_multi_omics_targets.tsv"
    )


def test_hyprcoloc_out():
    assert paths.hyprcoloc_out(PQTL, PHENO) == Path("results/hyprcoloc/ukb_ppp_AD_all_hyprcoloc.tsv")


def test_hyprcoloc_dataset_out():
    assert paths.hyprcoloc_dataset_out(PQTL, BULK, PHENO) == Path(
        "results/hyprcoloc/ukb_ppp/GTEx_v10/AD_hyprcoloc.tsv"
    )


def test_phewas_out():
    assert paths.phewas_out(PQTL, PHENO) == Path("results/PheWAS/ukb_ppp/AD/ukb_ppp_AD_PheWAS.tsv")


def test_phewas_ukbb_out():
    assert paths.phewas_ukbb_out(PQTL, PHENO) == Path(
        "results/PheWAS_UKBB/ukb_ppp/AD/ukb_ppp_AD_PheWAS.tsv"
    )


def test_smr_raw_prefix_matches_smr_binary_out_prefix():
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
        "runs/AD_ukb_ppp_20260811_a3a2aa1/results/cis-MR/ukb_ppp_AD_all_MR.tsv"
    )


def test_paths_match_real_files_on_disk():
    """Sanity check against actual files already on disk from real prior runs -
    not just the literals, but that paths.py resolves to files that exist."""
    repo_root = Path(__file__).resolve().parents[1]
    for p in [
        paths.mr_out(PQTL, PHENO),
        paths.coloc_out(PQTL, PHENO),
        paths.hyprcoloc_out(PQTL, PHENO),
        paths.smr_final_targets_out(PQTL, PHENO),
        paths.phewas_out(PQTL, PHENO),
        paths.phewas_ukbb_out(PQTL, PHENO),
        paths.target_stats_out(PQTL, PHENO),
    ]:
        assert (repo_root / p).exists(), f"{p} should exist from a real prior run"
