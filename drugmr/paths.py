#!/usr/bin/env python3
"""
Central output-path resolver for the drugMR pipeline.

Every function here returns a path *relative to the project root* (a
plain pathlib.Path, never made absolute) so both `drugmr/local.py`
(which joins it onto `project_root`) and `drugmr/hpc.py` (which uses it
directly as a string relative to the remote repo checkout, since every
`ssh(...)` call already does `cd "{remote}"` first) can share the exact
same path logic without this module knowing anything about SSH.

Migration status: this module currently reproduces the exact literals
that were previously hardcoded independently in `drugmr/local.py`,
`drugmr/hpc.py`, and `bin/*.py` - byte-for-byte, including the
`out_dir` argument being ignored by every caller that didn't already
thread it through (now fixed here: every function takes `out_dir` and
uses it, defaulting to "results" - the same effective value every
existing config falls back to, so no run today changes location).
Phase 3 (runs/<run_id>/ + registry): Tier-2 functions' `out_dir` is now
computed by callers as `run_results_dir(run_id)` instead of a bare
"results" default - see `drugmr/registry.py` for run_id computation and
the registry.json/manifest.json read/write logic (deliberately kept out
of this module, which stays pure path arithmetic with no I/O).

Phase 5 (network_mr + synthesis/): the standalone `network_mr_out()` gate
(which never matched what assort_network_mr.py actually wrote) is retired -
`drugmr/local.py`/`drugmr/hpc.py` now gate NetworkMR directly on
`network_mr_mediation_estimates_out()`, the function that's actually written
to. M_Y (mediator -> outcome, pheno-scoped only) moved to Tier 1
(`network_mr_m_y_*`); X_M and mediation_estimates (per pqtl_dataset) stay
Tier 2 (`network_mr_x_m_*`, `network_mr_mediation_estimates_out`).

Known pre-existing mismatch preserved as-is (not fixed by this module):
  - `dat/bulk-eQTL` (singular) vs `dat/bulk-eQTLs` (plural, only in
    `scripts/GTEx_v10/eqtl_gtex_eqtl.py`) is a separate, untouched
    input-path mismatch, unrelated to this module.
"""
from pathlib import Path


def qc_out(pheno_id: str) -> Path:
    """Tier-1 shared preprocessing output (dat/derived/) - not scoped to a
    run's out_dir, since it's reused across every pqtl_dataset run for this
    pheno_id rather than regenerated per run."""
    return Path("dat") / "derived" / pheno_id / "qc_gwas" / f"{pheno_id}.tsv"


def qc_mediator_dir() -> Path:
    """Tier-1 shared mediator QC directory - flat, one file per mediator,
    not scoped to pheno_id or pqtl_dataset (mediators are dataset-level)."""
    return Path("dat") / "derived" / "mediators" / "qc_gwas"


def qc_mediator_out(mediator_id: str) -> Path:
    return qc_mediator_dir() / f"{mediator_id}.tsv"


def mr_out(pqtl_dataset: str, pheno_id: str, out_dir: str = "results") -> Path:
    return Path(out_dir) / "cis-MR" / f"{pqtl_dataset}_{pheno_id}_all_MR.tsv"


def mr_instruments_out(pqtl_dataset: str, pheno_id: str, out_dir: str = "results") -> Path:
    return (
        Path(out_dir)
        / "cis-MR"
        / "instruments"
        / f"{pqtl_dataset}_{pheno_id}_all_MR_instruments.tsv"
    )


def coloc_out(pqtl_dataset: str, pheno_id: str, out_dir: str = "results") -> Path:
    return Path(out_dir) / "coloc" / pqtl_dataset / f"{pqtl_dataset}_{pheno_id}_all_coloc.tsv"


def coloc_susie_out(pqtl_dataset: str, pheno_id: str, out_dir: str = "results") -> Path:
    """SuSiE-based sensitivity check on top of coloc_out() - informational only,
    nothing downstream gates on this file (see bin/coloc_susie.py)."""
    return Path(out_dir) / "coloc_susie" / pqtl_dataset / f"{pqtl_dataset}_{pheno_id}_all_coloc_susie.tsv"


def coloc_sensitivity_out(pqtl_dataset: str, pheno_id: str, out_dir: str = "results") -> Path:
    """Prior sensitivity check on top of coloc_out() - reruns coloc.abf() per
    protein under the named prior scenarios in bin/coloc.R's sensitivity_priors
    list. Informational only, nothing downstream gates on this file. Lives
    alongside coloc_out() (not its own top-level dir) since bin/coloc.R writes
    both the primary and sensitivity TSV for a protein into the same out_dir."""
    return Path(out_dir) / "coloc" / pqtl_dataset / f"{pqtl_dataset}_{pheno_id}_all_coloc_sensitivity.tsv"


def coloc_mediator_out(pqtl_dataset: str, pheno_id: str, out_dir: str = "results") -> Path:
    """Protein x mediator COLOC results from coloc_with_mediators() (bin/coloc_targets.py) -
    the M-side pairs used to build the moloc.json candidate list, kept separate
    from coloc_out()'s protein x outcome (Y) pairs since it's a different pairing."""
    return Path(out_dir) / "coloc" / pqtl_dataset / f"{pqtl_dataset}_{pheno_id}_all_coloc_mediators.tsv"


def coloc_mediator_sensitivity_out(pqtl_dataset: str, pheno_id: str, out_dir: str = "results") -> Path:
    """Prior sensitivity check on top of coloc_mediator_out() - see coloc_sensitivity_out()."""
    return Path(out_dir) / "coloc" / pqtl_dataset / f"{pqtl_dataset}_{pheno_id}_all_coloc_mediators_sensitivity.tsv"


def pwcoco_raw_dir(pqtl_dataset: str, pheno_id: str, out_dir: str = "results") -> Path:
    """Per-protein PWCoCo --out prefix lives inside this dir. PWCoCo appends
    rather than overwrites its .coloc/.cojo output, so each protein needs its
    own unique prefix - this dir gives bin/pwcoco_targets.py somewhere to put
    per-protein raw PWCoCo output before it's parsed/aggregated into
    pwcoco_out() below."""
    return Path(out_dir) / "pwcoco" / pqtl_dataset / pheno_id


def pwcoco_raw_prefix(pqtl_dataset: str, pheno_id: str, protein: str, out_dir: str = "results") -> Path:
    return pwcoco_raw_dir(pqtl_dataset, pheno_id, out_dir) / protein


def pwcoco_out(pqtl_dataset: str, pheno_id: str, out_dir: str = "results") -> Path:
    """Aggregated PWCoCo results across all proteins for this pqtl_dataset x
    pheno_id run - what bin/pwcoco_targets.py parses pwcoco_raw_prefix()'s
    per-protein .coloc files into. Same shape/naming as coloc_out() so the
    coloc_support concordance join between the two can be a plain merge on
    (pqtl_dataset, pheno_id, protein)."""
    return Path(out_dir) / "pwcoco" / pqtl_dataset / f"{pqtl_dataset}_{pheno_id}_all_pwcoco.tsv"


def network_mr_m_y_dir(pheno_id: str) -> Path:
    """Tier-1: mediator -> outcome genome-wide MR. Shared across every
    pqtl_dataset run for this pheno_id (run_genomewide_mr() takes no
    pqtl_dataset at all), so this lives under dat/derived/ like qc_out."""
    return Path("dat") / "derived" / pheno_id / "network_mr" / "M_Y"


def network_mr_m_y_out(pheno_id: str) -> Path:
    return network_mr_m_y_dir(pheno_id) / f"{pheno_id}_mediator_genomewide_MR.tsv"


def network_mr_x_m_dir(pqtl_dataset: str, out_dir: str = "results") -> Path:
    """Tier-2: protein -> mediator cis-MR, one directory per pqtl_dataset run."""
    return Path(out_dir) / "network_mr" / "X_M" / pqtl_dataset


def network_mr_x_m_out(pqtl_dataset: str, mediator_id: str, out_dir: str = "results") -> Path:
    return network_mr_x_m_dir(pqtl_dataset, out_dir) / f"{pqtl_dataset}_{mediator_id}_all_MR.tsv"


def network_mr_mediation_estimates_out(pqtl_dataset: str, pheno_id: str, out_dir: str = "results") -> Path:
    """What bin/assort_network_mr.py's perform_network_mr() and
    bin/coloc_targets.py's coloc_with_mediators() both read/write - the final
    NetworkMR-package output, also what the orchestrator gates NetworkMR on."""
    return (
        Path(out_dir)
        / "network_mr"
        / "mediation_estimates"
        / pqtl_dataset
        / f"{pqtl_dataset}_{pheno_id}_networkMR.tsv"
    )


def target_stats_out(pqtl_dataset: str, pheno_id: str, out_dir: str = "results") -> Path:
    return (
        Path(out_dir)
        / "target_stats"
        / pqtl_dataset
        / pheno_id
        / f"{pqtl_dataset}_{pheno_id}_top_cis_hits.tsv"
    )


def smr_bulk_out(pqtl_dataset: str, pheno_id: str, bulk_dataset: str, out_dir: str = "results") -> Path:
    return (
        Path(out_dir)
        / "SMR"
        / "bulk"
        / bulk_dataset
        / pheno_id
        / f"{pqtl_dataset}_{pheno_id}_promising_targets_SMR.tsv"
    )


def smr_sc_out(pqtl_dataset: str, pheno_id: str, sc_eqtl_dataset: str, out_dir: str = "results") -> Path:
    return (
        Path(out_dir)
        / "SMR"
        / "sc"
        / sc_eqtl_dataset
        / pheno_id
        / f"{pqtl_dataset}_{pheno_id}_promising_targets_SMR.tsv"
    )


def smr_final_targets_out(pqtl_dataset: str, pheno_id: str, out_dir: str = "results") -> Path:
    return Path(out_dir) / "SMR" / f"{pqtl_dataset}_{pheno_id}_final_multi_omics_targets.tsv"


def smr_bulk_dir(bulk_dataset: str, out_dir: str = "results") -> Path:
    """Dataset-level bulk SMR dir (not yet scoped to pheno_id) - used to rglob for
    pre-existing .smr files regardless of their exact sub-nesting convention."""
    return Path(out_dir) / "SMR" / "bulk" / bulk_dataset


def hyprcoloc_out(pqtl_dataset: str, pheno_id: str, out_dir: str = "results") -> Path:
    return Path(out_dir) / "hyprcoloc" / f"{pqtl_dataset}_{pheno_id}_all_hyprcoloc.tsv"


def hyprcoloc_dataset_out(pqtl_dataset: str, hc_dataset: str, pheno_id: str, out_dir: str = "results") -> Path:
    return Path(out_dir) / "hyprcoloc" / pqtl_dataset / hc_dataset / f"{pheno_id}_hyprcoloc.tsv"


def phewas_out(pqtl_dataset: str, pheno_id: str, out_dir: str = "results") -> Path:
    return Path(out_dir) / "PheWAS-FinnGen" / pqtl_dataset / pheno_id / f"{pqtl_dataset}_{pheno_id}_PheWAS-FinnGen.tsv"


def phewas_ukbb_out(pqtl_dataset: str, pheno_id: str, out_dir: str = "results") -> Path:
    return Path(out_dir) / "PheWAS_UKBB" / pqtl_dataset / pheno_id / f"{pqtl_dataset}_{pheno_id}_PheWAS.tsv"


def phewas_finngen_coverage_out(pqtl_dataset: str, pheno_id: str, out_dir: str = "results") -> Path:
    """Per-target manifest of whether >=1 cis-MR instrument was found in FinnGen.

    Read by bin/ukb_phewas.py to restrict UKB PheWAS to the fallback set only
    (targets with zero retained instruments in FinnGen), matching the paper's
    "phenome-wide MR was instead performed across UKB" fallback design.
    """
    return Path(out_dir) / "PheWAS-FinnGen" / pqtl_dataset / pheno_id / f"{pqtl_dataset}_{pheno_id}_PheWAS-FinnGen_coverage.tsv"

def smr_raw_dir(eqtl_dataset, pheno_id: str, out_dir: str = "results") -> Path:
    """Directory for the raw `smr` binary's own output (drugmr.smr.SMR()'s
    --out prefix lives inside this dir - see smr_raw_prefix). eqtl_dataset
    may be a compound relative path (e.g. "bulk_raw/GTEx_v10/label/chr1")."""
    return Path(out_dir) / "SMR" / Path(eqtl_dataset) / pheno_id


def smr_raw_prefix(eqtl_dataset, pheno_id: str, out_dir: str = "results") -> Path:
    eqtl_name = Path(eqtl_dataset).name
    return smr_raw_dir(eqtl_dataset, pheno_id, out_dir) / f"{pheno_id}_{eqtl_name}"


def make_run_id(pheno_id: str, pqtl_dataset: str, date_str: str, git_sha7: str) -> str:
    return f"{pheno_id}_{pqtl_dataset}_{date_str}_{git_sha7}"


def run_dir(run_id: str, root: str = "runs") -> Path:
    return Path(root) / run_id


def run_results_dir(run_id: str, root: str = "runs") -> Path:
    return run_dir(run_id, root) / "results"


def run_work_dir(run_id: str, root: str = "runs") -> Path:
    return run_dir(run_id, root) / "work"


def run_logs_dir(run_id: str, root: str = "runs") -> Path:
    return run_dir(run_id, root) / "logs"


def run_manifest_path(run_id: str, root: str = "runs") -> Path:
    return run_dir(run_id, root) / "manifest.json"


def run_params_lock_path(run_id: str, root: str = "runs") -> Path:
    return run_dir(run_id, root) / "params.lock.yaml"


def registry_path(root: str = "runs") -> Path:
    return Path(root) / "registry.json"


def synthesis_dir(pheno_id: str, root: str = "synthesis") -> Path:
    """Tier-3: genuinely cross-dataset aggregation for a pheno_id. No automated
    producer exists yet as of Phase 5 - see synthesis_target_stats_out()."""
    return Path(root) / pheno_id


def synthesis_target_stats_out(pheno_id: str, root: str = "synthesis") -> Path:
    return synthesis_dir(pheno_id, root) / "target_stats" / "all_datasets_mined_targets.tsv"


def synthesis_manifest_path(pheno_id: str, root: str = "synthesis") -> Path:
    return synthesis_dir(pheno_id, root) / "manifest.json"


def pwcoco_qtl_raw_dir(combo: str, pqtl_dataset: str, out_dir: str = "results") -> Path:
    """Per-(protein, eqtl_source) PWCoCo --out prefix dir for the eQTL-informed combos
    (combo: "eqtl_pqtl" or "eqtl_gwas") bin/pwcoco_qtl_wrapper.py runs on SMR-passing
    targets - complements pwcoco_raw_dir() above (the pQTL-GWAS PWCoCo)."""
    return Path(out_dir) / "pwcoco_qtl" / combo / pqtl_dataset


def pwcoco_qtl_raw_prefix(combo: str, pqtl_dataset: str, protein: str, eqtl_source: str, out_dir: str = "results") -> Path:
    return pwcoco_qtl_raw_dir(combo, pqtl_dataset, out_dir) / f"{protein}_{eqtl_source}"


def pwcoco_eqtl_pqtl_out(pqtl_dataset: str, pheno_id: str, out_dir: str = "results") -> Path:
    """Aggregated eQTL-pQTL PWCoCo results (1 row per PWCoCo output row - unconditioned
    plus any conditioned rows - across every SMR-passing protein x eqtl_source pair)."""
    return Path(out_dir) / "pwcoco_qtl" / pqtl_dataset / f"{pqtl_dataset}_{pheno_id}_all_pwcoco_eqtl_pqtl.tsv"


def pwcoco_eqtl_gwas_out(pqtl_dataset: str, pheno_id: str, out_dir: str = "results") -> Path:
    """Same shape as pwcoco_eqtl_pqtl_out() above, for the eQTL-GWAS combo."""
    return Path(out_dir) / "pwcoco_qtl" / pqtl_dataset / f"{pqtl_dataset}_{pheno_id}_all_pwcoco_eqtl_gwas.tsv"


def pwcoco_qtl_shared_out(pqtl_dataset: str, pheno_id: str, out_dir: str = "results") -> Path:
    """Per-protein shared-SNP table across the 3 PWCoCo combos (pQTL-GWAS from
    pwcoco_out() above, eQTL-pQTL, eQTL-GWAS) - which SNPs colocalise (H4 above
    threshold) in ALL 3 combos simultaneously, with each combo's own H4 for that
    SNP. Co-equal to HyPrColoc (same "1 causal variant across pQTL/eQTL/GWAS"
    question, tested via conditioning instead of HyPrColoc's single-cluster
    assumption), not a downstream refinement of it."""
    return Path(out_dir) / "pwcoco_qtl" / pqtl_dataset / f"{pqtl_dataset}_{pheno_id}_pwcoco_shared_snps.tsv"
