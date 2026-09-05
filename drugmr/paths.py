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

NetworkMR (mediator/biomarker mediation analysis - protein -> mediator ->
outcome, plus its MOLOC/HyPrColoc-mediator colocalisation branches) was
removed from the active pipeline and archived to analysis/networkmr/ - real
bugs found in it (cis_mr.R silently ignoring the mediator outcome and
re-running protein -> AD instead, coloc.R hardcoding the mediator's trait
type as case-control, several stale-glob/idempotency issues) made it not
worth keeping wired in until revisited properly. See analysis/networkmr/
for the archived code if this is ever picked back up.

Phase 4 (2026-09-04): stage dirs below are normalised to one casing
convention (snake_case) and consistently named for what they are, not what
dataset/trait produced them - `run_id` (the run's own directory name) already
encodes both `pqtl_dataset` and `pheno_id`, so repeating them inside
results/ or work/ is redundant (a run is always exactly 1 dataset x 1
trait). `pqtl_dataset`/`pheno_id` stay as parameters on every function below
for call-site stability - many callers already have them in scope for other
reasons (filtering, logging, upstream lookups) - they're just no longer
embedded in the returned Path. `pwcoco` and `pwcoco_qtl` are consolidated
into one `pwcoco/` stage with `cis_pqtl/`, `eqtl_pqtl/`, `eqtl_gwas/`,
`summary/` children - they're the same tool run against different trait
pairs, not different stages. Per-locus PWCoCo output is grouped one
directory per protein (or protein_eqtlsource) instead of a flat dump of
hundreds of files. `bin/coloc.R`, `bin/hyprcoloc.R`, `bin/cis_mr.R` compute
matching `out_dir` values independently (see their own comments) and must
stay in lock-step with the functions below. `work_dir_for_results_dir()`
gives every stage's scratch/intermediate output (deleted before the
pipeline finishes, never a final artifact) the same run-scoping, without
needing `run_id` itself threaded through every bin/*.py call site.
`smr_raw_dir`/`smr_raw_prefix`/`smr_bulk_dir` and everything under
`synthesis/` are UNCHANGED - they're a dataset-independent raw-SMR cache
shared across every pqtl_dataset, not part of a run's results/ tree.
`runs/*/results/network_mr/` and `runs/*/results/colocboost/` from before
those stages were removed were archived to `runs/<run_id>/_archived_<stage>/`
by `analysis/migrate_results_schema.py`, not migrated into this schema.

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


def mr_out(pqtl_dataset: str, pheno_id: str, out_dir: str = "results") -> Path:
    return Path(out_dir) / "cis_mr" / "mr.tsv"


def mr_instruments_out(pqtl_dataset: str, pheno_id: str, out_dir: str = "results") -> Path:
    return Path(out_dir) / "cis_mr" / "instruments" / "mr_instruments.tsv"


def coloc_out(pqtl_dataset: str, pheno_id: str, out_dir: str = "results") -> Path:
    return Path(out_dir) / "coloc" / "coloc.tsv"


def coloc_susie_out(pqtl_dataset: str, pheno_id: str, out_dir: str = "results") -> Path:
    """SuSiE-based sensitivity check on top of coloc_out() - informational only,
    nothing downstream gates on this file. Currently dead code: no producer
    exists (there is no bin/coloc_susie.py/R) - kept for schema consistency."""
    return Path(out_dir) / "coloc" / "coloc_susie.tsv"


def coloc_sensitivity_out(pqtl_dataset: str, pheno_id: str, out_dir: str = "results") -> Path:
    """Prior sensitivity check on top of coloc_out() - reruns coloc.abf() per
    protein under the named prior scenarios in bin/coloc.R's sensitivity_priors
    list. Informational only, nothing downstream gates on this file. Lives
    alongside coloc_out() (not its own top-level dir) since bin/coloc.R writes
    both the primary and sensitivity TSV for a protein into the same out_dir."""
    return Path(out_dir) / "coloc" / "coloc_sensitivity.tsv"


def pwcoco_raw_dir(pqtl_dataset: str, pheno_id: str, out_dir: str = "results") -> Path:
    """Per-protein PWCoCo --out prefix lives inside this dir. PWCoCo appends
    rather than overwrites its .coloc/.cojo output, so each protein needs its
    own unique prefix - this dir gives bin/pwcoco_targets.py somewhere to put
    per-protein raw PWCoCo output before it's parsed/aggregated into
    pwcoco_out() below."""
    return Path(out_dir) / "pwcoco" / "cis_pqtl"


def pwcoco_raw_prefix(pqtl_dataset: str, pheno_id: str, protein: str, out_dir: str = "results") -> Path:
    # nested 1 extra level (protein/protein) so every file PWCoCo appends to
    # this prefix (.coloc, .sumstats*.included/.badfreq, per-rsID .cojo) lands
    # in its own directory instead of a flat dump shared by every protein
    return pwcoco_raw_dir(pqtl_dataset, pheno_id, out_dir) / protein / protein


def pwcoco_out(pqtl_dataset: str, pheno_id: str, out_dir: str = "results") -> Path:
    """Aggregated PWCoCo results across all proteins for this pqtl_dataset x
    pheno_id run - what bin/pwcoco_targets.py parses pwcoco_raw_prefix()'s
    per-protein .coloc files into. Same shape/naming as coloc_out() so the
    coloc_support concordance join between the two can be a plain merge on
    (pqtl_dataset, pheno_id, protein)."""
    return Path(out_dir) / "pwcoco" / "summary" / "pwcoco.tsv"


def target_stats_out(pqtl_dataset: str, pheno_id: str, out_dir: str = "results") -> Path:
    return Path(out_dir) / "target_stats" / "top_cis_hits.tsv"


def smr_bulk_out(pqtl_dataset: str, pheno_id: str, bulk_dataset: str, out_dir: str = "results") -> Path:
    return Path(out_dir) / "smr" / "bulk" / bulk_dataset / "promising_targets.tsv"


def smr_sc_out(pqtl_dataset: str, pheno_id: str, sc_eqtl_dataset: str, out_dir: str = "results") -> Path:
    return Path(out_dir) / "smr" / "sc" / sc_eqtl_dataset / "promising_targets.tsv"


def smr_final_targets_out(pqtl_dataset: str, pheno_id: str, out_dir: str = "results") -> Path:
    return Path(out_dir) / "smr" / "final_multi_omics_targets.tsv"


def smr_bulk_dir(bulk_dataset: str, out_dir: str = "results") -> Path:
    """Dataset-level bulk SMR dir (not yet scoped to pheno_id) - used to rglob for
    pre-existing .smr files regardless of their exact sub-nesting convention.
    UNCHANGED - synthesis/ tier, dataset-independent (see module docstring)."""
    return Path(out_dir) / "SMR" / "bulk" / bulk_dataset


def hyprcoloc_out(pqtl_dataset: str, pheno_id: str, out_dir: str = "results") -> Path:
    return Path(out_dir) / "hyprcoloc" / "hyprcoloc.tsv"


def hyprcoloc_dataset_out(pqtl_dataset: str, hc_dataset: str, pheno_id: str, out_dir: str = "results") -> Path:
    return Path(out_dir) / "hyprcoloc" / "by_eqtl_source" / hc_dataset / "hyprcoloc.tsv"


def phewas_out(pqtl_dataset: str, pheno_id: str, out_dir: str = "results") -> Path:
    return Path(out_dir) / "phewas" / "finngen" / "phewas.tsv"


def phewas_ukbb_out(pqtl_dataset: str, pheno_id: str, out_dir: str = "results") -> Path:
    return Path(out_dir) / "phewas" / "ukbb" / "phewas.tsv"


def phewas_finngen_coverage_out(pqtl_dataset: str, pheno_id: str, out_dir: str = "results") -> Path:
    """Per-target manifest of whether >=1 cis-MR instrument was found in FinnGen.

    Read by bin/ukb_phewas.py to restrict UKB PheWAS to the fallback set only
    (targets with zero retained instruments in FinnGen), matching the paper's
    "phenome-wide MR was instead performed across UKB" fallback design.
    """
    return Path(out_dir) / "phewas" / "finngen" / "phewas_coverage.tsv"


def smr_raw_dir(eqtl_dataset, pheno_id: str, out_dir: str = "results") -> Path:
    """Directory for the raw `smr` binary's own output (drugmr.smr.SMR()'s
    --out prefix lives inside this dir - see smr_raw_prefix). eqtl_dataset
    may be a compound relative path (e.g. "bulk_raw/GTEx_v10/label/chr1").
    UNCHANGED - synthesis/ tier, dataset-independent (see module docstring)."""
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


def work_dir_for_results_dir(local_results_dir) -> Path:
    """Sibling run-scoped work/ dir for a given run's results dir. Every
    bin/*.py step already receives `local_results_dir` (== run_results_dir(run_id)
    == runs/<run_id>/results), so this derives runs/<run_id>/work from it without
    needing run_id threaded through as a separate CLI arg. Scratch/intermediate
    output that gets deleted before the step finishes (trio parquets, raw PheWAS
    hit dumps, the SMR .ma GWAS reformat) belongs under here, one subdir per
    stage - never under results/, which is final-artifacts-only."""
    return Path(local_results_dir).parent / "work"


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
    targets - complements pwcoco_raw_dir() above (the pQTL-GWAS PWCoCo). Sibling of
    cis_pqtl/ under the same pwcoco/ stage - same tool, different trait pairs."""
    return Path(out_dir) / "pwcoco" / combo


def pwcoco_qtl_raw_prefix(combo: str, pqtl_dataset: str, protein: str, eqtl_source: str, out_dir: str = "results") -> Path:
    # nested 1 extra level, same reasoning as pwcoco_raw_prefix() above
    name = f"{protein}_{eqtl_source}"
    return pwcoco_qtl_raw_dir(combo, pqtl_dataset, out_dir) / name / name


def pwcoco_eqtl_pqtl_out(pqtl_dataset: str, pheno_id: str, out_dir: str = "results") -> Path:
    """Aggregated eQTL-pQTL PWCoCo results (1 row per PWCoCo output row - unconditioned
    plus any conditioned rows - across every SMR-passing protein x eqtl_source pair)."""
    return Path(out_dir) / "pwcoco" / "summary" / "pwcoco_eqtl_pqtl.tsv"


def pwcoco_eqtl_gwas_out(pqtl_dataset: str, pheno_id: str, out_dir: str = "results") -> Path:
    """Same shape as pwcoco_eqtl_pqtl_out() above, for the eQTL-GWAS combo."""
    return Path(out_dir) / "pwcoco" / "summary" / "pwcoco_eqtl_gwas.tsv"


def pwcoco_qtl_shared_out(pqtl_dataset: str, pheno_id: str, out_dir: str = "results") -> Path:
    """Per-protein shared-SNP table across the 3 PWCoCo combos (pQTL-GWAS from
    pwcoco_out() above, eQTL-pQTL, eQTL-GWAS) - which SNPs colocalise (H4 above
    threshold) in ALL 3 combos simultaneously, with each combo's own H4 for that
    SNP. Co-equal to HyPrColoc (same "1 causal variant across pQTL/eQTL/GWAS"
    question, tested via conditioning instead of HyPrColoc's single-cluster
    assumption), not a downstream refinement of it."""
    return Path(out_dir) / "pwcoco" / "summary" / "pwcoco_shared_snps.tsv"
