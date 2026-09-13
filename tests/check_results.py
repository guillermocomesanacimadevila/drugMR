from pathlib import Path


runs = sorted(Path("tests/output/runs").glob("TEST_test_pqtl_*"))
if not runs:
    raise SystemExit("The pipeline did not create a test run")

results = runs[-1] / "results"
required = [
    "cis_mr/mr.tsv",
    "cis_mr/instruments/mr_instruments.tsv",
    "coloc/test_pqtl_TEST_all_coloc.tsv",
    "pwcoco/summary/pwcoco.tsv",
    "target_stats/top_cis_hits.tsv",
    "locus_data/index.tsv",
    "smr/bulk/test_eqtl/promising_targets.tsv",
    "smr/final_multi_omics_targets.tsv",
    "hyprcoloc/by_eqtl_source/test_eqtl/hyprcoloc.tsv",
    "hyprcoloc/hyprcoloc.tsv",
    "phewas/finngen/phewas_coverage.tsv",
]

missing = [name for name in required if not (results / name).exists()]
if missing:
    raise SystemExit("Missing results: " + ", ".join(missing))

print("End-to-end pipeline test passed")
