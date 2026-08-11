#!/usr/bin/env python3
"""
Registry of pipeline runs.

`runs/registry.json` maps `{pheno_id}__{pqtl_dataset}` to `{"latest": run_id,
"history": [run_id, ...]}`. It is written ONLY after every step of a run has
succeeded (drugmr/local.py, drugmr/hpc.py call `record_successful_run` as
their last line). A failed or partial run's `runs/<run_id>/` directory is
never referenced here, so a consumer reading "latest" can never be pointed
at a broken run.

Kept separate from drugmr/paths.py, which stays pure path arithmetic with no
I/O, this module is where the actual JSON reads/writes and the git/date
lookups that feed a run_id live.
"""
import json
import subprocess
from pathlib import Path

from drugmr import paths


def _registry_key(pheno_id: str, pqtl_dataset: str) -> str:
    return f"{pheno_id}__{pqtl_dataset}"


def current_git_sha7(cwd=None) -> str:
    result = subprocess.run(
        ["git", "rev-parse", "--short=7", "HEAD"],
        cwd=cwd,
        check=True,
        capture_output=True,
        text=True,
    )
    return result.stdout.strip()


def load_registry(root: str = "runs") -> dict:
    path = paths.registry_path(root)
    if not path.exists():
        return {}
    with open(path, "r") as f:
        return json.load(f)


def save_registry(registry: dict, root: str = "runs") -> None:
    path = paths.registry_path(root)
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w") as f:
        json.dump(registry, f, indent=2, sort_keys=True)


def record_successful_run(pheno_id: str, pqtl_dataset: str, run_id: str, root: str = "runs") -> None:
    """Call this ONLY once every step of a run has succeeded."""
    registry = load_registry(root)
    key = _registry_key(pheno_id, pqtl_dataset)
    entry = registry.setdefault(key, {"latest": None, "history": []})
    if run_id not in entry["history"]:
        entry["history"].append(run_id)
    entry["latest"] = run_id
    save_registry(registry, root)


def get_latest_run_id(pheno_id: str, pqtl_dataset: str, root: str = "runs"):
    registry = load_registry(root)
    entry = registry.get(_registry_key(pheno_id, pqtl_dataset))
    return entry["latest"] if entry else None


def write_manifest(run_id: str, manifest: dict, root: str = "runs") -> None:
    path = paths.run_manifest_path(run_id, root)
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w") as f:
        json.dump(manifest, f, indent=2, sort_keys=True, default=str)
