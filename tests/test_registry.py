import subprocess

from drugmr import registry


def test_current_git_sha7_matches_git_rev_parse():
    expected = subprocess.run(
        ["git", "rev-parse", "--short=7", "HEAD"], check=True, capture_output=True, text=True
    ).stdout.strip()
    assert registry.current_git_sha7() == expected


def test_load_registry_missing_file_returns_empty_dict(tmp_path):
    assert registry.load_registry(root=str(tmp_path / "runs")) == {}


def test_record_and_get_latest_run_id(tmp_path):
    root = str(tmp_path / "runs")
    registry.record_successful_run("AD", "ukb_ppp", "AD_ukb_ppp_20260811_a3a2aa1", root=root)
    assert registry.get_latest_run_id("AD", "ukb_ppp", root=root) == "AD_ukb_ppp_20260811_a3a2aa1"


def test_history_accumulates_and_latest_advances(tmp_path):
    root = str(tmp_path / "runs")
    registry.record_successful_run("AD", "ukb_ppp", "run_1", root=root)
    registry.record_successful_run("AD", "ukb_ppp", "run_2", root=root)
    reg = registry.load_registry(root=root)
    entry = reg["AD__ukb_ppp"]
    assert entry["latest"] == "run_2"
    assert entry["history"] == ["run_1", "run_2"]


def test_recording_same_run_id_twice_does_not_duplicate_history(tmp_path):
    root = str(tmp_path / "runs")
    registry.record_successful_run("AD", "ukb_ppp", "run_1", root=root)
    registry.record_successful_run("AD", "ukb_ppp", "run_1", root=root)
    entry = registry.load_registry(root=root)["AD__ukb_ppp"]
    assert entry["history"] == ["run_1"]


def test_different_pqtl_dataset_gets_its_own_key(tmp_path):
    root = str(tmp_path / "runs")
    registry.record_successful_run("AD", "ukb_ppp", "run_ukb", root=root)
    registry.record_successful_run("AD", "wingo_brain", "run_wingo", root=root)
    reg = registry.load_registry(root=root)
    assert set(reg.keys()) == {"AD__ukb_ppp", "AD__wingo_brain"}


def test_get_latest_run_id_unknown_key_returns_none(tmp_path):
    assert registry.get_latest_run_id("AD", "nonexistent", root=str(tmp_path / "runs")) is None


def test_write_manifest(tmp_path):
    root = str(tmp_path / "runs")
    registry.write_manifest("run_1", {"git_sha7": "abc1234", "backfilled": True}, root=root)
    import json
    from pathlib import Path
    manifest = json.loads((Path(root) / "run_1" / "manifest.json").read_text())
    assert manifest == {"git_sha7": "abc1234", "backfilled": True}
