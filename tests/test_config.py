from pathlib import Path

import pytest
import yaml
from jsonschema import ValidationError

from drugmr.config import Config

VALID = {
    "pheno_id": "AD", "sumstats": "dat/gwas/AD.tsv", "n_cases": 1, "n_controls": 1,
    "genome_build": "GRCh38", "target_build": "GRCh38",
    "snp_col": "SNP", "a1_col": "A1", "a2_col": "A2", "beta_col": "BETA",
    "se_col": "SE", "p_col": "P", "pos_col": "BP", "chr_col": "CHR", "af_col": "FRQ",
    "pqtl_dataset": "wingo_brain", "pqtl_dir": "dat/pqtl/wingo_brain", "ref_bfile": "dat/ref/x",
}


def _write(tmp_path, data, name="config.yaml"):
    path = tmp_path / name
    path.write_text(yaml.dump(data))
    return path


def test_valid_config_loads(tmp_path):
    path = _write(tmp_path, VALID)
    cfg = Config(path)
    assert cfg.pheno_id == "AD"
    assert cfg.pqtl_dataset == "wingo_brain"


def test_real_param_files_validate():
    for path in Path("params").glob("*.yaml"):
        cfg = Config(path)
        assert cfg.pheno_id
        assert cfg.pqtl_dataset


def test_missing_required_field_rejected(tmp_path):
    bad = dict(VALID)
    del bad["ref_bfile"]
    path = _write(tmp_path, bad)
    with pytest.raises(ValidationError):
        Config(path)


def test_invalid_pqtl_dataset_enum_rejected(tmp_path):
    bad = dict(VALID)
    bad["pqtl_dataset"] = "not_a_real_dataset"
    path = _write(tmp_path, bad)
    with pytest.raises(ValidationError):
        Config(path)


def test_gates_absent_falls_back_to_default(tmp_path):
    path = _write(tmp_path, VALID)
    cfg = Config(path)
    assert cfg.gate("coloc", "pp4_threshold", 0.7) == 0.7
    assert cfg.gate("cis_mr", "ivw_fdr_q", 0.05) == 0.05


def test_gates_present_overrides_default(tmp_path):
    with_gates = dict(VALID)
    with_gates["gates"] = {"coloc": {"pp4_threshold": 0.8}}
    path = _write(tmp_path, with_gates)
    cfg = Config(path)
    assert cfg.gate("coloc", "pp4_threshold", 0.7) == 0.8
    # a gate not set under a present step still falls back
    assert cfg.gate("cis_mr", "ivw_fdr_q", 0.05) == 0.05


def test_unknown_gate_key_rejected(tmp_path):
    bad = dict(VALID)
    bad["gates"] = {"coloc": {"made_up_threshold": 1}}
    path = _write(tmp_path, bad)
    with pytest.raises(ValidationError):
        Config(path)
