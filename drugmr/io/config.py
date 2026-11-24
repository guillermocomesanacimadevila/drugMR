#!/usr/bin/env python3
import json
from pathlib import Path

import jsonschema
import yaml

_SCHEMA_PATH = Path(__file__).resolve().parents[1] / "params" / "schema.json"


class Config:

    def __init__(self, file):
        file = Path(file)
        with open(file, "r") as f:
            data = yaml.safe_load(f)

        schema = json.loads(_SCHEMA_PATH.read_text())
        jsonschema.validate(instance=data, schema=schema)

        self.__dict__.update(data)
        reference_data = data.get("reference_data", {}) or {}
        if reference_data.get("ref_bfile"):
            self.ref_bfile = reference_data["ref_bfile"]
        self.liftover_dir = reference_data.get(
            "liftover_dir", data.get("liftover_dir", "dat/ref/liftover")
        )
        self.ncbi_ref_path = reference_data.get(
            "gene_annotation", data.get("gene_annotation")
        )

    def gate(self, step: str, name: str, default=None):
        gates = getattr(self, "gates", {}) or {}
        return gates.get(step, {}).get(name, default)
