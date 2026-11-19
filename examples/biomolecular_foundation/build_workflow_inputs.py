"""Offline demonstration of agent planning and native input generation."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

from delfin.tools import run_step
from delfin.tools.bio.workflows import backend_catalog, validate_plan, workflow_template


def build_example(output: Path) -> dict:
    output = Path(output)
    output.mkdir(parents=True, exist_ok=False)
    catalog = backend_catalog()
    (output / "catalog.json").write_text(json.dumps(catalog, indent=2), encoding="utf-8")
    for name in ("docking", "medchem", "md", "alphafold3", "boltz2", "chai1"):
        plan = workflow_template(name)
        report = validate_plan(plan)
        if not report["valid"]:
            raise ValueError(report["errors"])
        (output / f"{name}_plan.json").write_text(json.dumps({"plan": plan, "validation": report}, indent=2), encoding="utf-8")
    request = {"schema_version": 1, "name": "synthetic_example", "seeds": [7],
               "msa_policy": "single_sequence", "entities": [
                   {"id": "A", "kind": "protein", "sequence": "MKWVTFISLL"},
                   {"id": "L", "kind": "ligand", "smiles": "CCO"}]}
    for backend in ("alphafold3", "boltz2", "chai1"):
        run_step("bio_prediction_input", backend=backend, request=request,
                 work_dir=output / backend).require()
    return {"backends": len(catalog["backends"]), "planning_templates_written": 6,
            "native_inputs_written": 3, "engines_executed": 0}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    print(json.dumps(build_example(parser.parse_args().output), indent=2))
