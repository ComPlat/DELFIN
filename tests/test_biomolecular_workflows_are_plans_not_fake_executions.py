import copy
import json

import pytest

from delfin.tools.bio.workflows import BACKENDS, backend_catalog, validate_plan, workflow_template
from delfin.tools.bio.prediction_inputs import build_prediction_input


def request(*entities, policy="single_sequence"):
    return {"schema_version": 1, "name": "test_job", "seeds": [7],
            "msa_policy": policy, "entities": list(entities) or [{"id": "A", "kind": "protein", "sequence": "MKWVTFISLL"}]}


def test_catalog_covers_the_domains_and_never_claims_engine_readiness():
    report = backend_catalog()
    rows = {r["id"]: r for r in report["backends"]}
    assert {"alphafold2", "alphafold_multimer", "alphafold3", "colabfold", "boltz2", "chai1", "esmfold",
            "vina", "gnina", "diffdock", "openmm", "gromacs", "proteinmpnn", "medchem_properties"} <= rows.keys()
    assert len(rows) == len(BACKENDS)
    assert all(not r["execution_ready"] and r["command_paths"] is None for r in rows.values())
    assert rows["rdkit_prepare"]["adapter"] == "ligand_prepare"
    assert rows["alphafold3"]["integration"] == "planned"


def test_probe_uses_only_supplied_path_and_does_not_execute(tmp_path):
    program = tmp_path / "vina"
    program.write_text("#!/bin/sh\nexit 99\n")
    program.chmod(0o755)
    rows = backend_catalog(probe=True, env={"PATH": str(tmp_path)})["backends"]
    vina = next(r for r in rows if r["id"] == "vina")
    assert vina["command_paths"]["vina"] == str(program)
    assert vina["execution_ready"] is False
    empty = backend_catalog(probe=True, env={})["backends"]
    assert all(not any(r["command_paths"].values()) for r in empty)


@pytest.mark.parametrize("name", ["docking", "medchem", "md"] + [b.id for b in BACKENDS if b.task == "prediction"])
def test_templates_have_valid_wiring_but_explicit_execution_blockers(name):
    plan = workflow_template(name)
    original = copy.deepcopy(plan)
    report = validate_plan(plan)
    assert report["valid"], report
    assert not report["execution_ready"]
    if name in {"docking", "medchem", "esmfold"}:
        assert all(step["integration"] == "adapter" for step in report["steps"])
        assert report["blockers"]  # Files and settings still need execution preflight.
    else:
        assert any("not implemented" in message for message in report["blockers"])
    assert plan == original


def test_a_raw_structure_cannot_be_silently_used_as_a_docking_receptor():
    plan = workflow_template("docking")
    plan["nodes"][4]["inputs"]["receptor"] = "inputs.structure"
    result = validate_plan(plan)
    assert not result["valid"]
    assert "expected receptor_pdbqt, received macromolecule" in result["errors"][0]


@pytest.mark.parametrize("mutation", ["cycle", "unknown", "duplicate", "wrong_input_dialect"])
def test_invalid_edges_are_actionable(mutation):
    plan = workflow_template("alphafold3")
    if mutation == "cycle":
        plan["nodes"][0]["inputs"]["request"] = "select.structure"
    elif mutation == "unknown":
        plan["nodes"][0]["backend"] = "imaginary_model"
    elif mutation == "duplicate":
        plan["nodes"][1]["id"] = "predict"
    else:
        plan["inputs"]["request"]["type"] = "boltz2_input"
    report = validate_plan(plan)
    assert not report["valid"] and report["errors"]


@pytest.mark.parametrize("value", [None, [], {}, {"schema_version": 1, "inputs": {"x": {"type": []}}, "nodes": [{}]}])
def test_malformed_plans_return_diagnostics(value):
    assert not validate_plan(value)["valid"]


def test_af3_input_keeps_typed_entities_and_explicit_msa_policy():
    result = build_prediction_input("alphafold3", request(
        {"id": "A", "kind": "protein", "sequence": "MKWV"},
        {"id": "B", "kind": "rna", "sequence": "ACGU"},
        {"id": "L", "kind": "ligand", "smiles": "C[NH3+]"}))
    native = json.loads(result["content"])
    assert native["dialect"] == "alphafold3"
    assert native["modelSeeds"] == [7]
    assert native["sequences"][0]["protein"]["unpairedMsa"] == ""
    assert native["sequences"][1]["rna"]["unpairedMsa"] == ""
    assert native["sequences"][2]["ligand"]["smiles"] == "C[NH3+]"
    assert result["execution_ready"] is False


def test_boltz_and_chai_get_their_own_dialects():
    boltz = build_prediction_input("boltz2", request())
    assert boltz["filename"].endswith(".yaml")
    assert json.loads(boltz["content"])["sequences"][0]["protein"]["msa"] == "empty"
    chai = build_prediction_input("chai1", request())
    assert chai["content"].startswith(">protein|name=A\n")


@pytest.mark.parametrize("backend,policy", [("alphafold2", "engine_pipeline"), ("alphafold_multimer", "engine_pipeline"),
                                          ("colabfold", "single_sequence"), ("esmfold", "single_sequence")])
def test_protein_exporters_write_fasta_and_record_runner_requirements(backend, policy):
    result = build_prediction_input(backend, request(policy=policy))
    assert result["content"] == ">A\nMKWVTFISLL\n"
    assert not result["execution_performed"]


@pytest.mark.parametrize("mutation", ["unknown", "duplicate", "newline", "seed", "wrong_kind", "missing_msa"])
def test_request_conversion_refuses_unsupported_or_ambiguous_input(mutation):
    value = request()
    if mutation == "unknown":
        value["covalent_bonds"] = []
    elif mutation == "duplicate":
        value["entities"] *= 2
    elif mutation == "newline":
        value["entities"][0]["sequence"] += "\n>injected"
    elif mutation == "seed":
        value["seeds"] = [True]
    elif mutation == "wrong_kind":
        value["entities"][0]["kind"] = "guess"
    else:
        value.pop("msa_policy")
    with pytest.raises(ValueError):
        build_prediction_input("alphafold3", value)


def test_multimer_is_not_silently_exported_as_a_batch_of_monomers():
    value = request({"id": "A", "kind": "protein", "sequence": "MKWV"},
                    {"id": "B", "kind": "protein", "sequence": "MKWV"})
    with pytest.raises(ValueError, match="multimer/batch"):
        build_prediction_input("colabfold", value)


def test_actual_registered_input_adapter_writes_manifest_and_preserves_source(tmp_path):
    from delfin.tools import run_step
    value = request()
    result = run_step("bio_prediction_input", backend="alphafold3", request=value, work_dir=tmp_path)
    assert result.ok
    assert result.artifacts["prediction_input"].is_file()
    manifest = json.loads(result.artifacts["json_report"].read_text())
    assert manifest["request"] == value
    assert len(manifest["input_sha256"]) == 64
    rerun = run_step("bio_prediction_input", backend="alphafold3", request=value, work_dir=tmp_path)
    assert not rerun.ok


def test_agent_handlers_return_the_same_contracts_and_blockers():
    from delfin.tools.mcp_server import h_bio_backends, h_bio_plan
    from delfin.tools.manifest import build_guide
    assert json.loads(h_bio_backends("prediction"))["backends"]
    result = json.loads(h_bio_plan(template="docking"))
    assert result["validation"]["valid"] and not result["validation"]["execution_ready"]
    assert "error" in json.loads(h_bio_plan())
    assert "bio_backends" in build_guide()["biomolecular_planning"]["tools"]
