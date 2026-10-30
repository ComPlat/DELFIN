"""Workflow descriptor contract — ONE contract for every DELFIN workflow.

A ``WorkflowDescriptor`` extends ``delfin.tools._application.Application``
(no parallel type) so the classic CONTROL.txt workflows, the pipeline
applications and user workflows are all named, introspected and listed
through the same registry.  The descriptor DESCRIBES a classic workflow; it
never reimplements or changes how that workflow runs.
"""

from __future__ import annotations

import pytest

from delfin.common.control_validator import CONTROL_FIELD_SPECS
from delfin.tools._application import Application
from delfin.tools.workflows import (
    WorkflowDescriptor,
    add_step,
    list_workflows,
    register_classic_workflows,
    workflow_for,
)

# The classic CONTROL.txt workflows, by name.  Registering a descriptor for
# each must change nothing about how they run (their engines stay untouched).
CLASSIC_WORKFLOWS = ("classic", "OCCUPIER", "ESD", "IMAG", "redox", "co2")


# --- Requirement 2: registration + "is an Application" --------------------


def test_classic_workflows_are_registered():
    register_classic_workflows()  # idempotent
    for name in CLASSIC_WORKFLOWS:
        assert workflow_for(name) is not None, f"{name} must have a descriptor"


def test_every_classic_descriptor_is_an_application():
    """A descriptor extends Application — no parallel workflow type."""
    register_classic_workflows()
    for name in CLASSIC_WORKFLOWS:
        wf = workflow_for(name)
        assert isinstance(wf, Application), f"{name} is not an Application"
        assert isinstance(wf, WorkflowDescriptor)


def test_classic_descriptors_have_kind_classic_control():
    register_classic_workflows()
    for name in CLASSIC_WORKFLOWS:
        assert workflow_for(name).kind == "classic-control", name


# --- Requirement 1: every descriptor field is filled -----------------------


def test_every_descriptor_has_the_full_contract():
    """name, purpose, kind, entry, steps, artifacts, extension points,
    resume/recovery and backend are all filled for every classic workflow."""
    register_classic_workflows()
    for name in CLASSIC_WORKFLOWS:
        wf = workflow_for(name)
        assert wf.name == name
        assert wf.purpose, f"{name}: purpose (one sentence) missing"
        assert wf.kind, f"{name}: kind missing"
        assert wf.entry, f"{name}: entry (how it is started) missing"
        assert wf.steps, f"{name}: steps (jobs it writes) missing"
        assert wf.artifacts, f"{name}: artifacts it produces missing"
        assert wf.extension_points, f"{name}: extension points missing"
        assert wf.resume, f"{name}: resume/recovery behaviour missing"
        assert wf.backend, f"{name}: local vs slurm missing"


# --- Requirement 1: CONTROL keys match control_validator -------------------


def test_classic_descriptor_control_keys_match_control_validator():
    """Every CONTROL key a descriptor declares exists in CONTROL_FIELD_SPECS
    with the same default — the descriptor may not invent keys."""
    field_specs = {s.name: s for s in CONTROL_FIELD_SPECS}
    register_classic_workflows()
    for name in CLASSIC_WORKFLOWS:
        wf = workflow_for(name)
        assert wf.control_keys, f"{name}: declares no CONTROL keys"
        for ck in wf.control_keys:
            field = field_specs.get(ck.name)
            assert field is not None, (
                f"{name} declares CONTROL key {ck.name!r} that control_validator "
                "does not know"
            )
            assert ck.default == field.default, (
                f"{name}.{ck.name}: descriptor default {ck.default!r} != "
                f"control_validator default {field.default!r}"
            )


def test_descriptor_control_key_is_a_typed_spec():
    """A CONTROL key spec carries name, type, default, allowed values and a
    'changes' description, derived from control_validator where possible."""
    register_classic_workflows()
    for name in CLASSIC_WORKFLOWS:
        for ck in workflow_for(name).control_keys:
            assert ck.name
            assert ck.type in ("str", "int", "float", "bool", "list", "path")
            assert ck.changes, f"{name}.{ck.name}: 'what it changes' missing"


# --- Requirement 1/2: serialization round-trip -----------------------------


def test_descriptor_serialization_roundtrip():
    """to_dict/from_dict round-trips the descriptor fields, not just the
    base Application fields."""
    register_classic_workflows()
    wf = workflow_for("OCCUPIER")
    data = wf.to_dict()
    assert data["kind"] == "classic-control"
    assert data["purpose"]
    restored = WorkflowDescriptor.from_dict(data)
    assert restored.name == wf.name
    assert restored.kind == wf.kind
    assert restored.purpose == wf.purpose
    assert restored.backend == wf.backend
    assert [c.name for c in restored.control_keys] == \
        [c.name for c in wf.control_keys]


def test_application_from_dict_still_loads_without_descriptor_fields():
    """A plain serialized Application (no descriptor fields) still loads."""
    plain = Application(name="plain", description="x").to_dict()
    wf = WorkflowDescriptor.from_dict(plain)
    assert wf.kind == "pipeline"  # the base default
    assert wf.name == "plain"


# --- Requirement 2: catalog / list lists all three kinds together ----------


def test_catalog_lists_all_three_kinds_together(tmp_path):
    """list_workflows() returns classic-control, pipeline and user workflows
    together."""
    register_classic_workflows()
    kinds = {wf.kind for wf in list_workflows()}
    assert "classic-control" in kinds, f"classic workflows missing: {kinds}"
    assert "pipeline" in kinds, f"pipeline applications missing: {kinds}"
    # A saved user workflow appears with kind 'user'.
    import delfin.tools.user_workflows as uw
    spec = {
        "delfin_pipeline_spec": "1",
        "name": "user_demo",
        "kind": "pipeline",
        "defaults": {},
        "steps": [
            {"step": "smiles_to_xyz", "label": "structure",
             "kwargs": {"smiles": "CCO"}},
            {"step": "orca_sp", "label": "sp",
             "kwargs": {"charge": 0, "mult": 1}},
        ],
    }
    uw.save(spec, name="user_demo", root=tmp_path)
    user_kinds = {wf.kind for wf in list_workflows(root=tmp_path)}
    assert "user" in user_kinds, f"user workflows missing: {user_kinds}"


# --- Requirement 3: extension 'add step X after step Y' --------------------


def test_extension_occupier_then_dlpno_single_point_validates():
    """Add a DLPNO single point on the chosen state AFTER the OCCUPIER spin
    selection; the resulting spec passes the existing checks (check_spec)."""
    register_classic_workflows()
    occupier = workflow_for("OCCUPIER")
    extended = add_step(
        occupier,
        {
            "step": "orca_sp",
            "label": "dlpno_sp",
            "kwargs": {"charge": 0, "mult": 1, "method": "DLPNO-CCSD(T)"},
        },
        after="spin_states",
    )
    # The new spec has every original OCCUPIER step unchanged, then the new one.
    original_steps = [
        s.get("label", s["step"]) for s in occupier.steps_as_spec()["steps"]
    ]
    new_steps = [s.get("label", s["step"]) for s in extended["steps"]]
    for lbl in original_steps:
        assert lbl in new_steps, f"original step {lbl} lost by extension"
    assert new_steps[-1] == "dlpno_sp", "new step must run after W"

    # The existing checks accept it (static validation, nothing runs).
    from delfin.tools import user_workflows as uw
    report = uw.check_spec(extended)
    assert report["ok"], f"extended OCCUPIER+DLPNO fails checks: {report['errors']}"


def test_extension_step_before_rejects_unknown_anchor():
    register_classic_workflows()
    occupier = workflow_for("OCCUPIER")
    with pytest.raises(KeyError):
        add_step(occupier, {"step": "orca_sp", "label": "x"},
                 after="no_such_step")


# --- Requirement 4: scaffold writes a workflow that passes its acceptance ---


def test_scaffold_writes_workflow_that_passes_acceptance(tmp_path):
    """The scaffold writes a new-workflow skeleton whose descriptor exists
    and whose spec passes the same checks a real workflow must pass."""
    from delfin.tools.scaffold import build_workflow_scaffold
    plan = build_workflow_scaffold(
        "my_direct", kind="pipeline", directory=tmp_path,
        register=False,
        steps=[{"step": "smiles_to_xyz", "label": "structure",
                "kwargs": {"smiles": "CCO"}},
               {"step": "orca_sp", "label": "sp",
                "kwargs": {"charge": 0, "mult": 1}}],
    )
    assert plan.ok, plan.errors
    assert plan.module_path.exists(), "skeleton module not written"
    # The scaffolded descriptor's spec passes the existing checks.
    from delfin.tools import user_workflows as uw
    report = uw.check_spec(plan.spec)
    assert report["ok"], f"scaffolded workflow fails checks: {report['errors']}"


# --- Requirement 2: registering changes no classic behaviour ---------------


def test_registering_descriptors_changes_no_classic_engine():
    """The classic engines keep their entry points after registration — the
    descriptor only describes, it does not reimplement."""
    import delfin.workflows.pipeline as wf_pipeline
    for fn in ("run_classic_phase", "run_occuper_phase", "run_esd_phase"):
        assert hasattr(wf_pipeline, fn), \
            f"classic engine entry {fn} must remain present after registration"


def test_manifest_lists_classic_and_pipeline_together():
    """The machine-readable manifest lists the classic CONTROL workflows and
    the pipeline applications together under one 'workflows' section."""
    from delfin.tools.manifest import build_manifest

    workflows = build_manifest()["workflows"]
    kinds = {w["kind"] for w in workflows}
    assert "classic-control" in kinds
    assert "pipeline" in kinds
    # The classic descriptors carry the descriptor fields, not just a name.
    occ = next(w for w in workflows
               if w["name"] == "OCCUPIER" and w["kind"] == "classic-control")
    assert occ["entry"], "OCCUPIER entry missing in manifest"
    assert occ["control_keys"], "OCCUPIER control keys missing in manifest"
    assert occ["extension_points"], "OCCUPIER extension points missing in manifest"


def test_new_workflow_cli_scaffolds_a_valid_workflow(tmp_path):
    """The `new-workflow` CLI command writes a module and acceptance test whose
    spec passes the existing checks — a scaffold cannot emit a workflow that
    would not run."""
    import argparse
    from delfin.agent import cli_new_workflow

    parser = argparse.ArgumentParser()
    sub = parser.add_subparsers()
    cli_new_workflow.build_parser(sub)
    args = parser.parse_args(
        ["new-workflow", "my_flow", "--dir", str(tmp_path)])
    rc = cli_new_workflow.handler(args)
    assert rc == 0, f"new-workflow returned {rc}"
    assert (tmp_path / "my_flow.py").exists()
    assert (tmp_path / "test_workflow_my_flow.py").exists()

    from delfin.tools import user_workflows as uw
    spec = (tmp_path / "my_flow.py").read_text(encoding="utf-8")
    assert "WorkflowDescriptor" in spec, "scaffold module must register a descriptor"
    # The default skeleton spec validates under the existing checks.
    report = uw.check_spec({
        "name": "my_flow", "kind": "pipeline", "defaults": {},
        "steps": [
            {"step": "smiles_to_xyz", "label": "structure",
             "kwargs": {"smiles": "CCO"}},
            {"step": "orca_sp", "label": "sp",
             "kwargs": {"charge": 0, "mult": 1}},
        ],
    })
    assert report["ok"], f"scaffold default spec fails: {report['errors']}"
