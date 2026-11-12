"""WB2 phase 3 — workflow_builder: request text -> checked draft spec.

Red control: every case failed on the previous commit with ImportError
(delfin.agent.workflow_builder did not exist).
"""

from __future__ import annotations

import pytest

from delfin.agent import workflow_builder as wb


REQUEST = "Run an opt+freq workflow on ethanol (SMILES CCO) at B3LYP/def2-SVP."


def test_request_yields_a_spec_that_passes_all_checks():
    result = wb.build_workflow(REQUEST)
    assert result["spec"]["steps"], "a request must produce steps"
    steps = [s["step"] for s in result["spec"]["steps"]]
    assert steps[:4] == ["smiles_to_xyz", "xtb_opt", "orca_opt", "orca_freq"]
    # every check ran, in order, and passed
    assert result["checks"]["ok"], result["checks"]["errors"]
    order = result["checks"]["order"]
    assert order == ["resolve", "validate", "lint", "dry_run", "cost"]
    # every step carries a one-line reason
    for entry in result["steps_reasons"]:
        assert entry["reason"].strip()
    # cost estimate per step came from DELFIN's own advise()
    assert result["cost"], "a cost estimate is mandatory before submit"
    assert all("PAL" in c for c in result["cost"])


def test_a_spec_failing_the_lint_is_not_submitted():
    bad = wb.build_workflow(
        REQUEST,
        base_spec={
            "name": "bad",
            "steps": [
                {"step": "orca_sp", "label": "orcatcharge",
                 "kwargs": {"charge": 0, "mult": -1,
                            "method": "B3LYP", "basis": "def2-SVP"}},
            ],
        },
    )
    assert not bad["checks"]["ok"]
    assert any("invalid-multiplicity" in e for e in bad["checks"]["errors"])
    # the fix is named in the returned material, never applied silently
    assert any("multiplicity" in f["message"] for f in bad["checks"]["lint"])
    assert bad.get("submitted") is not True


def test_extend_template_path_returns_new_spec_original_untouched():
    import delfin.tools.templates as templates

    before = templates.classic_opt_freq.to_dict()
    result = wb.build_workflow(
        "extend classic_opt_freq with an orca_sp single point",
        base_template="classic_opt_freq",
        edits={"add": [{"step": "orca_sp", "kwargs": {"charge": 0}}]},
    )
    assert templates.classic_opt_freq.to_dict() == before
    steps = [s["step"] for s in result["spec"]["steps"]]
    assert steps[-1] == "orca_sp"
    assert result["checks"]["ok"]


def test_draft_from_the_catalog_is_explicit():
    result = wb.build_workflow(REQUEST, base_template="classic_opt_freq")
    assert result["source"] == "template:classic_opt_freq"
    # extend/modify of a saved user workflow names it as source too
    saved = wb.build_workflow(REQUEST, base_spec={
        "name": "from_saved", "steps": [
            {"step": "xtb_opt", "kwargs": {"charge": 0, "method": "XTB2"}},
        ]})
    assert saved["source"] == "provided spec"
