"""WB2 adversarial review (nacht-s4) — attacks on user_workflows (a8fba967).

Each test names the candidate defect it probes. A red test here is a
FINDING to nacht-s3; a green test is committed evidence the code holds.
"""

from __future__ import annotations

import json

import pytest

from delfin.tools import user_workflows as uw

GOOD_SPEC = {
    "name": "my_opt_freq",
    "steps": [
        {"step": "smiles_to_xyz", "kwargs": {"smiles": "CCO"},
         "label": "smiles"},
    ],
}


def test_name_with_trailing_newline_is_rejected(tmp_path):
    """A1: _NAME_RE uses .match, and '$' also matches before a final
    newline — so 'my_flow\\n' passes as a name. The file name derived
    from it then contains a newline. skill_share._check_name guards
    exactly this with fullmatch (delfin/agent/skill_share.py:56-57)."""
    with pytest.raises(ValueError):
        uw._dir(None, "my_flow\n")  # must reject, not accept


def test_list_workflows_survives_a_corrupt_record(tmp_path):
    """A2: list_workflows reads every dir's newest vNNN.json and json.loads
    it (user_workflows.py:157). One corrupt file kills the WHOLE listing
    instead of skipping the broken workflow."""
    uw.save(GOOD_SPEC, name="good_flow", root=tmp_path)
    bad = tmp_path / "bad_flow"
    bad.mkdir()
    (bad / "v001.json").write_text("{not json", encoding="utf-8")
    names = [w["name"] for w in uw.list_workflows(root=tmp_path)]
    assert "good_flow" in names
    assert "bad_flow" not in names


def test_check_spec_fails_on_a_lint_error(tmp_path):
    """A3 (negative control for A4): check_spec must refuse a spec whose
    lint has an error level finding. Shape: a spec dict fed through the
    real platform.resolve_spec."""
    bad = dict(GOOD_SPEC, steps=[
        {"step": "orca_opt", "kwargs": {"method": "CCSD(T)",
                                        "basis": "def2-SVP"},
         "label": "x"},
    ])
    result = uw.check_spec(bad, geometry=True)
    assert result["ok"] is False, result["errors"]


def test_check_spec_order_chain_present():
    """A4: the mandatory chain resolve -> validate -> lint -> dry run
    must be visible in check_spec's result and must NOT execute anything.
    R1 note: pipeline.dry_run (pipeline.py:1648) is an alias of
    validate (pipeline.py:1629) — the builder must say what the dry-run
    stage adds beyond re-validating."""
    result = uw.check_spec(GOOD_SPEC, geometry=True)
    assert result["ok"] is True
    assert "resolved" in result and "report" in result and "lint" in result
    assert "dry_run" in result


def test_extend_template_rejects_bad_new_name(tmp_path):
    """A5: extend_template with a new_name that is path-unsafe must raise,
    never write a file with an unsafe name."""
    with pytest.raises(ValueError):
        uw.extend_template("classic_opt_freq", new_name="../escape")
    with pytest.raises(ValueError):
        uw.extend_template("classic_opt_freq", new_name="x\n")
