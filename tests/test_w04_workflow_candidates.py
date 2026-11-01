"""wave14/w04 RULE 2: existing workflows first.

``workflow_candidates(task_text, available)`` ranks the platform's
workflows against a chemist's request, deterministically (no model call).
Every candidate carries the reason it matched and the CONTROL keys to set.
When nothing fits, it says so, so the agent knows to extend or compose.
"""
from __future__ import annotations

import pytest

from delfin.tools.workflow_candidates import workflow_candidates


def _names(cands):
    return [c["workflow"] for c in cands]


# 12 chemist-worded requests -> the workflow that must win.
EXPECTED = [
    ("redox potential of this Fe complex", "redox"),
    ("compute the oxidation potential for the Cu pair", "redox"),
    ("which spin state is lowest for my iron complex", "OCCUPIER"),
    ("compare multiplicities and pick the ground state", "OCCUPIER"),
    ("excited states and emission of this dye", "ESD"),
    ("phosphorescence rates and ISC for the Pt complex", "ESD"),
    ("CO2 binding to my catalyst", "co2"),
    ("adduct formation energies for CO2 on the metal", "co2"),
    ("rotational barrier of butane", "classic"),
    ("TS for this step and its imaginary frequency", "classic"),
    ("I have an imaginary frequency after optimisation", "IMAG"),
    ("drive the saddle point down to a minimum", "IMAG"),
]


@pytest.mark.parametrize("text,expected", EXPECTED,
                         ids=[t[0][:38] for t in EXPECTED])
def test_expected_workflow_wins(text, expected):
    out = workflow_candidates(text)
    assert out, "always returns candidates or an explicit 'no fit'"
    assert _names(out)[0] == expected, (text, _names(out))


def test_every_candidate_carries_reason_and_control_keys():
    out = workflow_candidates("redox potential of this Fe complex")
    top = out[0]
    assert top["reason"]
    keys = top["control_keys"]
    assert isinstance(keys, list) and keys
    for key in keys:
        assert "name" in key and "type" in key
    # the keys are the workflow's validated CONTROL keys
    from delfin.tools.workflows import workflow_for
    valid = {c.name for c in workflow_for("redox").control_keys}
    assert {k["name"] for k in keys} <= valid


def test_ranking_is_deterministic():
    a = workflow_candidates("redox potential of this Fe complex")
    b = workflow_candidates("redox potential of this Fe complex")
    assert a == b


def test_no_fit_is_said_plainly():
    out = workflow_candidates("please reorder my pdf pages")
    assert not _names(out) or out[-1].get("workflow") == ""
    assert any(c.get("no_fit") for c in out) or out == []


def test_available_narrows_the_pool():
    out = workflow_candidates("redox potential of this Fe complex",
                              available=["ESD"])
    assert "redox" not in _names(out)
    assert out == [] or out[-1].get("no_fit")


def test_candidate_lists_the_steps_it_would_run():
    out = workflow_candidates("which spin state is lowest", available=None)
    top = next(c for c in out if c["workflow"] == "OCCUPIER")
    assert top.get("steps"), "the model must see what the workflow runs"


def test_tool_exposed_in_api_catalog():
    from delfin import api

    names = {e["name"] for e in api._TOOL_CATALOG}
    assert "workflow_candidates" in names
    assert "workflow_candidates" in api.__all__
    out = api.workflow_candidates("redox potential of this Fe complex")
    assert out[0]["workflow"] == "redox"


def test_schema_stays_short():
    from delfin.agent.api_client import _DOC_TOOLS_OPENAI, estimate_schema_tokens

    schema = next(t["function"] for t in _DOC_TOOLS_OPENAI
                  if t["function"]["name"] == "workflow_candidates")
    tokens = estimate_schema_tokens({"type": "function", "function": schema})
    assert tokens <= 130, f"workflow_candidates schema costs {tokens} tokens"
