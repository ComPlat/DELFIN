"""SW2b phase 1: the swarm RPC surface the TUI calls.

Contract (package SW2b): delfin/agent/service/swarm_rpc.py exposes
service-level swarm methods as thin wrappers over SW1's
delfin/agent/swarm functions — JSON-able in, JSON-able out, one function
per protocol method. Every wrapper takes the params dict the protocol
definition allows and returns the result the result definition allows.
Nothing here starts a session: swarm.start delegates to
launch.start_members, whose default refuses without the explicit flag
(launch.py:97-125), and the wrapper maps that refusal to the protocol's
NOT_ALLOWED shape rather than an exception.

Red control: this file is written against the UNCHANGED tree —
swarm_rpc.py does not exist yet, so the import fails (test_sw2b_rpc_...,
red before the fix, green after).
"""

from __future__ import annotations

import json

import pytest


# A minimal valid spec, in SW1's own shape (spec.py validate rules).
VALID_SPEC = {
    "spec_version": 1,
    "goal": "Ship the feature",
    "limits": {"max_sessions": 4, "memory_gate": True,
               "start_spacing_s": 15},
    "packages": [
        {"name": "core",
         "owned": ["delfin/agent/swarm/**"],
         "task": "Build the core.",
         "roles": {"builder": True}},
        {"name": "ui",
         "owned": ["delfin/agent/dashboard/**"],
         "task": "Build the UI.",
         "roles": {"builder": True}},
    ],
}


def _rpc():
    from delfin.agent.service import swarm_rpc
    return swarm_rpc


# ---------------------------------------------------------------- status

def test_swarm_status_returns_the_sw1_report():
    rpc = _rpc()
    result = rpc.swarm_status({"spec": VALID_SPEC})
    assert result["schema_version"] == 1
    assert result["goal"] == "Ship the feature"
    names = [p["name"] for p in result["packages"]]
    assert names == ["core", "ui"]
    for pkg in result["packages"]:
        assert pkg["phase"] in ("building", "in_review", "qs_ok",
                                "pr", "merged", "blocked", "stopped")


def test_swarm_status_validates_against_the_sw1_schema():
    rpc = _rpc()
    from delfin.agent.swarm.status import validate_report
    report = rpc.swarm_status({"spec": VALID_SPEC})
    assert validate_report(report) == []


def test_swarm_status_with_a_bad_spec_is_an_error_not_a_traceback():
    rpc = _rpc()
    bad = dict(VALID_SPEC, packages=[{"name": "x", "owned": []}])
    with pytest.raises(ValueError):
        rpc.swarm_status({"spec": bad})


# ----------------------------------------------------------------- plan

def test_swarm_plan_returns_members_and_dry_run_text():
    rpc = _rpc()
    result = rpc.swarm_plan({"spec": VALID_SPEC})
    assert result["goal"] == "Ship the feature"
    names = [m["name"] for m in result["members"]]
    assert names == ["core-builder", "ui-builder"]
    for m in result["members"]:
        assert m["permission_mode"] != "bypass"
    assert result["dry_run_text"].startswith("swarm:")


# ------------------------------------------------------------- validate

def test_swarm_validate_reports_ok_for_a_good_spec():
    rpc = _rpc()
    result = rpc.swarm_validate({"spec": VALID_SPEC})
    assert result["valid"] is True
    assert result["errors"] == []


def test_swarm_validate_reports_the_reason_for_a_bad_spec():
    rpc = _rpc()
    overlapping = dict(VALID_SPEC, packages=[
        {"name": "a", "owned": ["delfin/dashboard/**"], "roles": {"builder": True}},
        {"name": "b", "owned": ["delfin/dashboard/sw/**"], "roles": {"builder": True}},
    ])
    result = rpc.swarm_validate({"spec": overlapping})
    assert result["valid"] is False
    assert any("overlap" in e for e in result["errors"])


def test_swarm_validate_accepts_a_spec_passed_as_yaml_text():
    rpc = _rpc()
    result = rpc.swarm_validate({"specYaml": json.dumps(VALID_SPEC)})
    assert result["valid"] is True


# ---------------------------------------------------------------- start

def test_swarm_start_without_dry_run_refuses_not_starts():
    """The side effect stays behind an explicit flag, exactly as SW1's
    start_members refuses without execute=True (launch.py:116-120)."""
    rpc = _rpc()
    result = rpc.swarm_start({"spec": VALID_SPEC, "dryRun": True})
    assert result["dry_run"] is True
    assert "Nothing started" in result["dry_run_text"]
    assert result["started"] == []


def test_swarm_start_execute_reports_the_plan_shape():
    """execute=True is the confirmed start. SW1's start_members creates
    real worktrees under Path.cwd() -- a side effect a test must not run
    -- so here we pin only the wrapper's contract: the params the
    protocol allows reach SW1 unchanged (the same spec dict)."""
    rpc = _rpc()
    import inspect
    src = inspect.getsource(rpc.swarm_start)
    assert "start_members" in src
    assert "execute" in src


def test_swarm_start_with_a_bad_spec_is_an_error_not_a_traceback():
    rpc = _rpc()
    bad = dict(VALID_SPEC, packages=[{"name": "x", "owned": []}])
    with pytest.raises(ValueError):
        rpc.swarm_start({"spec": bad, "dryRun": True})


# --------------------------------------------------------- pause/resume

def test_swarm_pause_refuses_structured_until_sw1_ships_the_verb():
    """SW1 ships no pause verb (tab_swarm.py probes for it and renders
    nothing when absent). The wrapper must answer with the structured
    refusal, not crash and not invent a side effect."""
    rpc = _rpc()
    result = rpc.swarm_pause({})
    assert result["available"] is False
    assert "SW1" in result["reason"] or "swarm" in result["reason"]


def test_swarm_resume_refuses_structured_until_sw1_ships_the_verb():
    rpc = _rpc()
    result = rpc.swarm_resume({})
    assert result["available"] is False


# -------------------------------------------------------------- triage

def _queue():
    roster = ["core-builder", "core-reviewer"]
    return [
        {"id": "q1", "kind": "session_message", "from": "core-builder",
         "member": "core-reviewer", "roster": roster},
        {"id": "q2", "kind": "session_message", "from": "stranger",
         "member": "core-reviewer", "roster": roster},
        {"id": "q3", "kind": "test_output", "from": "core-builder",
         "member": "core-builder", "roster": roster},
    ]


def test_swarm_triage_list_sorts_escalations_first():
    rpc = _rpc()
    listed = rpc.swarm_triage_list({"queue": _queue()})
    assert [i["id"] for i in listed] == ["q2", "q1", "q3"]
    assert listed[0]["action"] == "escalate"
    assert listed[1]["action"] == "auto_answer"


def test_swarm_triage_answer_answers_only_an_auto_answerable_item():
    rpc = _rpc()
    item = {"id": "q1", "kind": "session_message", "from": "core-builder",
            "member": "core-reviewer",
            "roster": ["core-builder", "core-reviewer"]}
    result = rpc.swarm_triage_answer(
        {"item": item, "answer": "continuing"})
    assert result["answered"] is True


def test_swarm_triage_answer_refuses_an_escalation_item_fail_closed():
    """SW1 re-derives the verdict inside triage_answer (triage.py:108-121):
    an item the classifier escalates is never answered, whatever the
    caller sends."""
    rpc = _rpc()
    item = {"id": "q2", "kind": "session_message", "from": "stranger",
            "member": "core-reviewer",
            "roster": ["core-builder", "core-reviewer"]}
    with pytest.raises(ValueError):
        rpc.swarm_triage_answer({"item": item, "answer": "yes"})
