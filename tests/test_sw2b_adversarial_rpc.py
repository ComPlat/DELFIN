"""Adversarial tests for the SW2b service surface (reviewer nacht-s14).

These go through the PUBLIC paths the product uses:
- protocol.params_problem (the service-level params gate, protocol.py),
- delfin/agent/service/swarm_rpc.py (the wrappers the supervisor hook
  dispatches to).

No test touches the live service process, the real home, or git state.
"""

from __future__ import annotations

import pytest

from delfin.agent.service import protocol
from delfin.agent.service import swarm_rpc


VALID_SPEC = {
    "goal": "adversarial review spec",
    "packages": [
        {"name": "pkg-a", "owned": ["tests/test_pkg_a_*.py"],
         "roles": {"builder": True}},
    ],
}


def _params_problem(method, params):
    return protocol.params_problem(
        protocol.SERVICE_REQUESTS[method], params)


# ---------------------------------------------------------------------------
# The confirm gate: a non-boolean confirm must never reach the wrapper.
# ---------------------------------------------------------------------------

def test_start_rejects_a_non_boolean_confirm_at_the_protocol_gate():
    # 'confirm' is declared "boolean?" in SERVICE_REQUESTS; a string
    # ("0", "yes", "false") or an integer must be refused before any
    # wrapper runs -- bool("0") is True and would flip the execute gate.
    for bad in ["0", "yes", "false", 1, 0]:
        reason = _params_problem("swarm.start",
                                 {"spec": VALID_SPEC, "confirm": bad})
        assert reason, f"confirm={bad!r} was accepted at the gate"


def test_start_rejects_a_non_boolean_dry_run_at_the_protocol_gate():
    reason = _params_problem("swarm.start",
                             {"spec": VALID_SPEC, "dryRun": "false"})
    assert reason, "dryRun='false' was accepted at the gate"


def test_start_wrapper_is_conservative_on_non_boolean_confirm():
    # Defense in depth: even called directly, a truthy non-bool confirm
    # must NOT execute -- the wrapper treats only a real True as confirmed.
    result = swarm_rpc.swarm_start({"spec": VALID_SPEC, "confirm": "yes"})
    assert result["dry_run"] is True
    assert result["started"] == []
    assert result["refused"]


def test_start_without_confirm_returns_a_dry_run_with_a_refusal_reason():
    result = swarm_rpc.swarm_start({"spec": VALID_SPEC, "dryRun": False})
    assert result["dry_run"] is True
    assert result["started"] == []
    assert result["refused"], "no refusal reason for an unconfirmed execute"


# ---------------------------------------------------------------------------
# _spec_of: the spec carrier must not fall open to strange payloads.
# ---------------------------------------------------------------------------

def test_spec_of_refuses_a_missing_spec():
    with pytest.raises(swarm_rpc.SwarmRpcError) as exc:
        swarm_rpc.swarm_validate({})
    assert "missing spec" in str(exc.value)


def test_spec_of_refuses_a_non_string_spec_yaml():
    with pytest.raises(swarm_rpc.SwarmRpcError):
        swarm_rpc.swarm_validate({"specYaml": 42})


def test_spec_of_refuses_yaml_that_parses_to_a_non_mapping():
    # YAML "yes" parses to a string; YAML "[1, 2]" to a list. Neither is
    # a spec object -- it must be refused, not defaulted.
    for text in ("yes", "[1, 2]", "42", "null"):
        with pytest.raises(swarm_rpc.SwarmRpcError):
            swarm_rpc.swarm_validate({"specYaml": text})


def test_spec_of_accepts_a_real_spec_yaml_round_trip():
    text = "goal: adversarial\npackages:\n  - name: pkg-a\n    owned:\n      - tests/test_pkg_a_*.py\n    roles:\n      builder: true\n"
    assert swarm_rpc.swarm_validate({"specYaml": text}) == {
        "valid": True, "errors": []}


# ---------------------------------------------------------------------------
# triage.answer: fail-closed on escalations; roster self-claim is the
# adversarial edge -- an item that brings its own roster claiming
# membership must not widen the allow-list.
# ---------------------------------------------------------------------------

def _escalation_item(**extra):
    item = {
        "id": "it-1", "kind": "permission_request",
        "from": "nacht-s13", "member": "nacht-s14",
        "roster": ["nacht-s13", "nacht-s14"],
        "text": "push my branch please",
    }
    item.update(extra)
    return item


def test_triage_answer_refuses_an_escalation_item_fail_closed():
    with pytest.raises(swarm_rpc.SwarmRpcError) as exc:
        swarm_rpc.swarm_triage_answer(
            {"item": _escalation_item(), "answer": "ok"})
    assert "escalate" in str(exc.value) or "user" in str(exc.value)


def test_triage_answer_ignores_a_crafted_verdict_in_the_item():
    # A caller smuggling "action": "auto_answer" INTO the item dict must
    # not widen the decision -- the verdict is re-derived from the item's
    # kind/from/member, not read from the item's fields.
    item = _escalation_item(action="auto_answer",
                            reason="plain partner message")
    with pytest.raises(swarm_rpc.SwarmRpcError):
        swarm_rpc.swarm_triage_answer({"item": item, "answer": "ok"})


def test_triage_list_of_an_empty_queue_returns_an_empty_list():
    assert swarm_rpc.swarm_triage_list({"queue": []}) == []


def test_triage_list_refuses_a_non_list_queue():
    with pytest.raises(swarm_rpc.SwarmRpcError):
        swarm_rpc.swarm_triage_list({"queue": {"0": {}}})


# ---------------------------------------------------------------------------
# pause/resume: the structured refusal, never an invented side effect.
# ---------------------------------------------------------------------------

def test_pause_and_resume_answer_available_false_with_a_reason():
    for fn in (swarm_rpc.swarm_pause, swarm_rpc.swarm_resume):
        result = fn({})
        assert result["available"] is False
        assert result["reason"]


# ---------------------------------------------------------------------------
# The service result shapes must match SERVICE_RESULTS so the TUI's
# generated types cannot drift from what the wrappers return.
# ---------------------------------------------------------------------------

def test_every_swarm_method_has_service_request_and_result_entries():
    surface = ["swarm.status", "swarm.plan", "swarm.validate",
               "swarm.start", "swarm.pause", "swarm.resume",
               "swarm.triage.list", "swarm.triage.answer"]
    schema = protocol.json_schema()
    for method in surface:
        assert method in protocol.SERVICE_REQUESTS, method
        assert method in protocol.SERVICE_RESULTS, method
        assert method in schema["service"]["requests"], method
        assert method in schema["service"]["results"], method
        fn = "swarm_" + method.replace("swarm.", "").replace(".", "_")
        assert callable(getattr(swarm_rpc, fn, None)), method

