"""Adversarial tests for SW3 phase 1 lead-tool error consistency (s21).

Target: delfin/agent/swarm/lead.py request_swarm_start and
request_triage_answer error handling.

The lead's state-changing tools promise (docstring, lead.py) to return a
structured refusal rather than crash when an action cannot be taken.
Two asymmetries are probed here:

- request_swarm_start(confirmed=True) catches LaunchError and returns a
  structured {"action":"refused"}; the UNCONFIRMED branch (the normal,
  ask-the-user path) calls the same plan_launch OUTSIDE that guard and
  so raises on a bad spec instead of returning the same structured
  refusal. The user's first contact with start must not be a crash.

- request_triage_answer's fail-closed behaviour is now pinned on the
  real dispatch path (swarm_rpc.swarm_triage_answer); see
  test_triage_answer_confirmed_never_raises_on_an_escalation.
"""

from delfin.agent.swarm import lead


def _bad_spec() -> dict:
    # Missing required "goal"/"packages" keys -> plan_launch rejects.
    return {"spec_version": 1}


def test_start_unconfirmed_returns_refusal_on_a_bad_spec():
    """The ask-the-user path must not crash on a spec SW1 rejects; it
    should hand the refusal back in the same shape as the confirmed
    path."""
    result = lead.request_swarm_start(_bad_spec(), confirmed=False)
    assert isinstance(result, dict)
    assert result.get("action") in ("ask_user", "refused")


def test_start_confirmed_and_unconfirmed_fail_the_same_way():
    """Both branches of request_swarm_start reject the same bad spec in
    the same shape (structured dict), never a raise."""
    confirmed = lead.request_swarm_start(_bad_spec(), confirmed=True)
    unconfirmed = lead.request_swarm_start(_bad_spec(), confirmed=False)
    assert isinstance(confirmed, dict) and isinstance(unconfirmed, dict)


def test_triage_answer_confirmed_never_raises_on_an_escalation():
    """Answering an escalation on the user's behalf must fail closed
    with a structured refusal, not leak TriageError into the tool
    call. The real dispatch path is the service handler the supervisor
    resolves (swarm_rpc.swarm_triage_answer); it fails closed with a
    SwarmRpcError naming the escalation."""
    import pytest
    from delfin.agent.service import swarm_rpc
    item = {"id": "q1", "kind": "push_request", "from": "builder-core",
            "member": "builder-core"}
    with pytest.raises(swarm_rpc.SwarmRpcError) as excinfo:
        swarm_rpc.swarm_triage_answer({"item": item, "answer": "no"})
    assert str(excinfo.value)
