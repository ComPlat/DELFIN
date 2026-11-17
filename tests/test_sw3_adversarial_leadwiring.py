"""Adversarial tests, round 5: every registered lead tool must be
reachable through the REAL dispatch path.

Reviewer nacht-s21, package SW3. The nine lead tool names are registered
in api_client._DELFIN_ONLY_TOOL_NAMES (operator build b93f39465), but a
registered name that no service method can dispatch is a promise the
agent cannot keep: the model sees the tool, the user confirms, and the
call dies with METHOD_NOT_FOUND (supervisor.py:920-921) -- or worse,
after confirmation, with NotImplementedError.

Real path (verified this session): a lead tool call maps to a
``swarm.*`` service method; protocol.SERVICE_REQUESTS must list it
(protocol.py:317-325) and supervisor.py:930-934 resolves it with
``getattr(swarm_rpc, method.replace(".", "_"))``, so a handler of that
name must exist in delfin/agent/service/swarm_rpc.py.
"""
from delfin.agent.swarm import lead as _lead
from delfin.agent.service import protocol
from delfin.agent.service import swarm_rpc

# The wire contract: registered tool name -> service method the
# supervisor would dispatch. Pinned explicitly because the mapping is
# not mechanical (swarm_triage_list -> swarm.triage.list).
_TOOL_TO_METHOD = {
    "swarm_status": "swarm.status",
    "swarm_plan": "swarm.plan",
    "swarm_message": "swarm.message",
    "swarm_triage_list": "swarm.triage.list",
    "swarm_start": "swarm.start",
    "swarm_pause": "swarm.pause",
    "swarm_resume": "swarm.resume",
    "swarm_triage_answer": "swarm.triage.answer",
    "swarm_prcheck": "swarm.prcheck",
}


def test_every_registered_lead_tool_is_dispatchable():
    """Each registered name must be listed in SERVICE_REQUESTS and have
    the handler the supervisor resolves by getattr."""
    for tool, method in _TOOL_TO_METHOD.items():
        assert method in protocol.SERVICE_REQUESTS, (
            f"{tool}: {method} missing from protocol.SERVICE_REQUESTS -- "
            "the registered lead tool dies with METHOD_NOT_FOUND "
            "(supervisor.py:920-921) before any dispatch")
        handler_name = method.replace(".", "_")
        assert hasattr(swarm_rpc, handler_name), (
            f"{tool}: {method} has no swarm_rpc.{handler_name} -- "
            "supervisor.py:933 getattr would raise AttributeError")


def test_prcheck_confirmed_call_is_a_structured_result_not_a_crash():
    """request_prcheck behind the user's confirmation must return a
    structured dict (the package's own refusal pattern, lead.py:105-116),
    never raise NotImplementedError."""
    lead = _lead
    try:
        result = _lead.request_prcheck(642, confirmed=True)
    except NotImplementedError as exc:  # pragma: no cover - the defect
        raise AssertionError(
            "request_prcheck raised NotImplementedError after "
            f"confirmation: {exc}") from exc
    assert isinstance(result, dict), result
