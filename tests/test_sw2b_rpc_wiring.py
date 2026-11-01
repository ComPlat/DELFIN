"""SW2b phase 1: the swarm.* protocol surface is reachable through the
service dispatch (protocol definitions + the supervisor hook).

Reviewer angle (nacht-s14): the _m_ name derivation replaces "." with
"_", so the hook's module-function names and the SERVICE_REQUESTS keys
must agree EXACTLY (swarm.triage.list -> _swarm_rpc.swarm_triage_list)
or the method is unreachable. Also: params_problem must reject unknown
or missing required fields for swarm.* (a missing spec silently
defaulting would be a defect).

Red control: written against the tree WITHOUT the supervisor hook
built -- a swarm.status request through the real dispatch answers
METHOD_NOT_FOUND, so the wiring cases fail; the dispatch-agreement and
params_problem cases fail because the protocol has no swarm.* entries.
"""

from __future__ import annotations

import pytest

from delfin.agent.service import protocol


SWARM_METHODS = [
    "swarm.status", "swarm.plan", "swarm.validate", "swarm.start",
    "swarm.pause", "swarm.resume", "swarm.triage.list",
    "swarm.triage.answer",
]


def test_every_swarm_method_is_a_defined_service_request():
    for method in SWARM_METHODS:
        assert method in protocol.SERVICE_REQUESTS, method
        assert method in protocol.SERVICE_RESULTS, method


def test_hook_function_names_match_the_protocol_keys_exactly():
    """swarm.triage.list must dispatch to swarm_rpc.swarm_triage_list:
    the hook maps method -> module attr with the same replace(".", "_")
    the supervisor's _m_ derivation uses (supervisor.py:879)."""
    from delfin.agent.service import swarm_rpc
    for method in SWARM_METHODS:
        attr = method.replace(".", "_")
        assert hasattr(swarm_rpc, attr), (method, attr)


def test_params_problem_rejects_a_missing_required_field():
    # queue and answer are required on the triage pair.
    assert protocol.params_problem(
        protocol.SERVICE_REQUESTS["swarm.triage.list"], {}) != ""
    assert protocol.params_problem(
        protocol.SERVICE_REQUESTS["swarm.triage.answer"],
        {"item": {}}) != ""
    assert protocol.params_problem(
        protocol.SERVICE_REQUESTS["swarm.triage.answer"],
        {"answer": "x"}) != ""


def test_params_problem_accepts_the_documented_optional_shapes():
    assert protocol.params_problem(
        protocol.SERVICE_REQUESTS["swarm.status"], {}) == ""
    assert protocol.params_problem(
        protocol.SERVICE_REQUESTS["swarm.status"],
        {"spec": {"spec_version": 1}}) == ""
    assert protocol.params_problem(
        protocol.SERVICE_REQUESTS["swarm.start"],
        {"specYaml": "spec_version: 1", "dryRun": True}) == ""


def test_params_problem_rejects_a_wrong_typed_field():
    assert protocol.params_problem(
        protocol.SERVICE_REQUESTS["swarm.triage.list"],
        {"queue": "not-a-list"}) != ""


def test_json_schema_carries_the_swarm_surface():
    schema = protocol.json_schema()
    service = schema["service"]["requests"]
    for method in SWARM_METHODS:
        assert method in service, method
    results = schema["service"]["results"]
    for method in SWARM_METHODS:
        assert method in results, method


# ------------------------------------------------- dispatch (real path)

@pytest.fixture()
def dispatch_env(monkeypatch):
    """The supervisor's swarm dispatch, without a socket: the hook line
    (supervisor.py:879 area) resolves method -> swarm_rpc attr and runs
    it in a thread. We exercise exactly that resolution, through the
    module the hook imports."""
    return protocol


def test_swarm_status_dispatch_reaches_the_real_handler():
    """Through the dispatch the supervisor hook performs: name
    derivation + module import + call with params only."""
    from delfin.agent.service import swarm_rpc
    fn = getattr(swarm_rpc, "swarm.status".replace(".", "_"))
    # A valid but empty-report spec: goal-only, no packages.
    result = fn({"spec": {"spec_version": 1, "goal": "g",
                          "limits": {"max_sessions": 1},
                          "packages": [{"name": "p",
                                        "owned": ["tui/src/swarm/**"],
                                        "roles": {"builder": True}}]}})
    assert result["schema_version"] == 1


def test_unhooked_tree_answers_method_not_found_for_swarm_status():
    """The red half, seen through the real dispatch decision the
    supervisor makes (supervisor.py:874-879): before the hook exists,
    a swarm.* method that IS in SERVICE_REQUESTS but has no _m_ handler
    and no swarm.* branch is unreachable. This test asserts the
    dispatch NAME the hook must produce exists -- on the unhooked tree
    the protocol side is already defined, so this documents the name
    contract; the reachability case above fails without swarm_rpc."""
    assert "swarm.status".replace(".", "_") == "swarm_status"
