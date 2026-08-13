"""The ops server's destructive tools were gated by the caller's own word.

``delete_calc_folder``, ``kill_all_user_jobs``, ``pipeline_run``,
``move_to_archive`` and eleven others took ``allow_mutate: bool = False`` as
an MCP tool PARAMETER. An MCP tool parameter travels inside the model's
message, so the consent for a destructive act was supplied by the party the
consent exists to restrain — self-attested, and invisible to every gate
upstream, which saw only "an argument". The refusal even said so out loud:
"Pass allow_mutate=True to execute."

The parameter is gone from the MCP surface. The value comes from the process
that started the server. The Python functions keep the keyword so an
in-process caller (dashboard, tests) can still pass a decision it actually
made.
"""

from __future__ import annotations

import inspect

import pytest

from delfin.ops_server import server as ops


_MUTATING = [
    "submit_calculation", "cancel_calculation", "rename_calc_folder",
    "create_calc_folder", "move_calc_folder", "move_to_archive",
    "delete_calc_folder", "kill_all_user_jobs", "prepare_recalc",
    "run_calc_option", "cleanup", "stop", "pipeline_prepare",
    "pipeline_run", "run_orca_input", "co2", "tadf_xtb", "hyperpol",
]


def test_every_mutating_tool_still_declares_the_keyword_internally():
    """The premise: these are the tools whose consent was a parameter."""
    for name in _MUTATING:
        fn = getattr(ops, f"tool_{name}")
        assert "allow_mutate" in inspect.signature(fn).parameters, name


@pytest.mark.parametrize("name", _MUTATING)
def test_the_model_facing_wrapper_has_no_allow_mutate_parameter(name):
    exposed = ops._host_gated(getattr(ops, f"tool_{name}"))
    assert "allow_mutate" not in inspect.signature(exposed).parameters, (
        f"{name} still lets its caller pass its own permission")


@pytest.mark.parametrize("name", _MUTATING)
def test_the_description_no_longer_tells_the_model_to_pass_it(name):
    exposed = ops._host_gated(getattr(ops, f"tool_{name}"))
    assert "allow_mutate" not in (exposed.__doc__ or ""), name


def test_a_model_supplied_flag_is_discarded(monkeypatch):
    monkeypatch.delenv(ops._MUTATION_ENV, raising=False)
    seen: dict = {}

    def _fake(folder: str = "", allow_mutate: bool = False) -> str:
        seen["allow_mutate"] = allow_mutate
        return "called"

    gated = ops._host_gated(_fake)
    assert gated(folder="x", allow_mutate=True) == "called"
    assert seen["allow_mutate"] is False, (
        "the model's own allow_mutate reached the implementation")


def test_the_host_grant_is_what_lets_it_through(monkeypatch):
    seen: dict = {}

    def _fake(folder: str = "", allow_mutate: bool = False) -> str:
        seen["allow_mutate"] = allow_mutate
        return "called"

    gated = ops._host_gated(_fake)
    monkeypatch.setenv(ops._MUTATION_ENV, "1")
    gated(folder="x")
    assert seen["allow_mutate"] is True
    monkeypatch.setenv(ops._MUTATION_ENV, "0")
    gated(folder="x")
    assert seen["allow_mutate"] is False


@pytest.mark.parametrize("raw,expected", [
    ("1", True), ("true", True), ("TRUE", True), ("yes", True), ("on", True),
    ("0", False), ("false", False), ("", False), ("maybe", False),
])
def test_the_grant_is_read_strictly(monkeypatch, raw, expected):
    monkeypatch.setenv(ops._MUTATION_ENV, raw)
    assert ops.host_grants_mutation() is expected


def test_the_refusal_does_not_teach_the_model_to_pass_a_flag():
    import json
    payload = json.loads(ops._refuse_mutation("delete_calc_folder"))
    assert payload["error"] == "mutation_blocked"
    assert "allow_mutate" not in payload["message"]
    assert "USER" in payload["message"]


@pytest.mark.parametrize("name,kwargs", [
    ("pipeline_run", {"control_file": "CONTROL.txt"}),
    ("run_orca_input", {"input_file": "x.inp"}),
    ("tadf_xtb", {}),
    ("hyperpol", {}),
    ("co2", {}),
    ("pipeline_prepare", {}),
])
def test_a_mutating_tool_refuses_without_the_grant(monkeypatch, name, kwargs):
    import json
    monkeypatch.delenv(ops._MUTATION_ENV, raising=False)
    gated = ops._host_gated(getattr(ops, f"tool_{name}"))
    payload = json.loads(gated(allow_mutate=True, **kwargs))
    assert payload["error"] == "mutation_blocked", (
        f"{name} ran on the caller's own say-so")


def test_the_destructive_folder_delete_never_sees_a_model_grant(monkeypatch):
    """delete_calc_folder delegates its three locks to the API layer, so
    what matters here is which value reaches it."""
    seen: dict = {}
    monkeypatch.delenv(ops._MUTATION_ENV, raising=False)
    monkeypatch.setattr(
        ops.delfin_api, "delete_calc_folder",
        lambda folder, **kw: seen.update(folder=folder, **kw) or {"ok": False})
    gated = ops._host_gated(ops.tool_delete_calc_folder)
    gated(folder="x", confirm_token="x", allow_mutate=True)
    assert seen["allow_mutate"] is False


def test_an_in_process_caller_can_still_decide(monkeypatch):
    """The dashboard asks the user and then calls the function directly.
    That decision is real and must keep working."""
    monkeypatch.delenv(ops._MUTATION_ENV, raising=False)
    seen: dict = {}
    monkeypatch.setattr(
        ops.delfin_api, "delete_calc_folder",
        lambda folder, **kw: seen.update(folder=folder, **kw) or {"ok": True})
    ops.tool_delete_calc_folder(folder="x", confirm_token="x",
                                allow_mutate=True)
    assert seen["allow_mutate"] is True


def test_the_read_only_tools_are_untouched():
    """Nothing was wrapped that did not need it."""
    for name in ("qm_check", "parse_orca_output", "list_active_calculations"):
        fn = getattr(ops, f"tool_{name}")
        assert "allow_mutate" not in inspect.signature(fn).parameters, name


def test_the_module_docstring_no_longer_advertises_the_parameter():
    doc = ops.__doc__ or ""
    assert "require explicit ``allow_mutate=True``" not in doc
    assert "DELFIN_OPS_ALLOW_MUTATE" in doc
