"""A read-only delegate can use the lookups its preset offers it.

Input: the permissions ``run_subagent`` builds for a preset (``_derive_perms``
with the preset's mode and tool list, depth 1), and a call through the real
executor -- ``_DocToolExecutor.execute`` for the native name and
``_gate_mcp_tool`` for the ``mcp__kit-coding__`` spelling the
OpenAI-compatible backends send.

Semantics pinned here:

* every tool a plan-mode preset lists passes the plan-mode gate. The list
  is the delegate's advertised surface, so a listed tool the mode refuses
  is an offer the gate withdraws on the first call. Recorded failure: an
  explore delegate asked to survey the literature called web_search, was
  refused by plan mode on every call, and spent its 300 s budget on the
  refusals;
* web_search runs for a delegate whose preset names it. It sends a query
  to a fixed set of search providers, never to a host the model chooses,
  and the query passes the outbound secret check first;
* web_fetch stays refused in plan mode and is no longer listed: its URL
  names an arbitrary host and carries arbitrary data, which is an outbound
  channel rather than a lookup;
* project_introspect stays refused in plan mode (it runs the workspace's
  own interpreter) and is no longer listed;
* subagent_message runs for a delegate: it is how the delegate reports to
  the session running it, like subagent_result;
* nothing changes for a session at depth 0, and a delegate of a plan-mode
  parent does not get the search its parent could not run.
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from delfin.agent import api_client as A
from delfin.agent import subagents as S


@pytest.fixture(autouse=True)
def _no_network_no_model(monkeypatch):
    # The executors are replaced so a call that passes the gate returns a
    # marker instead of leaving the machine; the classifier is off so no
    # second opinion is requested from a model endpoint.
    monkeypatch.setenv("DELFIN_CLASSIFIER", "off")
    monkeypatch.setattr(A._DocToolExecutor, "_execute_web_search",
                        lambda self, args: json.dumps({"ran": "web_search"}))
    monkeypatch.setattr(A._DocToolExecutor, "_execute_web_fetch",
                        lambda self, args: json.dumps({"ran": "web_fetch"}))
    monkeypatch.setattr(S, "post_reply",
                        lambda me, text: {"ran": "subagent_message"})


def _parent(tmp_path: Path, mode: str = "default") -> A.KitToolPermissions:
    return A.KitToolPermissions(workspace=tmp_path, mode=mode)


def _child(tmp_path: Path, preset: str, parent_mode: str = "default"):
    p = S.SUBAGENT_PRESETS[preset]
    child = S._derive_perms(_parent(tmp_path, parent_mode), p.mode,
                            preset_tools=p.tools)
    child.subagent_id = "sa-test"
    return child


def _run(name, args, perms) -> str:
    ex = A._DocToolExecutor()
    return ex.execute(name, args, perms)


def _plan_refused(out: str) -> bool:
    try:
        return "plan mode" in str((json.loads(out) or {}).get("error", ""))
    except Exception:
        return False


def test_an_explore_child_runs_web_search(tmp_path):
    child = _child(tmp_path, "explore")
    assert child.mode == "plan" and child.subagent_depth == 1
    out = _run("web_search", {"query": "DLPNO-CCSD(T) benchmark"}, child)
    assert json.loads(out) == {"ran": "web_search"}, out


def test_the_namespaced_search_passes_the_mcp_gate(tmp_path):
    child = _child(tmp_path, "explore")
    ex = A._DocToolExecutor()
    assert ex._gate_mcp_tool("mcp__kit-coding__web_search",
                             {"query": "r2SCAN-3c"}, child) is None


def test_a_research_preset_runs_web_search(tmp_path):
    child = _child(tmp_path, "method-researcher")
    out = _run("web_search", {"query": "implicit solvation pKa"}, child)
    assert json.loads(out) == {"ran": "web_search"}, out


def test_a_child_reports_to_its_session_in_plan_mode(tmp_path):
    child = _child(tmp_path, "explore")
    out = _run("subagent_message", {"message": "need a decision"}, child)
    assert json.loads(out) == {"ran": "subagent_message"}, out


def test_web_fetch_stays_refused_and_is_not_offered(tmp_path):
    for name, preset in S.SUBAGENT_PRESETS.items():
        if preset.mode == "plan":
            assert "web_fetch" not in preset.tools, name
    child = _child(tmp_path, "explore")
    out = _run("web_fetch", {"url": "https://example.org"}, child)
    assert "ran" not in json.loads(out), out
    # A preset that lists it anyway (a markdown file on disk) is still
    # refused by the mode, not only by the list.
    listed = S._derive_perms(_parent(tmp_path), "plan",
                             preset_tools=("web_fetch", "web_search"))
    listed.subagent_id = "sa-test"
    assert _plan_refused(_run("web_fetch", {"url": "https://example.org"},
                              listed))
    assert json.loads(_run("web_search", {"query": "x"}, listed)) == {
        "ran": "web_search"}


def test_a_session_in_plan_mode_still_cannot_search(tmp_path):
    perms = _parent(tmp_path, "plan")
    assert _plan_refused(_run("web_search", {"query": "x"}, perms))
    ex = A._DocToolExecutor()
    assert ex._gate_mcp_tool("mcp__kit-coding__web_search",
                             {"query": "x"}, perms)


def test_a_plan_parent_does_not_search_through_its_child(tmp_path):
    child = _child(tmp_path, "explore", parent_mode="plan")
    assert "web_search" not in (child.session_allowed_tools or frozenset())
    out = _run("web_search", {"query": "x"}, child)
    assert "ran" not in json.loads(out), out


def test_a_child_without_a_list_does_not_gain_the_search(tmp_path):
    # general-purpose declares no list; clamped to plan by its parent it
    # must not be admitted to a lookup no preset granted.
    child = _child(tmp_path, "general-purpose", parent_mode="plan")
    assert child.mode == "plan"
    assert _plan_refused(_run("web_search", {"query": "x"}, child))


def test_a_locked_scope_still_asks_before_a_child_searches(tmp_path):
    asked = []
    parent = A.KitToolPermissions(
        workspace=tmp_path, mode="default", lock_workspace=True,
        confirm_callback=lambda n, a, p: asked.append(n) or False)
    p = S.SUBAGENT_PRESETS["explore"]
    child = S._derive_perms(parent, p.mode, preset_tools=p.tools)
    out = _run("web_search", {"query": "x"}, child)
    assert asked == ["web_search"]
    assert "denied" in json.loads(out)["error"]


@pytest.mark.parametrize("preset", sorted(
    n for n, p in S.SUBAGENT_PRESETS.items() if p.mode == "plan"))
def test_no_plan_preset_offers_a_tool_its_mode_refuses(tmp_path, preset):
    # A listed ``mcp__<server>__<tool>`` name is judged where the tool loop
    # sends it: ``_gate_mcp_tool``, not the native ``execute``, which knows
    # no server tool and would refuse every one of them for being unknown.
    child = _child(tmp_path, preset)
    ex = A._DocToolExecutor()

    def refused_for(t: str) -> bool:
        if t.startswith("mcp__"):
            return ex._gate_mcp_tool(t, {}, child) is not None
        return _plan_refused(_run(t, {}, child))

    refused = [t for t in S.SUBAGENT_PRESETS[preset].tools if refused_for(t)]
    assert refused == [], (preset, refused)
