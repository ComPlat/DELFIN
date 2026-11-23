"""A research role could set a lab device, because prose is not a gate.

MEASURED. Five modes are selectable from the manifest, and the execution
allow-list resolves for only four of them:

    dashboard          dashboard_agent     18 tools
    office             office_agent        38 tools
    workflow           workflow-operator   40 tools
    solo               solo_agent          None  <- the generalist
    research           research_agent      None  <- nobody decided this

``allow is None`` means ALLOW in ``_tool_denied_for_role``, so the
research route was offered the whole surface. Measured at the gate:

    research_agent  bash             denied   (another mechanism)
    research_agent  write_file       denied   (another mechanism)
    research_agent  run_application  ALLOWED
    research_agent  cluster_run      ALLOWED
    research_agent  lab_set          ALLOWED

A role whose job is reading literature could submit to the cluster and
move a lab device.

THE ROOT IS ONE MISSING LINE, not a missing rule.
``pack/agents/research_agent.md`` has no YAML front matter at all -- it
opens with ``# Research Agent`` -- so there is no ``tools:`` key for
``_role_allowlist`` to derive from, and the derivation that works for
workflow-operator (40 tools, bash and write_file correctly denied)
returns None here. Its own prose names a bounded set: search_docs,
read_section, list_docs, list_sections, WebSearch, WebFetch, reading the
codebase, and bash for git history. Nothing enforced it.

AND ``None`` STAYS LEGAL FOR ONE ROLE, which is why this is not simply a
new rule. The generalist is meant to hold the whole surface. The fix is
that ``None`` must be a DECISION with a name beside it rather than the
value a missing key happens to produce -- so a mode added next month
cannot arrive unrestricted by omission.
"""
from __future__ import annotations

import pathlib

import pytest
import yaml

from delfin.agent import api_client as ac

_MANIFEST = (pathlib.Path(ac.__file__).parent
             / "pack_lite" / "manifest.yaml")


def _route_roles() -> dict[str, list[str]]:
    data = yaml.safe_load(_MANIFEST.read_text()) or {}
    out = {}
    for mode in (data.get("modes") or []):
        name = str(mode.get("name") or mode.get("id") or "")
        out[name] = [str(r) for r in (mode.get("route") or [])]
    return out


# --- the gate, which is what counts -----------------------------------

@pytest.mark.parametrize("tool", ["lab_set", "lab_estop", "cluster_run",
                                  "cluster_submit", "run_application"])
def test_the_research_route_cannot_touch_hardware(tool):
    """It reads literature. It does not move a lab device or submit to a
    cluster, and until this change nothing said so except prose."""
    assert ac._tool_denied_for_role("research_agent", tool), (
        f"the research route may still call {tool}")


def test_the_research_route_keeps_what_its_own_prompt_names():
    """Restricting it must not break it: these are the tools its own
    prose describes."""
    for tool in ("search_docs", "read_section", "list_docs",
                 "list_sections", "read_file", "grep_file", "list_files",
                 "web_search", "web_fetch"):
        assert not ac._tool_denied_for_role("research_agent", tool), tool


# --- and the root: None is a decision, not an omission ----------------

def test_every_route_resolves_to_a_list_or_is_named_unrestricted():
    """The rule that makes the above permanent. A mode added later
    cannot arrive unrestricted because somebody forgot a front-matter
    key -- it either derives a list or is named here with a reason."""
    unexplained = []
    for mode, roles in _route_roles().items():
        for role in roles:
            if ac._role_allowlist(role) is None and (
                    role not in ac._ROLES_WITHOUT_A_TOOL_LIMIT):
                unexplained.append(f"{mode} -> {role}")
    assert not unexplained, (
        "these routes resolve to allow=None, which means ALLOW: "
        f"{unexplained}. Give the role a `tools:` front-matter key, or "
        "name it in _ROLES_WITHOUT_A_TOOL_LIMIT with the reason.")


def test_the_unrestricted_roles_are_few_and_reasoned():
    assert ac._ROLES_WITHOUT_A_TOOL_LIMIT
    for role, reason in ac._ROLES_WITHOUT_A_TOOL_LIMIT.items():
        assert len(reason) > 20, f"{role} has no stated reason"
    assert "solo_agent" in ac._ROLES_WITHOUT_A_TOOL_LIMIT


def test_a_new_mode_cannot_be_unrestricted_by_omission(tmp_path):
    """The failure mode this prevents, driven rather than described: a
    role whose prompt has no front matter resolves to None, and the
    check above is what turns that into a test failure instead of a
    silent full surface."""
    assert ac._role_allowlist("a-role-nobody-defined") is None
    assert "a-role-nobody-defined" not in ac._ROLES_WITHOUT_A_TOOL_LIMIT


# --- the escape hatch has to be verifiable, not just reasoned --------

_GENERALIST_ROLES = {"solo_agent", ""}


def test_a_named_role_is_either_the_generalist_or_really_bounded():
    """Otherwise _ROLES_WITHOUT_A_TOOL_LIMIT becomes the place omissions
    are filed. A role named there without an allow-list must carry a
    deny-list that covers everything acting outside the tree -- which
    is checked here against the derived set, not against a comment.
    """
    external = ac._external_effect_tool_names()
    assert external, "the derived set is empty; this test measures nothing"
    for role in ac._ROLES_WITHOUT_A_TOOL_LIMIT:
        if role in _GENERALIST_ROLES:
            continue
        deny = ac._ROLE_EXEC_DENYLIST.get(role) or frozenset()
        missing = sorted(external - deny)
        assert not missing, (
            f"{role} is named as unrestricted but its deny-list does not "
            f"cover: {missing}")


def test_every_reading_role_refuses_what_acts_outside_the_tree():
    """The six roles whose stated purpose is reading. Measured before
    this change: each was allowed all 22 acting tools."""
    reading = ["research_agent", "critic_agent", "reviewer_agent",
               "chief_agent", "session_manager", "runtime_agent"]
    for role in reading:
        for tool in sorted(ac._external_effect_tool_names()):
            assert ac._tool_denied_for_role(role, tool), f"{role} / {tool}"


def test_a_namespaced_call_is_judged_by_the_same_names():
    """_tool_denied_for_role strips mcp__server__, which is why the
    runner names can be listed without their namespace -- and why a
    rule written for `lab_set` is not escaped by calling
    mcp__delfin-tools__lab_set."""
    assert ac._tool_denied_for_role(
        "reviewer_agent", "mcp__delfin-tools__lab_set")
    assert ac._tool_denied_for_role(
        "reviewer_agent", "mcp__anything__run_application")


def test_the_reporting_tools_stay_reachable():
    """No regression: reading state moves nothing, and these six are
    what the reading roles have a use for."""
    for tool in sorted(ac._HARDWARE_REPORTING_TOOLS):
        assert not ac._tool_denied_for_role("reviewer_agent", tool), tool


def test_a_new_tool_in_a_family_is_denied_without_an_edit():
    """The direction that matters. The previous arrangement was a
    deny-list of what somebody had thought of; this one denies a
    lab_/cluster_ tool added later from the moment it exists."""
    import json

    from delfin.agent import tool_catalog

    entry = {"type": "function", "function": {
        "name": "lab_move_arm",
        "description": "move the arm",
        "parameters": {"type": "object", "properties": {}}}}
    original = list(tool_catalog._DOC_TOOLS_OPENAI)
    try:
        tool_catalog._DOC_TOOLS_OPENAI.append(entry)
        assert "lab_move_arm" in ac._external_effect_tool_names()
    finally:
        tool_catalog._DOC_TOOLS_OPENAI[:] = original
    assert "lab_move_arm" not in ac._external_effect_tool_names()
