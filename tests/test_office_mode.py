"""Office mode: administrative work, with the chemistry surface removed.

The separation is the point of the mode. A mode that merely *says* it is
about documents while the calc search and the ORCA manual stay reachable
would let an administrative question be answered with methodology, so
the exclusion is enforced where tools are authorised, not in prose.

Also covers the mode selector remembering the user's choice: it reset to
Dashboard on every dashboard start, which made every other mode a thing
you had to re-pick each session.
"""

from __future__ import annotations

import inspect
from pathlib import Path

import pytest

from delfin.agent import api_client as A
from delfin.agent.api_client import (
    _DELFIN_ONLY_TOOL_NAMES,
    _DOC_TOOLS_OPENAI,
    _ROLE_EXEC_DENYLIST,
    ToolSurfaceContext,
    _tool_denied_for_role,
    advertisable_tools,
    role_tool_surface_report,
    tool_unavailable_reason,
)
from delfin.agent.prompt_loader import PromptLoader


# ---------------------------------------------------------------------------
# The mode is registered and routes to its own role
# ---------------------------------------------------------------------------

def test_office_mode_is_available_and_routes_to_the_office_agent():
    loader = PromptLoader()
    assert "office" in loader.available_modes()
    mode = loader.load_mode("office")
    assert mode["route"] == ["office_agent"]
    assert mode["description"].strip().startswith("# Office Mode")


def test_the_office_role_prompt_exists_and_is_about_documents():
    prompt = PromptLoader().load_role_prompt("office_agent")
    assert prompt.strip().startswith("# Office Agent")
    for tool in ("read_document", "edit_sheet", "fill_pdf_form",
                 "fill_docx_template", "create_docx"):
        assert tool in prompt, tool


def test_the_office_prompt_does_not_send_the_model_to_chemistry_tools():
    """It may mention that chemistry is out of scope; it must not offer
    the tools this role cannot execute."""
    prompt = PromptLoader().load_role_prompt("office_agent")
    for tool in _DELFIN_ONLY_TOOL_NAMES:
        assert tool not in prompt, tool


def test_the_office_agent_is_declared_in_the_pack_manifest():
    import delfin.agent as agent_pkg

    manifest = (Path(agent_pkg.__file__).parent / "pack" / "manifest.yaml"
                ).read_text(encoding="utf-8")
    assert "id: office_agent" in manifest
    assert "agents/office_agent.md" in manifest


# ---------------------------------------------------------------------------
# The chemistry surface is denied, not merely hidden
# ---------------------------------------------------------------------------

def test_chemistry_tools_are_denied_for_the_office_role():
    for tool in _DELFIN_ONLY_TOOL_NAMES:
        assert _tool_denied_for_role("office_agent", tool), tool
        # …including through the MCP namespace, which is how a
        # namespaced call would otherwise slip past a bare-name check.
        assert _tool_denied_for_role(
            "office_agent", f"mcp__delfin-docs__{tool}"), tool


def test_the_office_role_keeps_the_tools_the_work_needs():
    for tool in ("read_document", "edit_sheet", "fill_pdf_form",
                 "fill_docx_template", "create_docx", "bash", "read_file",
                 "write_file", "task_create", "ask_user_question"):
        assert not _tool_denied_for_role("office_agent", tool), tool


def test_other_roles_are_unaffected_by_the_deny_list():
    for tool in _DELFIN_ONLY_TOOL_NAMES:
        assert not _tool_denied_for_role("solo_agent", tool), tool
        assert not _tool_denied_for_role("", tool), tool


def test_chemistry_tools_leave_the_advertised_surface_too():
    ctx = ToolSurfaceContext(role="office_agent")
    advertised = {t["function"]["name"]
                  for t in advertisable_tools(_DOC_TOOLS_OPENAI, ctx)}
    assert not (_DELFIN_ONLY_TOOL_NAMES & advertised)
    for tool in _DELFIN_ONLY_TOOL_NAMES:
        assert tool_unavailable_reason(tool, ctx) is not None


def test_advertising_stays_a_subset_of_what_may_execute():
    """The invariant the whole surface layer rests on, for this role too."""
    ctx = ToolSurfaceContext(role="office_agent")
    for tool in advertisable_tools(_DOC_TOOLS_OPENAI, ctx):
        name = tool["function"]["name"]
        assert not _tool_denied_for_role("office_agent", name), name


def test_the_office_role_appears_in_the_surface_report():
    report = role_tool_surface_report()
    assert "office_agent" in report
    assert report["office_agent"]["count"] > 0


def test_the_deny_list_only_names_real_tools():
    catalogue = {t["function"]["name"] for t in _DOC_TOOLS_OPENAI}
    for role, denied in _ROLE_EXEC_DENYLIST.items():
        assert denied <= catalogue, role


# ---------------------------------------------------------------------------
# The mode selector remembers the user's choice
# ---------------------------------------------------------------------------

def _tab_agent_source() -> str:
    from delfin.dashboard import tab_agent

    return inspect.getsource(tab_agent)


def test_office_is_offered_in_the_mode_selector():
    assert '("Office", "office")' in _tab_agent_source()


def test_a_user_mode_change_is_persisted():
    source = _tab_agent_source()
    handler = source.split("def _on_mode_change(", 1)[1][:2000]
    assert '_s["agent"]["mode"] = new_mode' in handler
    # Only a user-driven switch: an internal one is not a preference.
    assert handler.index('if not state.get("_mode_change_internal")') < \
        handler.index('_s["agent"]["mode"] = new_mode')


def test_the_saved_mode_is_restored_after_the_observer_is_wired():
    """Order matters: restoring before the observer exists would set the
    label without applying the mode's UI rules."""
    source = _tab_agent_source()
    observe_at = source.index('mode_dropdown.observe(_on_mode_change')
    restore_at = source.index('_set_mode_programmatically(_saved_mode)')
    assert observe_at < restore_at


def test_the_restore_goes_through_the_internal_setter():
    """_set_mode_programmatically clamps a retired/unknown saved mode, so a
    stale settings file cannot raise a TraitError on startup."""
    source = _tab_agent_source()
    assert "_set_mode_programmatically(_saved_mode)" in source


def test_mode_round_trips_through_user_settings(tmp_path, monkeypatch):
    import delfin.user_settings as us

    monkeypatch.setattr(Path, "home", lambda: tmp_path)
    settings = us.load_settings()
    settings.setdefault("agent", {})
    settings["agent"]["mode"] = "office"
    us.save_settings(settings)
    assert us.load_settings()["agent"]["mode"] == "office"
