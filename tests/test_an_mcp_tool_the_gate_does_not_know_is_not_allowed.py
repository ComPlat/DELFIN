"""The MCP gate allowed everything it did not recognise.

``_gate_mcp_tool`` classified side effects with a five-name write map
(``write_file``, ``edit_file``, ``multi_edit``, ``notebook_edit``,
``apply_patch``) plus the shell and network families, and returned ``None``
for every other name. "Every other name" is the normal case: a server names
its own actions, and no map written here can contain them.

Measured before the fix, on the always-registered tools server:

    mcp__delfin-tools__register_module   default            -> None (allowed)
    mcp__delfin-tools__register_module   acceptEdits        -> None (allowed)
    mcp__delfin-tools__register_module   bypassPermissions  -> None (allowed)
    mcp__kit-coding__write_file (~/.delfin/adapters/x.py)   -> BLOCKED

``register_module(name, code)`` writes ``<adapters_dir>/<name>.py`` and
imports it in the same call, so it is a write plus an execution — the one
thing the write gate exists to stop — and it went through the gate that
blocked the plain write to the same directory. ``name`` was never
validated: ``"../escaped_outside_adapters"`` resolves one directory ABOVE
the adapters directory, an absolute name resolves wherever it points, and
the cleanup that unlinks the file on failure runs only after the import has
already executed it.

The same hole covered the ops server: delete_calc_folder,
kill_all_user_jobs, pipeline_run and move_to_archive were gated by an
``allow_mutate`` boolean the model itself passed.

Nothing here executes a payload. The assertions are on the DECISION the
gate returns and on the PATH a name resolves to.
"""

from __future__ import annotations

import os
from pathlib import Path

import pytest

from delfin.agent import api_client as A
from delfin.agent.api_client import KitToolPermissions, _doc_executor


MODES = ("default", "acceptEdits", "bypassPermissions")


def _perms(tmp_path, mode="default", **kw):
    return KitToolPermissions(workspace=tmp_path, mode=mode, **kw)


def _gate(name, args, perms):
    return _doc_executor._gate_mcp_tool(name, args, perms)


# ---------------------------------------------------------------------------
# 1. The measured case: writing and importing Python through the tools server
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("mode", ["default", "acceptEdits"])
def test_register_module_is_not_allowed_by_silence(tmp_path, mode):
    """The call the gate used to return None for, head-less."""
    out = _gate("mcp__delfin-tools__register_module",
                {"name": "my_new_block", "code": "print(1)"},
                _perms(tmp_path, mode))
    assert out is not None, (
        "a tool that writes source and imports it passed the gate unread")
    assert "executes it" in out


@pytest.mark.parametrize("mode", MODES)
def test_a_traversal_module_name_is_refused_in_every_mode(tmp_path, mode):
    """Refused before any approval, so nobody is asked to approve a
    traversal — including under bypassPermissions, which opts out of being
    ASKED, not out of the adapters directory having a boundary."""
    out = _gate("mcp__delfin-tools__register_module",
                {"name": "../escaped_outside_adapters", "code": "print(1)"},
                _perms(tmp_path, mode))
    assert out is not None and "invalid module name" in out


@pytest.mark.parametrize("bad", [
    "../escaped_outside_adapters",
    "/tmp/absolute_escape",
    "sub/nested",
    "dotted.name",
    "with-dash",
    "",
])
def test_the_writer_refuses_the_same_names_the_gate_refuses(bad):
    from delfin.tools.platform import module_name_error
    assert module_name_error(bad) is not None, bad


@pytest.mark.parametrize("good", ["my_block", "_private", "b2plyp_parser", "x1"])
def test_a_legitimate_module_name_still_works(good):
    from delfin.tools.platform import module_name_error
    assert module_name_error(good) is None, good


def test_a_traversal_name_resolves_outside_the_adapters_dir(tmp_path,
                                                            monkeypatch):
    """The premise, asserted on the resolved PATH and nothing else.

    No file is written and no code is imported here — this only shows that
    the name the gate now refuses is a name that leaves the directory.
    """
    adapters = tmp_path / "adapters"
    adapters.mkdir()
    monkeypatch.setenv("DELFIN_ADAPTERS_DIR", str(adapters))
    from delfin.tools import _registry
    target = _registry._user_adapter_dirs()[0]
    assert target.resolve() == adapters.resolve()

    escaped = (target / "../escaped_outside_adapters.py").resolve()
    assert not str(escaped).startswith(str(adapters.resolve()) + os.sep)
    assert (target / "ok_name.py").resolve().parent == adapters.resolve()


def test_the_writer_refuses_a_traversal_without_writing_anything(tmp_path):
    from delfin.tools import platform
    adapters = tmp_path / "adapters"
    adapters.mkdir()
    out = platform.register_module("../escaped", "print(1)", directory=adapters)
    assert out["ok"] is False
    assert "invalid module name" in out["diagnostics"][0]
    assert not (tmp_path / "escaped.py").exists()
    assert list(adapters.iterdir()) == []


# ---------------------------------------------------------------------------
# 2. Deny-by-default in general
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("tool", [
    "mcp__delfin-ops__delete_calc_folder",
    "mcp__delfin-ops__kill_all_user_jobs",
    "mcp__delfin-ops__pipeline_run",
    "mcp__delfin-ops__move_to_archive",
    "mcp__github__create_pull_request",
    "mcp__jira__create_issue",
    "mcp__db__delete_row",
    "mcp__slack__post_message",
])
def test_an_unrecognised_mcp_tool_is_refused_head_less(tmp_path, tool):
    out = _gate(tool, {"anything": "at all"}, _perms(tmp_path, "default"))
    assert out is not None, f"{tool} passed the gate unread"
    assert "read-only" in out


def test_a_self_supplied_allow_mutate_does_not_change_the_verdict(tmp_path):
    """The ops server's consent model, seen from the gate: the flag is the
    model's own word about the model, so it must not move the decision."""
    without = _gate("mcp__delfin-ops__delete_calc_folder",
                    {"folder": "x"}, _perms(tmp_path, "default"))
    with_flag = _gate("mcp__delfin-ops__delete_calc_folder",
                      {"folder": "x", "allow_mutate": True},
                      _perms(tmp_path, "default"))
    assert without is not None and with_flag is not None


def test_the_user_is_asked_when_someone_can_be_asked(tmp_path):
    asked: list[tuple] = []

    def _confirm(name, args, preview):
        asked.append((name, preview))
        return True

    perms = _perms(tmp_path, "default", confirm_callback=_confirm)
    assert _gate("mcp__delfin-ops__pipeline_run",
                 {"control_file": "CONTROL.txt"}, perms) is None
    assert asked and asked[0][0] == "mcp__delfin-ops__pipeline_run"
    # The user is judging a tool name that means nothing outside its
    # server, so the arguments have to be in front of them.
    assert "CONTROL.txt" in asked[0][1]


def test_a_refusal_is_a_refusal(tmp_path):
    perms = _perms(tmp_path, "default",
                   confirm_callback=lambda n, a, p: False)
    out = _gate("mcp__delfin-ops__delete_calc_folder", {"folder": "x"}, perms)
    assert out is not None and "denied" in out
    # ...and it is remembered, so the same call through any server stays
    # refused without asking again.
    again = _gate("mcp__delfin-ops__delete_calc_folder", {"folder": "x"},
                  perms)
    assert again is not None and "already refused" in again


def test_bypass_permissions_still_asks_nothing(tmp_path):
    asked: list[str] = []
    perms = _perms(tmp_path, "bypassPermissions",
                   confirm_callback=lambda n, a, p: asked.append(n) or True)
    assert _gate("mcp__delfin-ops__pipeline_run", {}, perms) is None
    assert asked == []


def test_a_dry_run_needs_no_decision(tmp_path):
    assert _gate("mcp__github__create_pull_request", {"check_only": True},
                 _perms(tmp_path, "default")) is None


# ---------------------------------------------------------------------------
# 3. The legitimate case: read-only tools are untouched
# ---------------------------------------------------------------------------

READ_ONLY_CALLS = [
    ("mcp__delfin-tools__list_capabilities", {}),
    ("mcp__delfin-tools__get_guide", {}),
    ("mcp__delfin-tools__describe_capability", {"name": "orca_opt"}),
    ("mcp__delfin-tools__compatible_successors", {"name": "orca_freq"}),
    ("mcp__delfin-ops__parse_orca_output", {"folder": "calc/x"}),
    ("mcp__delfin-ops__extract_energy_table", {"folders": "a,b"}),
    ("mcp__delfin-ops__qm_check", {}),
    ("mcp__delfin-ops__list_tools", {"category": "parsing"}),
    ("mcp__delfin-docs__search_docs", {"query": "RIJCOSX"}),
    ("mcp__coding__read_file", {"path": "a.py"}),
    ("mcp__fs__list_directory", {"path": "."}),
]


@pytest.mark.parametrize("mode", list(MODES) + ["plan"])
@pytest.mark.parametrize("name,args", READ_ONLY_CALLS)
def test_a_read_only_tool_is_never_gated(tmp_path, mode, name, args):
    assert _gate(name, args, _perms(tmp_path, mode)) is None, name


def test_plan_mode_can_still_discover(tmp_path):
    """Plan mode refused every MCP name outside three families, so a
    planning turn could not call the discovery tools the role prompt sends
    it to first. Investigating is what the mode is for."""
    perms = _perms(tmp_path, "plan")
    for name, args in READ_ONLY_CALLS:
        assert _gate(name, args, perms) is None, name
    assert _gate("mcp__delfin-tools__run_application", {}, perms) is not None


# ---------------------------------------------------------------------------
# 4. The allow-list is a snapshot, so a NEW server tool defaults to deny
# ---------------------------------------------------------------------------

def _registered_tool_names() -> dict[str, set[str]]:
    """Tool names each built-in DELFIN server actually registers."""
    import re
    root = Path(__file__).resolve().parents[1] / "delfin"
    out: dict[str, set[str]] = {}
    tools_src = (root / "tools" / "mcp_server.py").read_text(encoding="utf-8")
    out["delfin-tools"] = set(re.findall(
        r"@mcp\.tool\(\)\s*\n\s*def\s+([a-z_][a-z0-9_]*)", tools_src))
    ops_src = (root / "ops_server" / "server.py").read_text(encoding="utf-8")
    out["delfin-ops"] = set(re.findall(
        r"mcp\.tool\(name=\"([a-z_][a-z0-9_]*)\"", ops_src))
    docs_src = (root / "doc_server" / "server.py").read_text(encoding="utf-8")
    out["delfin-docs"] = set(re.findall(
        r"@mcp\.tool\(\)\s*\n\s*def\s+([a-z_][a-z0-9_]*)", docs_src))
    return out


# The DELFIN server tools the gate lets through untouched. A snapshot on
# purpose: a tool added to a server is NOT on it, so it starts out denied
# and someone has to decide it is read-only and say so here. That is the
# whole property — the previous arrangement made a new tool allowed by
# default and nothing noticed for the seventy-odd that already were.
_EXPECTED_READ_ONLY = {
    # delfin-tools
    "get_manifest", "get_guide", "resolve_spec", "scientific_lint",
    "list_capabilities", "describe_capability", "catalog",
    "compatible_successors", "new_capability_template", "list_keys",
    "describe_key", "list_applications", "describe_application",
    "validate_application", "validate_spec", "run_diagnostics",
    "run_status", "list_runs", "run_metrics", "probe", "install_plan",
    # delfin-docs
    "search_docs", "read_section", "list_docs", "list_sections",
    "search_calcs", "get_calc_info", "calc_summary",
    # delfin-ops
    "qm_check", "csp_check", "mlp_check", "analysis_check", "stop_dry_run",
    "parse_orca_output", "find_orca_errors", "extract_thermochem",
    "extract_energy_table", "find_calculation_extreme",
    "extract_imaginary_frequencies", "extract_orbital_energies",
    "extract_excited_states", "extract_dipole",
    "extract_optimization_trajectory", "extract_scf_convergence",
    "extract_mulliken_charges", "extract_loewdin_charges",
    "extract_vibrational_modes", "extract_delfin_json",
    "extract_calc_summary_table", "compare_calculations",
    "compare_across_functionals", "list_tools", "describe_tool",
    "list_dashboard_patterns", "get_dashboard_pattern",
    "list_dashboard_widgets", "get_widget_options", "validate_orca_input",
    "list_delfin_features", "explain_delfin_feature",
    "list_active_calculations", "list_ssh_transfer_jobs",
    "list_calc_options", "list_literature_files",
    "check_orca_manual_indexed", "read_pdf", "search_pdf_local",
    "extract_pdf_section",
}


def test_the_read_only_allow_list_is_a_deliberate_snapshot():
    registered = set().union(*_registered_tool_names().values())
    allowed = registered & A._MCP_READONLY_TOOL_BASES
    assert allowed == _EXPECTED_READ_ONLY, (
        "the set of DELFIN server tools that pass the gate untouched "
        "changed. A NEW tool defaults to deny — if it really changes "
        "nothing, add it to _MCP_READONLY_TOOL_BASES and to this snapshot; "
        "otherwise leave it denied.\n"
        f"newly allowed: {sorted(allowed - _EXPECTED_READ_ONLY)}\n"
        f"no longer allowed: {sorted(_EXPECTED_READ_ONLY - allowed)}")


def test_every_mutating_server_tool_is_denied(tmp_path):
    registered = set().union(*_registered_tool_names().values())
    perms = _perms(tmp_path, "default")
    for base in sorted(registered - A._MCP_READONLY_TOOL_BASES):
        out = _gate(f"mcp__delfin-ops__{base}", {}, perms)
        assert out is not None, f"{base} passes the gate without a decision"
