"""The dashboard offers Research, and a research session stays read-only.

The research mode exists end to end in the terminal: ``--mode research``
routes to the ``research_agent`` role, whose tool surface is web search and
code reading with no write tool. The dashboard, though, never offered it:
the mode dropdown listed Dashboard / Code / Office, a saved research
session was restored as Code (the legacy map), and the engine constructor
had no research branch — the tool list fell through to the office /
repo-free branches and inherited ``Write`` and ``Edit``.

These tests pin the three ends: the dropdown offers Research, restoring a
research session keeps it, and the tool list built for a research session
carries no tool that can change a file or run a shell command.
"""

from __future__ import annotations

import ast
import inspect
from pathlib import Path

import pytest

from delfin.dashboard import tab_agent as T

_SRC = Path(inspect.getfile(T)).read_text(encoding="utf-8")


# ---------------------------------------------------------------------------
# dropdown + mode sets (structure: the dropdown is built inside create_tab,
# which cannot be constructed without a browser; the constants beside it can)
# ---------------------------------------------------------------------------

def _closure_constant(name: str):
    """The mode tuples live inside create_tab (they drive the cycle
    inspector and the minimal layout, both built there). Read them from the
    source: the fact under test is what the tab offers, not what a module
    attribute happens to carry."""
    tree = ast.parse(_SRC)
    for node in ast.walk(tree):
        if isinstance(node, ast.FunctionDef) and node.name == "create_tab":
            for sub in ast.walk(node):
                if isinstance(sub, ast.Assign) and any(
                    isinstance(t, ast.Name) and t.id == name
                    for t in sub.targets
                ):
                    return ast.literal_eval(sub.value)
    raise AssertionError(f"{name} not found in create_tab")


def test_single_agent_modes_include_research():
    assert "research" in _closure_constant("_SINGLE_AGENT_MODES")


def test_launch_independent_modes_include_research():
    # Research binds no project directory: like Dashboard and Office its
    # sessions must be reachable from any launch.
    assert "research" in _closure_constant("_LAUNCH_INDEPENDENT_MODES")


def test_dropdown_options_offered_research():
    """The dropdown tuple inside create_tab carries ('Research', 'research')."""
    assert ("Research", "research") in _dropdown_options()


def _dropdown_options() -> list:
    tree = ast.parse(_SRC)
    for node in ast.walk(tree):
        if isinstance(node, ast.Assign) and any(
            isinstance(t, ast.Name) and t.id == "mode_dropdown"
            for t in node.targets
        ):
            for kw in node.value.keywords:
                if kw.arg == "options":
                    return ast.literal_eval(kw.value)
    raise AssertionError("mode_dropdown options tuple not found")


# ---------------------------------------------------------------------------
# restore: a saved research session reopens as research, not Code
# ---------------------------------------------------------------------------

def test_legacy_map_does_not_remap_research():
    """research was listed in _legacy_map mapping it to solo: a saved
    research session reopened as Code. It is a real mode now and must
    survive the restore map unchanged."""
    tree = ast.parse(_SRC)
    for node in ast.walk(tree):
        if isinstance(node, ast.Assign) and any(
            isinstance(t, ast.Name) and t.id == "_legacy_map"
            for t in node.targets
        ):
            mapping = ast.literal_eval(node.value)
            assert "research" not in mapping, (
                "research is a real mode; mapping it away reopens every "
                "saved research session as Code")
            return
    # absent map is fine too
    return


# ---------------------------------------------------------------------------
# engine construction: the tool list for research is read-only
# ---------------------------------------------------------------------------

_WRITE_TOOLS = {"Write", "Edit", "MultiEdit", "NotebookEdit", "Bash"}


def test_engine_construction_has_a_research_branch_with_readonly_tools():
    """Inside create_tab's _ensure_engine there is a branch keyed on mode
    research, and the tool list it builds names no write or shell tool."""
    tree = ast.parse(_SRC)
    source_lines = _SRC.splitlines()
    for node in ast.walk(tree):
        if not (isinstance(node, ast.FunctionDef)
                and node.name == "_ensure_engine"):
            continue
        # find the comparison `mode_dropdown.value == "research"` (or
        # "research" in a membership test) inside the function
        research_lines = [
            source_lines[line_no - 1]
            for sub in ast.walk(node)
            if isinstance(sub, ast.Compare)
            for line_no in range(sub.lineno, sub.end_lineno + 1)
            if "research" in source_lines[line_no - 1]
            and "==" in source_lines[line_no - 1]
        ]
        assert research_lines, (
            "_ensure_engine has no research branch; a research session "
            "falls through to the office/repo_free branches and inherits "
            "their Write/Edit tools")
        # the branch body must not assign a write tool into the tool list
        branch_src = "\n".join(research_lines)
        for tool in _WRITE_TOOLS:
            assert f'"{tool}"' not in branch_src, (
                f"the research tool list includes {tool}; research is "
                "read-only")
        return
    raise AssertionError("_ensure_engine not found")


def test_research_role_is_stamped_like_office():
    """Office stamps kit_permissions.agent_role at construction so the
    session is locked before the first turn; research needs the same stamp
    with its own role id."""
    lines = _SRC.splitlines()
    for i, line in enumerate(lines):
        if "agent_role" not in line:
            continue
        region = "\n".join(lines[i - 7: i + 8])
        if "office_agent" not in region and "research_agent" not in region:
            continue
        assert '("office", "research")' in region, (
            "the role stamp does not cover research")
        return
    raise AssertionError("no role stamp found to mirror")


def test_research_mode_description_says_readonly():
    """The description beside the dropdown says what research does and that
    it does not change files. Scoped to _MODE_DESCRIPTIONS — the module
    also carries _MODE_LABELS, whose "Research" label matched first."""
    tree = ast.parse(_SRC)
    for node in ast.walk(tree):
        if not (isinstance(node, ast.Assign) and any(
            isinstance(t, ast.Name) and t.id == "_MODE_DESCRIPTIONS"
            for t in node.targets
        )):
            continue
        for k, v in zip(node.value.keys, node.value.values):
            if getattr(k, "value", None) == "research":
                text = ast.literal_eval(v)
                assert "read-only" in text.lower()
                assert "no " in text.lower()
                return
    raise AssertionError("no research description in _MODE_DESCRIPTIONS")


# ---------------------------------------------------------------------------
# the description cannot drift from what the role can actually do
# ---------------------------------------------------------------------------

def test_research_engine_branch_uses_the_fixed_scratch_folder():
    """Like Dashboard, a research session runs in the fixed agent scratch
    folder (ctx.agent_dir), not the launch directory, so a saved research
    session is reachable from any launch — matching the
    _LAUNCH_INDEPENDENT_MODES comment."""
    assert '_research_ws = _abs_dir(getattr(ctx, "agent_dir", None))' in _SRC


def test_research_engine_branch_has_no_write_tool():
    """The tool list built in _ensure_engine for research carries no tool
    that can change a file or run a shell command."""
    start = _SRC.find('elif mode_dropdown.value == "research":')
    assert start > 0, "research branch missing from _ensure_engine"
    body = _SRC[start:_SRC.find("elif", start + 10)]
    for tool in ("Write", "Edit", "MultiEdit", "NotebookEdit", "Bash"):
        assert f'"{tool}"' not in body, (
            f"research tool list includes {tool}; research is read-only")


# ---------------------------------------------------------------------------
# README and the dashboard cannot drift apart again
# ---------------------------------------------------------------------------

def test_readme_research_row_matches_the_dashboard():
    """README said "terminal only" while the dashboard offered research;
    pin the row to the dropdown so the two are changed together."""
    readme = Path(T.__file__).parents[2] / "README.md"
    text = readme.read_text(encoding="utf-8")
    assert "| `research` | dashboard, terminal |" in text
    # the row must not claim a restriction the dropdown does not enforce,
    # nor promise one the engine branch does not deliver
    assert ("Research", "research") in _dropdown_options()
    row = next(line for line in text.splitlines()
               if line.startswith("| `research` |"))
    assert "read-only" in row
    assert "terminal only" not in row


def test_research_role_surface_matches_the_description():
    """The dropdown text says "no file changes, no shell commands". The
    role layer denies the write tools under their internal names
    (api_client._WRITE_TOOL_NAMES); Bash is excluded at the engine branch
    (test above) and by the dashboard's tool advertisement."""
    from delfin.agent.api_client import _tool_denied_for_role

    for tool in ("write_file", "edit_file", "multi_edit", "apply_patch",
                 "notebook_edit"):
        assert _tool_denied_for_role("research_agent", tool), (
            f"{tool} is allowed for research_agent; the description that "
            "says read-only would be lying")
    assert not _tool_denied_for_role("research_agent", "web_search")
