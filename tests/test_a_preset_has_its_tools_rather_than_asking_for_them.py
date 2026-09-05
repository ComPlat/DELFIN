"""A sub-agent preset's tool surface, and the gate that actually holds it.

The ``explore`` preset used to say, in its system prompt: "You may use
read_file, grep_file, list_files, web_search, web_fetch, notebook_read. You
may NOT edit, write, or run bash." That was a sentence. Nothing below the
model read it, ``SubagentPreset`` carried no tool list, and the only thing
that bound was ``mode="plan"`` refusing writes at execution time — a weaker
and different statement, and for a preset whose mode is not ``plan``, no
statement at all.

``SubagentPreset.tools`` replaces the sentence with a list that
``_derive_perms`` intersects into the child's ``session_allowed_tools`` —
the field the EXECUTOR is handed on every call. So the tests below call the
executor directly, the way an MCP backend, a replayed tool call or a nested
sub-agent would reach it, with no advertised list anywhere in the path. A
preset that only filtered the schema sent to the model would pass a test
that checked the schema and still run the tool. That is the defect, not the
fix for it.

Three semantics are silent when they are backwards, so each is pinned here:
an empty list is NO restriction (the inverse disables every preset that
declares nothing), a preset can only ever narrow its parent (the inverse is
a privilege escalation via a markdown file on disk), and an intersection
that comes out empty leaves the child able to call nothing (the inverse
hands the most freedom to the child that earned the least).
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from delfin.agent import subagents as S
from delfin.agent.api_client import (
    KitToolPermissions,
    ToolSurfaceContext,
    _DOC_TOOLS_OPENAI,
    _DocToolExecutor,
    tool_unavailable_reason,
    unknown_tool_names,
)
from delfin.agent.subagents import (
    SUBAGENT_PRESETS,
    SubagentPreset,
    _derive_perms,
    _narrow_allowed_tools,
    _parse_preset_tools,
)

# The built-ins that declare a list. ``general-purpose`` declares none on
# purpose and is covered by its own tests below.
RESTRICTED_PRESETS = ("explore", "plan", "code-reviewer")

# Outside every read-only preset's list, and the three the deleted prose
# named: "You may NOT edit, write, or run bash".
OUTSIDE = (
    ("bash", {"command": "echo escaped"}),
    ("write_file", {"path": "new.txt", "content": "escaped"}),
    ("edit_file", {"path": "a.txt", "old_string": "hello", "new_string": "x"}),
)


@pytest.fixture
def workspace(tmp_path) -> Path:
    (tmp_path / "a.txt").write_text("hello\n", encoding="utf-8")
    return tmp_path


def _parent(workspace, **over) -> KitToolPermissions:
    base = dict(workspace=workspace, mode="default")
    base.update(over)
    return KitToolPermissions(**base)


def _child_of(preset: SubagentPreset, parent) -> KitToolPermissions:
    """The permissions ``run_subagent`` builds for *preset* — same call."""
    return _derive_perms(parent, preset.mode, preset_tools=preset.tools)


def _error_of(result: str) -> str:
    try:
        return json.loads(result).get("error", "")
    except Exception:
        return ""


# ---------------------------------------------------------------------------
# A tool outside the list is REFUSED. By the executor, when called.
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("preset_name", RESTRICTED_PRESETS)
@pytest.mark.parametrize("tool,args", OUTSIDE, ids=[t for t, _ in OUTSIDE])
def test_a_tool_outside_the_preset_list_is_refused_by_the_executor(
    preset_name, tool, args, workspace
):
    """Called directly, with no advertised list in the path.

    Hiding a tool from the model is not refusing it: an MCP backend, a
    replayed call or a nested sub-agent reaches this same executor without
    ever consulting the schema. If this passed only because the tool was
    unadvertised, the tool would still run here.
    """
    preset = SUBAGENT_PRESETS[preset_name]
    child = _child_of(preset, _parent(workspace))

    result = _DocToolExecutor().execute(tool, args, child)
    error = _error_of(result)
    assert error, f"{preset_name}: {tool} was not refused: {result[:200]}"
    assert "--allowed-tools" in error, (
        f"{preset_name}: {tool} was refused, but not by the preset's tool "
        f"list: {error}"
    )
    # And the side effect did not happen.
    assert not (workspace / "new.txt").exists()
    assert (workspace / "a.txt").read_text(encoding="utf-8") == "hello\n"


@pytest.mark.parametrize("preset_name", RESTRICTED_PRESETS)
def test_the_list_binds_where_the_mode_does_not(preset_name, workspace):
    """The reason a list is worth having at all.

    ``mode="plan"`` was the only thing that used to bind, and it binds only
    in plan mode. Here the child is derived at ``bypassPermissions`` — every
    mode-based refusal is off — and the tool is still refused, because the
    refusal comes from the role's list rather than from the mode.
    """
    preset = SUBAGENT_PRESETS[preset_name]
    child = _derive_perms(
        _parent(workspace, mode="bypassPermissions"),
        "bypassPermissions",
        preset_tools=preset.tools,
    )
    assert child.mode == "bypassPermissions"

    result = _DocToolExecutor().execute("bash", {"command": "echo escaped"}, child)
    assert "--allowed-tools" in _error_of(result), result[:200]


# ---------------------------------------------------------------------------
# A tool inside the list still works. A list that broke the preset's own
# work would be a worse defect than the one it fixes.
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("preset_name", RESTRICTED_PRESETS)
def test_a_tool_inside_the_preset_list_still_works(preset_name, workspace):
    preset = SUBAGENT_PRESETS[preset_name]
    child = _child_of(preset, _parent(workspace))
    assert "hello" in _DocToolExecutor().execute(
        "read_file", {"path": "a.txt"}, child)


@pytest.mark.parametrize("preset_name", RESTRICTED_PRESETS)
@pytest.mark.parametrize(
    "tool",
    # Measured, not guessed: the tools 111 recorded ``explore`` runs actually
    # called (~/.delfin/subagent_sessions). ``find_references``,
    # ``find_definition``, ``read_document`` and ``exit_plan_mode`` are the
    # four the old prose never named — a list copied from that sentence
    # would have refused 67 calls that worked.
    ["read_file", "grep_file", "list_files", "find_references",
     "find_definition", "read_document", "exit_plan_mode"],
)
def test_every_tool_the_recorded_runs_called_is_still_permitted(
    preset_name, tool, workspace
):
    preset = SUBAGENT_PRESETS[preset_name]
    child = _child_of(preset, _parent(workspace))
    from delfin.agent.api_client import _tool_denied_for_session
    assert _tool_denied_for_session(child, tool) is None
    # The MCP-namespaced form of the same tool is how the recorded runs
    # actually called it, and it must not be a route around the list either.
    assert _tool_denied_for_session(child, f"mcp__delfin-docs__{tool}") is None


# ---------------------------------------------------------------------------
# A preset may only ever NARROW its parent.
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("preset_name", RESTRICTED_PRESETS)
def test_a_preset_cannot_give_back_a_tool_the_parent_denied(
    preset_name, workspace
):
    """The escalation this guards against.

    Presets are loaded from markdown files on disk. If naming a tool in one
    could re-admit it, dropping a file into ``~/.delfin/subagents/`` would
    be a way around ``--disallowed-tools``.
    """
    preset = SUBAGENT_PRESETS[preset_name]
    assert "read_file" in preset.tools, "the premise: the preset asks for it"

    parent = _parent(workspace, session_denied_tools=frozenset({"read_file"}))
    child = _child_of(preset, parent)

    result = _DocToolExecutor().execute("read_file", {"path": "a.txt"}, child)
    assert "hello" not in result
    assert "--disallowed-tools" in _error_of(result)


@pytest.mark.parametrize("preset_name", RESTRICTED_PRESETS)
def test_a_preset_cannot_widen_the_parents_own_allow_list(
    preset_name, workspace
):
    """The parent allowed exactly one tool. The preset names dozens.

    The child ends up with the intersection — not the union, and not the
    preset's list on its own.
    """
    preset = SUBAGENT_PRESETS[preset_name]
    parent = _parent(workspace, session_allowed_tools=frozenset({"read_file"}))
    child = _child_of(preset, parent)

    assert set(child.session_allowed_tools) == {"read_file"}
    ex = _DocToolExecutor()
    assert "hello" in ex.execute("read_file", {"path": "a.txt"}, child)
    # grep_file is on the preset's list but was never the parent's to give.
    assert "--allowed-tools" in _error_of(
        ex.execute("grep_file", {"path": "a.txt", "pattern": "hello"}, child))


def test_an_empty_intersection_permits_nothing_rather_than_everything(
    workspace
):
    """The one case where "empty means unrestricted" would invert the rule.

    A parent allowing only ``bash`` and a read-only preset share no tool. If
    the empty intersection were stored as an empty frozenset, every gate
    downstream would read it as "no restriction" and the most-constrained
    child in the system would become the least.
    """
    preset = SUBAGENT_PRESETS["explore"]
    parent = _parent(workspace, session_allowed_tools=frozenset({"bash"}))
    child = _child_of(preset, parent)

    assert child.session_allowed_tools, "an empty set here would mean 'anything'"
    ex = _DocToolExecutor()
    for tool, args in (("read_file", {"path": "a.txt"}),
                       ("bash", {"command": "echo escaped"})):
        assert "--allowed-tools" in _error_of(ex.execute(tool, args, child))


# ---------------------------------------------------------------------------
# An empty list is NO restriction — not "no tools".
# ---------------------------------------------------------------------------

def test_general_purpose_declares_no_list_and_is_therefore_unrestricted(
    workspace
):
    """Getting this backwards would disable the preset entirely.

    ``general-purpose`` IS the full tool surface, so any list on it could
    only be a stale copy of the catalogue that silently drops each newly
    added tool.
    """
    preset = SUBAGENT_PRESETS["general-purpose"]
    assert preset.tools == ()

    child = _child_of(preset, _parent(workspace))
    assert not child.session_allowed_tools

    ex = _DocToolExecutor()
    assert "hello" in ex.execute("read_file", {"path": "a.txt"}, child)
    assert not _error_of(ex.execute("bash", {"command": "echo ok"}, child))


def test_a_preset_without_a_list_leaves_the_parents_lists_untouched(workspace):
    """No list means the preset abstains, not that it overwrites."""
    parent = _parent(
        workspace,
        session_allowed_tools=frozenset({"read_file", "bash"}),
        session_denied_tools=frozenset({"bash"}),
    )
    child = _derive_perms(parent, "default", preset_tools=())
    assert set(child.session_allowed_tools) == {"read_file", "bash"}
    assert set(child.session_denied_tools) == {"bash"}
    assert _narrow_allowed_tools(parent, ()) is None


def test_derive_perms_keeps_its_old_signature_for_callers_without_a_preset(
    workspace
):
    """``preset_tools`` is optional: the other call sites predate it."""
    child = _derive_perms(_parent(workspace), "plan")
    assert child.mode == "plan"
    assert not child.session_allowed_tools


# ---------------------------------------------------------------------------
# The list is the only place the rule is stated.
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("preset_name", sorted(SUBAGENT_PRESETS))
def test_no_preset_prompt_still_enumerates_its_tools_in_prose(preset_name):
    """A rule written in two places drifts, and prose is the copy that
    cannot be enforced.

    The measurement is what makes this concrete rather than tidy: the
    ``explore`` prose named six tools, four of which the recorded runs never
    called, and omitted four it called 67 times between them. It had been
    wrong for a while, and nothing could tell.
    """
    prompt = SUBAGENT_PRESETS[preset_name].system_prompt.lower()
    for phrase in ("you may use", "you may not", "available tools:",
                   "do not edit", "do not use"):
        assert phrase not in prompt, (
            f"{preset_name}: '{phrase}' restates in prose what "
            f"SubagentPreset.tools now enforces"
        )


@pytest.mark.parametrize("preset_name", sorted(SUBAGENT_PRESETS))
def test_every_name_in_a_preset_list_is_a_real_tool(preset_name):
    """A misspelt name narrows the preset further than anyone meant, and
    does it silently — the same failure ``unknown_tool_names`` exists to
    report for the CLI flags."""
    tools = SUBAGENT_PRESETS[preset_name].tools
    assert unknown_tool_names(tools, _DOC_TOOLS_OPENAI) == []


# ---------------------------------------------------------------------------
# Markdown presets can declare a list too.
# ---------------------------------------------------------------------------

def test_the_packaged_markdown_preset_declares_its_tools_in_frontmatter():
    preset = SUBAGENT_PRESETS.get("chemistry-reviewer")
    assert preset is not None
    assert "read_file" in preset.tools and "search_calcs" in preset.tools
    assert "bash" not in preset.tools and "write_file" not in preset.tools


def test_a_user_markdown_preset_can_declare_a_list(tmp_path, monkeypatch):
    d = tmp_path / ".delfin" / "subagents"
    d.mkdir(parents=True)
    (d / "narrow_subagent.md").write_text(
        "---\nname: narrow\nmode: plan\ntools: read_file, grep_file\n---\n"
        "Look at one file and report.\n",
        encoding="utf-8",
    )
    monkeypatch.setattr(Path, "home", classmethod(lambda cls: tmp_path))
    try:
        registry = S.reload_subagent_presets()
        assert registry["narrow"].tools == ("read_file", "grep_file")
    finally:
        monkeypatch.undo()
        S.reload_subagent_presets()


@pytest.mark.parametrize(
    "raw,expected",
    [
        ("read_file, grep_file", ("read_file", "grep_file")),
        ("[read_file, grep_file]", ("read_file", "grep_file")),
        ("['read_file', \"grep_file\"]", ("read_file", "grep_file")),
        ("mcp__kit-coding__find_references",
         ("mcp__kit-coding__find_references",)),
        ("", ()),
        (None, ()),
        (["read_file"], ()),
    ],
)
def test_the_frontmatter_tools_value_is_read_in_the_shapes_people_write(
    raw, expected
):
    """An unreadable value reads as "no list", i.e. no restriction.

    The opposite default would let one bad frontmatter line reduce a preset
    to nothing, and a sub-agent that can call no tool fails in a way nobody
    traces back to a markdown file.
    """
    assert _parse_preset_tools(raw) == expected


# ---------------------------------------------------------------------------
# The advertised surface follows the refusal — it is not a second policy.
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("preset_name", RESTRICTED_PRESETS)
def test_the_advertised_surface_narrows_from_the_same_field(
    preset_name, workspace
):
    """Advertising is derived from the executor's decision, never instead of
    it: the model is not offered a tool whose every call would come back a
    refusal."""
    preset = SUBAGENT_PRESETS[preset_name]
    child = _child_of(preset, _parent(workspace))
    ctx = ToolSurfaceContext(
        session_allowed=child.session_allowed_tools,
        session_denied=child.session_denied_tools,
    )
    assert tool_unavailable_reason("bash", ctx) is not None
    assert tool_unavailable_reason("write_file", ctx) is not None
    assert tool_unavailable_reason("read_file", ctx) is None
