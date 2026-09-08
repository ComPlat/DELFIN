"""The approval dialog shows the diff, or it asks you to approve a name.

Before a file-mutating tool runs, the dashboard renders a unified diff
under the approval prompt so the decision is about the CHANGE and not
about the word "write_file". The renderer recognised the tool by name and
knew three: ``Edit``, ``Write`` and ``multi_edit`` — two of the CLI
backend's spellings plus one of the OpenAI-compatible backend's, which is
the shape of a list somebody extended once and did not finish.

On the KIT and Ollama backends the two most common write tools are
``write_file`` and ``edit_file``. Neither matched, so on the backend
these models are served from the user was asked to approve a file change
with nothing shown. Approving blind is the failure this preview exists to
prevent.

The argument names were never the problem — ``path``, ``old_string``,
``new_string``, ``content`` and ``edits`` are the same on both sides, and
the reader already accepted ``path`` as well as ``file_path``. Only the
tool name was.
"""

from __future__ import annotations

import pytest

from delfin.dashboard.tab_agent import _compute_approval_diff


@pytest.fixture
def target(tmp_path):
    p = tmp_path / "bookmark_store.py"
    p.write_text("def load():\n    return []\n", encoding="utf-8")
    return p


def _detail(tool: str, **inp) -> str:
    return repr({"tool_name": tool, "tool_input": inp})


@pytest.mark.parametrize("tool", ["Write", "write_file",
                                  "mcp__kit-coding__write_file"])
def test_a_whole_file_write_is_shown(target, tool):
    out = _compute_approval_diff(_detail(
        tool, path=str(target), content="def load():\n    return [1]\n"))
    assert "```diff" in out, f"{tool}: nothing was shown"
    assert "+    return [1]" in out


@pytest.mark.parametrize("tool", ["Edit", "edit_file",
                                  "mcp__kit-coding__edit_file"])
def test_a_substring_edit_is_shown(target, tool):
    out = _compute_approval_diff(_detail(
        tool, path=str(target), old_string="return []",
        new_string="return [1]"))
    assert "```diff" in out, f"{tool}: nothing was shown"
    assert "+    return [1]" in out


@pytest.mark.parametrize("tool", ["multi_edit", "MultiEdit",
                                  "mcp__kit-coding__multi_edit"])
def test_several_edits_at_once_are_shown(target, tool):
    out = _compute_approval_diff(_detail(
        tool, path=str(target),
        edits=[{"old_string": "def load", "new_string": "def load_all"},
               {"old_string": "return []", "new_string": "return [1]"}]))
    assert "```diff" in out, f"{tool}: nothing was shown"
    assert "+def load_all" in out


def test_a_patch_is_its_own_preview(target):
    """apply_patch carries the diff in its argument. Recomputing one from
    the file would be work to arrive back where we started."""
    patch = ("--- a/bookmark_store.py\n+++ b/bookmark_store.py\n"
             "@@ -1,2 +1,2 @@\n def load():\n-    return []\n+    return [1]\n")
    out = _compute_approval_diff(_detail("apply_patch", diff=patch))
    assert "```diff" in out
    assert "+    return [1]" in out


def test_a_tool_that_changes_no_file_shows_nothing(target):
    for tool in ("bash", "read_file", "Read", "search_docs"):
        assert _compute_approval_diff(
            _detail(tool, path=str(target))) == "", tool


def test_an_edit_whose_anchor_is_gone_says_so_rather_than_lying(target):
    out = _compute_approval_diff(_detail(
        "edit_file", path=str(target), old_string="nicht vorhanden",
        new_string="x"))
    assert "unavailable" in out


def test_a_write_that_changes_nothing_says_that(target):
    out = _compute_approval_diff(_detail(
        "write_file", path=str(target), content=target.read_text()))
    assert "no change" in out


def test_a_malformed_detail_is_not_an_exception(target):
    for detail in ("", "not a dict", "{unclosed", repr({"tool_name": "Write"})):
        assert _compute_approval_diff(detail) == ""


# ---------------------------------------------------------------------------
# ...and the line above it says what the tool would do
# ---------------------------------------------------------------------------

def _readable(tool: str, **inp) -> str:
    from delfin.dashboard.tab_agent import _format_tool_description
    return _format_tool_description(repr({"tool_name": tool,
                                          "tool_input": inp}))


@pytest.mark.parametrize("tool", ["Write", "write_file",
                                  "mcp__kit-coding__write_file"])
def test_the_approval_line_names_the_file_being_written(tool):
    """The generic fallback rendered `write_file: path=… content=…` with
    the content cut at 80 characters — a wall of the file's own text
    where the sentence should say what is about to happen."""
    assert _readable(tool, path="/a/b/bookmark_store.py",
                     content="x" * 500) == "Write file: bookmark_store.py"


@pytest.mark.parametrize("tool", ["Edit", "edit_file"])
def test_the_approval_line_shows_both_sides_of_an_edit(tool):
    out = _readable(tool, path="x.py", old_string="return []",
                    new_string="return [1]")
    assert "return []" in out and "return [1]" in out


def test_a_multi_edit_says_how_many_changes():
    out = _readable("multi_edit", path="x.py",
                    edits=[{"old_string": "a"}, {"old_string": "b"}])
    assert "2 change" in out


@pytest.mark.parametrize("tool", ["Bash", "bash"])
def test_a_shell_call_shows_the_command(tool):
    assert _readable(tool, command="pytest -q tests/") == "Run: pytest -q tests/"


def test_an_unknown_tool_still_says_something():
    out = _readable("some_new_tool", subject="x")
    assert "some_new_tool" in out and "x" in out


# ---------------------------------------------------------------------------
# ...and the spinner says what is happening while it happens
# ---------------------------------------------------------------------------

def _label(tool: str, **parsed) -> str:
    from delfin.dashboard.tab_agent import _tool_activity_label
    return _tool_activity_label(tool, parsed)


@pytest.mark.parametrize("tool,expected", [
    ("write_file", "Writing"), ("Write", "Writing"),
    ("edit_file", "Editing"), ("Edit", "Editing"),
    ("multi_edit", "Editing"), ("apply_patch", "Patching"),
    ("read_file", "Reading"), ("Read", "Reading"),
])
def test_the_spinner_says_what_is_happening_on_either_backend(tool, expected):
    """The one place a long turn tells the user what it is doing. The
    table held the CLI backend's seven names, so on KIT every line read
    "Running mcp__kit-coding__write_file…" — the transport's vocabulary
    where a sentence belongs."""
    assert _label(tool, path="/a/b/bookmark_store.py") == (
        f"{expected} bookmark_store.py...")


def test_the_mcp_namespace_does_not_hide_the_verb():
    assert _label("mcp__kit-coding__write_file", path="x.py") == "Writing x.py..."


@pytest.mark.parametrize("tool", ["bash", "Bash"])
def test_a_shell_line_shows_the_command(tool):
    assert _label(tool, command="pytest -q").startswith("$ pytest -q")


def test_a_tool_nobody_mapped_still_reads_as_a_sentence():
    assert _label("some_new_tool") == "Running some_new_tool..."


def test_a_write_without_a_path_still_names_the_verb():
    assert _label("write_file") == "Writing..."
