"""A long prompt can be written in an editor.

Nothing in this codebase read $VISUAL or $EDITOR, so a multi-paragraph
prompt had to be typed into a single-line reader or pasted. `/editor`
opens the configured editor on a scratch file and sends what comes back
as the next prompt.

Two things this file is careful about, both of which have cost a green
run before.

NO EDITOR IS INSTALLED ON A RUNNER. Every test here writes its own
"editor" -- a small python script -- and points $EDITOR at it through
monkeypatch. A test that assumed `vi` exists passes on a developer box
and fails where it matters.

A CANCEL MUST NOT SEND. A non-zero exit, an empty file and an untouched
template are all "no", and each gets its own test: sending a draft the
user backed out of puts words in their mouth, and the model cannot tell
that from a decision.
"""

from __future__ import annotations

import json
import shlex
import subprocess
import sys
from pathlib import Path

import pytest

from delfin.agent import repl_commands as rc


# ---------------------------------------------------------------------------
# An editor that is a python script, so no binary has to exist
# ---------------------------------------------------------------------------

#: Recorded by every fake editor before it does anything else: the path it
#: was handed, that path's mode, and what was already in it. One test
#: asserts each, and they are facts only the editor process can observe --
#: the file is gone by the time the handler returns.
_PREAMBLE = (
    "import json, os, pathlib, sys\n"
    "target = pathlib.Path(sys.argv[1])\n"
    "pathlib.Path(%r).write_text(json.dumps({\n"
    "    'target': str(target),\n"
    "    'mode': oct(os.stat(target).st_mode & 0o777),\n"
    "    'seen': target.read_text(encoding='utf-8'),\n"
    "}), encoding='utf-8')\n"
)


def _fake_editor(tmp_path: Path, monkeypatch, body: str = "") -> Path:
    """Point $EDITOR at a script that does *body*, and return its note file.

    The interpreter running the tests is the one that runs the script, so
    this works wherever pytest works. Both paths are quoted because the
    handler splits the variable the way a shell would.
    """
    note = tmp_path / "editor_note.json"
    script = tmp_path / "fake_editor.py"
    script.write_text(_PREAMBLE % str(note) + body, encoding="utf-8")
    monkeypatch.setenv("EDITOR", "%s %s" % (shlex.quote(sys.executable),
                                            shlex.quote(str(script))))
    return note


@pytest.fixture(autouse=True)
def _no_inherited_editor(monkeypatch):
    """Whatever the developer has set must not decide these tests.

    A machine with VISUAL exported would otherwise run the real editor
    for the tests that only set EDITOR, and hang the suite on a terminal
    nobody is watching.
    """
    monkeypatch.delenv("VISUAL", raising=False)
    monkeypatch.delenv("EDITOR", raising=False)


def _run(tmp_path: Path, args: str = "") -> rc.CommandResult:
    return rc.BUILTINS["/editor"].handler(
        rc.ReplContext(engine=None, workspace=tmp_path), args)


# ---------------------------------------------------------------------------
# The happy path
# ---------------------------------------------------------------------------

def test_what_was_written_becomes_the_next_prompt(tmp_path, monkeypatch):
    _fake_editor(tmp_path, monkeypatch,
                 "target.write_text('first paragraph\\n\\n"
                 "second paragraph\\n', encoding='utf-8')\n")

    result = _run(tmp_path)

    assert result.handled is True
    assert result.prompt == "first paragraph\n\nsecond paragraph"
    assert result.quit is False and result.clear is False


def test_the_paragraphs_survive_the_trip(tmp_path, monkeypatch):
    """The blank line is the whole point of writing this in an editor.

    `expand_at_references` re-joins a line on single spaces, so routing
    the draft through it would flatten exactly what /editor exists to
    make possible. This asserts it is not.
    """
    _fake_editor(tmp_path, monkeypatch,
                 "target.write_text('a\\n\\nb\\n', encoding='utf-8')\n")

    assert _run(tmp_path).prompt.count("\n\n") == 1


def test_the_argument_pre_fills_the_draft(tmp_path, monkeypatch):
    """A line already half typed is carried into the editor, not lost."""
    note = _fake_editor(tmp_path, monkeypatch,
                        "target.write_text(target.read_text(encoding='utf-8')"
                        " + 'and the rest\\n', encoding='utf-8')\n")

    result = _run(tmp_path, "half a sentence")

    assert json.loads(note.read_text())["seen"] == "half a sentence\n"
    assert result.prompt == "half a sentence\nand the rest"


def test_the_editor_that_opened_is_named_on_screen(tmp_path, monkeypatch):
    _fake_editor(tmp_path, monkeypatch,
                 "target.write_text('one\\n', encoding='utf-8')\n")

    assert sys.executable in _run(tmp_path).output


def test_visual_is_preferred_over_editor(tmp_path, monkeypatch):
    """The convention's order: VISUAL is the full-screen one."""
    _fake_editor(tmp_path, monkeypatch,
                 "target.write_text('from EDITOR\\n', encoding='utf-8')\n")
    visual = tmp_path / "visual_editor.py"
    visual.write_text("import pathlib, sys\n"
                      "pathlib.Path(sys.argv[1]).write_text("
                      "'from VISUAL\\n', encoding='utf-8')\n",
                      encoding="utf-8")
    monkeypatch.setenv("VISUAL", "%s %s" % (shlex.quote(sys.executable),
                                            shlex.quote(str(visual))))

    assert _run(tmp_path).prompt == "from VISUAL"


# ---------------------------------------------------------------------------
# Every way of saying no
# ---------------------------------------------------------------------------

def test_no_editor_configured_says_which_variable_to_set(tmp_path):
    result = _run(tmp_path)

    assert result.prompt == ""
    assert "VISUAL" in result.output and "EDITOR" in result.output


def test_no_editor_configured_does_not_guess_one(tmp_path):
    """No silent fallback: an editor nobody named is one nobody can leave."""
    out = _run(tmp_path).output
    assert "vi" not in out.split(), "a fallback editor must be named, not implied"


def test_a_non_zero_exit_sends_nothing(tmp_path, monkeypatch):
    """The draft is written AND the editor fails -- the draft still loses."""
    _fake_editor(tmp_path, monkeypatch,
                 "target.write_text('do not send this\\n', encoding='utf-8')\n"
                 "sys.exit(3)\n")

    result = _run(tmp_path)

    assert result.prompt == ""
    assert "3" in result.output


def test_an_untouched_draft_sends_nothing(tmp_path, monkeypatch):
    """Quitting without editing is a cancel, not a decision to send."""
    _fake_editor(tmp_path, monkeypatch)          # writes the note, edits nothing

    result = _run(tmp_path, "half a sentence")

    assert result.prompt == ""
    assert result.output.strip()


def test_an_emptied_draft_sends_nothing(tmp_path, monkeypatch):
    """Deleting everything is how a person cancels in an editor."""
    _fake_editor(tmp_path, monkeypatch,
                 "target.write_text('', encoding='utf-8')\n")

    result = _run(tmp_path)

    assert result.prompt == ""
    assert result.output.strip()


def test_a_draft_of_only_whitespace_sends_nothing(tmp_path, monkeypatch):
    _fake_editor(tmp_path, monkeypatch,
                 "target.write_text('  \\n\\n', encoding='utf-8')\n")

    assert _run(tmp_path).prompt == ""


def test_an_editor_binary_that_is_not_installed_is_one_line(tmp_path,
                                                            monkeypatch):
    """The common server case, and it must not be a traceback."""
    monkeypatch.setenv("EDITOR", "delfin-no-such-editor-9d3f")

    result = _run(tmp_path)

    assert result.prompt == ""
    assert "delfin-no-such-editor-9d3f" in result.output
    assert len(result.output.splitlines()) == 1


def test_an_editor_that_raises_is_one_line_too(tmp_path, monkeypatch):
    """Whatever `subprocess` throws, a handler in this table answers text."""
    _fake_editor(tmp_path, monkeypatch)

    def _boom(*_a, **_kw):
        raise RuntimeError("the terminal went away")

    monkeypatch.setattr(subprocess, "call", _boom)

    result = _run(tmp_path)

    assert result.prompt == ""
    assert "the terminal went away" in result.output


def test_a_variable_set_to_whitespace_counts_as_unset(tmp_path, monkeypatch):
    monkeypatch.setenv("EDITOR", "   ")

    assert "VISUAL" in _run(tmp_path).output


# ---------------------------------------------------------------------------
# The scratch file
# ---------------------------------------------------------------------------

def test_the_draft_is_readable_only_by_its_owner(tmp_path, monkeypatch):
    """It holds what the user is about to say to a model, in a shared /tmp."""
    note = _fake_editor(tmp_path, monkeypatch,
                        "target.write_text('secret\\n', encoding='utf-8')\n")

    _run(tmp_path)

    assert json.loads(note.read_text())["mode"] == "0o600"


def test_the_draft_is_a_markdown_file(tmp_path, monkeypatch):
    """So the editor picks prose wrapping instead of treating it as a blob."""
    note = _fake_editor(tmp_path, monkeypatch,
                        "target.write_text('x\\n', encoding='utf-8')\n")

    _run(tmp_path)

    assert json.loads(note.read_text())["target"].endswith(".md")


def test_the_draft_is_removed_afterwards(tmp_path, monkeypatch):
    note = _fake_editor(tmp_path, monkeypatch,
                        "target.write_text('sent\\n', encoding='utf-8')\n")

    _run(tmp_path)

    assert not Path(json.loads(note.read_text())["target"]).exists()


def test_the_draft_is_removed_when_the_editor_raises(tmp_path, monkeypatch):
    """The failure path is where an unsent prompt would be left behind."""
    _fake_editor(tmp_path, monkeypatch)
    seen: list[str] = []

    def _boom(argv, *_a, **_kw):
        seen.append(argv[-1])
        raise RuntimeError("no terminal")

    monkeypatch.setattr(subprocess, "call", _boom)

    _run(tmp_path)

    assert seen and not Path(seen[0]).exists()


def test_the_draft_is_removed_when_the_editor_fails(tmp_path, monkeypatch):
    note = _fake_editor(tmp_path, monkeypatch, "sys.exit(1)\n")

    _run(tmp_path)

    assert not Path(json.loads(note.read_text())["target"]).exists()


# ---------------------------------------------------------------------------
# The table and the help page
# ---------------------------------------------------------------------------

def test_the_command_is_in_the_table():
    cmd = rc.BUILTINS.get("/editor")
    assert cmd is not None, "a handler nobody registered is dead code"
    assert callable(cmd.handler)
    assert cmd.takes_args is True, "/editor <text> pre-fills the draft"
    assert cmd.summary.strip(), "a row with no summary renders a blank line"


def test_help_lists_the_command():
    """Help is generated from the table; an unlisted command is unfindable."""
    text = rc.dispatch("/help", rc.ReplContext()).output
    assert "/editor" in text
    assert rc.BUILTINS["/editor"].summary in text


def test_the_summary_stands_on_its_own():
    summary = rc.BUILTINS["/editor"].summary
    assert "EDITOR" in summary or "editor" in summary


def test_it_routes_as_a_command_and_not_as_a_path():
    assert rc.looks_like_command("/editor") is True
    assert rc.looks_like_command("/editor write a long thing") is True


def test_it_is_reached_through_the_table_not_a_special_case():
    import inspect

    assert "/editor" not in inspect.getsource(rc.dispatch), (
        "the table is the router; a name spelled out in dispatch is a "
        "second place to keep in step")
