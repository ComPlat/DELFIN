"""A ``read_file(around=...)`` pattern anchored with ``^`` matches a line start.

Both halves of this file come from real runs of the product on 2026-11-26
(w60 package r2), where the model reached for the two pattern spellings
that ``grep`` accepts and got an answer that blamed the wrong thing:

* ``around="^def check"`` -- the agent's own summary says it read the
  definition and then *"re-read that region by line offset instead of
  retrying the pattern"*. ``re.compile`` without ``re.MULTILINE`` makes
  ``^`` anchor at the start of the *file*, so a pattern that plainly
  exists matched nothing and the tool answered "pattern ... not found --
  no window to show", which reads as "there is no such definition".
* ``around="cage"`` on a file where it exists -- answered correctly, and
  is the control here.

The question a caller asks with ``around=`` is always "where is this in
the file", never "is the file exactly this", so every line of the file is
a legitimate anchor position. MULTILINE is what makes ``^``/``$`` mean
that; it changes nothing for a pattern without anchors.

Red on the parent commit for the first test only, with the reason in its
name: the window is not found, so ``read_file`` returns an error string
instead of a numbered line. The control case passes on both.
"""

import pytest

from delfin.agent.api_client import _DocToolExecutor

SOURCE = "\n".join([
    "import os",
    "",
    "",
    "def helper(value):",
    "    return value",
    "",
    "",
    "def check(schema, arguments):",
    "    return True",
])


def _perms(ws):
    from delfin.agent.api_client import KitToolPermissions
    perms = KitToolPermissions(workspace=str(ws))
    perms.mode = "acceptEdits"
    perms.task_session_id = "line-anchor-test"
    return perms


def _read(ws, **arguments):
    """Through ``_dispatch`` -- the entry a model's tool call reaches."""
    return _DocToolExecutor()._dispatch("read_file", arguments, _perms(ws))


@pytest.fixture
def ws(tmp_path):
    d = tmp_path / "ws"
    d.mkdir()
    (d / "mod.py").write_text(SOURCE, encoding="utf-8")
    return d


def test_a_caret_anchored_pattern_finds_the_definition_line(ws):
    """``^def check`` is the spelling grep accepts; it must find line 8."""
    out = _read(ws, path="mod.py", around="^def check", context=1)
    assert "not found" not in out, out
    # The window is reported 1-based and grep-compatible.
    assert out.splitlines()[0].startswith("7  "), out
    assert "def check(schema, arguments):" in out


def test_a_dollar_anchored_pattern_finds_the_end_of_a_line(ws):
    """``return value$`` must match the line, not only end-of-file."""
    out = _read(ws, path="mod.py", around="return value$", context=0)
    assert "not found" not in out, out
    assert "return value" in out


def test_an_unanchored_pattern_still_behaves_as_before(ws):
    """Control: no anchors, so MULTILINE must change nothing at all."""
    out = _read(ws, path="mod.py", around="def helper", context=1)
    assert "not found" not in out, out
    assert "def helper(value):" in out
