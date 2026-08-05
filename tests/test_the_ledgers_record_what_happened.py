"""Evidence ledgers must record outcomes, not intentions.

The functional-claim guard asks: was this artifact ever exercised in this
session? Its answer comes from the exec ledger, which was written in the
`tool_use` branch -- at the moment the model DECIDED to run something,
before anything had happened. Everything that can go wrong in between was
therefore invisible to it: a permission gate refusing the command, a hook
blocking it, a missing interpreter, a traceback, a non-zero exit. The
agent runs `python app.py`, it dies with ModuleNotFoundError, and the
guard that exists to stop "the script works now" finds its evidence.

A ledger of attempts cannot ground a claim about results.
"""

from __future__ import annotations

from unittest.mock import MagicMock, patch

import pytest

from delfin.agent import engine as E
from delfin.agent.api_client import StreamEvent


@pytest.fixture
def eng(tmp_path):
    with patch("delfin.agent.engine.create_client", return_value=MagicMock()):
        return E.AgentEngine(repo_dir=tmp_path, backend="api", provider="kit",
                             model="kit.qwen3.5-397b-A17b", mode="solo")


def _run(engine, tool_name, tool_input, tool_output):
    def stream(**kw):
        yield StreamEvent(type="tool_use", tool_name=tool_name,
                          tool_input=tool_input)
        yield StreamEvent(type="tool_result", tool_name=tool_name,
                          tool_output=tool_output)
        yield StreamEvent(type="text_delta", text="done")

    engine.client.stream_message = MagicMock(side_effect=stream)
    engine.stream_response("do it")
    return engine._exec_commands_session


# ---------------------------------------------------------------------------
# What counts as having run
# ---------------------------------------------------------------------------

def test_a_command_that_succeeded_is_recorded(eng):
    ledger = _run(eng, "mcp__kit-coding__bash", '{"command": "python app.py"}',
                  "Hello from app\n")
    assert any("app.py" in c for c in ledger), ledger


def test_a_command_the_gate_refused_is_not_recorded(eng):
    ledger = _run(eng, "mcp__kit-coding__bash", '{"command": "python app.py"}',
                  '{"error": "command rejected by deny-pattern"}')
    assert ledger == [], ledger


def test_a_command_that_crashed_is_not_recorded(eng):
    ledger = _run(
        eng, "bash", '{"command": "python app.py"}',
        "Traceback (most recent call last):\n"
        "  File \"app.py\", line 1\nModuleNotFoundError: no module named x")
    assert ledger == [], ledger


def test_a_missing_interpreter_is_not_recorded(eng):
    ledger = _run(eng, "mcp__kit-coding__bash", '{"command": "python3.13 app.py"}',
                  "bash: python3.13: command not found")
    assert ledger == [], ledger


def test_a_role_denial_is_not_recorded(eng):
    ledger = _run(eng, "mcp__kit-coding__bash", '{"command": "python app.py"}',
                  "Tool 'bash' is not available to the 'dashboard_agent' role")
    assert ledger == [], ledger


def test_a_nonzero_exit_is_not_recorded(eng):
    ledger = _run(eng, "mcp__kit-coding__bash", '{"command": "pytest tests/"}',
                  "5 failed, 2 passed\nexit code 1")
    assert ledger == [], ledger


# ---------------------------------------------------------------------------
# The mechanism cannot leak or grow
# ---------------------------------------------------------------------------

def test_a_call_with_no_result_never_reaches_the_ledger(eng):
    """A turn that ends mid-loop must not leave an attempt behind."""
    def stream(**kw):
        yield StreamEvent(type="tool_use", tool_name="mcp__kit-coding__bash",
                          tool_input='{"command": "python app.py"}')
        yield StreamEvent(type="text_delta", text="stopped")

    eng.client.stream_message = MagicMock(side_effect=stream)
    eng.stream_response("do it")
    assert eng._exec_commands_session == []


def test_the_pending_list_is_bounded(eng):
    def stream(**kw):
        for i in range(300):
            yield StreamEvent(type="tool_use", tool_name="mcp__kit-coding__bash",
                              tool_input='{"command": "echo %d"}' % i)
        yield StreamEvent(type="text_delta", text="done")

    eng.client.stream_message = MagicMock(side_effect=stream)
    eng.stream_response("do it")
    assert len(eng._exec_pending) <= 64


def test_results_are_paired_with_their_own_call(eng):
    """Two commands, the second fails: only the first is evidence."""
    def stream(**kw):
        yield StreamEvent(type="tool_use", tool_name="mcp__kit-coding__bash",
                          tool_input='{"command": "python good.py"}')
        yield StreamEvent(type="tool_result", tool_name="mcp__kit-coding__bash",
                          tool_output="ok\n")
        yield StreamEvent(type="tool_use", tool_name="mcp__kit-coding__bash",
                          tool_input='{"command": "python bad.py"}')
        yield StreamEvent(type="tool_result", tool_name="mcp__kit-coding__bash",
                          tool_output="boom\nexit code 1")
        yield StreamEvent(type="text_delta", text="done")

    eng.client.stream_message = MagicMock(side_effect=stream)
    eng.stream_response("do it")
    joined = " ".join(eng._exec_commands_session)
    assert "good.py" in joined
    assert "bad.py" not in joined


def test_a_broken_ledger_write_does_not_end_the_turn(eng, monkeypatch):
    monkeypatch.setattr(
        E.AgentEngine, "_note_exec_command",
        lambda self, *a: (_ for _ in ()).throw(RuntimeError("boom")))
    assert _run(eng, "mcp__kit-coding__bash", '{"command": "echo hi"}', "hi") == []


# ---------------------------------------------------------------------------
# The observed-files ledger needs the read to have shown something
# ---------------------------------------------------------------------------

def _observed(fn_name, result, path="delfin/agent/engine.py"):
    from delfin.agent.api_client import _observe_read_files
    seen = set()
    _observe_read_files(seen, fn_name, {"path": path}, result)
    return seen


def test_a_grep_with_no_matches_does_not_ground_the_file():
    """One deliberately non-matching grep used to disarm the location
    guard for a whole file."""
    assert _observed("grep_file", "No matches found.") == set()


def test_a_read_that_returned_nothing_does_not_ground_the_file():
    assert _observed("read_file", "") == set()


def test_a_plain_text_refusal_does_not_ground_the_file():
    """Only JSON-shaped errors were excluded; a bare refusal was not."""
    assert _observed("read_file", "Permission denied") == set()


def test_a_read_that_returned_content_still_grounds_it():
    assert _observed("read_file", "1  import os\n2  import sys") == {
        "delfin/agent/engine.py"}


def test_a_grep_with_a_hit_still_grounds_it():
    assert _observed("grep_file", "delfin/agent/engine.py:12: import os") == {
        "delfin/agent/engine.py"}


def test_a_write_grounds_the_file_whatever_it_returns():
    """The agent produced the bytes, so it knows them."""
    assert _observed("write_file", "") == {"delfin/agent/engine.py"}
    assert _observed("edit_file", "ok") == {"delfin/agent/engine.py"}


def test_a_json_error_is_still_excluded():
    assert _observed("read_file", '{"error": "not found"}') == set()


# ---------------------------------------------------------------------------
# What the first version of these fixes got wrong
# ---------------------------------------------------------------------------

def test_a_failing_test_run_still_counts_as_having_run(eng):
    """The first marker list matched "traceback" anywhere in the output.

    A pytest run with one failure contains a traceback, so the run that
    demonstrably happened was discarded and the guard then told an
    accurate, appropriately-scoped report that the file had never been
    exercised. A guard that fires on correct answers teaches the model to
    write around it.
    """
    ledger = _run(
        eng, "mcp__kit-coding__bash", '{"command": "pytest tests/test_x.py"}',
        "tests/test_x.py ...F\nTraceback (most recent call last):\n"
        "AssertionError\n3 passed, 1 failed")
    assert any("pytest" in c for c in ledger), ledger


def test_ordinary_output_mentioning_a_denial_still_counts(eng):
    ledger = _run(
        eng, "mcp__kit-coding__bash", '{"command": "python scan.py"}',
        "scanning... 3 files skipped (Permission denied), 41 processed. Done.")
    assert any("scan.py" in c for c in ledger), ledger


def test_a_reported_exit_code_still_disqualifies(eng):
    ledger = _run(eng, "mcp__kit-coding__bash", '{"command": "python app.py"}',
                  "some output\nexit code 1")
    assert ledger == [], ledger


def test_a_program_that_merely_prints_about_exit_codes_counts(eng):
    ledger = _run(
        eng, "mcp__kit-coding__bash", '{"command": "python doc.py"}',
        "usage: the tool returns exit code 1 on failure\n"
        "and exit code 2 on a usage error\n"
        "Generating documentation...\nWrote doc.txt\nAll done.")
    assert any("doc.py" in c for c in ledger), ledger


def test_two_calls_in_one_round_are_paired_in_order(eng):
    """Taking the NEWEST pending entry was guaranteed wrong for a batch:
    the first command's failure popped the second command's entry, so the
    crashed one entered the ledger and the successful one was dropped."""
    def stream(**kw):
        yield StreamEvent(type="tool_use", tool_name="mcp__kit-coding__bash",
                          tool_input='{"command": "python crashed.py"}')
        yield StreamEvent(type="tool_use", tool_name="mcp__kit-coding__bash",
                          tool_input='{"command": "python worked.py"}')
        yield StreamEvent(type="tool_result", tool_name="mcp__kit-coding__bash",
                          tool_output="boom\nexit code 1")
        yield StreamEvent(type="tool_result", tool_name="mcp__kit-coding__bash",
                          tool_output="fine\n")
        yield StreamEvent(type="text_delta", text="done")

    eng.client.stream_message = MagicMock(side_effect=stream)
    eng.stream_response("do it")
    joined = " ".join(eng._exec_commands_session)
    assert "worked.py" in joined, joined
    assert "crashed.py" not in joined, joined


def test_a_read_under_an_argument_alias_grounds_the_file():
    """The executor accepts file_path/filename/file/target; the ledger
    read only "path", so a successful read under an alias grounded
    nothing -- and alias-using models are the weak ones."""
    from delfin.agent.api_client import _observe_read_files
    for key in ("path", "file_path", "filename", "file", "target"):
        seen = set()
        _observe_read_files(seen, "read_file", {key: "a/b.py"}, "1  import os")
        assert seen == {"a/b.py"}, key


def test_a_repo_wide_grep_grounds_the_files_it_showed():
    """grep with no path argument is the normal way to search a repo."""
    from delfin.agent.api_client import _observe_read_files
    seen = set()
    _observe_read_files(
        seen, "grep_file", {"pattern": "_note_exec"},
        "delfin/agent/engine.py:1654: self._note_exec_command(\n"
        "delfin/agent/engine.py:2014:     def _note_exec_command(")
    assert seen == {"delfin/agent/engine.py"}


def test_a_colon_in_matched_text_does_not_invent_a_path():
    from delfin.agent.api_client import _observe_read_files
    seen = set()
    _observe_read_files(seen, "grep_file", {"pattern": "x"},
                        "a/b.py:12: url = \"http://x/y:9\"")
    assert seen == {"a/b.py"}
