"""The trace said every tool call succeeded, including the refusals.

`ok=True` was hardcoded on the tool_result branch, and the only writer of
`ok=False` was the `permission_denied` event -- which only the CLI backend
emits. On the API and KIT backends a gate refusal comes back as an
ordinary tool_result carrying `{"error": ...}`, so it was recorded green.

Five surfaces read that flag, and all five were wrong together: `/trace`
printed a checkmark beside the blocked call, `aggregate_tools` reported an
error rate of zero permanently, `/agents tools` showed a `err%` column
that could only ever say 0%, the live panel rendered a green tick, and the
bug report embedded all of it. A user asking "why did it refuse to do
anything for an hour" received a clean list of successes.

Whether a call worked is a property of its result. It is now read, not
asserted.
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


# ---------------------------------------------------------------------------
# The predicate
# ---------------------------------------------------------------------------

def test_a_json_error_is_a_failure():
    failed, reason = E.AgentEngine._tool_result_failed(
        '{"error": "command rejected by deny-pattern"}')
    assert failed
    assert "deny-pattern" in reason


def test_ordinary_output_is_not_a_failure():
    for text in ("ok\n", "", "1  import os", "No matches found.",
                 '{"status": "ok"}'):
        assert E.AgentEngine._tool_result_failed(text) == (False, ""), text


def test_a_command_that_ran_and_exited_non_zero_is_not_a_tool_failure():
    """A successful call reporting a failed command. Conflating the two
    would make the error rate measure the user's code, not the tooling."""
    failed, _ = E.AgentEngine._tool_result_failed(
        "5 failed, 2 passed\nexit code 1")
    assert failed is False


def test_a_malformed_error_object_is_still_a_failure():
    failed, reason = E.AgentEngine._tool_result_failed('{"error": ')
    assert failed and reason


def test_the_reason_is_bounded():
    failed, reason = E.AgentEngine._tool_result_failed(
        '{"error": "%s"}' % ("x" * 5000))
    assert failed and len(reason) <= 300


# ---------------------------------------------------------------------------
# ...and it reaches the trace
# ---------------------------------------------------------------------------

def _traced(engine, tool_output):
    seen = []

    def spy(name, output="", *, ok=True, error=""):
        seen.append({"name": name, "ok": ok, "error": error})

    engine._record_tool_trace = spy

    def stream(**kw):
        yield StreamEvent(type="tool_use", tool_name="mcp__kit-coding__bash",
                          tool_input='{"command": "rm -rf build"}')
        yield StreamEvent(type="tool_result", tool_name="mcp__kit-coding__bash",
                          tool_output=tool_output)
        yield StreamEvent(type="text_delta", text="done")

    engine.client.stream_message = MagicMock(side_effect=stream)
    engine.stream_response("clean up")
    return seen


def test_a_refused_call_is_traced_as_failed(eng):
    seen = _traced(eng, '{"error": "command rejected by deny-pattern"}')
    assert seen and seen[0]["ok"] is False
    assert "deny-pattern" in seen[0]["error"]


def test_a_working_call_is_still_traced_as_ok(eng):
    seen = _traced(eng, "removed 3 files\n")
    assert seen and seen[0]["ok"] is True
    assert seen[0]["error"] == ""


def test_the_permission_denied_event_still_works(eng):
    """The CLI backend's own path must not regress."""
    seen = []
    eng._record_tool_trace = lambda name, output="", *, ok=True, error="": \
        seen.append({"ok": ok, "error": error})

    def stream(**kw):
        yield StreamEvent(type="permission_denied", tool_name="bash")
        yield StreamEvent(type="text_delta", text="done")

    eng.client.stream_message = MagicMock(side_effect=stream)
    eng.stream_response("do it")
    assert seen and seen[0]["ok"] is False
