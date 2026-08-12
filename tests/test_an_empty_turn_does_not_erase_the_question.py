"""A turn that answers nothing must not take the question with it.

WHAT THE USER SAW. A long, carefully written task went in. The agent
thought for a while and produced nothing — no text, no error, just an
empty bubble. They typed "hm?" — and from then on the agent behaved as if
"hm?" were the whole assignment. Scrolling up, their task was still there
on screen; the model simply never saw it again.

WHAT ACTUALLY HAPPENED. ``stream_response`` appends the user message
before streaming and commits the answer after. The commit was written as
``if full_response: ... elif self._stop_requested: ...`` with no ``else``,
so a turn that produced no text and raised nothing fell through both
branches and left the user message in the history with nothing after it.
On the next turn ``_sanitize_messages`` restores alternation by resolving
two consecutive user messages in favour of the NEWEST — the follow-up
overwrote the task, verbatim and silently.

It takes no crash to get there. A reasoning model that streams its answer
on the reasoning channel emits thinking deltas only, and those are handed
to ``on_thinking`` and never enter ``chunks``. And the CLI backend's
``_read_turn`` fell off the end of its generator with neither an event nor
an exception in three shapes: the process still alive with its stdout
closed, a clean exit with no result event, and a non-zero exit with EMPTY
stderr — which is exactly what a SIGKILL from the out-of-memory killer
looks like.

WHY THE OBVIOUS FIX IS NOT THE FIX. Making ``_sanitize_messages`` keep the
OLDER message instead would trade this bug for its mirror image (every
genuine correction would then be ignored), and it would leave the caller
still receiving "" for a turn that failed — indistinguishable from a model
that legitimately had nothing to say. The empty turn has to be named where
it happens: take the unanswered message back out, say what was observed,
and record it so a turn log can tell it apart from a normal finish.

The exception path had already been repaired for exactly this hazard
(partial text is committed, marked as ended-in-error). The quiet path was
the same hazard with no crash to make it visible.
"""

from __future__ import annotations

import json
from unittest.mock import MagicMock, patch

import pytest

from delfin.agent import engine as E
from delfin.agent.api_client import CLIClient, StreamEvent


@pytest.fixture
def eng(tmp_path):
    with patch("delfin.agent.engine.create_client", return_value=MagicMock()):
        return E.AgentEngine(
            repo_dir=tmp_path, backend="api", provider="kit",
            model="kit.qwen3.5-397b-A17b", mode="solo")


# ---------------------------------------------------------------------------
# The engine: an answerless turn leaves the history clean
# ---------------------------------------------------------------------------

def test_a_turn_that_answers_only_in_reasoning_keeps_the_question(eng):
    """The incident: everything came back on the reasoning channel."""
    def all_reasoning(**kw):
        yield StreamEvent(type="thinking_delta", text="thinking hard " * 20)
        yield StreamEvent(type="message_delta", input_tokens=900,
                          output_tokens=400, cost_usd=0.03)

    def fine(**kw):
        yield StreamEvent(type="text_delta", text="ok")

    eng.client.stream_message = MagicMock(side_effect=all_reasoning)
    eng.stream_response("rewrite the ORCA input parser to accept blocks")

    assert eng.messages == [], (
        "the unanswered question was left in the history, where the next "
        "message overwrites it")

    eng.client.stream_message = MagicMock(side_effect=fine)
    eng.stream_response("hm?")
    joined = " ".join(str(m.get("content")) for m in eng.messages)
    assert "hm?" in joined
    assert "ORCA input parser" not in joined or joined.count("hm?") == 1


def test_the_orphan_is_gone_before_the_sanitiser_can_merge_it(eng):
    """The sanitiser is what destroys the question, and it runs at the
    START of the next turn — so there must be nothing left for it to
    merge by the time the next message arrives."""
    def nothing(**kw):
        yield StreamEvent(type="thinking_delta", text="...")

    def fine(**kw):
        yield StreamEvent(type="text_delta", text="ok")

    eng.client.stream_message = MagicMock(side_effect=nothing)
    eng.stream_response("refactor the ORCA parser")
    assert not any(m.get("role") == "user" for m in eng.messages), (
        "an unanswered question is sitting where the sanitiser will "
        "overwrite it with whatever the user types next")

    eng.client.stream_message = MagicMock(side_effect=fine)
    eng.stream_response("hm?")
    users = [m for m in eng.messages if m.get("role") == "user"]
    assert len(users) == 1 and users[0]["content"] == "hm?"


def test_the_caller_is_told_what_happened_not_handed_an_empty_string(eng):
    def nothing(**kw):
        yield StreamEvent(type="thinking_delta", text="x" * 1234)

    eng.client.stream_message = MagicMock(side_effect=nothing)
    out = eng.stream_response("do the thing")
    assert "empty turn" in out.lower()
    assert "1234" in out, "the diagnostic does not say how much reasoning ran"


def test_the_empty_turn_is_recorded_so_a_log_can_see_it(eng):
    def nothing(**kw):
        yield StreamEvent(type="thinking_delta", text="abc")

    eng.client.stream_message = MagicMock(side_effect=nothing)
    eng.stream_response("do the thing")
    rec = eng.last_empty_turn
    assert rec and rec["thinking_chars"] == 3
    assert rec["duration_s"] >= 0.0


def test_a_normal_turn_is_not_marked_empty(eng):
    def fine(**kw):
        yield StreamEvent(type="text_delta", text="the answer")

    eng.client.stream_message = MagicMock(side_effect=fine)
    assert eng.stream_response("q") == "the answer"
    assert eng.last_empty_turn is None
    assert [m["role"] for m in eng.messages] == ["user", "assistant"]


def test_the_diagnostic_is_never_stored_as_the_models_answer(eng):
    """It reports on the turn; feeding it back as context would teach the
    model that "[empty turn] ..." is a thing it says."""
    def nothing(**kw):
        yield StreamEvent(type="thinking_delta", text="...")

    eng.client.stream_message = MagicMock(side_effect=nothing)
    eng.stream_response("q")
    assert not any(m.get("role") == "assistant" for m in eng.messages)


def test_a_stopped_empty_turn_still_reads_as_stopped(eng):
    """The stop branch must keep its own wording — an empty turn the user
    ended is not a backend that produced nothing."""
    def stopped(**kw):
        eng.request_stop()
        yield StreamEvent(type="thinking_delta", text="...")

    eng.client.stream_message = MagicMock(side_effect=stopped)
    out = eng.stream_response("q")
    assert "empty turn" not in out.lower()
    assert eng.messages == []
    assert eng.last_empty_turn is None


# ---------------------------------------------------------------------------
# The CLI backend: no silent end of turn
# ---------------------------------------------------------------------------

def _client(tmp_path):
    c = CLIClient.__new__(CLIClient)
    c.claude_path = "cli"
    c.model = "m"
    c.permission_mode = ""
    c.cwd = str(tmp_path)
    c.mcp_config = ""
    c.allowed_tools = None
    c.extra_dirs = None
    c._proc = None
    c._session_id = ""
    return c


def _proc(lines, rc, stderr=""):
    p = MagicMock()
    p.stdout = iter(lines)
    p.poll.return_value = rc
    p.stderr.read.return_value = stderr
    return p


def test_a_killed_process_with_no_stderr_raises(tmp_path):
    """SIGKILL (OOM) is the case that produced no diagnostic at all: a
    non-zero exit and not one byte of stderr."""
    c = _client(tmp_path)
    with pytest.raises(RuntimeError) as exc:
        list(c._read_turn(_proc([], -9, stderr="")))
    assert "signal 9" in str(exc.value)


def test_a_clean_exit_without_a_result_event_raises(tmp_path):
    c = _client(tmp_path)
    with pytest.raises(RuntimeError) as exc:
        list(c._read_turn(_proc([], 0)))
    assert "without completing the turn" in str(exc.value)


def test_a_live_process_whose_output_ended_raises(tmp_path):
    c = _client(tmp_path)
    p = _proc([], None)
    with pytest.raises(RuntimeError) as exc:
        list(c._read_turn(p))
    assert "stopped streaming" in str(exc.value)
    assert c._proc is None, "the dead pipe would be handed to the next send"


def test_a_stderr_diagnostic_is_still_preferred(tmp_path):
    c = _client(tmp_path)
    with pytest.raises(RuntimeError) as exc:
        list(c._read_turn(_proc([], 1, stderr="config parse error")))
    assert "config parse error" in str(exc.value)


def test_a_complete_turn_still_ends_with_the_accounting(tmp_path):
    c = _client(tmp_path)
    lines = [
        json.dumps({"type": "result", "usage": {"input_tokens": 10,
                                                "output_tokens": 3},
                    "total_cost_usd": 0.01, "result": "hello",
                    "session_id": "s1"}),
    ]
    events = list(c._read_turn(_proc(lines, None)))
    assert [e.type for e in events] == ["text_delta", "message_delta"]
    assert events[-1].cost_usd == 0.01


def test_the_engine_turns_that_raise_into_a_visible_error(tmp_path):
    """The two halves together: the CLI raises instead of returning
    nothing, and the engine's exception path takes the question back out
    rather than leaving it to be overwritten."""
    with patch("delfin.agent.engine.create_client", return_value=MagicMock()):
        engine = E.AgentEngine(
            repo_dir=tmp_path, backend="api", provider="kit",
            model="kit.qwen3.5-397b-A17b", mode="solo")

    c = _client(tmp_path)

    def killed(**kw):
        yield from c._read_turn(_proc([], -9))

    engine.client.stream_message = MagicMock(side_effect=killed)
    with pytest.raises(RuntimeError):
        engine.stream_response("the long task")
    assert engine.messages == []
