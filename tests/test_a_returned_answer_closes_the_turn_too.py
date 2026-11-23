"""A turn ends on its last output, however that output was delivered.

A headless run continues a turn whose final message announces an action
it did not take. The decision stands down for a turn that ended ON a
tool round -- that turn is mid-work, not stalled -- and that is read off
``engine._tools_after_last_text``, which ``_run_once`` maintains from two
callbacks: the token callback clears it, the tool-use callback sets it.

``_run_once`` accepts the answer in two shapes. One streams through
``on_token``; the other is the RETURN VALUE of ``stream_response``, which
the result falls back to (``full_text or "".join(chunks)``). Only the
first shape reached the flag. A backend that reports its tool calls and
then returns its final answer therefore left the flag standing from the
last tool call, so a turn that called tools and closed on an
announcement read as "ended on a tool round" and was never continued --
the exact failure the continuation was built for, unreachable for every
backend that does not stream its answer.

Measured 2026-11-23 with one scripted turn delivered both ways: two
``bash`` calls, then "The cgroup read is blocked. Let me get the CPU
count a different way and write the ORCA input files." Streamed, the
flag read False and the continuation fired; returned, it read True and
the continuation was blocked.

Input: a turn's events. Output: the flag the continuation reads.

The existing contract test for the continuation substitutes ``_run_once``
and sets the flag by hand, so it asserts what the loop does with the
flag and never what the flag says. This module drives the real
``_run_once``.
"""
from __future__ import annotations

import pytest

from delfin.agent import cli as agent_cli

ANNOUNCING = ("The cgroup read is blocked. Let me get the CPU count a "
              "different way and write the ORCA input files.")
FINISHED = ("The excitation energies are 3.41 eV (PBE0) and 3.88 eV "
            "(wB97X-D); the report is in report.md.")


class _Engine:
    """Replays one turn: some tool calls, then a final answer.

    ``stream_text`` picks the delivery shape. ``trailing_tool`` puts a
    tool call AFTER the answer, which is the turn that must keep its
    stand-down.
    """

    session_id = "sess-returned-1"
    token_usage = {"input": 1, "output": 1}
    run_budget_turns = 0

    def __init__(self, text=ANNOUNCING, *, stream_text=True, tools=2,
                 trailing_tool=False):
        self.text = text
        self.stream_text = stream_text
        self.tools = tools
        self.trailing_tool = trailing_tool
        self.messages = []

    def stream_response(self, **kw):
        on_tool = kw.get("on_tool_use")
        on_token = kw.get("on_token")
        for _ in range(self.tools):
            if on_tool:
                on_tool("bash", '{"command": "nproc"}')
        if self.stream_text and on_token:
            on_token(self.text)
            if self.trailing_tool and on_tool:
                on_tool("bash", '{"command": "ls"}')
            return ""
        if self.trailing_tool and on_tool:
            on_tool("bash", '{"command": "ls"}')
        return self.text

    def export_state(self):
        return {}

    def record_cycle_outcome(self, *a, **k):
        return None


def _flag(engine) -> bool:
    agent_cli._run_once(engine, "do it")
    return bool(getattr(engine, "_tools_after_last_text", False))


# --- the defect ---------------------------------------------------------

@pytest.mark.parametrize("streamed", [True, False])
def test_a_turn_that_closes_on_its_answer_is_not_mid_work(streamed):
    """Both deliveries describe the same turn, so both read the same."""
    assert _flag(_Engine(stream_text=streamed)) is False


def test_the_answer_is_the_same_whichever_way_it_arrives():
    """The flag is the only thing that differed; the answer never did."""
    streamed, returned = _Engine(), _Engine(stream_text=False)
    assert (agent_cli._run_once(streamed, "do it")["text"]
            == agent_cli._run_once(returned, "do it")["text"] == ANNOUNCING)


# --- what must not move with it -----------------------------------------

def test_a_turn_that_ends_on_a_tool_round_still_stands_down():
    """A tool call after the answer means the turn is mid-work. That is
    the case the stand-down exists for, and it keeps it.

    Streamed only, and not by omission: an engine that returns its answer
    cannot report a tool call after it, because the answer IS the return.
    Anything it reports, it reports before returning, so for that shape
    "after the answer" is not a state the caller can observe and the
    returned text is the turn's last output by construction.
    """
    assert _flag(_Engine(trailing_tool=True)) is True


def test_an_engine_that_hands_back_its_own_accumulation_decides_nothing():
    """Returning the text it already streamed is an engine repeating
    itself, not a second answer -- the callbacks keep the verdict.

    Without this the fix would read every accumulating engine as having
    closed on text, and the stand-down above would never fire for one.
    """
    class _Accumulating(_Engine):
        def stream_response(self, **kw):
            on_tool, on_token = kw.get("on_tool_use"), kw.get("on_token")
            if on_token:
                on_token(self.text)
            if on_tool:
                on_tool("bash", '{"command": "ls"}')
            return self.text            # the same text, handed back

    assert _flag(_Accumulating()) is True


def test_a_turn_with_no_answer_at_all_keeps_the_tool_verdict():
    """An empty answer is not a close: there is no final message to
    announce anything, and the turn's last output was the tool round."""
    assert _flag(_Engine(text="", stream_text=False)) is True


def test_a_finished_answer_needs_no_continuation_either_way():
    """The flag says the turn closed on text; whether that text asks for
    a continuation is the announcement detector's question, not this
    one. Pinned so the two are not confused."""
    from delfin.agent.api_client import _announces_an_action
    for streamed in (True, False):
        eng = _Engine(text=FINISHED, stream_text=streamed)
        assert _flag(eng) is False
    assert _announces_an_action(FINISHED) is False
    assert _announces_an_action(ANNOUNCING) is True
