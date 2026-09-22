"""A turn cut at the token ceiling continues itself instead of stalling.

The terminal chat's ``--max-tokens`` default (4096) overrode the role
budget, so a solo session -- whose role asks for 32768 -- ran on an
eighth of it. On a model that thinks invisibly (kit.glm-5.3), a large
edit spends the ceiling on thinking and the visible answer ends at
``finish_reason == "length"`` with no tool call: the turn just stops,
half-way through the edit, and stands at the prompt (observed with
sessions 4 and 7, 2026-09-22, operator's point d).

Two rules, judged here:

  the terminal budget      the role's, unless ``--max-tokens`` was
                           typed: the flag wins, the default does not
                           halve the role's budget silently

  length continues        a turn that ends at ``length`` WITHOUT a
                           tool call continues itself, at most twice;
                           the third time the operator hears it, over
                           the channel that reaches a person who is
                           not watching six panes

What is judged here, and the instrument:

  budget            a parser-level check: the chat default passes 0
                    (the role decides), an explicit --max-tokens
                    passes through unchanged; and the engine's
                    effective budget for a solo session is the role's
                    32768 when no flag was given

  continuation      a scripted turn that ends length-without-tool:
                    run() re-sends "continue" up to twice, and the
                    third length end leaves a note in the operator
                    inbox when an operator session is open

  the counter-test  a turn that ends length WITH a tool call, and a
                    turn that ends normally, continue nothing: the
                    tool loop owns the first case, and normality
                    needs no help
"""

from __future__ import annotations

from pathlib import Path

import pytest


# -----------------------------------------------------------------------
# The budget
# -----------------------------------------------------------------------

def test_the_chat_default_leaves_the_budget_to_the_role():
    """--max-tokens was never required to type; its 4096 default must
    not silently cut a solo session's 32768."""
    from delfin.agent import cli as C
    parser = C.build_parser()
    args = parser.parse_args(["chat"])
    assert args.max_tokens == 0, (
        "the default must hand the decision to the role (0); 4096 "
        "would cut a solo session to an eighth of its budget")
    args = parser.parse_args(["chat", "--max-tokens", "6000"])
    assert args.max_tokens == 6000, "an explicit flag wins, verbatim"


def test_the_repl_passes_the_role_deciding_zero_through():
    """ReplOptions carries the parser's 0 into run_turn, and run_turn
    passes 0 to stream_response -- which resolves it to the role's
    budget (`max_tokens or max_tokens_for_role`)."""
    from delfin.agent import repl as R

    seen = {}

    class _Engine:
        session_id = "budget-check-0001"
        token_usage = {"input": 0, "output": 0}
        messages = []

        def stream_response(self, user_message="", max_tokens=0, **kw):
            seen["max_tokens"] = max_tokens
            self.messages.append({"role": "user",
                                  "content": user_message})
            self.messages.append({"role": "assistant", "content": "ok"})
            return "ok"

        def get_status(self):
            return {}

        def export_state(self):
            return {"engine_messages": list(self.messages),
                    "token_usage": dict(self.token_usage)}

    engine = _Engine()
    import io
    out, err = io.StringIO(), io.StringIO()
    out.isatty = lambda: False
    err.isatty = lambda: False
    lines = ["hello"]                    # one turn, then the prompt

    def _read_line(_prompt=""):
        if lines:
            return lines.pop(0)
        raise KeyboardInterrupt

    agent = R.TerminalAgent(
        engine, opts=R.ReplOptions(cwd=Path.cwd(), max_tokens=0),
        out=out, err=err, read_line=_read_line)
    agent._stdin = type("S", (), {"isatty": lambda self: False})()
    agent._idle_interrupts = 1           # the next interrupt leaves (130)
    agent.run()
    assert seen["max_tokens"] == 0, (
        "the zero that means 'the role decides' must reach the engine "
        "uncut; run_turn's own default of 0 does the same, but an "
        "explicit opts.max_tokens of 0 must not be replaced")


# -----------------------------------------------------------------------
# The continuation
# -----------------------------------------------------------------------

class _LengthEngine:
    """A turn that ends at length; configurable per turn."""

    session_id = "length-check-0001"
    token_usage = {"input": 0, "output": 0}
    last_turn_stop_reason = ""

    def __init__(self, script):
        self.script = list(script)        # ("length"|"normal", had_tools)
        self.messages = []
        self.turns = 0

    def stream_response(self, user_message="", max_tokens=0, **kw):
        self.turns += 1
        kind, had_tools = self.script.pop(0) if self.script \
            else ("normal", False)
        self.messages.append({"role": "user", "content": user_message})
        self.messages.append({"role": "assistant", "content": "half"})
        if had_tools:
            # The real route a tool takes: the callback run_turn passes
            # in, which is what populates result.tool_calls.
            cb = kw.get("on_tool_use")
            if cb:
                cb("bash", "{}")
        self.last_turn_stop_reason = "length" if kind == "length" \
            else "end_turn"
        return "half"

    def _notify_operator_of_stall(self, what):
        # The engine's method the real object carries; routed to the
        # same inbox the real one writes.
        from delfin.agent import session_messages as _m
        _m.send("operator", what, from_key="length-check",
                from_title="a stalled session")

    def get_status(self):
        return {}

    def export_state(self):
        return {"engine_messages": list(self.messages),
                "token_usage": dict(self.token_usage)}


def _run_and_count_turns(script, tmp_path):
    from delfin.agent import repl as R
    import io
    engine = _LengthEngine(script)
    out, err = io.StringIO(), io.StringIO()
    out.isatty = lambda: False
    err.isatty = lambda: False
    lines = ["do the big edit"]

    def _read_line(_prompt=""):
        if lines:
            return lines.pop(0)
        raise KeyboardInterrupt

    agent = R.TerminalAgent(
        engine, opts=R.ReplOptions(cwd=tmp_path, max_tokens=0),
        out=out, err=err, read_line=_read_line)
    agent._stdin = type("S", (), {"isatty": lambda self: False})()
    agent._idle_interrupts = 1           # one interrupt leaves (130)
    try:
        agent.run()
    except KeyboardInterrupt:
        pass
    return engine


def test_a_length_end_without_a_tool_call_continues_itself(tmp_path):
    engine = _run_and_count_turns(
        [("length", False), ("length", False), ("length", False)],
        tmp_path)
    assert engine.turns == 3, (
        f"a length end continued {engine.turns - 1} times; the rule is "
        "up to two self-continuations (three turns in all)")


def test_the_third_length_end_tells_the_operator(tmp_path, monkeypatch):
    from delfin.agent import session_messages as msgs
    from delfin.agent import session_presence as pres
    monkeypatch.setattr(pres, "_DIR", tmp_path / "presence")
    monkeypatch.setattr(msgs, "_DIR", tmp_path / "inbox")
    pres.announce("operator", session_id="op", title="operator",
                  workspace=str(tmp_path))
    engine = _run_and_count_turns(
        [("length", False), ("length", False), ("length", False)],
        tmp_path)
    assert engine.turns == 3
    note = msgs.take("operator")
    assert note and "token ceiling" in note[0]["text"], (
        "after the second self-continuation the standstill goes to the "
        "operator -- the channel that reaches a person not watching "
        "six panes")


def test_a_length_end_with_a_tool_call_is_not_continued(tmp_path):
    """The tool loop owns that case: the model meant to go on and the
    next round carries the work. An extra 'continue' on top would
    double-drive it."""
    engine = _run_and_count_turns([("length", True)], tmp_path)
    assert engine.turns == 1


def test_a_normal_end_is_not_continued(tmp_path):
    engine = _run_and_count_turns([("normal", False)], tmp_path)
    assert engine.turns == 1
