"""Ctrl+C has to end the turn without destroying what was asked.

THE HAZARD, which predates the terminal. `stream_response` catches
`Exception`, not `BaseException`. A KeyboardInterrupt raised inside it
therefore skips the cleanup that takes the unanswered user message back
out of the history — and the alternation sanitiser resolves two
consecutive user messages by keeping the NEWEST. So an interrupt at the
wrong moment silently replaces the question with whatever is typed next,
and the model answers something nobody asked. `cli run` has this today.

THE FIX IS STRUCTURAL, not a try/except. Python delivers signals to the
main thread; the turn runs on a worker. So the interrupt cannot reach the
turn body at all — it arrives as `request_stop()`, and the turn unwinds
through the stop path it already has.

The rest of this file pins the ordering rules that make the loop safe to
return from: join before clearing the stop, and never return to the
prompt while a worker is still unwinding.
"""

from __future__ import annotations

import io
import signal
import threading
from pathlib import Path

import pytest

from delfin.agent import repl, repl_render as rr


PLAIN = rr.Theme(enabled=False)


class FakeEngine:
    """An engine whose turn blocks until stopped, like a long tool call."""

    def __init__(self):
        self.client = _FakeClient()
        self.token_usage = {"input": 0, "output": 0}
        self.stop_requested = threading.Event()
        self.entered = threading.Event()
        self.cleared = 0
        self.busy = False
        self.concurrent = False
        self.finished_naturally = False

    def stream_response(self, **kwargs):
        if self.busy:
            self.concurrent = True
            return "A turn is already running on this session."
        self.busy = True
        try:
            self.entered.set()
            if self.stop_requested.wait(timeout=5):
                return "[stopped] Stopped before any answer text."
            self.finished_naturally = True
            return "done"
        finally:
            self.busy = False

    def request_stop(self):
        self.stop_requested.set()

    def clear_stop(self):
        self.cleared += 1


class _FakeClient:
    def __init__(self):
        self.signalled = 0

    def signal_stop(self):
        self.signalled += 1


def _agent(engine, **kw):
    out, err = io.StringIO(), io.StringIO()
    agent = repl.TerminalAgent(engine, out=out, err=err, **kw)
    agent.transcript.theme = PLAIN
    return agent, out, err


# ---------------------------------------------------------------------------
# The interrupt
# ---------------------------------------------------------------------------

def test_the_turn_never_sees_the_interrupt():
    """The structural guarantee.

    If the turn ran on the main thread, this signal would land inside
    stream_response and skip the cleanup that protects the question.
    """
    engine = FakeEngine()
    agent, out, err = _agent(engine)

    def _interrupt_once_running():
        engine.entered.wait(timeout=5)
        agent._on_sigint(signal.SIGINT, None)

    threading.Thread(target=_interrupt_once_running, daemon=True).start()
    result = agent.turn("do the long thing")

    assert engine.stop_requested.is_set()
    assert engine.client.signalled == 1
    assert engine.finished_naturally is False
    assert "stopped" in result.text.lower()


def test_the_first_interrupt_does_not_raise():
    engine = FakeEngine()
    agent, _out, _err = _agent(engine)
    agent._turn_active.set()
    agent._on_sigint(signal.SIGINT, None)      # must not raise
    agent._on_sigint(signal.SIGINT, None)      # nor the second
    with pytest.raises(KeyboardInterrupt):
        agent._on_sigint(signal.SIGINT, None)  # the third abandons


def test_an_interrupt_at_the_prompt_raises_straight_away():
    engine = FakeEngine()
    agent, _out, _err = _agent(engine)
    assert not agent._turn_active.is_set()
    with pytest.raises(KeyboardInterrupt):
        agent._on_sigint(signal.SIGINT, None)


def test_the_handler_never_writes_to_the_terminal():
    """It runs between two bytecodes, possibly inside a write.

    Printing from there interleaves with the pump's output once every few
    interrupts — reproducible only by luck.
    """
    engine = FakeEngine()
    agent, out, err = _agent(engine)
    agent._turn_active.set()
    agent._on_sigint(signal.SIGINT, None)
    assert out.getvalue() == ""
    assert err.getvalue() == ""
    assert not agent._q.empty(), "the notice has to be queued, not printed"


# ---------------------------------------------------------------------------
# Ordering
# ---------------------------------------------------------------------------

def test_a_second_turn_never_starts_while_the_first_unwinds():
    """The turn gate refuses by RETURNING a sentence, not by raising.

    Return to the prompt too early and the next turn renders machinery
    speech as the model's answer.
    """
    engine = FakeEngine()
    agent, out, err = _agent(engine)

    threading.Thread(
        target=lambda: (engine.entered.wait(timeout=5),
                        agent._on_sigint(signal.SIGINT, None)),
        daemon=True).start()
    agent.turn("first")

    engine.stop_requested.clear()
    engine.entered.clear()
    threading.Thread(
        target=lambda: (engine.entered.wait(timeout=5), engine.request_stop()),
        daemon=True).start()
    second = agent.turn("second")

    assert engine.concurrent is False
    assert "already running" not in second.text
    assert "already running" not in out.getvalue()


def test_the_stop_is_cleared_only_after_the_worker_is_gone():
    engine = FakeEngine()
    agent, _out, _err = _agent(engine)
    order: list[str] = []

    real_clear = engine.clear_stop

    def _clear():
        order.append("cleared while alive" if engine.busy else "cleared after")
        real_clear()

    engine.clear_stop = _clear
    threading.Thread(
        target=lambda: (engine.entered.wait(timeout=5), engine.request_stop()),
        daemon=True).start()
    agent.turn("go")

    assert order == ["cleared after"]
    assert engine.cleared == 1


def test_a_worker_that_dies_is_reported_not_swallowed():
    class Exploding(FakeEngine):
        def stream_response(self, **kwargs):
            raise SystemExit("backend vanished")

    agent, out, err = _agent(Exploding())
    result = agent.turn("go")
    assert "backend vanished" in err.getvalue()
    assert result.error


# ---------------------------------------------------------------------------
# Input
# ---------------------------------------------------------------------------

def test_a_pasted_block_is_one_message():
    lines = iter(['"""', "line a", "line b", '"""'])
    assert repl.read_block(lambda _p: next(lines)) == "line a\nline b"


def test_a_trailing_backslash_continues_the_line():
    lines = iter(["do this \\", "and that"])
    assert repl.read_block(lambda _p: next(lines)) == "do this and that"


def test_an_ordinary_line_is_itself():
    assert repl.read_block(lambda _p: "just this") == "just this"


def test_the_loop_leaves_on_ctrl_d():
    engine = FakeEngine()

    def _eof(_prompt):
        raise EOFError

    agent, _out, _err = _agent(engine, read_line=_eof)
    assert agent.run() == 0


def test_two_interrupts_at_the_prompt_leave():
    engine = FakeEngine()

    def _interrupt(_prompt):
        raise KeyboardInterrupt

    agent, _out, err = _agent(engine, read_line=_interrupt)
    assert agent.run() == 130
    assert "once more to leave" in err.getvalue()


def test_the_exit_command_leaves_without_a_turn():
    engine = FakeEngine()
    replies = iter(["/exit"])
    agent, _out, _err = _agent(engine, read_line=lambda _p: next(replies))
    assert agent.run() == 0
    assert engine.entered.is_set() is False


# ---------------------------------------------------------------------------
# History
# ---------------------------------------------------------------------------

def test_the_history_file_is_the_users_alone(monkeypatch, tmp_path):
    pytest.importorskip("readline")
    monkeypatch.setattr(Path, "home", classmethod(lambda cls: tmp_path))
    engine = FakeEngine()
    replies = iter(["/exit"])
    agent, _out, _err = _agent(engine, read_line=lambda _p: next(replies))
    agent.run()

    path = tmp_path / ".delfin" / repl.HISTORY_NAME
    assert path.exists()
    assert oct(path.stat().st_mode & 0o777) == "0o600", (
        "prompts describe the user's project")


# ---------------------------------------------------------------------------
# The subprocesses a terminal Ctrl+C would otherwise take with it
# ---------------------------------------------------------------------------

def test_the_long_lived_subprocesses_are_out_of_the_terminals_process_group():
    """A stdio MCP server and the CLI backend both outlive a turn.

    Under a terminal front-end they would sit in the foreground process
    group and take every Ctrl+C aimed at the agent, so interrupting one
    turn would tear down every configured server for the rest of the
    session. The dashboard never had a controlling terminal, which is why
    nothing surfaced this before.
    """
    import inspect
    from delfin.agent import mcp_client
    from delfin.agent import api_client

    src = inspect.getsource(mcp_client.MCPServer.start)
    assert "subprocess.Popen(" in src, "the spawn moved; re-point this test"
    assert "start_new_session=True" in src

    src = inspect.getsource(api_client.CLIClient._ensure_proc)
    assert "self._proc = subprocess.Popen(" in src, (
        "the spawn moved; re-point this test")
    assert "start_new_session=True" in src
