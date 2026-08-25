"""Keys during a turn: stop it, say the next thing, change the posture.

Two properties matter more than the key map itself.

THE TERMINAL IS ALWAYS GIVEN BACK. Raw mode is entered per turn and
restored on the way out AND from an atexit hook. A terminal left without
echo is the classic way this kind of code ruins someone's afternoon, and
the failure is invisible until the process is already gone.

CTRL+C STAYS A SIGNAL. `tty.setraw` clears ISIG, which would deliver
Ctrl+C as the byte 0x03 and silently disable the entire interrupt ladder.
`setcbreak` keeps it. The test for that reads the code, because the
alternative is a test that has to press Ctrl+C for real.
"""

from __future__ import annotations

import inspect
import io

import pytest

from delfin.agent import repl, repl_keys as rk, repl_render as rr


PLAIN = rr.Theme(enabled=False)


# ---------------------------------------------------------------------------
# The map
# ---------------------------------------------------------------------------

def test_a_lone_escape_ends_the_turn():
    assert rk.KeyDecoder().feed("\x1b") == [rk.KeyEvent(rk.INTERRUPT)]


def test_an_arrow_key_is_not_an_escape():
    """Esc and the start of a sequence are the same byte.

    Telling them apart by what arrives alongside is the usual heuristic;
    getting it wrong would make every arrow key abort the turn.
    """
    decoder = rk.KeyDecoder()
    assert decoder.feed("\x1b[A") == []
    assert decoder.buffer == "", "the sequence must not land in the line"


def test_shift_tab_asks_for_the_next_posture():
    assert rk.KeyDecoder().feed("\x1b[Z") == [rk.KeyEvent(rk.CYCLE_MODE)]


def test_typing_builds_a_line_and_enter_submits_it():
    decoder = rk.KeyDecoder()
    for ch in "stop and check":
        decoder.feed(ch)
    assert decoder.buffer == "stop and check"
    assert decoder.feed("\r") == [rk.KeyEvent(rk.SUBMIT, text="stop and check")]
    assert decoder.buffer == ""


def test_backspace_and_ctrl_u_edit_the_line():
    decoder = rk.KeyDecoder()
    decoder.feed("abc")
    decoder.feed("\x7f")
    assert decoder.buffer == "ab"
    decoder.feed("\x15")
    assert decoder.buffer == ""


def test_ctrl_o_and_ctrl_l_are_their_own_intentions():
    assert rk.KeyDecoder().feed("\x0f") == [rk.KeyEvent(rk.EXPAND)]
    assert rk.KeyDecoder().feed("\x0c") == [rk.KeyEvent(rk.REDRAW)]


def test_a_control_byte_the_map_does_not_claim_is_left_alone():
    decoder = rk.KeyDecoder()
    assert decoder.feed("\x03") == []          # Ctrl+C belongs to the signal
    assert decoder.buffer == ""


def test_a_pasted_line_arrives_as_one_message():
    decoder = rk.KeyDecoder()
    events = decoder.feed("do the thing\r")
    assert events[-1] == rk.KeyEvent(rk.SUBMIT, text="do the thing")


# ---------------------------------------------------------------------------
# The terminal is given back
# ---------------------------------------------------------------------------

def test_ctrl_c_is_left_as_a_signal():
    """setraw would clear ISIG and disable the interrupt ladder silently."""
    src = inspect.getsource(rk.RawMode.__enter__)
    assert "setcbreak" in src
    assert "setraw" not in src, (
        "setraw clears ISIG, so Ctrl+C would arrive as a byte and the "
        "interrupt ladder in repl.py would stop working")


def test_raw_mode_is_a_no_op_where_it_cannot_work():
    with rk.RawMode(io.StringIO()) as raw:
        assert raw.active is False
        assert raw.read_ready(0.0) == ""


def test_the_terminal_is_restored_even_when_the_turn_explodes():
    calls: list[str] = []

    class _FakeRaw(rk.RawMode):
        def __enter__(self):
            calls.append("enter")
            return self

        def restore(self):
            calls.append("restore")

    with pytest.raises(RuntimeError):
        with _FakeRaw(io.StringIO()):
            raise RuntimeError("boom")
    assert calls == ["enter", "restore"]


def test_restoring_twice_is_harmless():
    raw = rk.RawMode(io.StringIO())
    raw.restore()
    raw.restore()


def test_an_atexit_hook_covers_the_gap(monkeypatch):
    """__exit__ is the normal path; the hook is what covers a crash."""
    src = inspect.getsource(rk.RawMode)
    assert "atexit.register" in src
    assert "atexit.unregister" in src, (
        "a hook that outlives the object would restore a terminal that has "
        "since been handed to something else")


# ---------------------------------------------------------------------------
# What the loop does with the events
# ---------------------------------------------------------------------------

class _Perms:
    def __init__(self, mode="plan"):
        self.mode = mode


class _Engine:
    def __init__(self, mode="plan"):
        self.kit_permissions = _Perms(mode)
        self.token_usage = {"input": 0, "output": 0}
        self.client = None

    def set_kit_permission_mode(self, mode):
        self.kit_permissions.mode = mode

    def request_stop(self):
        self.stopped = True


class _Tty(io.StringIO):
    """A stream that claims to be a terminal, so cursor control is allowed."""

    def isatty(self):
        return True


def _agent(engine, *, tty_err=True):
    out = io.StringIO()
    err = _Tty() if tty_err else io.StringIO()
    agent = repl.TerminalAgent(engine, out=out, err=err)
    agent.transcript.theme = PLAIN
    return agent, out, err


def test_a_line_typed_during_a_turn_is_queued_not_injected():
    """Queued, deliberately.

    engine.steer exists and works on this backend, but a queued message
    cannot be lost and cannot land in a context the user could not see.
    """
    agent, _out, err = _agent(_Engine())
    agent._on_key(rk.KeyEvent(rk.SUBMIT, text="also check the docs"),
                  rk.KeyDecoder())
    assert agent.queued == ["also check the docs"]
    assert "queued" in err.getvalue()


def test_shift_tab_walks_the_ladder_and_never_reaches_bypass():
    engine = _Engine("plan")
    agent, _out, _err = _agent(engine)
    seen = []
    for _ in range(6):
        agent._cycle_mode()
        seen.append(engine.kit_permissions.mode)
    assert seen == ["default", "acceptEdits", "plan",
                    "default", "acceptEdits", "plan"]
    assert "bypassPermissions" not in seen, (
        "unattended execution must stay something typed on purpose")


def test_a_backend_without_a_gate_says_so_instead_of_crashing():
    class _NoPerms:
        token_usage = {"input": 0, "output": 0}

        @property
        def kit_permissions(self):
            raise RuntimeError("no permissions on this provider")

    agent, _out, err = _agent(_NoPerms())
    agent._cycle_mode()
    assert "no permission gate" in err.getvalue()


def test_the_typed_line_survives_a_tool_line_landing_on_top_of_it():
    """The bottom-line problem, and the reason typing is worth having."""
    agent, _out, err = _agent(_Engine())
    agent._draw_input_line("half a sen")
    agent._render_around_input(
        repl.RenderItem("tool_use", name="bash", text='{"command": "ls"}'))

    text = err.getvalue()
    assert "ls" in text
    assert text.rstrip().endswith("half a sen"), (
        "the line being typed has to be put back after the transcript line")


def test_ctrl_o_expands_the_last_result_and_says_when_there_is_none():
    agent, _out, err = _agent(_Engine())
    agent._expand_last_result()
    assert "nothing to expand" in err.getvalue()

    agent._render_around_input(
        repl.RenderItem("tool_result", name="read_file", text="alpha\nbeta"))
    agent._expand_last_result()
    assert "alpha" in err.getvalue() and "beta" in err.getvalue()


def test_an_expanded_result_cannot_repaint_the_terminal():
    agent, _out, err = _agent(_Engine())
    agent._render_around_input(
        repl.RenderItem("tool_result", name="read_file",
                        text="ok\x1b[2J\x1b]0;pwned\x07"))
    agent._expand_last_result()
    assert "\x1b" not in err.getvalue()


def test_a_queued_line_is_sent_before_the_next_prompt_is_read():
    engine = _Engine()
    asked: list[str] = []

    def _read(_prompt):
        asked.append("read")
        return "/exit"

    agent, _out, _err = _agent(engine)
    agent._read_line = _read
    agent.queued.append("the queued one")
    turns: list[str] = []
    agent.turn = lambda prompt: turns.append(prompt) or repl.TurnResult()

    agent.run()
    assert turns == ["the queued one"], (
        "a message typed during a turn must go out before the prompt asks "
        "for another")


def test_cursor_control_never_reaches_a_redirected_stderr():
    """The redraw needs \\r and an erase sequence.

    Those are precisely the characters this codebase strips out of tool
    output before printing it, so writing them into a log file for a line
    nobody can see would be the same defect from the other side.
    """
    agent, _out, err = _agent(_Engine(), tty_err=False)
    agent._draw_input_line("typing away")
    agent._clear_input_line()
    assert err.getvalue() == ""
    assert agent._input_line == "", "the buffer state still has to be tracked"
