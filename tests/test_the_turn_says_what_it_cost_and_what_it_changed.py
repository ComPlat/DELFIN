"""A running turn is not a black box, and a finished one owes an account.

Three things a terminal has to say that the dashboard said with widgets.

WHAT IT IS COSTING, while it runs. The engine's counters are cumulative
for the SESSION, so a live line that reads them raw shows the whole
conversation's cost as this turn's — right on the first turn and wrong on
every one after it. Differencing against a baseline taken at turn start
is the whole mechanism, and it is the kind of thing that looks correct
until the second turn.

WHAT IT DISCARDED. Auto-compaction drops the early history to stay under
the window. Saying nothing makes a long session quietly lose what it was
told first.

WHAT IT LEFT OPEN. Tasks the agent created and did not finish.
"""

from __future__ import annotations

import io


from delfin.agent import repl, repl_render as rr


PLAIN = rr.Theme(enabled=False)


class _Tty(io.StringIO):
    def isatty(self):
        return True


class _Engine:
    """Counters that only ever go up, like the real ones."""

    def __init__(self):
        self.client = type("C", (), {"model": "kit.qwen3.5-397b-A17b"})()
        self.token_usage = {"input": 0, "output": 0}
        self.kit_permissions = type("P", (), {"mode": "plan"})()
        self.session_id = "sid"
        self.last_compaction_info = None
        self._status = {"input_tokens": 0, "output_tokens": 0,
                        "cost_usd": 0.0, "mode": "solo"}

    def get_status(self):
        return dict(self._status)

    def spend(self, tin, tout, cost):
        self._status["input_tokens"] += tin
        self._status["output_tokens"] += tout
        self._status["cost_usd"] += cost

    def _compaction_status_line(self):
        return "- Last compaction: 12 msg(s) compacted, ~4000 tokens saved"


def _agent(engine=None, *, tty=True):
    engine = engine or _Engine()
    out = io.StringIO()
    err = _Tty() if tty else io.StringIO()
    agent = repl.TerminalAgent(engine, out=out, err=err)
    agent.transcript.theme = PLAIN
    return agent, engine, err


# ---------------------------------------------------------------------------
# The live line
# ---------------------------------------------------------------------------

def test_the_live_line_reports_this_turn_and_not_the_session():
    """The defect that only shows up on the second turn."""
    agent, engine, _err = _agent()
    engine.spend(5000, 900, 0.40)          # an earlier turn already happened

    agent._turn_base = (5000, 900, 0.40)   # baseline as turn() takes it
    engine.spend(120, 30, 0.01)            # what THIS turn has cost so far

    line = agent._status_line()
    assert "↑120 ↓30" in line, f"reported the session, not the turn: {line}"
    assert "$0.0100" in line


def test_the_live_line_names_the_model_the_posture_and_the_way_out():
    agent, _engine, _err = _agent()
    line = agent._status_line()
    assert "kit.qwen3.5-397b-A17b" in line
    assert "plan" in line
    assert "esc to interrupt" in line


def test_the_live_line_never_grows_past_the_terminal():
    agent, _engine, _err = _agent()
    agent.transcript.width = 40
    assert len(agent._status_line()) <= 40


def test_a_backend_that_cannot_be_asked_still_paints_a_line():
    class _Mute:
        client = None
        token_usage = {}

        def get_status(self):
            raise RuntimeError("backend is gone")

    agent, _engine, _err = _agent(_Mute())
    assert agent._status_line()          # must not raise


# ---------------------------------------------------------------------------
# One row, one owner
# ---------------------------------------------------------------------------

def test_what_is_being_typed_wins_the_row():
    agent, _engine, err = _agent()
    agent._turn_active.set()
    agent._draw_input_line("half a sentence")
    assert agent._bottom == "» half a sentence"

    agent._repaint_bottom(force=True)
    assert agent._bottom == "» half a sentence", (
        "the status line must not overwrite what the user is typing")


def test_the_status_takes_the_row_back_when_the_line_is_cleared():
    agent, _engine, _err = _agent()
    agent._turn_active.set()
    agent._draw_input_line("typing")
    agent._clear_input_line()
    agent._repaint_bottom(force=True)
    assert "esc to interrupt" in agent._bottom


def test_nothing_is_painted_between_turns():
    agent, _engine, _err = _agent()
    agent._repaint_bottom(force=True)
    assert agent._bottom == ""


def test_a_redirected_stderr_gets_no_cursor_control():
    agent, _engine, err = _agent(tty=False)
    agent._turn_active.set()
    agent._repaint_bottom(force=True)
    assert err.getvalue() == ""


# ---------------------------------------------------------------------------
# The account a finished turn owes
# ---------------------------------------------------------------------------

def test_a_compacted_turn_says_so():
    agent, engine, err = _agent()
    before = engine.last_compaction_info
    engine.last_compaction_info = {"kind": "compaction", "messages_compacted": 12}
    agent._report_compaction(before)
    assert "12 msg(s) compacted" in err.getvalue()


def test_a_turn_that_compacted_nothing_says_nothing():
    agent, engine, err = _agent()
    before = engine.last_compaction_info
    agent._report_compaction(before)
    assert err.getvalue() == ""


def test_the_task_list_is_off_until_it_is_asked_for(monkeypatch):
    agent, _engine, err = _agent()
    monkeypatch.setattr("delfin.agent.task_ticker.render_text",
                        lambda ws, **kw: "[ ] finish the migration")
    agent._report_tasks()
    assert err.getvalue() == ""

    agent._show_tasks = True
    agent._report_tasks()
    assert "finish the migration" in err.getvalue()


def test_an_empty_task_store_prints_nothing(monkeypatch):
    agent, _engine, err = _agent()
    agent._show_tasks = True
    monkeypatch.setattr("delfin.agent.task_ticker.render_text",
                        lambda ws, **kw: "(no tasks)")
    agent._report_tasks()
    assert err.getvalue() == ""


def test_a_broken_task_store_does_not_end_the_session(monkeypatch):
    agent, _engine, err = _agent()
    agent._show_tasks = True

    def _boom(ws, **kw):
        raise RuntimeError("task store is corrupt")

    monkeypatch.setattr("delfin.agent.task_ticker.render_text", _boom)
    agent._report_tasks()          # must not raise


def test_the_users_own_status_line_is_reused_not_re_derived(monkeypatch):
    """It already never raises and already refuses a workspace command."""
    seen = {}

    def _render(ctx):
        seen["ctx"] = ctx
        return "42 tokens | mode=plan"

    monkeypatch.setattr("delfin.agent.status_line.render_status_line", _render)
    agent, engine, err = _agent()
    engine.spend(30, 12, 0.02)
    agent._report_status()

    assert "42 tokens" in err.getvalue()
    assert seen["ctx"].tokens == 42
    assert seen["ctx"].mode == "plan"


def test_a_status_line_that_fails_is_simply_absent(monkeypatch):
    def _boom(ctx):
        raise RuntimeError("bad template")

    monkeypatch.setattr("delfin.agent.status_line.render_status_line", _boom)
    agent, _engine, err = _agent()
    agent._report_status()
    assert err.getvalue() == ""


# ---------------------------------------------------------------------------
# The key
# ---------------------------------------------------------------------------

def test_ctrl_t_toggles_the_task_list():
    from delfin.agent import repl_keys as rk

    assert rk.KeyDecoder().feed("\x14") == [rk.KeyEvent(rk.TASKS)]

    agent, _engine, err = _agent()
    agent._on_key(rk.KeyEvent(rk.TASKS), rk.KeyDecoder())
    assert agent._show_tasks is True
    assert "task list on" in err.getvalue()
    agent._on_key(rk.KeyEvent(rk.TASKS), rk.KeyDecoder())
    assert agent._show_tasks is False


def test_the_status_line_stands_down_while_the_answer_streams():
    """Found by watching a real turn, not by reading the code.

    Both streams share one cursor. The status repaint erases the current
    line, and during streaming that line holds the answer — so the spinner
    was rubbing out the sentence as the model wrote it, four times a
    second.
    """
    agent, _engine, err = _agent()
    agent._turn_active.set()
    agent.transcript.answer("Der Test schlaegt fehl")   # no trailing newline
    err.truncate(0), err.seek(0)

    agent._repaint_bottom(force=True)
    assert agent._bottom == "", "nothing may touch the answer's own line"
    assert err.getvalue() == ""

    agent.transcript.answer("\n")                       # answer closes its line
    agent._repaint_bottom(force=True)
    assert "esc to interrupt" in agent._bottom


def test_the_row_is_given_up_when_the_answer_starts_mid_turn():
    agent, _engine, _err = _agent()
    agent._turn_active.set()
    agent._repaint_bottom(force=True)
    assert agent._bottom, "the status owns the row before any answer text"

    agent._render_around_bottom(repl.RenderItem("text", text="Antwort"))
    assert agent._bottom == "", "the answer took the row and keeps it"
