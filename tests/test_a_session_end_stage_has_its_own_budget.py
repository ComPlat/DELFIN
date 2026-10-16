"""``repl._LeaveStages``: the session-end stages of an attended leave.

The pty tests in test_every_leave_is_bounded_and_keeps_its_save.py drive
the whole leave; these pin the runner's rules without a terminal: a
stage out of time does not take the later stages with it, a skip skips
everything after it and arms the short deadline once, the model stage's
limit follows the request timeout, and the total comes from the setting.
"""

from __future__ import annotations

import io
import threading
import time

from delfin.agent import repl


class _Agent:
    def __init__(self, budgeted: bool = True):
        self.budgeted = budgeted
        self.skips = 0

    def _stages_begin(self, runner):
        return self.budgeted

    def _stages_skipped(self):
        self.skips += 1


def _runner(monkeypatch, *, total=300.0, per_stage=None, budgeted=True):
    monkeypatch.setattr(repl, "_leave_stage_budget_s", lambda: total)
    if per_stage is not None:
        monkeypatch.setattr(repl, "_stage_budget_s", lambda name: per_stage)
    monkeypatch.setattr(repl, "_STAGE_SAY_AFTER_S", 0.2)
    agent = _Agent(budgeted)
    out = io.StringIO()
    return repl._LeaveStages(agent, out=out), agent, out


def test_a_stage_out_of_time_does_not_take_the_next_one(monkeypatch):
    run, agent, out = _runner(monkeypatch, per_stage=1.0)
    release = threading.Event()
    ran: list[str] = []
    t0 = time.monotonic()
    run("skill learning", lambda: release.wait(30))
    run("indexing", lambda: ran.append("indexing"))
    release.set()
    assert ran == ["indexing"]
    assert time.monotonic() - t0 < 5
    text = out.getvalue()
    assert "Finishing the session: skill learning" in text
    assert "Ctrl+C skips what is left" in text
    assert "skill learning did not finish within 1 s" in text
    assert agent.skips == 0


def test_a_skip_skips_the_rest_and_arms_the_short_deadline(monkeypatch):
    run, agent, out = _runner(monkeypatch)
    release = threading.Event()
    ran: list[str] = []

    def _press_ctrl_c():
        time.sleep(0.3)
        run.skipped_by = 2          # what the SIGINT handler sets

    threading.Thread(target=_press_ctrl_c, daemon=True).start()
    t0 = time.monotonic()
    run("skill learning", lambda: release.wait(30))
    run("indexing", lambda: ran.append("indexing"))
    run("eval loop", lambda: ran.append("eval loop"))
    release.set()
    assert ran == []
    assert time.monotonic() - t0 < 5
    assert agent.skips >= 1
    assert out.getvalue().count("Skipped the rest of the session end") == 1


def test_without_a_budget_the_stages_run_inline(monkeypatch):
    run, _agent, out = _runner(monkeypatch, budgeted=False)
    seen: list[str] = []
    run("report", lambda: seen.append(threading.current_thread().name))
    assert seen == [threading.main_thread().name]
    assert out.getvalue() == ""


def test_the_model_stage_follows_the_request_timeout(monkeypatch):
    monkeypatch.setenv("DELFIN_REQUEST_TIMEOUT_S", "45")
    assert repl._stage_budget_s("skill learning") == 75.0
    assert repl._stage_budget_s("indexing") == repl._STAGE_DEFAULT_S


def test_the_total_comes_from_the_setting(monkeypatch):
    import delfin.user_settings as us
    monkeypatch.setattr(us, "load_settings",
                        lambda *a, **k: {"agent": {"leave_stage_budget_s": 7}})
    assert repl._leave_stage_budget_s() == 7.0
    monkeypatch.setattr(us, "load_settings", lambda *a, **k: {})
    assert repl._leave_stage_budget_s() == repl._LEAVE_STAGE_BUDGET_S
