"""Auto mode's state is changed by one thread at a time, and a slow check
does not hold the others.

Input: the dashboard's three threads on one AutoController -- the tick
(next_prompt), the turn's worker (finish_iteration, which runs the loop's
check) and the UI (stop, the toggle, a typed message) -- with a check that
takes a second. Output: the UI's calls return while the check runs; a stop
that lands during the check ends the loop cleanly (the iteration's commit
stays, its log row says why it was not judged); the tick starts nothing
while an iteration is being judged; the dashboard makes one controller
however many threads ask for it first.
"""

from __future__ import annotations

import json
import subprocess
import threading
import time
from pathlib import Path
from types import SimpleNamespace

import pytest

from delfin.agent import auto_mode as am
from delfin.agent import loops
from delfin.agent.harness_rows import HarnessText


def _git(ws, *a):
    return subprocess.run(["git", "-C", str(ws), *a], capture_output=True,
                          text=True, check=True).stdout.strip()


@pytest.fixture()
def repo(tmp_path, monkeypatch):
    home = tmp_path / "home"
    home.mkdir()
    monkeypatch.setattr(Path, "home", classmethod(lambda cls: home))
    ws = tmp_path / "ws"
    ws.mkdir()
    _git(ws, "init", "-q", "-b", "main")
    _git(ws, "config", "user.name", "t")
    _git(ws, "config", "user.email", "t@example.org")
    (ws / "r.txt").write_text("0\n")
    _git(ws, "add", "-A")
    _git(ws, "commit", "-q", "-m", "i")
    d = ws / ".delfin" / "loops"
    d.mkdir(parents=True)
    (d / "slow.md").write_text(
        "---\ngoal: g\ncheck: python check.py\nmetric: 'n=(\\d+)'\n---\n1. a\n")
    loops.confirm(loops.load("slow", ws))
    return ws


class Engine:
    session_id = "s"

    def __init__(self):
        self.kit_permissions = SimpleNamespace(mode="acceptEdits", auto_veto=None)
        self.token_usage = {}
        self.messages = []
        self.pending_turn_continuation = ""

    def request_stop(self):
        pass


def _slow_checks():
    calls = {"n": 0}

    def run(engine, command):
        calls["n"] += 1
        if calls["n"] > 1:            # the baseline is quick
            time.sleep(1.0)
        return 1, f"n={10 - calls['n']}"
    return run


def test_a_stop_during_the_check_returns_at_once_and_ends_cleanly(repo):
    ctl = am.AutoController(repo, run_check=_slow_checks(),
                            notify=lambda *a: None,
                            settings={"loop_gap_seconds": 0})
    eng = Engine()
    ctl.start_loop(loops.load("slow", repo), eng)
    assert isinstance(ctl.next_prompt(eng), HarnessText)
    (repo / "r.txt").write_text("1\n")
    done = {}
    worker = threading.Thread(
        target=lambda: done.setdefault("digest",
                                       ctl.finish_iteration(eng, "LOOP: progress")))
    worker.start()
    time.sleep(0.3)                    # the check is running now
    assert ctl.next_prompt(eng) == ""  # the tick starts nothing meanwhile
    t0 = time.monotonic()
    ctl.note_user_input(eng)
    msg = ctl.disable("turned off in the dashboard")
    assert time.monotonic() - t0 < 0.5, "the UI waited for the check"
    assert "stopped" in msg
    worker.join(5)
    assert not worker.is_alive()
    assert ctl.loop is None and not ctl.enabled
    rows = [json.loads(l) for l in ctl.log_path.read_text().splitlines()]
    assert "stopped while" in rows[-1]["why"]
    assert _git(repo, "log", "--format=%s", "-1").startswith("slow iteration 1")


def test_threads_hammering_the_controller_leave_it_consistent(repo):
    ctl = am.AutoController(repo, run_check=lambda e, c: (1, "n=5"),
                            notify=lambda *a: None,
                            settings={"loop_gap_seconds": 0, "idle_seconds": 0})
    eng = Engine()
    ctl.start_loop(loops.load("slow", repo), eng)
    errors = []

    def tick():
        for _ in range(40):
            try:
                if ctl.next_prompt(eng):
                    ctl.finish_iteration(eng, "LOOP: progress")
            except Exception as exc:          # pragma: no cover - the failure
                errors.append(exc)

    def ui():
        for _ in range(40):
            try:
                ctl.status_line()
                ctl.state()
                ctl.pop_notices()
                ctl._queue("outward", f"x {time.monotonic()}")
            except Exception as exc:          # pragma: no cover
                errors.append(exc)

    threads = [threading.Thread(target=tick) for _ in range(3)] + \
              [threading.Thread(target=ui) for _ in range(2)]
    for t in threads:
        t.start()
    for t in threads:
        t.join(30)
    assert errors == []
    assert ctl.iteration is None
    ids = [n["id"] for n in ctl.needs_you]
    assert len(ids) == len(set(ids))


def test_the_dashboard_makes_one_controller(tmp_path, monkeypatch):
    pytest.importorskip("ipywidgets")
    from delfin.agent import scheduler as S
    from delfin.dashboard import tab_agent
    from delfin.dashboard.context import DashboardContext
    monkeypatch.setattr(S, "_GLOBAL", S.Scheduler(path=tmp_path / "cron.json"))
    for name in ("calc", "archive", "office", "ws"):
        (tmp_path / name).mkdir()
    ctx = DashboardContext(calc_dir=tmp_path / "calc",
                           archive_dir=tmp_path / "archive",
                           office_dir=tmp_path / "office")
    ctx.run_js = lambda script, **kw: None
    ctx.agent_workspace = str(tmp_path / "ws")
    _tab, refs = tab_agent.create_tab(ctx)
    try:
        made = []
        barrier = threading.Barrier(8)

        def ask():
            barrier.wait()
            made.append(refs["auto_ctl"]())

        threads = [threading.Thread(target=ask) for _ in range(8)]
        for t in threads:
            t.start()
        for t in threads:
            t.join(10)
        assert len({id(c) for c in made}) == 1
    finally:
        refs["shutdown"]()
