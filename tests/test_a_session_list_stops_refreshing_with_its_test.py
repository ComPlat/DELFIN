"""A session list's refresh loop ends with the test that built it.

``agent_sessions.create_tab`` starts a refresh loop: a ``threading.Timer``
that runs ``refresh`` every ``_REFRESH_S`` and reschedules itself. Each
tick announces every open session through ``session_presence.announce``,
which writes to ``session_presence._DIR`` -- a module attribute the suite
points into each test's own home (``_isolate_user_state``). The loop had
no stop, so a list built by one test kept ticking through every later test
and wrote its presence record into THAT test's home (observed in CI:
``tests/test_v2_make_plot.py`` found ``user_home/.delfin/session_presence``,
a directory it never created).

Input: two tests in one file, run in order. The first builds a session
list. The second waits longer than one refresh period with the presence
heartbeat set to zero (so an unchanged record is rewritten on every tick)
and looks into its own presence directory. Output: that directory stays
absent. A third test drives the stop the panel exposes directly.

Rejected: resolving the presence directory once when the list is built.
It keeps the write inside the first test's directory but leaves the loop
running: it still calls ``_persist`` (``_OPEN_SESSIONS_PATH``, redirected
per test the same way), delivers messages and reads the roster of every
later test. A loop that outlives its owner is the defect; the write was
one symptom.
"""

from __future__ import annotations

import threading
import time

import pytest

pytest.importorskip("ipywidgets")

import ipywidgets as widgets                                # noqa: E402

from delfin.agent import scheduler as S                     # noqa: E402
from delfin.agent import session_presence as P              # noqa: E402
from delfin.dashboard import agent_sessions as AS           # noqa: E402
from delfin.dashboard.context import DashboardContext       # noqa: E402


class _Build:
    def __call__(self, ctx):
        state = {"active_session_id": ctx.initial_session_id or "",
                 "chat_messages": [], "streaming": False, "engine": None}
        ctx.agent_state = state
        return widgets.VBox(), {"state": state, "shutdown": lambda: None,
                                "workspace": lambda: ctx.agent_workspace}


@pytest.fixture
def home(tmp_path, monkeypatch):
    monkeypatch.setattr(AS, "_OPEN_SESSIONS_PATH", tmp_path / "open.json")
    monkeypatch.setattr(S, "_GLOBAL", S.Scheduler(path=tmp_path / "cron.json"))
    monkeypatch.delenv("DELFIN_RESUME_SESSION", raising=False)
    for name in ("launch", "calc", "agent_workspace", "archive", "office"):
        (tmp_path / name).mkdir()
    monkeypatch.setenv("DELFIN_LAUNCH_CWD", str(tmp_path / "launch"))
    monkeypatch.setattr(AS, "_saved_session", lambda sid: None)
    ctx = DashboardContext(calc_dir=tmp_path / "calc",
                           agent_dir=tmp_path / "agent_workspace",
                           archive_dir=tmp_path / "archive",
                           office_dir=tmp_path / "office")
    ctx.run_js = lambda script: None
    return ctx


def _refresh_timers() -> list[threading.Thread]:
    return [t for t in threading.enumerate()
            if isinstance(t, threading.Timer) and t.is_alive()
            and getattr(t, "function", None) is not None
            and getattr(t.function, "__name__", "") == "_tick"]


def test_1_a_session_list_is_built_and_announces_its_session(home):
    _widget, refs = AS.create_tab(home, build=_Build())
    assert len(refs["sessions"]()) == 1
    assert list(P._DIR.glob("*.json")), "the list announced nothing"
    assert _refresh_timers(), "no refresh loop was started"


def test_2_a_later_test_finds_no_presence_it_did_not_write(monkeypatch):
    # Every unchanged record is rewritten on each tick, as it is once a
    # minute in a real run (the heartbeat); one period is then enough.
    monkeypatch.setattr(P, "_HEARTBEAT_S", 0.0)
    mine = P._DIR
    deadline = time.monotonic() + AS._REFRESH_S + 1.5
    while time.monotonic() < deadline and not mine.exists():
        time.sleep(0.1)
    assert not mine.exists(), (
        f"a refresh loop of an earlier test wrote {sorted(mine.glob('*'))}")
    assert _refresh_timers() == [], "the earlier test's refresh loop still runs"


def test_3_the_stop_ends_the_loop_and_waits_for_it(home, monkeypatch):
    monkeypatch.setattr(AS, "_REFRESH_S", 0.05)
    _widget, refs = AS.create_tab(home, build=_Build())
    time.sleep(0.2)
    assert _refresh_timers(), "no refresh loop was started"
    refs["stop_refresh"]()
    assert _refresh_timers() == []
    monkeypatch.setattr(P, "_HEARTBEAT_S", 0.0)
    for f in P._DIR.glob("*.json"):
        f.unlink()
    time.sleep(0.3)
    assert list(P._DIR.glob("*.json")) == [], "a stopped loop announced again"
