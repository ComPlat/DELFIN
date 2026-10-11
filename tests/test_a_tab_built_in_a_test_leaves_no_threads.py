"""A tab built in a test leaves no thread behind.

Every Agent tab starts a watcher thread that refreshes the subagent panel.
In a kernel it ends at exit; in the suite it never ended: 8 were alive
after 37 tests of four files, and one of them was counted by a neighbouring
test that patches time.sleep process-wide (CI on PR #138). The watcher now
waits on an event the tab's shutdown sets, and the suite closes every tab
it built when the test ends (conftest).
"""

from __future__ import annotations

import threading
import time


def _tab(tmp_path):
    from delfin.agent import scheduler as S
    from delfin.dashboard import tab_agent
    from delfin.dashboard.context import DashboardContext

    S._GLOBAL = S.Scheduler(path=tmp_path / "cron.json")
    for name in ("calc", "archive", "office"):
        (tmp_path / name).mkdir(exist_ok=True)
    ctx = DashboardContext(calc_dir=tmp_path / "calc",
                           archive_dir=tmp_path / "archive",
                           office_dir=tmp_path / "office")
    ctx.run_js = lambda script, **kw: None
    return tab_agent.create_tab(ctx)


def _watchers():
    return [t for t in threading.enumerate()
            if t.name == "delfin-subagent-live" and t.is_alive()]


def test_closing_the_tab_ends_its_watcher_at_once(tmp_path):
    from delfin.dashboard import tab_agent

    tab_agent.close_open_tabs()
    before = len(_watchers())
    _tab(tmp_path)
    assert len(_watchers()) == before + 1
    began = time.monotonic()
    assert tab_agent.close_open_tabs() == 1
    while _watchers() and time.monotonic() - began < 1.0:
        time.sleep(0.01)
    # Woken, not left to finish its 1.5 s interval.
    assert len(_watchers()) == before
    assert time.monotonic() - began < 1.0


def test_the_watcher_never_calls_time_sleep(tmp_path, monkeypatch):
    from delfin.dashboard import tab_agent

    tab_agent.close_open_tabs()
    seen = []
    real = time.sleep

    def _sleep(s):
        if threading.current_thread().name == "delfin-subagent-live":
            seen.append(s)
        return real(s)

    monkeypatch.setattr(time, "sleep", _sleep)
    _tab(tmp_path)
    real(0.3)
    assert seen == []
