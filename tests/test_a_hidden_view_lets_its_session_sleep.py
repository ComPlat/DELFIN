"""A view the person does not look at keeps no session awake.

The service puts a session nobody watches to sleep (no window attached,
nothing to do) and wakes it when it is used. The dashboard attaches a
session's view only while that view is the one on screen: a hidden view
is detached, so its session may sleep, and a page that comes back (a
reload) wakes only the session it shows -- not every one it remembers.

Real path: a private service with a one-second sleep, the dashboard's
session list with two remembered service sessions. Control: the one on
screen is awake; switching views wakes the other and lets the first sleep.
"""

from __future__ import annotations

import time

import pytest

pytest.importorskip("ipywidgets")

from delfin.agent.service.client import ServiceClient           # noqa: E402
from delfin.dashboard import live_sessions as LS                 # noqa: E402
from tui_service import ServiceProcess                           # noqa: E402


def _until(cond, timeout: float = 40.0):
    deadline = time.monotonic() + timeout
    while time.monotonic() < deadline:
        got = cond()
        if got:
            return got
        time.sleep(0.1)
    return cond()


def _state(c, sid):
    return next(s for s in c.request("sessions.list")["sessions"] if s["sessionId"] == sid)["state"]


class _Build:
    def __call__(self, ctx):
        raise AssertionError("no session may start in the kernel here")


@pytest.fixture
def home(tmp_path, monkeypatch):
    from delfin.agent import scheduler as S
    from delfin.dashboard import agent_sessions as AS

    monkeypatch.setattr(AS, "_OPEN_SESSIONS_PATH", tmp_path / "open.json")
    monkeypatch.setattr(S, "_GLOBAL", S.Scheduler(path=tmp_path / "cron.json"))
    monkeypatch.delenv("DELFIN_RESUME_SESSION", raising=False)
    for name in ("launch", "calc", "agent_workspace", "archive", "office"):
        (tmp_path / name).mkdir(exist_ok=True)
    monkeypatch.setenv("DELFIN_LAUNCH_CWD", str(tmp_path / "launch"))
    monkeypatch.setattr(AS, "_saved_session", lambda sid: None)
    return tmp_path


def _ctx(tmp_path):
    from delfin.dashboard.context import DashboardContext

    ctx = DashboardContext(calc_dir=tmp_path / "calc", agent_dir=tmp_path / "agent_workspace",
                           archive_dir=tmp_path / "archive", office_dir=tmp_path / "office")
    ctx.run_js = lambda script: None
    return ctx


def test_a_reload_wakes_only_the_session_on_screen_and_a_hidden_one_sleeps(home, tmp_path):
    from delfin.dashboard import agent_sessions as AS

    (tmp_path / "svc").mkdir()
    svc = ServiceProcess(tmp_path / "svc", start=False)
    svc.env["DELFIN_AGENT_SLEEP_AFTER_S"] = "1"
    svc.start()
    dash = None
    try:
        tui = ServiceClient(svc.socket, client="delfin-tui")
        svc.clients.append(tui)
        one = svc.session(tui, attach=False, cwd=str(svc.workspace("one")))
        two = svc.session(tui, attach=False, cwd=str(svc.workspace("two")))
        assert _until(lambda: _state(tui, one) == "sleeping" and _state(tui, two) == "sleeping")
        # The page remembers both; the second was on screen last.
        AS.save_open_sessions([{"session_id": one, "workspace": "", "live": True},
                               {"session_id": two, "workspace": "", "live": True}])
        dash = LS.LiveClient(svc.socket)
        _widget, refs = AS.create_tab(_ctx(home), build=_Build(), service=lambda: dash,
                                      session_defaults=lambda: {"fake": True})
        shown = refs["active"]()
        assert shown["session_id"] == two
        # Control: the one on screen is woken and stays awake.
        assert _until(lambda: _state(tui, two) not in ("sleeping",))
        time.sleep(3)
        assert _state(tui, two) != "sleeping"
        assert _state(tui, one) == "sleeping"                 # not woken by the reload
        # Switching views: the other wakes, the first may sleep again.
        first = next(r for r in refs["sessions"]() if r["session_id"] == one)
        refs["activate"](first["key"])
        assert _until(lambda: _state(tui, one) != "sleeping")
        assert _until(lambda: _state(tui, two) == "sleeping")
        # The hidden view's row still says what its session does.
        refs["refresh"]()
        hidden = next(r for r in refs["sessions"]() if r["session_id"] == two)
        assert hidden["title_btn"].description == f"two · {two[:8]}"
        refs["stop_refresh"]()
    finally:
        if dash is not None:
            dash.close()
        svc.stop()
