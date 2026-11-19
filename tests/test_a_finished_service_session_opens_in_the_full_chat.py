"""A finished session of the agent service, clicked in the running list, continues in the full chat.

Dashboard sessions run in the kernel by default (#636). A session the service
holds that is done or asleep is taken over: the service saves and ends it and
the kernel resumes the same session id in the agent tab's full chat. A session
that is working or waits for an answer stays in the service and opens in its
view. Driven against a real private service with the scripted engine; the
kernel tab is a stand-in that records what it was asked to open.
"""
from __future__ import annotations

import time

import pytest

pytest.importorskip("ipywidgets")

import ipywidgets as widgets                                     # noqa: E402

from delfin.agent.service.client import ServiceClient           # noqa: E402
from delfin.dashboard import live_sessions as LS                 # noqa: E402
from tui_service import ServiceProcess, turn                     # noqa: E402


def _until(cond, timeout: float = 40.0):
    deadline = time.monotonic() + timeout
    while time.monotonic() < deadline:
        got = cond()
        if got:
            return got
        time.sleep(0.1)
    return cond()


@pytest.fixture
def setup(tmp_path, monkeypatch):
    from delfin.agent import scheduler as S
    from delfin.dashboard import agent_sessions as AS
    from delfin.dashboard import tab_agent
    from delfin.dashboard.context import DashboardContext

    monkeypatch.setattr(AS, "_OPEN_SESSIONS_PATH", tmp_path / "open.json")
    monkeypatch.setattr(S, "_GLOBAL", S.Scheduler(path=tmp_path / "cron.json"))
    monkeypatch.delenv("DELFIN_RESUME_SESSION", raising=False)
    monkeypatch.delenv(AS.SESSIONS_ENV, raising=False)
    for name in ("launch", "calc", "agent_workspace", "archive", "office"):
        (tmp_path / name).mkdir(exist_ok=True)
    monkeypatch.setenv("DELFIN_LAUNCH_CWD", str(tmp_path / "launch"))
    (tmp_path / "svc").mkdir()
    svc = ServiceProcess(tmp_path / "svc")
    opened = []

    def build(ctx):
        opened.append(getattr(ctx, "resume_session_id", None) or getattr(ctx, "session_id", None))
        state = {"active_session_id": "", "chat_messages": [], "streaming": False, "engine": None}
        return widgets.VBox(), {"state": state, "shutdown": lambda: None}

    monkeypatch.setattr(tab_agent, "create_tab", build)
    dash = LS.LiveClient(svc.socket)
    monkeypatch.setattr(LS, "connect", lambda **kw: dash)
    ctx = DashboardContext(calc_dir=tmp_path / "calc", agent_dir=tmp_path / "agent_workspace",
                           archive_dir=tmp_path / "archive", office_dir=tmp_path / "office")
    ctx.run_js = lambda script: None
    try:
        yield svc, dash, ctx, opened, AS, monkeypatch
    finally:
        dash.close()
        svc.stop()


def test_a_done_session_is_moved_into_the_kernel(setup):
    svc, dash, ctx, opened, AS, monkeypatch = setup
    tui = ServiceClient(svc.socket, client="delfin-tui"); svc.clients.append(tui)
    sid = svc.session(tui, attach=True)
    turn(tui, sid, "hello")
    # the scripted engine does not write a session file: the take-over reads one
    monkeypatch.setattr(AS, "_saved_session", lambda s: {"session_id": s} if s == sid else None)
    _w, refs = AS.create_tab(ctx)
    assert _until(lambda: refs["running_panel"].running.summary(sid))
    before = len(opened)
    rec = refs["open_from_list"](sid)
    assert rec is not None and not rec.get("live")             # the kernel's full chat
    assert len(opened) == before + 1
    states = {s["sessionId"]: s["state"] for s in tui.request("sessions.list")["sessions"]}
    assert sid not in states                                    # the service let it go


def test_a_working_session_stays_in_the_service(setup):
    svc, dash, ctx, opened, AS, monkeypatch = setup
    tui = ServiceClient(svc.socket, client="delfin-tui"); svc.clients.append(tui)
    sid = svc.session(tui, attach=True, fakeDelay=0.5)
    tui.request("session.prompt", {"sessionId": sid, "text": "a long one " * 40})
    _w, refs = AS.create_tab(ctx)
    assert _until(lambda: LS.state_word(refs["running_panel"].running.summary(sid) or {}) == "working")
    before = len(opened)
    rec = refs["open_from_list"](sid)
    assert rec is not None and rec.get("live")                  # the service's view
    assert len(opened) == before
    assert sid in {s["sessionId"] for s in tui.request("sessions.list")["sessions"]}
