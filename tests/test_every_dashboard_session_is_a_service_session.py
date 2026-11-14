"""Every session the agent tab starts is a session of the agent service.

It runs where it always ran -- in the dashboard's kernel, on the agent
tab's own turn code -- and is hosted in the person's agent service
(delfin/dashboard/hosted_sessions.py) as soon as it has an id: listed
once, under one name, by every dashboard and the TUI, joinable from each.
Another dashboard opens it as a view of the same chat (never taken over
while its kernel runs it), and what that view sends runs in the hosting
tab. Keep session is the service's switch: off, the session ends when its
last window lets go; on, it runs on.

Driven against a real private service (tui_service.ServiceProcess); the
agent tab is a stand-in with the real tab's host adapter shape, except in
the tests of the adapter itself, which build the real tab.
"""
from __future__ import annotations

import time

import pytest

pytest.importorskip("ipywidgets")

import ipywidgets as widgets                                     # noqa: E402

from delfin.agent.service.client import ServiceClient           # noqa: E402
from delfin.dashboard import agent_sessions as AS                # noqa: E402
from delfin.dashboard import kept_agent_sessions as K            # noqa: E402
from delfin.dashboard import live_sessions as LS                 # noqa: E402
from tui_service import ServiceProcess                           # noqa: E402


def _until(cond, timeout: float = 30.0):
    deadline = time.monotonic() + timeout
    while time.monotonic() < deadline:
        got = cond()
        if got:
            return got
        time.sleep(0.05)
    return cond()


class _Chat:
    def __init__(self):
        self.value = '<template class="delfin-chat-src"><div class="delfin-agent-chat">' \
                     '</div></template>'


class _Build:
    """A stand-in agent tab: its state, a control row with a Stop, and the
    host adapter tab_agent.create_tab gives (the tab's own callbacks)."""

    def __init__(self):
        self.contexts = []
        self.tabs = []

    def __call__(self, ctx):
        state = {"active_session_id": ctx.initial_session_id or "",
                 "chat_messages": [], "streaming": False, "engine": None}
        ctx.agent_state = state
        self.contexts.append(ctx)
        stop = widgets.Button(description="Stop")
        row = widgets.HBox([widgets.Dropdown(options=["solo"], description="Mode:"), stop])
        tab = {"state": state, "sent": [], "stops": 0, "chat": _Chat(), "closed": False}
        self.tabs.append(tab)

        def send(text, _tab=tab):
            _tab["sent"].append(text)
            _tab["state"]["chat_messages"].append({"role": "user", "content": text})
            return True

        def stop_turn(_tab=tab):
            _tab["stops"] += 1

        def adapter(_tab=tab):
            return {"state": _tab["state"], "send": send, "stop": stop_turn,
                    "answer": lambda text: False, "plan": lambda yes: None,
                    "approval": lambda yes: None, "question": lambda text: False,
                    "set_permission": lambda profile, confirmed: profile,
                    "choices": lambda what: {"current": "", "choices": []},
                    "info": lambda: {"provider": "kit", "model": "m1", "effort": "medium",
                                     "mode": "solo", "approval": "default"},
                    "chat_html": _tab["chat"], "status": lambda: {}}

        def shutdown(_tab=tab):
            _tab["closed"] = True

        return widgets.VBox([row]), {"state": state, "shutdown": shutdown,
                                     "control_row": row, "stop_button": stop,
                                     "host_adapter": adapter}


@pytest.fixture
def home(tmp_path, monkeypatch):
    from delfin.agent import scheduler as S

    monkeypatch.setattr(AS, "_OPEN_SESSIONS_PATH", tmp_path / "open.json")
    monkeypatch.setattr(S, "_GLOBAL", S.Scheduler(path=tmp_path / "cron.json"))
    monkeypatch.delenv("DELFIN_RESUME_SESSION", raising=False)
    monkeypatch.setenv(K.DIR_ENV, str(tmp_path / "kept"))
    for name in ("launch", "calc", "agent_workspace", "archive", "office"):
        (tmp_path / name).mkdir(exist_ok=True)
    monkeypatch.setenv("DELFIN_LAUNCH_CWD", str(tmp_path / "launch"))
    monkeypatch.setattr(AS, "_saved_session", lambda sid: None)
    monkeypatch.setattr(AS, "keep_session_refusal", lambda: "")
    return tmp_path


def _ctx(tmp_path):
    from delfin.dashboard.context import DashboardContext

    ctx = DashboardContext(calc_dir=tmp_path / "calc", agent_dir=tmp_path / "agent_workspace",
                           archive_dir=tmp_path / "archive", office_dir=tmp_path / "office")
    ctx.run_js = lambda script: None
    return ctx


@pytest.fixture
def service(tmp_path):
    svc = ServiceProcess(tmp_path)
    try:
        yield svc
    finally:
        svc.stop()


@pytest.fixture
def tui(service):
    c = ServiceClient(service.socket, client="delfin-tui")
    service.clients.append(c)
    return c


def _dashboard(home, service, build=None, **kw):
    """A dashboard whose sessions run in its kernel, hosted in ``service``."""
    client = LS.LiveClient(service.socket)
    service.clients.append(client)
    build = build or _Build()
    _w, refs = AS.create_tab(_ctx(home), build=build, service=lambda: client,
                             start_in="kernel", hosting={"fake": True}, **kw)
    return build, refs, client


def _summary(c, sid):
    return next((s for s in c.request("sessions.list")["sessions"]
                 if s["sessionId"] == sid), {})


def _give_id(refs, sid):
    rec = refs["sessions"]()[0]
    rec["refs"]["state"]["active_session_id"] = sid
    refs["refresh"]()
    return rec


SID_A = "aaaaaaaa11112222333344445555aaaa"


def test_a_session_with_an_id_is_hosted_and_listed_once(home, service, tui):
    build, refs, _client = _dashboard(home, service)
    rec = _give_id(refs, SID_A)
    summary = _until(lambda: _summary(tui, SID_A))
    assert summary["host"] == "dashboard" and summary["keep"] is False
    assert not rec.get("live"), "it runs in this kernel, as before"
    panel = refs["running_panel"]
    assert _until(lambda: panel.running.summary(SID_A))
    refs["refresh"]()
    assert SID_A not in {r["sid"] for r in panel.ordered()}, \
        "the session open here is not listed a second time"
    listing = tui.request("sessions.list")
    assert SID_A not in {e.get("sessionId") for e in listing["elsewhere"]}
    refs["stop_refresh"]()


def test_keep_session_is_the_services_switch(home, service, tui):
    build, refs, _client = _dashboard(home, service)
    rec = _give_id(refs, SID_A)
    assert _until(lambda: _summary(tui, SID_A))
    refs["set_keep"](rec["key"], True)
    assert _until(lambda: _summary(tui, SID_A).get("keep") is True)
    assert K.records() == [], "the service decides; no kept record besides"
    tui.request("session.setKeep", {"sessionId": SID_A, "keep": False})
    assert _until(lambda: (refs["refresh"]() or True) and rec["keep_switch"].value is False), \
        "switched off in the terminal, the switch here follows"
    refs["stop_refresh"]()


def test_another_dashboard_opens_it_as_a_view_and_its_prompt_runs_here(home, service, tui):
    build_a, refs_a, _a = _dashboard(home, service)
    rec_a = _give_id(refs_a, SID_A)
    assert _until(lambda: _summary(tui, SID_A))
    build_b, refs_b, _b = _dashboard(home, service)
    assert _until(lambda: refs_b["running"].summary(SID_A))
    view = refs_b["open_from_list"](SID_A)
    assert view is not None and view.get("live"), "a view of the same session"
    assert _summary(tui, SID_A).get("host") == "dashboard", "not taken over from A"
    assert not rec_a.get("live") and rec_a in refs_a["sessions"]()
    view["refs"]["input"].value = "sent from B"
    view["refs"]["send"]()
    assert _until(lambda: "sent from B" in build_a.tabs[0]["sent"])
    refs_a["stop_refresh"]()
    refs_b["stop_refresh"]()


def test_closed_with_keep_off_and_held_by_nobody_it_ends(home, service, tui):
    build, refs, _client = _dashboard(home, service)
    refs["open"]()
    rec = _give_id(refs, SID_A)
    assert _until(lambda: _summary(tui, SID_A))
    refs["close"](rec["key"])
    assert build.tabs[0]["closed"], "the tab saved and closed first"
    assert _until(lambda: not _summary(tui, SID_A))
    refs["stop_refresh"]()


def test_closed_with_keep_on_it_asks_and_keep_running_runs_on(home, service, tui):
    build, refs, _client = _dashboard(home, service)
    refs["open"]()
    rec = _give_id(refs, SID_A)
    assert _until(lambda: _summary(tui, SID_A))
    refs["set_keep"](rec["key"], True)
    refs["close"](rec["key"])
    assert rec["close_ask"].layout.display == "", "asked first"
    refs["close"](rec["key"], how="keep")
    moved = _until(lambda: (lambda s: s if s and not s.get("host") else None)(
        _summary(tui, SID_A)), 60)
    assert moved and moved["keep"] is True, "it runs on in the service"
    assert K.records() == [], "handed over by the service, not by a record"
    refs["stop_refresh"]()


def test_ended_from_the_terminal_it_is_saved_and_listed_as_not_running(home, service, tui):
    build, refs, _client = _dashboard(home, service)
    rec = _give_id(refs, SID_A)
    assert _until(lambda: _summary(tui, SID_A))
    tui.request("session.shutdown", {"sessionId": SID_A})
    assert _until(lambda: build.tabs[0]["closed"])
    assert _until(lambda: not _summary(tui, SID_A))
    rows = refs["sessions"]()
    assert rec not in rows
    assert [r.get("dormant") for r in rows if r.get("session_id") == SID_A] == [True]
    refs["stop_refresh"]()


def test_a_dashboard_session_without_an_id_is_no_elsewhere_row(home, service, tui):
    build, refs, _client = _dashboard(home, service)
    refs["refresh"]()
    from delfin.agent import session_presence
    assert _until(lambda: [r for r in session_presence.open_sessions()
                           if r.get("surface") == "dashboard"], 10), "it announces itself"
    listing = tui.request("sessions.list")
    assert not [e for e in listing["elsewhere"] if not e.get("sessionId")], \
        "a dashboard session with nothing to join is not listed elsewhere"
    refs["stop_refresh"]()


# -- the real agent tab's host adapter ------------------------------------------------

@pytest.fixture
def real_tab(tmp_path, monkeypatch):
    from delfin.agent import scheduler as S
    from delfin.dashboard import tab_agent
    from delfin.dashboard.context import DashboardContext

    monkeypatch.setattr(S, "_GLOBAL", S.Scheduler(path=tmp_path / "cron.json"))
    for name in ("calc", "archive", "office", "work"):
        (tmp_path / name).mkdir(exist_ok=True)
    ctx = DashboardContext(calc_dir=tmp_path / "calc", archive_dir=tmp_path / "archive",
                           office_dir=tmp_path / "office")
    ctx.run_js = lambda script: None
    ctx.agent_workspace = str(tmp_path / "work")
    _w, refs = tab_agent.create_tab(ctx)
    try:
        yield refs
    finally:
        refs["shutdown"]()


def test_the_real_tab_gives_its_own_callbacks(real_tab):
    adapter = real_tab["host_adapter"]()
    for key in ("state", "send", "stop", "answer", "plan", "approval", "question",
                "set_permission", "choices", "info", "chat_html", "status", "input"):
        assert key in adapter, key
    assert adapter["state"] is real_tab["state"]
    assert adapter["chat_html"] is real_tab["chat_html"]
    info = adapter["info"]()
    assert {"provider", "model", "effort", "mode", "approval"} <= set(info)


def test_a_line_from_another_window_runs_as_typed_and_keeps_the_draft(real_tab):
    adapter = real_tab["host_adapter"]()
    adapter["input"].value = "my half-written draft"
    assert adapter["send"]("/help") is True
    said = [m for m in real_tab["state"]["chat_messages"] if m.get("role") == "user"]
    assert said and said[-1]["content"] == "/help", "a slash command runs in the tab"
    assert adapter["input"].value == "my half-written draft"


def test_a_posture_from_another_window_moves_the_selector_bypass_only_confirmed(real_tab):
    adapter = real_tab["host_adapter"]()
    assert adapter["set_permission"]("plan", False) == "plan"
    assert adapter["info"]()["approval"] == "plan"
    assert adapter["set_permission"]("all_free", False) == "plan", \
        "Bypass is not switched on by one unconfirmed request"
    assert adapter["set_permission"]("all_free", True) == "all_free"


def test_the_pickers_values_are_offered(real_tab):
    adapter = real_tab["host_adapter"]()
    effort = adapter["choices"]("effort")
    assert effort["current"] and {c["value"] for c in effort["choices"]} >= {effort["current"]}
    assert adapter["choices"]("model")["choices"]
