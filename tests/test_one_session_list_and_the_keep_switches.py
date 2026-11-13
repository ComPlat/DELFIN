"""The agent tab's session sidebar: one list, the resume search, and the
two keep switches (reported 2026-11-13).

* Keep dashboard (header): keeps this dashboard's kernel; its tooltip
  names what holds a dashboard that is not kept (the server's own grace
  numbers, a running turn).
* Keep session (each session's control row, after Stop): on, the session
  runs on after the dashboard closes -- an in-kernel session is handed to
  the agent service when its kernel ends (kept_agent_sessions), a service
  session stays kept there; off, it ends with the dashboard.
* One list: every session once, running first, then "Not running", then
  "Background runs"; no "Running sessions" panel of its own.
* The resume field is a search whose list pops over the page; a row
  resumes, the field empties and the list closes; the keys are taken by
  _RESUME_KEYS_JS (Up/Down, Enter, Esc).

Driven against a real private agent service (tui_service.ServiceProcess)
where the service matters, with the scripted engine.
"""
from __future__ import annotations

import time

import pytest

pytest.importorskip("ipywidgets")

import ipywidgets as widgets                                     # noqa: E402

from delfin.dashboard import agent_sessions as AS                # noqa: E402
from delfin.dashboard import kept_agent_sessions as K            # noqa: E402
from delfin.dashboard import live_sessions as LS                 # noqa: E402
from tui_service import ServiceProcess, turn                     # noqa: E402

_FAKE = {"fake": True, "fakeDelay": 0.003}


def _until(cond, timeout: float = 30.0):
    deadline = time.monotonic() + timeout
    while time.monotonic() < deadline:
        got = cond()
        if got:
            return got
        time.sleep(0.05)
    return cond()


def _descendants(w):
    for child in getattr(w, "children", ()):
        yield child
        yield from _descendants(child)


def _classes(w) -> tuple:
    return tuple(getattr(w, "_dom_classes", ()))


class _Build:
    """A stand-in agent tab: its state, and a control row with a Stop."""

    def __init__(self):
        self.contexts = []
        self.rows = []

    def __call__(self, ctx):
        state = {"active_session_id": ctx.initial_session_id or "",
                 "chat_messages": [], "streaming": False, "engine": None}
        ctx.agent_state = state
        self.contexts.append(ctx)
        mode = widgets.Dropdown(options=["dashboard"], description="Mode:")
        stop = widgets.Button(description="Stop")
        git = widgets.HBox()
        row = widgets.HBox([mode, stop, git])
        self.rows.append(row)
        return widgets.VBox([row]), {"state": state, "shutdown": lambda: None,
                                     "control_row": row, "stop_button": stop}


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
    from delfin.agent.service.client import ServiceClient
    c = ServiceClient(service.socket, client="delfin-tui")
    service.clients.append(c)
    return c


@pytest.fixture
def dash(service):
    c = LS.LiveClient(service.socket)
    yield c
    c.close()


# -- Keep dashboard ---------------------------------------------------------------

def test_the_header_switch_is_keep_dashboard_and_says_what_holds_it(monkeypatch):
    from delfin.dashboard import session as S

    monkeypatch.setenv("DELFIN_SESSION_GRACE_SECONDS", "90")
    monkeypatch.setenv("DELFIN_SESSION_DROPPED_SECONDS", "900")
    strip = S.build_status_strip()
    (toggle,) = [w for w in strip.children if isinstance(w, widgets.ToggleButton)]
    assert toggle.description == "Keep dashboard"
    tip = toggle.tooltip
    for words in ("90 s", "15 min", "turn that is running finishes first",
                  "Keep session"):
        assert words in tip, (words, tip)


# -- Keep session in the kernel ---------------------------------------------------

def test_the_switch_sits_after_stop_and_arms_the_record(home):
    build = _Build()
    _w, refs = AS.create_tab(_ctx(home), build=build, service=False)
    (rec,) = refs["sessions"]()
    row = build.rows[0]
    names = [getattr(w, "description", "") for w in row.children]
    assert names[:3] == ["Mode:", "Stop", "Keep session"], names
    switch = row.children[2]
    assert "delfin-keep-session" in _classes(switch) and switch.value is False
    assert "ends when this dashboard closes" in switch.tooltip
    rec["refs"]["state"]["active_session_id"] = "abcdef0123456789"
    switch.value = True
    assert K.is_armed("abcdef0123456789")
    assert "moves to the agent service" in switch.tooltip
    assert "delfin-kept-pill" in rec["info"].value
    switch.value = False
    assert not K.is_armed("abcdef0123456789")
    assert "delfin-kept-pill" not in rec["info"].value
    refs["stop_refresh"]()


def test_a_session_without_an_id_is_armed_once_it_has_one(home):
    build = _Build()
    _w, refs = AS.create_tab(_ctx(home), build=build, service=False)
    (rec,) = refs["sessions"]()
    refs["set_keep"](rec["key"], True)
    assert K.records() == []                       # nothing to name yet
    rec["refs"]["state"]["active_session_id"] = "1234567890abcdef"
    refs["refresh"]()
    assert K.is_armed("1234567890abcdef")
    refs["stop_refresh"]()


def test_closing_a_kept_session_ends_it_here(home):
    build = _Build()
    _w, refs = AS.create_tab(_ctx(home), build=build, service=False)
    refs["open"]()                                 # a second one, so a close leaves one
    rec = refs["sessions"]()[0]
    rec["refs"]["state"]["active_session_id"] = "feedfacefeedface"
    refs["set_keep"](rec["key"], True)
    assert K.is_armed("feedfacefeedface")
    refs["close"](rec["key"])
    assert not K.is_armed("feedfacefeedface")
    refs["stop_refresh"]()


def test_where_the_service_may_not_run_the_switch_says_why(home, monkeypatch):
    monkeypatch.setattr(AS, "keep_session_refusal", lambda: "login node: keep the dashboard")
    build = _Build()
    _w, refs = AS.create_tab(_ctx(home), build=build, service=False)
    (rec,) = refs["sessions"]()
    switch = rec["keep_switch"]
    assert switch.disabled and "login node" in switch.tooltip
    assert refs["set_keep"](rec["key"], True) is False
    assert K.records() == []
    refs["stop_refresh"]()


def test_a_kept_service_session_taken_over_stays_kept(home, service, tui, dash, monkeypatch):
    """The take-over into the full chat (kernel mode, the default) keeps
    the switch: a session kept in the service goes back there when this
    dashboard closes."""
    from delfin.dashboard import tab_agent

    sid = service.session(tui, cwd=str(service.workspace("from-tui")), attach=True)
    turn(tui, sid, "started in the terminal")
    build = _Build()
    monkeypatch.setattr(tab_agent, "create_tab", build)
    monkeypatch.setattr(LS, "connect", lambda **kw: dash)
    monkeypatch.delenv(AS.SESSIONS_ENV, raising=False)
    monkeypatch.setattr(AS, "_saved_session", lambda s: {"session_id": s} if s == sid else None)
    _w, refs = AS.create_tab(_ctx(home))
    assert _until(lambda: refs["running_panel"].running.summary(sid))
    rec = refs["open_from_list"](sid)
    assert rec is not None and not rec.get("live")
    assert rec["keep_switch"].value is True and rec["kept"]
    refs["stop_refresh"]()


# -- Keep session in the service --------------------------------------------------

def test_a_service_session_started_here_ends_with_this_dashboard(home, service, tui, dash):
    build = _Build()
    _w, refs = AS.create_tab(_ctx(home), build=build, service=lambda: dash,
                             session_defaults=lambda: dict(_FAKE))
    (rec,) = refs["sessions"]()
    sid = rec["session_id"]
    assert rec["live"] and rec["keep_switch"].value is False
    assert _until(lambda: (refs["running"].summary(sid) or {}).get("keep") is False)
    refs["set_keep"](rec["key"], True)
    assert _until(lambda: (LS.RunningList.of(tui.request("sessions.list"))
                           .summary(sid) or {}).get("keep") is True)
    assert rec["keep_switch"].value is True
    refs["stop_refresh"]()


def test_closing_the_view_of_a_session_not_kept_ends_it(home, service, tui, dash):
    build = _Build()
    _w, refs = AS.create_tab(_ctx(home), build=build, service=lambda: dash,
                             session_defaults=lambda: dict(_FAKE))
    (rec,) = refs["sessions"]()
    sid = rec["session_id"]
    mark = len(tui.events)
    refs["close"](rec["key"])
    tui.wait_for(lambda e: e.get("type") == "session.removed" and e.get("session") == sid,
                 after=mark)
    refs["stop_refresh"]()


# -- one list --------------------------------------------------------------------

def test_the_sidebar_has_one_list_and_no_running_panel(home):
    widget, refs = AS.create_tab(_ctx(home), build=_Build(), service=False)
    texts = [str(getattr(w, "value", "")) for w in _descendants(widget)
             if isinstance(w, widgets.HTML)]
    assert not any("Running sessions" in t for t in texts)
    assert not [w for w in _descendants(widget) if "delfin-running-list" in _classes(w)]
    sidebar = next(w for w in _descendants(widget) if "delfin-sessions" in _classes(w))
    order = [c for c in sidebar.children
             if {"delfin-resume-combo", "delfin-session-list"} & set(_classes(c))]
    assert ["delfin-resume-combo" in _classes(order[0]),
            "delfin-session-list" in _classes(order[1])] == [True, True]
    refs["stop_refresh"]()


def test_every_session_once_running_first_background_apart(home, service, tui, dash):
    asleep_ws = service.workspace("dozing")
    bg_ws = service.workspace("proj")
    (bg_ws / ".delfin" / "worktrees" / "delfin-wt-1234abcd").mkdir(parents=True)
    awake = service.session(tui, cwd=str(service.workspace("awake")))
    turn(tui, awake, "hello")
    background = service.session(tui, cwd=str(bg_ws / ".delfin" / "worktrees"
                                                / "delfin-wt-1234abcd"))
    turn(tui, background, "nightly")
    dozing = service.session(tui, cwd=str(asleep_ws))
    turn(tui, dozing, "later")
    build = _Build()
    _w, refs = AS.create_tab(_ctx(home), build=build, service=lambda: dash)
    panel = refs["running_panel"]
    assert _until(lambda: {awake, background, dozing} <= set(panel.buttons))
    # Asleep, as the service says of a session without a worker.
    summary = dict(panel.running.summary(dozing))
    panel.running.apply({"type": "session.summary", "session": dozing,
                         "summary": {**summary, "state": "sleeping"}})
    panel.render()
    box = refs["list_box"]
    kids = list(box.children)
    open_row = refs["sessions"]()[0]["row"]

    def at(sid):
        return kids.index(panel._rows[sid]["box"])

    labels = {w.value: kids.index(w) for w in kids if isinstance(w, widgets.HTML)}
    assert kids.index(open_row) < at(awake) < labels["Not running"] < at(dozing) \
        < labels["Background runs"] < at(background)
    assert "delfin-state-running" in _classes(open_row)
    assert "delfin-state-running" in _classes(panel._rows[awake]["box"])
    assert "delfin-state-stopped" in _classes(panel._rows[dozing]["box"])
    # Open here, it is the open row -- never a second one.
    panel.buttons[awake].click()
    _until(lambda: any(r.get("session_id") == awake for r in refs["sessions"]()))
    kids = list(box.children)
    assert panel._rows[awake]["box"] not in kids
    assert sum(1 for r in refs["sessions"]() if r.get("session_id") == awake) == 1
    refs["stop_refresh"]()


# -- the resume search ------------------------------------------------------------

def test_resume_offers_match_every_word_own_first():
    rows = [{"session_id": "a", "title": "Fix the failing test", "where": "/proj",
             "background": False},
            {"session_id": "b", "title": "nightly test run", "where": "/proj/.delfin/worktrees/x",
             "background": True},
            {"session_id": "c", "title": "Plot spectra", "where": "/data", "background": False}]
    own, bg = AS.resume_offers(rows, "TEST")
    assert [r["session_id"] for r in own] == ["a"] and [r["session_id"] for r in bg] == ["b"]
    own, bg = AS.resume_offers(rows, "fix test")
    assert [r["session_id"] for r in own] == ["a"] and bg == []
    own, bg = AS.resume_offers(rows, "")
    assert len(own) == 2 and len(bg) == 1
    own, bg = AS.resume_offers(rows, "", limit=2)
    assert len(own) == 2 and bg == []


def test_the_resume_list_pops_over_the_page_and_its_keys_are_bound():
    css = AS._SIDEBAR_CSS
    assert ".delfin-resume-pop { display: none !important; position: fixed;" in css
    assert ".delfin-resume-open .delfin-resume-pop { display: flex !important; }" in css
    js = AS._RESUME_KEYS_JS
    for key in ("'ArrowDown'", "'ArrowUp'", "'Enter'", "'Escape'"):
        assert key in js
    assert "delfin-resume-open" in js and "getBoundingClientRect" in js
    import inspect
    assert "_RESUME_KEYS_JS" in inspect.getsource(AS.create_tab)


def test_a_chosen_row_resumes_and_the_field_empties(home, monkeypatch):
    import delfin.agent.session_store as store
    monkeypatch.setattr(store, "list_sessions", lambda limit=50: [
        {"session_id": "front", "title": "Fix the failing test", "workspace": "/proj",
         "updated_at": 70}])
    monkeypatch.setattr(AS, "_saved_session", lambda sid: {"session_id": sid})
    build = _Build()
    _w, refs = AS.create_tab(_ctx(home), build=build, service=False)
    refs["resume_search"].value = "fix"
    (btn,) = [b for b in refs["resume_rows"].children if isinstance(b, widgets.Button)]
    btn.click()
    assert build.contexts[-1].initial_session_id == "front"
    assert refs["resume_search"].value == ""
    refs["resume_search"].value = "nothing like it"
    (empty,) = refs["resume_rows"].children
    assert "No saved session matches" in empty.value
    refs["stop_refresh"]()
