"""The dashboard is a window of the agent service, like the TUI and the phone.

``delfin.dashboard.live_sessions`` is the dashboard's client of the
service: it lists every live session of the person's service (all
folders), attaches to one (its recent turns replayed, then live), sends
prompts, answers approvals and questions, renames, detaches, and stops a
session only when asked to. ``live_tab`` draws one such session;
``agent_sessions`` starts new sessions in the service and falls back to
the kernel, saying so, when there is none.

Driven against a real private service (``tui_service.ServiceProcess``)
with the scripted engine: no model, no network, no external agent CLI.
The "TUI" here is a plain service client named as the TUI names itself.
"""

from __future__ import annotations

import time

import pytest

pytest.importorskip("ipywidgets")

import ipywidgets as widgets                                     # noqa: E402

from delfin.agent.service import protocol                        # noqa: E402
from delfin.agent.service.client import ServiceClient, ServiceError  # noqa: E402
from delfin.dashboard import live_sessions as LS                 # noqa: E402
from tui_service import ServiceProcess, turn                     # noqa: E402


def _until(cond, timeout: float = 30.0):
    deadline = time.monotonic() + timeout
    while time.monotonic() < deadline:
        got = cond()
        if got:
            return got
        time.sleep(0.05)
    return cond()


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


@pytest.fixture
def dash(service):
    c = LS.LiveClient(service.socket)
    yield c
    c.close()


def _texts(transcript: LS.Transcript, role: str) -> list[str]:
    return [e["text"] for e in transcript.entries if e["role"] == role]


# -- the client ----------------------------------------------------------------------

def test_the_dashboard_lists_every_live_session_of_its_service(service, tui, dash):
    one = service.session(tui, cwd=str(service.workspace("scan-a")))
    two = service.session(tui, cwd=str(service.workspace("scan-b")), title="the B scan")
    rows = {r["sessionId"]: r for r in LS.RunningList.of(dash.sessions()).rows()}
    assert set(rows) == {one, two}
    assert rows[one]["name"] == f"scan-a · {one[:8]}"
    assert rows[two]["name"] == "the B scan"
    assert rows[two]["folder"] == "scan-b"
    assert rows[one]["state"] == "done" and rows[one]["model"] == rows[one]["model"]
    assert dash.hello["mayAct"] is True


def test_attaching_to_a_tui_session_replays_its_turns_then_follows_live(service, tui, dash):
    sid = service.session(tui)
    turn(tui, sid, "the first question")
    view = LS.Transcript(sid)
    dash.subscribe(view.apply)
    view.load(dash.attach(sid))
    assert _texts(view, "user") == ["the first question"]
    assert _texts(view, "agent"), view.entries
    turn(tui, sid, "a second one, typed in the terminal")
    assert _until(lambda: _texts(view, "user") == ["the first question",
                                                   "a second one, typed in the terminal"])
    assert _until(lambda: len(_texts(view, "agent")) == 2)
    seqs = [e["seq"] for e in view.entries if e.get("seq")]
    assert len(seqs) == len(set(seqs))                     # nothing twice


def test_a_session_the_dashboard_starts_opens_in_the_tui_like_its_own(service, tui, dash):
    started = dash.create(str(service.workspace("from-dash")), fake=True, fakeDelay=0.003)
    sid = started["sessionId"]
    listed = tui.request("sessions.list")
    assert sid in {s["sessionId"] for s in listed["sessions"]}
    assert sid not in {e.get("sessionId") for e in listed["elsewhere"]}
    tui.request("session.attach", {"sessionId": sid})
    mark = len(tui.events)
    dash.prompt(sid, "typed in the dashboard")
    tui.wait_for(lambda e: e.get("type") == "user.message" and e.get("session") == sid
                 and e.get("text") == "typed in the dashboard", after=mark)
    tui.wait_for(lambda e: e.get("type") == "turn.ended" and e.get("session") == sid, after=mark)


def test_either_window_types_and_both_see_the_same_stream(service, tui, dash):
    sid = service.session(tui)
    view = LS.Transcript(sid)
    dash.subscribe(view.apply)
    view.load(dash.attach(sid))
    mark = len(tui.events)
    dash.prompt(sid, "from the dashboard")
    tui.wait_for(lambda e: e.get("type") == "turn.ended" and e.get("session") == sid, after=mark)
    turn(tui, sid, "from the terminal")
    assert _until(lambda: _texts(view, "user") == ["from the dashboard", "from the terminal"])
    shown = [e for e in tui.events[mark:] if e.get("session") == sid and e.get("type") == "text.ended"]
    assert _until(lambda: len(_texts(view, "agent")) == 2)
    assert _texts(view, "agent")[0] == shown[0]["text"]


def test_the_first_approval_wins_and_the_other_window_is_told_who(service, tui, dash):
    sid = service.session(tui, permissionMode="default")
    view = LS.Transcript(sid)
    dash.subscribe(view.apply)
    view.load(dash.attach(sid))
    tui.request("session.attach", {"sessionId": sid})
    mark = len(tui.events)
    tui.request("session.prompt", {"sessionId": sid, "text": "edit it"})
    asked = _until(lambda: next(iter(view.pending.values()), None))
    assert asked and view.needs_you
    dash.reply(sid, asked["requestId"], choice="y")
    with pytest.raises(ServiceError) as late:
        tui.request("permission.reply", {"sessionId": sid, "requestId": asked["requestId"],
                                         "choice": "n"})
    assert late.value.code == protocol.UNKNOWN_REQUEST
    resolved, _ = tui.wait_for("permission.resolved", after=mark)
    assert resolved["client"] == "delfin-dashboard"
    assert _until(lambda: not view.pending)
    said = [e for e in view.entries if e["role"] == "resolved"]
    assert said and said[-1]["client"] == "delfin-dashboard" and "dashboard" in said[-1]["text"]


def test_an_approval_answered_in_the_tui_closes_in_the_dashboard(service, tui, dash):
    sid = service.session(tui, permissionMode="default")
    view = LS.Transcript(sid)
    dash.subscribe(view.apply)
    view.load(dash.attach(sid))
    mark = len(tui.events)
    tui.request("session.prompt", {"sessionId": sid, "text": "edit it"})
    asked, _ = tui.wait_for("permission.asked", after=mark)
    assert _until(lambda: asked["requestId"] in view.pending)
    tui.request("permission.reply", {"sessionId": sid, "requestId": asked["requestId"], "choice": "n"})
    assert _until(lambda: not view.pending)
    said = [e for e in view.entries if e["role"] == "resolved"]
    assert said[-1]["client"] == "delfin-tui" and "terminal" in said[-1]["text"]
    with pytest.raises(ServiceError):
        dash.reply(sid, asked["requestId"], choice="y")


def test_a_question_is_answered_from_the_dashboard(service, tui, dash):
    sid = service.session(tui)
    view = LS.Transcript(sid)
    dash.subscribe(view.apply)
    view.load(dash.attach(sid))
    dash.prompt(sid, "ask me")
    asked = _until(lambda: next(iter(view.pending.values()), None))
    assert asked["kind"] == "ask" and asked["choices"]
    dash.reply(sid, asked["requestId"], answers=[asked["choices"][0]])
    assert _until(lambda: not view.pending)


def test_closing_the_dashboard_view_never_stops_the_session(service, tui, dash):
    sid = service.session(tui, fakeDelay=0.02)
    dash.attach(sid)
    mark = len(tui.events)
    tui.request("session.prompt", {"sessionId": sid, "text": "slow please"})
    tui.wait_for("text.started", after=mark)
    dash.detach(sid)
    dash.close()
    tui.wait_for(lambda e: e.get("type") == "turn.ended" and e.get("session") == sid, after=mark)
    assert sid in {s["sessionId"] for s in tui.request("sessions.list")["sessions"]}


def test_stop_is_explicit_and_ends_the_session_for_every_window(service, tui, dash):
    sid = service.session(tui)
    mark = len(tui.events)
    dash.stop(sid)
    tui.wait_for(lambda e: e.get("type") == "session.removed" and e.get("session") == sid, after=mark)


def test_the_running_list_follows_the_service_without_asking_again(service, tui, dash):
    running = LS.RunningList.of(dash.sessions())
    dash.subscribe(running.apply)
    asked = []
    real = dash.request

    def counting(method, params=None, **kw):
        asked.append(method)
        return real(method, params, **kw)

    dash.request = counting
    sid = service.session(tui, cwd=str(service.workspace("late")), permissionMode="default")
    assert _until(lambda: sid in {r["sessionId"] for r in running.rows()})
    tui.request("session.setTitle", {"sessionId": sid, "title": "renamed in the TUI"})
    assert _until(lambda: any(r["name"] == "renamed in the TUI" for r in running.rows()))
    tui.request("session.prompt", {"sessionId": sid, "text": "edit it"})
    row = _until(lambda: next((r for r in running.rows() if r["sessionId"] == sid and r["needs"]), None))
    assert row["state"] == "needs you" and running.rows()[0]["sessionId"] == sid
    tui.request("session.shutdown", {"sessionId": sid})
    assert _until(lambda: sid not in {r["sessionId"] for r in running.rows()})
    assert "sessions.list" not in asked                    # events, not polling


def test_a_rename_from_the_dashboard_reaches_the_tui(service, tui, dash):
    sid = service.session(tui)
    mark = len(tui.events)
    assert dash.rename(sid, "  named\x1b[1m in the dashboard ") == "named in the dashboard"
    tui.wait_for(lambda e: e.get("type") == "session.summary" and e.get("session") == sid
                 and e["summary"]["name"] == "named in the dashboard", after=mark)
    assert dash.rename(sid, "") == ""
    tui.wait_for(lambda e: e.get("type") == "session.summary" and e.get("session") == sid
                 and e["summary"]["name"].endswith(sid[:8]), after=mark)


def test_states_read_as_four_words():
    assert LS.state_word({"state": "running", "pending": 0}) == "working"
    assert LS.state_word({"state": "starting", "pending": 0}) == "working"
    assert LS.state_word({"state": "idle", "pending": 1}) == "needs you"
    assert LS.state_word({"state": "needs_input", "pending": 0}) == "needs you"
    assert LS.state_word({"state": "idle", "pending": 0}) == "done"
    assert LS.state_word({"state": "sleeping", "pending": 0}) == "asleep"
    assert LS.age_text(5) == "now" and LS.age_text(125) == "2 min" and LS.age_text(7300) == "2 h"


# -- people mode ---------------------------------------------------------------------

def test_a_person_s_dashboard_reaches_only_their_own_sessions(tmp_path, monkeypatch):
    roots = [tmp_path / "anna", tmp_path / "ben"]
    for r in roots:
        r.mkdir()
    anna, ben = ServiceProcess(roots[0]), ServiceProcess(roots[1])
    try:
        mine = anna.session(anna.client(), title="anna's")
        theirs = ben.session(ben.client(), title="ben's")
        monkeypatch.setenv("DELFIN_RUNTIME_DIR", str(anna.runtime))
        client = LS.connect(start=False)
        try:
            assert {r["sessionId"] for r in LS.RunningList.of(client.sessions()).rows()} == {mine}
            with pytest.raises(ServiceError):
                client.attach(theirs)
            with pytest.raises(ServiceError):
                client.rename(theirs, "not yours")
        finally:
            client.close()
    finally:
        anna.stop()
        ben.stop()


# -- no service ----------------------------------------------------------------------

def test_without_a_service_connect_says_why(tmp_path, monkeypatch):
    runtime = tmp_path / "rt"
    runtime.mkdir()
    monkeypatch.setenv("DELFIN_RUNTIME_DIR", str(runtime))
    with pytest.raises(LS.ServiceUnavailable) as none:
        LS.connect(start=False)
    assert "agent service" in str(none.value)
    # A runtime directory that is not safely the user's: never started there.
    import os
    (runtime / f"delfin-{os.getuid()}").symlink_to(tmp_path)
    with pytest.raises(LS.ServiceUnavailable) as unsafe:
        LS.connect(start=True)
    assert "not a plain directory" in str(unsafe.value)


# -- the tab that shows one session -------------------------------------------------

def test_the_live_tab_sends_what_is_typed_and_answers_with_its_buttons(service, tui, dash):
    from delfin.dashboard import live_tab

    sid = service.session(tui, permissionMode="default")
    widget, refs = live_tab.build_live_tab(dash, sid)
    assert isinstance(widget, widgets.Widget)
    tui.request("session.attach", {"sessionId": sid})
    mark = len(tui.events)
    refs["input"].value = "edit it"
    refs["send_btn"].click()
    tui.wait_for(lambda e: e.get("type") == "user.message" and e.get("text") == "edit it", after=mark)
    assert refs["input"].value == ""
    buttons = _until(lambda: [b for b in refs["approval_box"].children[0].children[-1].children]
                     if refs["approval_box"].children else None)
    labels = [b.description for b in buttons]
    assert any("Yes" in label or "yes" in label for label in labels), labels
    yes = next(b for b in buttons if b.tooltip.startswith("y"))
    yes.click()
    resolved, _ = tui.wait_for("permission.resolved", after=mark)
    assert (resolved["decision"], resolved["client"]) == ("approved", "delfin-dashboard")
    assert _until(lambda: not refs["approval_box"].children)
    tui.wait_for(lambda e: e.get("type") == "turn.ended" and e.get("session") == sid, after=mark)
    refs["render_now"]()
    assert "edit it" in refs["chat"].value
    refs["close"]()


def test_the_live_tab_renames_and_stops_only_after_a_confirmation(service, tui, dash):
    from delfin.dashboard import live_tab

    sid = service.session(tui, cwd=str(service.workspace("rename-me")))
    _widget, refs = live_tab.build_live_tab(dash, sid)
    refs["rename_btn"].click()
    assert refs["rename_box"].layout.display != "none"
    refs["rename_input"].value = "a better name"
    refs["rename_ok"].click()
    assert _until(lambda: "a better name" in refs["name_html"].value)
    assert _until(lambda: next(s for s in tui.request("sessions.list")["sessions"]
                               if s["sessionId"] == sid)["name"] == "a better name")
    refs["stop_btn"].click()                               # asks first
    assert refs["stop_confirm"].layout.display != "none"
    assert sid in {s["sessionId"] for s in tui.request("sessions.list")["sessions"]}
    refs["stop_no"].click()
    assert refs["stop_confirm"].layout.display == "none"
    mark = len(tui.events)
    refs["stop_btn"].click()
    refs["stop_yes"].click()
    tui.wait_for(lambda e: e.get("type") == "session.removed" and e.get("session") == sid, after=mark)


def test_every_control_of_the_live_tab_reaches_the_session(service, tui, dash):
    from delfin.dashboard import live_tab

    sid = service.session(tui)
    _widget, refs = live_tab.build_live_tab(dash, sid)
    tui.request("session.attach", {"sessionId": sid})
    mark = len(tui.events)
    refs["perm"].value = "plan"
    tui.wait_for(lambda e: e.get("type") == "session.updated" and e.get("approval") == "plan", after=mark)
    assert _until(lambda: refs["transcript"].info.get("approval") == "plan")
    # The control shows what the session has, not what was clicked.
    assert _until(lambda: refs["perm"].value == "plan")
    # Bypass asks a second time before it is sent.
    mark = len(tui.events)
    refs["perm"].value = "bypassPermissions"
    assert "again" in refs["notice"].value and refs["perm"].value == "plan"
    time.sleep(0.3)
    assert not [e for e in tui.events[mark:] if e.get("type") == "session.updated"
                and e.get("approval") == "bypassPermissions"]
    mark = len(tui.events)
    refs["effort"].value = "high"
    tui.wait_for(lambda e: e.get("type") == "user.message" and e.get("text") == "/effort high", after=mark)


# -- the session list ----------------------------------------------------------------

class _Build:
    def __init__(self):
        self.contexts = []

    def __call__(self, ctx):
        state = {"active_session_id": ctx.initial_session_id or "",
                 "chat_messages": [], "streaming": False, "engine": None}
        ctx.agent_state = state
        self.contexts.append(ctx)
        return widgets.VBox(), {"state": state, "shutdown": lambda: None}


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


def test_without_the_service_sessions_run_in_the_kernel_and_the_page_says_so(home):
    from delfin.dashboard import agent_sessions as AS

    def no_service():
        raise LS.ServiceUnavailable("the agent service could not be started: test")

    build = _Build()
    _widget, refs = AS.create_tab(_ctx(home), build=build, service=no_service)
    assert len(build.contexts) == 1                         # the kernel's session
    assert "agent service" in refs["notice"].value and "kernel" in refs["notice"].value
    assert refs["service"]() is None


_FAKE = {"fake": True, "fakeDelay": 0.003}


def test_a_new_session_runs_in_the_service_and_any_window_opens_it(home, service, tui, dash):
    from delfin.dashboard import agent_sessions as AS

    build = _Build()
    _widget, refs = AS.create_tab(_ctx(home), build=build, service=lambda: dash,
                                  session_defaults=lambda: dict(_FAKE))
    assert build.contexts == []                             # nothing in the kernel
    (rec,) = refs["sessions"]()
    assert rec["live"] and rec["workspace"] == str(home / "launch")
    sid = rec["session_id"]
    listed = tui.request("sessions.list")
    assert sid in {s["sessionId"] for s in listed["sessions"]}
    # The running list names it and marks it open here.
    assert _until(lambda: any(r["sessionId"] == sid for r in refs["running"].rows()))
    # Closing the view leaves it running.
    refs["close"](rec["key"])
    time.sleep(0.3)
    assert sid in {s["sessionId"] for s in tui.request("sessions.list")["sessions"]}
    refs["stop_refresh"]()


def test_one_click_in_the_running_list_opens_a_tui_session_here(home, service, tui, dash):
    from delfin.dashboard import agent_sessions as AS

    other = service.session(tui, cwd=str(service.workspace("tui-made")))
    turn(tui, other, "started in the terminal")
    build = _Build()
    _widget, refs = AS.create_tab(_ctx(home), build=build, service=lambda: dash,
                                  session_defaults=lambda: dict(_FAKE))
    panel = refs["running_panel"]
    assert _until(lambda: other in panel.buttons)
    panel.buttons[other].click()
    rec = refs["active"]()
    assert rec["live"] and rec["session_id"] == other
    transcript = rec["refs"]["transcript"]
    assert [e["text"] for e in transcript.entries if e["role"] == "user"] == ["started in the terminal"]
    # A second click shows the same view, never a second one.
    panel.buttons[other].click()
    assert sum(1 for r in refs["sessions"]() if r.get("session_id") == other) == 1
    refs["stop_refresh"]()
