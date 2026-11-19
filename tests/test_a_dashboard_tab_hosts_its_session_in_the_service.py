"""The kernel side of a hosted dashboard session
(delfin/dashboard/hosted_sessions.py, delfin/agent/service/host.py).

The agent tab is not changed to be hosted: its chat, its streaming flag
and its open questions are read and said as protocol events
(ChatFollower); its chat HTML goes whole to the person's other dashboards
(dashboard.chat diffs); a request from any window is run by the tab's own
callbacks. Driven against a real private service with a stand-in tab.
"""
from __future__ import annotations

import random
import time

import pytest

from delfin.agent.service import protocol
from delfin.agent.service.client import ServiceClient, ServiceError
from delfin.dashboard import hosted_sessions as H
from tui_service import ServiceProcess

SID = "fedcba9876543210fedcba9876543210"


# -- what the tab says, as events --------------------------------------------------

class _Said:
    def __init__(self):
        self.events: list[tuple[str, dict]] = []

    def __call__(self, type_, **fields):
        self.events.append((type_, fields))

    def types(self):
        return [t for t, _ in self.events]

    def of(self, type_):
        return [f for t, f in self.events if t == type_]


def test_a_turn_is_said_in_order_and_its_answer_streams():
    said = _Said()
    f = H.ChatFollower(said)
    msgs = [{"role": "user", "content": "hello"}]
    f.step(msgs, True, [])
    answer = {"role": "assistant", "content": "Hel", "_streaming": True}
    msgs.append(answer)
    f.step(msgs, True, [])
    answer["content"] = "Hello there"
    f.step(msgs, True, [])
    answer["_streaming"] = False
    msgs.append({"role": "tool", "content": '<span class="tool-name">Read</span>  a.py'})
    msgs.append({"role": "assistant", "content": "Done."})
    f.step(msgs, False, [])
    assert said.types() == [
        "user.message", "turn.started", "session.state", "text.started", "text.delta",
        "text.delta", "text.ended", "tool.started", "tool.ended", "text.started",
        "text.ended", "turn.ended", "session.state"], said.types()
    assert said.of("user.message")[0] == {"text": "hello", "source": "user"}
    assert "".join(d["delta"] for d in said.of("text.delta")) == "Hello there"
    assert said.of("text.ended")[0]["text"] == "Hello there"
    assert said.of("tool.started")[0]["name"] == "Read"
    assert [s["state"] for s in said.of("session.state")] == ["running", "idle"]


def test_a_live_line_waits_until_something_follows_it():
    said = _Said()
    f = H.ChatFollower(said)
    msgs = [{"role": "user", "content": "go"},
            {"role": "system", "content": "waiting 3 s", "_live_key": "wait"}]
    f.step(msgs, True, [])
    assert "notice" not in said.types()
    msgs[1]["content"] = "waited 9 s"
    msgs.append({"role": "assistant", "content": "ok"})
    f.step(msgs, False, [])
    assert said.of("notice") == [{"text": "waited 9 s"}]


def test_a_turn_nobody_typed_is_the_harness():
    said = _Said()
    f = H.ChatFollower(said)
    f.step([{"role": "user", "content": "job 7 finished", "origin": "event"}], True, [])
    assert said.of("user.message")[0]["source"] == "harness"
    assert said.of("turn.started")[0]["source"] == "harness"


def test_open_questions_come_and_go_and_name_who_answered():
    said = _Said()
    f = H.ChatFollower(said)
    asked = {"requestId": "r1", "kind": "confirm", "tool": "dashboard", "preview": "rm x",
             "options": [{"key": "y", "label": "Approve"}]}
    f.step([], False, [asked])
    assert said.of("permission.asked")[0]["requestId"] == "r1"
    assert said.of("session.state")[-1]["state"] == "needs_input"
    f.repliers["r1"] = "delfin-tui"
    f.step([], False, [])
    assert said.of("permission.resolved") == [{"requestId": "r1", "decision": "answered",
                                               "by": "client", "client": "delfin-tui"}]
    assert said.of("session.state")[-1]["state"] == "idle"


def test_history_is_said_once():
    said = _Said()
    history = [{"role": "user", "content": "before"}, {"role": "assistant", "content": "yes"}]
    H.ChatFollower(said, published=len(history)).step(history, False, [])
    assert said.events == [], "a session whose log has its history says nothing again"
    H.ChatFollower(said).step(history, False, [])
    assert said.types() == ["user.message", "text.started", "text.ended"]


def test_a_cleared_chat_is_said():
    said = _Said()
    f = H.ChatFollower(said, published=3)
    f.step([], False, [])
    assert said.of("notice") and f.published == 0


def test_the_chat_diff_rebuilds_the_new_html_exactly():
    rng = random.Random(7)
    for _ in range(300):
        old = "".join(rng.choice("ab<>/ ") for _ in range(rng.randint(0, 60)))
        new = old[:rng.randint(0, len(old))] + "".join(
            rng.choice("abc") for _ in range(rng.randint(0, 20)))
        keep, tail = H.chat_diff(old, new)
        assert old[:keep] + tail == new
        assert keep == len(old) or keep == len(new) or old[keep] != new[keep]


def test_a_tool_line_gives_its_name_and_text():
    assert H.html_text('<span class="tool-name">$</span> ls -la') == ("$", "$ ls -la")
    assert H.html_text("<b>Plan</b> &amp; go") == ("", "Plan & go")


# -- hosted, against a real service -------------------------------------------------

@pytest.fixture
def service(tmp_path):
    svc = ServiceProcess(tmp_path)
    try:
        yield svc
    finally:
        svc.stop()


class _Chat:
    def __init__(self, value=""):
        self.value = value


class _Tab:
    """What tab_agent.create_tab's host_adapter gives: the state and the
    tab's own callbacks."""

    def __init__(self):
        self.state = {"active_session_id": SID, "chat_messages": [], "streaming": False,
                      "engine": None}
        self.sent: list[str] = []
        self.stops = 0
        self.approvals: list[bool] = []
        self.profiles: list[tuple] = []
        self.chat = _Chat('<template class="delfin-chat-src"><div>empty</div></template>')

    def adapter(self) -> dict:
        def send(text):
            self.sent.append(text)
            self.state["chat_messages"].append({"role": "user", "content": text})
            return True

        def approval(yes):
            self.approvals.append(yes)
            self.state.pop("_pending_approval", None)

        def set_permission(profile, confirmed):
            self.profiles.append((profile, confirmed))
            return profile

        return {"state": self.state, "send": send, "stop": self._stop,
                "answer": lambda text: False, "plan": lambda yes: None,
                "approval": approval, "question": lambda text: False,
                "set_permission": set_permission,
                "choices": lambda what: {"current": "m1", "choices": [{"value": "m1"}]},
                "info": lambda: {"provider": "kit", "model": "m1", "effort": "medium",
                                 "mode": "solo", "approval": "default"},
                "chat_html": self.chat, "status": lambda: {}}

    def _stop(self):
        self.stops += 1


def _until(probe, timeout=20.0):
    deadline = time.monotonic() + timeout
    while time.monotonic() < deadline:
        got = probe()
        if got:
            return got
        time.sleep(0.05)
    return probe()


def _summary(c, sid=SID):
    return next((s for s in c.request("sessions.list")["sessions"]
                 if s["sessionId"] == sid), {})


def _host(service, tab, **kw):
    kw.setdefault("title", "hosted tab")
    host = H.TabHost(tab.adapter(), sid=SID, cwd=str(service.workspace()),
                     path=service.socket, **kw)
    assert host.start(), host.why
    return host


def test_the_tab_is_listed_joined_and_driven_from_the_terminal(service):
    tab = _Tab()
    tab.state["chat_messages"].append({"role": "user", "content": "earlier"})
    host = _host(service, tab, keep=True)
    tui = ServiceClient(service.socket, client="delfin-tui")
    service.clients.append(tui)
    summary = _until(lambda: _summary(tui))
    assert summary["host"] == "dashboard" and summary["name"] == "hosted tab"
    got = tui.request("session.attach", {"sessionId": SID, "turns": 0})
    assert [e["text"] for e in got["events"] if e["type"] == "user.message"] == ["earlier"], \
        "a session new to the service brings its history once"
    start = len(tui.events)
    tui.request("session.prompt", {"sessionId": SID, "text": "/model m1"})
    assert tab.sent == ["/model m1"], "a slash command from a window runs in the tab"
    tui.wait_for(lambda e: e.get("type") == "user.message" and e.get("text") == "/model m1",
                 after=start)
    tui.request("session.prompt", {"sessionId": SID, "text": "/keep off"})
    assert tab.sent == ["/model m1"], "/keep is the service's switch, not a message"
    tui.request("session.cancel", {"sessionId": SID})
    assert tab.stops == 1
    tui.request("session.setMode", {"sessionId": SID, "mode": "acceptEdits"})
    assert tab.profiles == [("repo_free", False)]
    host.close()


def test_an_approval_in_the_tab_is_answered_from_another_window(service):
    tab = _Tab()
    host = _host(service, tab, keep=True)
    tui = ServiceClient(service.socket, client="delfin-tui")
    service.clients.append(tui)
    tui.request("session.attach", {"sessionId": SID})
    start = len(tui.events)
    tab.state["_pending_approval"] = {"tool": "Bash", "detail": "rm -rf build"}
    asked, _ = tui.wait_for(lambda e: e.get("type") == "permission.asked", after=start)
    tui.request("permission.reply", {"sessionId": SID, "requestId": asked["requestId"],
                                     "choice": "y"})
    assert tab.approvals == [True]
    resolved, _ = tui.wait_for(lambda e: e.get("type") == "permission.resolved", after=start)
    assert resolved["client"] == "delfin-tui"
    host.close()


def test_a_command_from_the_phone_is_not_run_in_the_tab(service):
    tab = _Tab()
    host = _host(service, tab, keep=True)
    bridge = ServiceClient(service.socket, client="delfin-messenger-bridge")
    service.clients.append(bridge)
    with pytest.raises(ServiceError) as refused:
        bridge.request("session.prompt", {"sessionId": SID, "text": "/clear"})
    assert refused.value.code == protocol.NOT_ALLOWED
    bridge.request("session.prompt", {"sessionId": SID, "text": "status please"})
    assert tab.sent == ["status please"]
    host.close()


def test_another_dashboard_gets_the_chat_whole_then_its_changes(service):
    tab = _Tab()
    host = _host(service, tab, keep=True)
    dash = ServiceClient(service.socket, client="delfin-dashboard")
    service.clients.append(dash)
    got = _until(lambda: [e for e in dash.request("session.attach", {"sessionId": SID})["events"]
                          if e.get("type") == "dashboard.chat"])
    shown = {"v": got[-1]["v"], "html": got[-1]["html"]}
    assert shown["html"] == tab.chat.value
    start = len(dash.events)
    tab.chat.value = tab.chat.value.replace("empty", "a 3D card <div class='mol'></div>")
    event, _ = dash.wait_for(lambda e: e.get("type") == "dashboard.chat", after=start)
    assert event["base"] == shown["v"]
    assert shown["html"][:event["keep"]] + event["tail"] == tab.chat.value
    host.close()


def test_ended_from_another_window_the_tab_is_told(service):
    tab = _Tab()
    ended = []
    host = _host(service, tab, keep=True, on_end=lambda: ended.append(True))
    tui = ServiceClient(service.socket, client="delfin-tui")
    service.clients.append(tui)
    assert tui.request("session.shutdown", {"sessionId": SID, "dryRun": True})["closing"] is False
    assert tui.request("session.shutdown", {"sessionId": SID})["closing"] is True
    assert _until(lambda: ended)
    host.close()


def test_released_off_and_alone_it_ends_and_a_rehost_says_no_history_twice(service):
    tab = _Tab()
    tab.state["chat_messages"].append({"role": "user", "content": "once"})
    host = _host(service, tab, keep=False)
    watcher = service.client()
    assert _until(lambda: _summary(watcher))
    host.release(end=False)
    assert _until(lambda: not _summary(watcher)), "keep off, nobody holds it: it ends"
    again = _host(service, tab, keep=True)
    got = watcher.request("session.attach", {"sessionId": SID, "turns": 0})
    assert [e["text"] for e in got["events"] if e["type"] == "user.message"] == ["once"]
    again.close()


def test_a_lost_link_is_said_to_the_tab(service):
    tab = _Tab()
    lost = []
    host = _host(service, tab, keep=True, on_lost=lambda: lost.append(True))
    with ServiceClient(service.socket) as c:
        c.request("service.shutdown", timeout=10)
    assert _until(lambda: lost, 30)
    assert not host.hosted
