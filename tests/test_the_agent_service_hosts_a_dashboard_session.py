"""A session that runs in a dashboard's kernel, hosted in the agent service.

The kernel stays the turn driver; the service holds the session's
identity, transcript and liveness. A kernel says ``session.host`` on a
connection of its own, and from then on that connection is the session's
pipe: the service lists the session like any other (``host:
"dashboard"``), logs and fans out its events, forwards every window's
requests to it, and decides what happens when the host goes away --
keep on, or held by another window: it runs on in a worker; keep off and
held by nobody: it ends. The chat itself goes, as ``dashboard.chat``, to
dashboard windows only.

Driven against a real private service (tui_service.ServiceProcess) with a
host spoken by hand (fake_host.FakeHost).
"""
from __future__ import annotations

import os
import subprocess
import sys
import time

import pytest

from delfin.agent.service import protocol, supervisor
from delfin.agent.service.client import ServiceClient, ServiceError
from fake_host import FakeHost
from tui_service import ServiceProcess, turn

SID = "0123456789abcdef0123456789abcdef"
WAIT_S = supervisor.UNUSED_GRACE_S + 40.0


@pytest.fixture
def service(tmp_path):
    svc = ServiceProcess(tmp_path)
    try:
        yield svc
    finally:
        svc.stop()


def _summary(c, sid: str = SID) -> dict:
    rows = c.request("sessions.list")["sessions"]
    return next((s for s in rows if s["sessionId"] == sid), {})


def _until(probe, timeout: float = WAIT_S):
    deadline = time.monotonic() + timeout
    while time.monotonic() < deadline:
        got = probe()
        if got:
            return got
        time.sleep(0.2)
    return probe()


def _host(service, **params) -> FakeHost:
    params.setdefault("info", {"provider": "fake", "model": "fake-model",
                               "approval": "default", "mode": "solo"})
    return FakeHost(service.socket, SID, str(service.workspace()), **params)


def test_the_protocol_names_hosting():
    spec = protocol.SERVICE_REQUESTS["session.host"]
    assert spec["sessionId"] == "string" and spec["cwd"] == "string"
    assert spec["keep"] == "boolean?" and spec["info"] == "object?"
    assert "session.host" in supervisor.ACTING
    assert protocol.SESSION_SUMMARY["host"] == "string?"
    assert protocol.SERVICE_REQUESTS["session.detach"]["release"] == "boolean?"
    assert supervisor.SERVICE_VERSION >= 8


def test_a_hosted_session_is_listed_once_named_and_marked(service):
    host = _host(service, title="scan a", keep=False)
    assert host.answer["hosted"] is True and host.answer["keep"] is False
    window = service.client()
    summary = _until(lambda: _summary(window))
    assert summary["host"] == "dashboard"
    assert summary["name"] == "scan a" and summary["keep"] is False
    assert summary["model"] == "fake-model" and summary["state"] == "idle"
    assert summary["pid"] == os.getpid(), "the host's process, not a worker's"
    listing = window.request("sessions.list")
    assert [s["sessionId"] for s in listing["sessions"]].count(SID) == 1
    assert SID not in {e.get("sessionId") for e in listing["elsewhere"]}, \
        "one of the service's own: never an elsewhere row"
    host.close()


def test_a_window_joins_it_and_its_requests_reach_the_host(service):
    host = _host(service, keep=True)
    host.emit("user.message", text="what is 2+2", source="user")
    host.emit("turn.started", turnId="t1", source="user")
    host.emit("text.ended", turnId="t1", partId="p1", text="four")
    host.emit("turn.ended", turnId="t1", durationMs=5, inputTokens=1, outputTokens=1,
              costUsd=0.0, toolCalls=0)
    tui = ServiceClient(service.socket, client="delfin-tui")
    service.clients.append(tui)
    got = _until(lambda: (lambda a: a if any(e.get("type") == "text.ended"
                                             for e in a["events"]) else None)(
        tui.request("session.attach", {"sessionId": SID, "turns": 0})), 10)
    assert got, "the host's turn is replayed to a window that joins"
    assert [e["text"] for e in got["events"] if e.get("type") == "user.message"] == [
        "what is 2+2"]
    start = len(tui.events)
    host.emit("notice", text="live from the dashboard")
    tui.wait_for(lambda e: e.get("type") == "notice"
                 and e.get("text") == "live from the dashboard", after=start)
    tui.request("session.prompt", {"sessionId": SID, "text": "and 3+3?"})
    assert host.wait_asked("session.prompt")["text"] == "and 3+3?"
    host.emit("permission.asked", requestId="r1", kind="confirm", tool="bash",
              preview="ls", options=[{"key": "y", "label": "Yes"}], flags={})
    assert _until(lambda: _summary(tui).get("pending") == 1, 10)
    tui.request("permission.reply", {"sessionId": SID, "requestId": "r1", "choice": "y"})
    reply = host.wait_asked("permission.reply")
    assert reply["requestId"] == "r1" and reply["client"] == "delfin-tui"
    tui.request("session.cancel", {"sessionId": SID})
    host.wait_asked("session.cancel")
    tui.request("session.setTitle", {"sessionId": SID, "title": "renamed"})
    assert host.wait_asked("session.setTitle")["title"] == "renamed"
    host.close()


def test_a_second_host_is_refused(service):
    first = _host(service)
    with pytest.raises((ServiceError, ConnectionError)):
        _host(service)
    assert first.wait_closed(0.5) is False, "the first host keeps it"
    first.close()


def test_off_and_held_by_nobody_it_ends_with_its_host(service):
    host = _host(service, keep=False)
    watcher = service.client()
    assert _until(lambda: _summary(watcher), 10)
    host.close()
    assert _until(lambda: not _summary(watcher)), "it ended with its host"
    status = watcher.request("service.status")
    assert SID not in {(r.get("summary") or r).get("sessionId")
                       for r in status.get("dormant") or []}, "not kept asleep either"


def test_off_but_held_by_a_window_it_runs_on_in_a_worker(service):
    host = _host(service, keep=False, fake=True)
    tui = ServiceClient(service.socket, client="delfin-tui")
    service.clients.append(tui)
    tui.request("session.attach", {"sessionId": SID, "turns": 0})
    host.emit("host.released", end=False)
    host.close()
    moved = _until(lambda: (lambda s: s if s and not s.get("host") and s.get("state") == "idle"
                            else None)(_summary(tui)))
    assert moved, "the window that holds it keeps it: it runs on in a worker"
    assert moved["pid"] != os.getpid() and moved["keep"] is False
    turn(tui, SID, "hello from the terminal")


def test_kept_it_runs_on_in_a_worker_once_its_host_is_gone(service):
    host = _host(service, keep=True, fake=True, title="kept one")
    watcher = service.client()
    host.emit("host.released", end=False)
    host.close()
    moved = _until(lambda: (lambda s: s if s and not s.get("host") and s.get("state") == "idle"
                            else None)(_summary(watcher)))
    assert moved and moved["keep"] is True and moved["name"] == "kept one"


def test_a_host_that_dies_without_a_word_is_waited_for_then_kept(service):
    code = ("import sys, time; sys.path.insert(0, sys.argv[1]);"
            "from fake_host import FakeHost;"
            f"h = FakeHost(sys.argv[2], {SID!r}, sys.argv[3], keep=True, fake=True);"
            "print('hosted', flush=True); time.sleep(600)")
    from conftest import child_env
    env = child_env(service.tmp)
    env["PYTHONPATH"] = service.env["PYTHONPATH"]
    child = subprocess.Popen([sys.executable, "-c", code, str(os.path.dirname(__file__)),
                              str(service.socket), str(service.workspace())],
                             stdout=subprocess.PIPE, env=env)
    try:
        assert child.stdout.readline().strip() == b"hosted"
        watcher = service.client()
        assert _until(lambda: _summary(watcher).get("host") == "dashboard", 10)
        child.kill()
        child.wait(10)
        moved = _until(lambda: (lambda s: s if s and not s.get("host") else None)(
            _summary(watcher)))
        assert moved and moved["keep"] is True
    finally:
        if child.poll() is None:
            child.kill()
            child.wait(10)


def test_ended_in_its_dashboard_it_ends_for_every_window(service):
    host = _host(service, keep=True, fake=True)
    tui = ServiceClient(service.socket, client="delfin-tui")
    service.clients.append(tui)
    tui.request("session.attach", {"sessionId": SID, "turns": 0})
    host.emit("host.released", end=True)
    host.close()
    assert _until(lambda: not _summary(tui))


def test_the_service_stopping_neither_ends_nor_keeps_a_hosted_session(tmp_path):
    svc = ServiceProcess(tmp_path)
    try:
        host = FakeHost(svc.socket, SID, str(svc.workspace()), keep=True)
        with ServiceClient(svc.socket) as c:
            c.request("service.shutdown", timeout=10)
        assert host.wait_closed(30), "the link closes with the service"
        assert host.asked_for("session.shutdown") == [], \
            "the dashboard's session is not ended by the service stopping"
        svc.proc.wait(60)
        svc.start()
        with ServiceClient(svc.socket) as c:
            assert not _summary(c), "a later service does not wake what a kernel runs"
    finally:
        svc.stop()


def test_an_unused_hosted_session_is_neither_ended_nor_put_to_sleep(tmp_path):
    svc = ServiceProcess(tmp_path, start=False)
    svc.env["DELFIN_AGENT_SLEEP_AFTER_S"] = "1"
    svc.start()
    try:
        host = FakeHost(svc.socket, SID, str(svc.workspace()), keep=True)
        time.sleep(supervisor.UNUSED_GRACE_S + 4.0)
        with ServiceClient(svc.socket) as c:
            assert _summary(c).get("host") == "dashboard"
        assert host.asked_for("session.shutdown") == []
        host.close()
    finally:
        svc.stop()


def test_a_process_under_the_host_may_watch_but_not_act(service):
    host = _host(service, keep=True)
    code = ("import sys; from delfin.agent.service.client import ServiceClient, ServiceError\n"
            "c = ServiceClient(sys.argv[1])\n"
            "try:\n"
            f"    c.request('session.prompt', {{'sessionId': {SID!r}, 'text': 'x'}})\n"
            "    print('acted')\n"
            "except ServiceError as e:\n"
            "    print('refused', e.code)\n")
    from conftest import child_env
    env = child_env(service.tmp)
    env["PYTHONPATH"] = service.env["PYTHONPATH"]
    out = subprocess.run([sys.executable, "-c", code, str(service.socket)], env=env,
                         capture_output=True, text=True, timeout=60)
    assert out.stdout.strip() == f"refused {protocol.NOT_ALLOWED}", out
    assert host.asked_for("session.prompt") == []
    own = service.client()
    own.request("session.prompt", {"sessionId": SID, "text": "from the kernel's own window"})
    assert host.wait_asked("session.prompt")["text"] == "from the kernel's own window"
    host.close()


def test_the_chat_goes_whole_to_a_dashboard_and_not_to_the_terminal(service):
    host = _host(service, keep=True)
    host.emit("dashboard.chat", v=1, base=0, keep=0, tail="<div>one</div>")
    dash = ServiceClient(service.socket, client="delfin-dashboard")
    tui = ServiceClient(service.socket, client="delfin-tui")
    service.clients += [dash, tui]

    def chat_in(answer):
        return [e for e in answer["events"] if e.get("type") == "dashboard.chat"]

    got = _until(lambda: chat_in(dash.request("session.attach", {"sessionId": SID})), 10)
    assert got and got[-1]["html"] == "<div>one</div>" and got[-1]["v"] == 1
    assert chat_in(tui.request("session.attach", {"sessionId": SID})) == []
    start_d, start_t = len(dash.events), len(tui.events)
    host.emit("dashboard.chat", v=2, base=1, keep=len("<div>one"), tail=" two</div>")
    event, _ = dash.wait_for(lambda e: e.get("type") == "dashboard.chat", after=start_d)
    assert (event["v"], event["base"], event["keep"], event["tail"]) == (
        2, 1, len("<div>one"), " two</div>")
    host.emit("notice", text="after the chat")
    tui.wait_for(lambda e: e.get("type") == "notice", after=start_t)
    assert not [e for e in tui.events[start_t:] if e.get("type") == "dashboard.chat"]
    host.close()


def test_a_window_that_lets_go_releases_its_hold(service):
    window = service.client()
    sid = service.session(window, keep=False)
    turn(window, sid, "hello")
    window.request("session.detach", {"sessionId": sid, "release": True})
    assert _until(lambda: not _summary(window, sid)), \
        "its only holder let go: a session switched off ends"


_CHILD = r"""
import ctypes, sys
# What process_guard.protect does in a dashboard kernel: its /proc
# environment becomes unreadable to the user's other processes.
assert ctypes.CDLL(None).prctl(4, 0, 0, 0, 0) == 0
from delfin.agent.service.client import ServiceClient, ServiceError
c = ServiceClient(sys.argv[1])
try:
    c.request("session.prompt", {"sessionId": "0" * 32, "text": "x"})
except ServiceError as e:
    print("code", e.code)
"""

_PARENT = r"""
import subprocess, sys
# argv names what started the child: a dashboard server is a Voila process.
out = subprocess.run([sys.executable, "-c", sys.argv[1], sys.argv[2], "ipykernel_launcher"],
                     capture_output=True, text=True, timeout=120)
print(out.stdout.strip() or out.stderr[-500:])
"""


@pytest.mark.skipif(not sys.platform.startswith("linux"), reason="PR_SET_DUMPABLE is Linux")
@pytest.mark.parametrize("server, acts", [("voila", True), ("some-tool", False)])
def test_a_protected_dashboard_kernel_is_a_window(service, server, acts):
    """A dashboard kernel makes itself unreadable (process_guard) before it
    connects for a session that gets its id later: unreadable is no reason
    to refuse it when it is a kernel a Voila server started. Control: the
    same protected process under another parent stays watch-only."""
    from conftest import child_env
    env = child_env(service.tmp)
    env["PYTHONPATH"] = service.env["PYTHONPATH"]
    env.pop("DELFIN_PROCESS_GUARD", None)
    out = subprocess.run([sys.executable, "-c", _PARENT, _CHILD, str(service.socket), server],
                         env=env, capture_output=True, text=True, timeout=180)
    code = int(out.stdout.strip().split()[-1]) if "code" in out.stdout else None
    if acts:
        assert code == protocol.UNKNOWN_SESSION, out
    else:
        assert code == protocol.NOT_ALLOWED, out


@pytest.mark.parametrize("after", [0.0, 0.05, 0.15, 0.3])
def test_off_and_let_go_while_a_worker_takes_it_over_it_ends(service, after):
    """Measured in the real run: the last window let go while the worker
    that takes the session over was still starting; the wake then brought
    it up held by nobody, listed as kept, and it ran on. Let go at several
    moments of the take-over."""
    host = _host(service, keep=False, fake=True)
    tui = ServiceClient(service.socket, client="delfin-tui")
    tui.request("session.attach", {"sessionId": SID, "turns": 0})
    host.emit("host.released", end=False)
    host.close()
    time.sleep(after)
    tui.close()
    watcher = service.client()
    assert _until(lambda: not _summary(watcher), 20)
    time.sleep(4.0)                    # and it does not come back with the wake
    assert not _summary(watcher), \
        "let go by every window, a session switched off ends: " + str(_summary(watcher))
