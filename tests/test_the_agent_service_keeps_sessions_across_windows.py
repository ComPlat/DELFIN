"""Sessions live in the agent service and outlive the windows that watch them.

The service (``delfin.agent.service.supervisor``) runs one worker per
session and fans its events out to every window attached to it. These
tests drive a private service with the scripted engine:

* two windows on one session see the same turn, and the first answer to
  a question wins while the other window is told;
* a window that closes leaves its session working, and a window that
  attaches later gets the session's last turns first and the rest live,
  nothing twice;
* three sessions survive a disconnect and are all still there after a
  reconnect (the phase's "done when");
* a stale socket and lock do not stop a new service, a second service
  leaves the first alone, and the service ends on its own when nothing
  is left, on shutdown, and on an emergency stop;
* a process an agent session started may watch, never act.
"""

from __future__ import annotations

import json
import os
import signal
import subprocess
import sys
import time

import pytest

from conftest import child_env
from delfin.agent.service import protocol
from delfin.agent.service.client import ServiceClient, ServiceError
from tui_service import REPO, ServiceProcess, of, turn


def _env(service: ServiceProcess, **extra: str) -> dict:
    """A child of this test, on the service's private state: child_env,
    plus where the service's socket is."""
    env = child_env(service.tmp)
    env.update(DELFIN_RUNTIME_DIR=service.env["DELFIN_RUNTIME_DIR"],
               PYTHONPATH=service.env["PYTHONPATH"], **extra)
    return env


@pytest.fixture
def service(tmp_path):
    svc = ServiceProcess(tmp_path)
    try:
        yield svc
    finally:
        svc.stop()


def test_two_windows_on_one_session_see_the_same_turn(service):
    a, b = service.client(), service.client()
    sid = service.session(a)
    turn(a, sid, "hello")
    assert [s["sessionId"] for s in b.request("sessions.list")["sessions"]] == [sid]
    attached = b.request("session.attach", {"sessionId": sid})
    replayed = [e["type"] for e in attached["events"]]
    assert replayed[0] == "session.info" and "turn.ended" in replayed
    mark_a, mark_b = len(a.events), len(b.events)
    turn(a, sid, "hello again")
    b.wait_for(lambda e: e.get("type") == "turn.ended" and e.get("session") == sid,
               after=mark_b)
    seen_a = [e.get("seq") for e in of(a, sid, mark_a) if "seq" in e]
    seen_b = [e.get("seq") for e in of(b, sid, mark_b) if "seq" in e]
    assert seen_b == seen_a and len(set(seen_b)) == len(seen_b)


def test_the_first_answer_wins_and_every_window_hears_it(service):
    a, b = service.client(), service.client()
    sid = service.session(a, permissionMode="default")
    b.request("session.attach", {"sessionId": sid})
    mark_a, mark_b = len(a.events), len(b.events)
    a.request("session.prompt", {"sessionId": sid, "text": "edit it"})
    asked, _ = b.wait_for("permission.asked", after=mark_b)
    a.wait_for("permission.asked", after=mark_a)
    a.request("permission.reply", {"sessionId": sid,
                                   "requestId": asked["requestId"], "choice": "y"})
    with pytest.raises(ServiceError) as late:
        b.request("permission.reply", {"sessionId": sid,
                                       "requestId": asked["requestId"], "choice": "n"})
    assert late.value.code == protocol.UNKNOWN_REQUEST
    for c, mark in ((a, mark_a), (b, mark_b)):
        resolved, _ = c.wait_for("permission.resolved", after=mark)
        assert (resolved["decision"], resolved["by"]) == ("approved", "client")
    summary, _ = b.wait_for(lambda e: e.get("type") == "session.summary"
                            and e["summary"]["pending"] == 0, after=mark_b)
    assert summary["session"] == sid


def test_a_window_that_closes_leaves_its_session_working(service):
    a = service.client()
    sid = service.session(a, fakeDelay=0.02)
    mark = len(a.events)
    a.request("session.prompt", {"sessionId": sid, "text": "slow please"})
    a.wait_for("text.started", after=mark)
    a.close()
    time.sleep(1.0)
    b = service.client()
    attached = b.request("session.attach", {"sessionId": sid})
    if not any(e["type"] == "turn.ended" for e in attached["events"]):
        b.wait_for(lambda e: e.get("type") == "turn.ended", timeout=30)
    full = len("Working through a long calculation step by step. " * 4)
    events = attached["events"] + b.events
    text = "".join(e.get("text", "") for e in events if e.get("type") == "text.ended")
    assert len(text) >= full                          # it was not cut short


def test_three_sessions_survive_a_disconnect_and_a_reconnect(service):
    a = service.client()
    sids = [service.session(a, cwd=str(service.workspace(f"w{n}"))) for n in range(3)]
    for n, sid in enumerate(sids):
        turn(a, sid, f"first message to session {n}")
    a.close()
    b = service.client()
    listed = {s["sessionId"]: s for s in b.request("sessions.list")["sessions"]}
    assert set(listed) == set(sids)
    assert all(s["state"] == "idle" for s in listed.values())
    for n, sid in enumerate(sids):
        replay = b.request("session.attach", {"sessionId": sid})["events"]
        assert any(e.get("type") == "user.message"
                   and e["text"] == f"first message to session {n}" for e in replay)
        start, end = turn(b, sid, f"second message to session {n}")
        assert b.events[end].get("error") is None


def test_a_reattach_gets_what_happened_meanwhile_once(service):
    a, b = service.client(), service.client()
    sid = service.session(a)
    b.request("session.attach", {"sessionId": sid})
    b.request("session.detach", {"sessionId": sid})
    mark_b = len(b.events)
    turn(a, sid, "while you were away")
    assert not [e for e in of(b, sid, mark_b) if e.get("type") == "turn.ended"]
    # Attach while a turn is running: the replay and the live events meet,
    # and no durable event arrives twice.
    mark_a = len(a.events)
    a.request("session.prompt", {"sessionId": sid, "text": "slow please"})
    a.wait_for("text.started", after=mark_a)
    attached = b.request("session.attach", {"sessionId": sid, "turns": 2})
    b.wait_for(lambda e: e.get("type") == "turn.ended" and e.get("session") == sid,
               after=mark_b)
    seqs = [e["seq"] for e in attached["events"] if e.get("seq")]
    seqs += [e["seq"] for e in of(b, sid, mark_b) if "seq" in e]
    assert len(seqs) == len(set(seqs))
    assert any(e.get("type") == "user.message" and e["text"] == "while you were away"
               for e in attached["events"])


def test_a_session_that_ends_is_removed_for_every_window(service):
    a, b = service.client(), service.client()
    sid = service.session(a)
    mark = len(b.events)
    assert a.request("session.shutdown", {"sessionId": sid}) == {"closing": True}
    removed, _ = b.wait_for(lambda e: e.get("type") == "session.removed"
                            and e.get("session") == sid, after=mark)
    assert removed["reason"] == "shutdown"
    assert b.request("sessions.list")["sessions"] == []
    with pytest.raises(ServiceError) as gone:
        b.request("session.prompt", {"sessionId": sid, "text": "hello?"})
    assert gone.value.code == -32004


def test_a_session_that_cannot_start_says_why(service):
    a = service.client()
    with pytest.raises(ServiceError) as refused:
        a.request("session.create", {"cwd": "/", "fake": False}, timeout=60)
    assert refused.value.message
    assert a.request("sessions.list")["sessions"] == []


def test_a_process_an_agent_started_may_watch_but_not_act(service):
    a = service.client()
    sid = service.session(a)
    probe = r"""
import json, sys
sys.path.insert(0, sys.argv[2])
from delfin.agent.service.client import ServiceClient, ServiceError
c = ServiceClient(sys.argv[1])
out = {"mayAct": c.hello["mayAct"], "listed": len(c.request("sessions.list")["sessions"])}
for method, params in (("session.prompt", {"sessionId": sys.argv[3], "text": "rm -rf"}),
                       ("permission.reply", {"sessionId": sys.argv[3], "requestId": "x",
                                             "choice": "y"}),
                       ("session.create", {"cwd": "/tmp", "fake": True})):
    try:
        c.request(method, params)
        out[method] = "done"
    except ServiceError as exc:
        out[method] = exc.code
print(json.dumps(out))
"""
    run = subprocess.run([sys.executable, "-c", probe, str(service.socket), str(REPO), sid],
                         env=_env(service, DELFIN_AGENT_SESSION="12345"),
                         capture_output=True, text=True, timeout=60)
    got = json.loads(run.stdout.strip().splitlines()[-1])
    assert got["mayAct"] is False and got["listed"] == 1
    for method in ("session.prompt", "permission.reply", "session.create"):
        assert got[method] == protocol.NOT_ALLOWED, got
    assert a.hello["mayAct"] is True


def test_nothing_is_accepted_before_hello(service):
    raw = ServiceClient(service.socket, hello=False)
    try:
        with pytest.raises(ServiceError) as early:
            raw.request("sessions.list")
        assert early.value.code == protocol.NOT_ALLOWED
        with pytest.raises(ServiceError) as old:
            raw.request("hello", {"protocolVersion": 0})
        assert old.value.code == protocol.VERSION_MISMATCH
    finally:
        raw.close()


# -- the service's own life --------------------------------------------------------

def test_a_second_service_leaves_the_first_alone(service):
    run = subprocess.run([sys.executable, "-m", "delfin.agent.service.supervisor"],
                         env=_env(service), capture_output=True, text=True, timeout=60)
    assert run.returncode == 0 and "already runs" in run.stderr
    service.client().request("sessions.list")         # the first one still answers


def test_a_stale_socket_and_lock_do_not_stop_a_new_service(tmp_path):
    svc = ServiceProcess(tmp_path)
    try:
        svc.proc.send_signal(signal.SIGKILL)
        svc.proc.wait(timeout=10)
        assert svc.socket.exists()                    # left behind
        svc.start()
        svc.client().request("sessions.list")
    finally:
        svc.stop()


def test_the_service_ends_its_sessions_when_it_stops(tmp_path):
    svc = ServiceProcess(tmp_path)
    try:
        a = svc.client()
        sid = svc.session(a)
        worker = next(s["pid"] for s in a.request("sessions.list")["sessions"])
        assert sid
        a.request("service.shutdown")
        assert svc.proc.wait(timeout=60) == 0
        deadline = time.monotonic() + 10
        while time.monotonic() < deadline and _alive(worker):
            time.sleep(0.1)
        assert not _alive(worker)
        assert not svc.socket.exists()
    finally:
        svc.stop()


def test_the_service_ends_on_its_own_when_nothing_is_left(tmp_path):
    svc = ServiceProcess(tmp_path, idle_exit=2.0)
    try:
        a = svc.client()
        sid = svc.session(a)
        a.request("session.shutdown", {"sessionId": sid})
        a.wait_for("session.removed")
        a.close()
        assert svc.proc.wait(timeout=30) == 0
    finally:
        svc.stop()


def test_an_emergency_stop_ends_the_service_and_its_sessions(tmp_path):
    # The stop is a file every agent process polls (stop_all): written here
    # in the service's private home, it reaches this service and nothing
    # else -- `delfin-agent stop-all` itself would also end the user's
    # own sessions on this machine.
    svc = ServiceProcess(tmp_path)
    try:
        a = svc.client()
        svc.session(a)
        worker = next(s["pid"] for s in a.request("sessions.list")["sessions"])
        stop = svc.home / ".delfin" / "stop_all.json"
        stop.parent.mkdir(parents=True, exist_ok=True)
        stop.write_text(json.dumps({"at": time.time() + 1, "host": "test",
                                    "pid": os.getpid(), "reason": "test", "count": 1}))
        assert svc.proc.wait(timeout=60) == 0
        deadline = time.monotonic() + 20
        while time.monotonic() < deadline and _alive(worker):
            time.sleep(0.2)
        assert not _alive(worker)
    finally:
        svc.stop()


def _alive(pid: int) -> bool:
    try:
        os.kill(pid, 0)
    except ProcessLookupError:
        return False
    except PermissionError:
        return True
    # A zombie is gone for every purpose here.
    try:
        with open(f"/proc/{pid}/stat") as fh:
            return fh.read().rsplit(")", 1)[1].split()[0] != "Z"
    except OSError:
        return False


# -- the command line --------------------------------------------------------------

def test_service_status_and_stop_from_the_command_line(service):
    a = service.client()
    service.session(a, name="parser-fix")
    cli = [sys.executable, "-m", "delfin.agent.cli", "service"]
    status = subprocess.run(cli + ["status"], env=_env(service), capture_output=True,
                            text=True, timeout=60)
    assert status.returncode == 0, status.stderr
    assert "1 session" in status.stdout and "parser-fix" in status.stdout
    stopped = subprocess.run(cli + ["stop"], env=_env(service), capture_output=True,
                             text=True, timeout=120)
    assert stopped.returncode == 0 and "stopped" in stopped.stdout
    assert service.proc.wait(timeout=30) == 0
    after = subprocess.run(cli + ["status"], env=_env(service), capture_output=True,
                           text=True, timeout=60)
    assert "no agent service runs" in after.stdout


# -- other sessions; the guard; the launcher ---------------------------------------

def test_other_open_sessions_are_listed_and_can_be_messaged(service):
    # A terminal session elsewhere: a process that announces itself the way
    # the REPL does, on the service's private state.
    other = subprocess.Popen(
        [sys.executable, "-c",
         "import time; from delfin.agent import state_paths as sp; "
         "sp.scratch_state_from_environment(); "
         "from delfin.agent import session_presence as p; "
         "p.announce('dashboard-7', session_id='abc123', title='redox table', "
         "workspace='/tmp'); time.sleep(60)"],
        env=_env(service))
    try:
        a = service.client()
        mine = service.session(a)
        deadline = time.monotonic() + 15
        elsewhere: list = []
        while time.monotonic() < deadline:
            elsewhere = a.request("sessions.list")["elsewhere"]
            if any(e["key"] == "dashboard-7" for e in elsewhere):
                break
            time.sleep(0.2)
        entry = next(e for e in elsewhere if e["key"] == "dashboard-7")
        assert entry["title"] == "redox table"
        assert all(e.get("sessionId") != mine for e in elsewhere)
        assert a.request("session.message", {"key": "dashboard-7",
                                             "text": "please re-check row 4"}) == {"sent": True}
        inbox = subprocess.run(
            [sys.executable, "-c",
             "import json; from delfin.agent import state_paths as sp; "
             "sp.scratch_state_from_environment(); "
             "from delfin.agent import session_messages as m; "
             "print(json.dumps(m.take('dashboard-7')))"],
            env=_env(service), capture_output=True, text=True, timeout=60)
        taken = json.loads(inbox.stdout.strip().splitlines()[-1])
        assert [t["text"] for t in taken] == ["please re-check row 4"]
    finally:
        other.kill()
        other.wait(timeout=10)


def test_the_gate_refuses_a_command_that_names_the_services_socket():
    from delfin.agent import approval_guard
    from delfin.agent.service import rundir

    rooms = approval_guard.answer_rooms()
    assert rundir.runtime_dir().resolve() in rooms
    assert rundir.run_file().parent.resolve() in rooms
    assert approval_guard.is_answer_room_path(rundir.socket_path())
    for command in ("socat - UNIX-CONNECT:/tmp/delfin-1000/agent.sock",
                    "cat ~/.delfin/agent_service/run/node1.json",
                    'python -c "os.path.join(\'/tmp/delfin-1000\', \'agent.sock\')"'):
        assert approval_guard.answer_command(command), command
    for command in ("ls ~/.delfin/agent_service/sessions", "grep agent.sock README.md"):
        assert approval_guard.answer_command(command) is None, command


def test_the_launcher_starts_the_service_once_and_finds_it_after(tmp_path):
    svc = ServiceProcess(tmp_path, start=False)
    probe = r"""
import json, sys
from delfin.agent import state_paths as sp
sp.scratch_state_from_environment()
from delfin.agent.service import launcher
first = launcher.ensure_service(idle_exit_s=60)
second = launcher.ensure_service()
print(json.dumps({"first": first["pid"], "second": second["pid"]}))
"""
    try:
        run = subprocess.run([sys.executable, "-c", probe], env=_env(svc),
                             capture_output=True, text=True, timeout=120)
        assert run.returncode == 0, run.stderr[-2000:]
        got = json.loads(run.stdout.strip().splitlines()[-1])
        assert got["first"] == got["second"]
        with ServiceClient(svc.socket) as c:
            assert c.hello["pid"] == got["first"]
            c.request("service.shutdown")
    finally:
        svc.stop()


def test_tui_starts_the_service_and_says_when_the_window_is_not_built(tmp_path, monkeypatch):
    # Until the TypeScript window exists in a checkout, `delfin-agent tui`
    # still brings the service up, and says what is missing.
    from delfin.agent.service import launcher

    svc = ServiceProcess(tmp_path, start=False)
    try:
        if launcher.tui_command() is not None:
            pytest.skip("the TUI is built in this checkout; covered by its own tests")
        run = subprocess.run([sys.executable, "-m", "delfin.agent.cli", "tui"],
                             env=_env(svc), capture_output=True, text=True,
                             timeout=120, stdin=subprocess.DEVNULL)
        assert run.returncode == 2, run.stderr[-2000:]
        assert "not built" in run.stderr
        with ServiceClient(svc.socket) as c:
            c.request("service.shutdown")
    finally:
        svc.stop()
