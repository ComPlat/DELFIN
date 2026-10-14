"""Sessions without a worker that wake when used: dormant sessions.

* A session no window shows, with nothing to do, sleeps after a while: its
  worker ends (the session is saved), its summary says "sleeping", and it
  wakes when used -- sent a prompt, attached to. A peek reads its log and
  does not wake it. A session a window shows never sleeps.
* What a service ran when it stopped -- or died -- comes back in the next
  one, to be woken. Stopping a dormant session forgets it.
* A roster another live service writes is not taken over.
* A service that stops on purpose says so to its windows first (they leave
  instead of starting another); `service restart` stops it, starts the
  next one, and that one offers the sessions again. The launcher, run as a
  module, starts a service when none runs (what a window does after a
  crash).
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
from delfin.agent.service import protocol, rundir, supervisor
from delfin.agent.service.client import ServiceError
from tui_service import ServiceProcess, of, turn


def _sleepy(tmp_path, after: float = 1.0) -> ServiceProcess:
    svc = ServiceProcess(tmp_path, start=False)
    svc.env["DELFIN_AGENT_SLEEP_AFTER_S"] = str(after)
    svc.start()
    return svc


def _listed(c, sid: str) -> dict | None:
    return next((s for s in c.request("sessions.list")["sessions"] if s["sessionId"] == sid), None)


def _wait_state(c, sid: str, state: str, timeout: float = 30) -> dict:
    deadline = time.monotonic() + timeout
    while time.monotonic() < deadline:
        row = _listed(c, sid)
        if row and row["state"] == state:
            return row
        time.sleep(0.2)
    raise AssertionError(f"{sid} never became {state}: {_listed(c, sid)}")


def _gone(pid: int, timeout: float = 30) -> bool:
    deadline = time.monotonic() + timeout
    while time.monotonic() < deadline:
        try:
            os.kill(pid, 0)
        except ProcessLookupError:
            return True
        try:
            with open(f"/proc/{pid}/stat") as fh:
                if fh.read().rsplit(")", 1)[1].split()[0] == "Z":
                    return True
        except OSError:
            return True
        time.sleep(0.1)
    return False


def test_an_unused_session_sleeps_and_wakes_when_used(tmp_path):
    svc = _sleepy(tmp_path)
    try:
        a = svc.client()
        sid = svc.session(a)
        turn(a, sid, "hello")
        worker = _listed(a, sid)["pid"]
        a.request("session.detach", {"sessionId": sid})     # no window shows it now
        mark = len(a.events)
        row = _wait_state(a, sid, "sleeping")
        assert row["pid"] == 0 and row["pending"] == 0
        assert _gone(worker), "its worker should have ended"
        a.wait_for(lambda e: e.get("type") == "session.summary" and e.get("session") == sid
                   and e["summary"]["state"] == "sleeping", after=mark)
        # A peek reads the log; it does not wake it.
        replay = a.request("session.replay", {"sessionId": sid, "turns": 1})["events"]
        assert [e["text"] for e in replay if e.get("type") == "user.message"] == ["hello"]
        assert a.request("session.status", {"sessionId": sid})["state"] == "sleeping"
        assert _listed(a, sid)["state"] == "sleeping"
        env = child_env(svc.tmp)
        env.update(DELFIN_RUNTIME_DIR=svc.env["DELFIN_RUNTIME_DIR"], PYTHONPATH=svc.env["PYTHONPATH"])
        status = subprocess.run([sys.executable, "-m", "delfin.agent.cli", "service", "status"],
                                env=env, capture_output=True, text=True, timeout=60)
        assert "asleep" in status.stdout and sid[:8] in status.stdout, status.stdout + status.stderr
        # A prompt wakes it, under its own id, and is answered.
        b = svc.client()
        mark = len(b.events)
        b.request("session.attach", {"sessionId": sid, "turns": 0}, timeout=120)
        b.request("session.prompt", {"sessionId": sid, "text": "again"})
        b.wait_for(lambda e: e.get("type") == "turn.ended" and e.get("session") == sid,
                   after=mark, timeout=60)
        assert _listed(b, sid)["state"] == "idle"
        said = [e["text"] for e in b.request("session.replay", {"sessionId": sid})["events"]
                if e.get("type") == "user.message"]
        assert said == ["hello", "again"]
    finally:
        svc.stop()


def test_a_session_a_window_shows_never_sleeps(tmp_path):
    svc = _sleepy(tmp_path)
    try:
        a = svc.client()
        sid = svc.session(a)
        turn(a, sid, "hello")
        time.sleep(3.5)
        assert _listed(a, sid)["state"] == "idle"
    finally:
        svc.stop()


def test_a_prompt_wakes_a_sleeping_session(tmp_path):
    svc = _sleepy(tmp_path)
    try:
        a = svc.client()
        sid = svc.session(a, attach=False)
        _wait_state(a, sid, "sleeping")
        answer = a.request("session.prompt", {"sessionId": sid, "text": "wake up"}, timeout=120)
        assert answer["accepted"] is True
        a.request("session.attach", {"sessionId": sid, "turns": 0})
        deadline = time.monotonic() + 60
        while time.monotonic() < deadline and _listed(a, sid)["state"] != "idle":
            time.sleep(0.2)
        said = [e["text"] for e in a.request("session.replay", {"sessionId": sid})["events"]
                if e.get("type") == "user.message"]
        assert said == ["wake up"]
    finally:
        svc.stop()


def test_what_a_stopped_service_ran_comes_back_to_be_woken(tmp_path):
    svc = ServiceProcess(tmp_path)
    try:
        a = svc.client()
        sid = svc.session(a)
        turn(a, sid, "before the update")
        a.request("service.shutdown")
        assert svc.proc.wait(timeout=60) == 0
        svc.start()                                   # the next service
        b = svc.client()
        row = next(s for s in b.hello["sessions"] if s["sessionId"] == sid)
        assert row["state"] == "sleeping"
        attached = b.request("session.attach", {"sessionId": sid, "turns": 2}, timeout=120)
        said = [e["text"] for e in attached["events"] if e.get("type") == "user.message"]
        assert said == ["before the update"]
        assert attached["state"] == "idle" and _listed(b, sid)["state"] == "idle"
        assert not [e for e in attached["events"] if e.get("type") == "session.ended"]
    finally:
        svc.stop()


def test_what_a_service_ran_when_it_died_comes_back(tmp_path):
    svc = ServiceProcess(tmp_path)
    try:
        a = svc.client()
        sid = svc.session(a)
        turn(a, sid, "before the crash")
        worker = _listed(a, sid)["pid"]
        os.kill(svc.proc.pid, signal.SIGKILL)
        svc.proc.wait(timeout=30)
        assert _gone(worker), "a worker ends (and saves) when its service's pipe closes"
        svc.start()
        b = svc.client()
        row = next(s for s in b.hello["sessions"] if s["sessionId"] == sid)
        assert row["state"] == "sleeping"
        b.request("session.attach", {"sessionId": sid, "turns": 0}, timeout=120)
        mark = len(b.events)
        b.request("session.prompt", {"sessionId": sid, "text": "after the crash"})
        b.wait_for(lambda e: e.get("type") == "turn.ended" and e.get("session") == sid,
                   after=mark, timeout=60)
        said = [e["text"] for e in b.request("session.replay", {"sessionId": sid})["events"]
                if e.get("type") == "user.message"]
        assert said == ["before the crash", "after the crash"]
    finally:
        svc.stop()


def test_stopping_a_dormant_session_forgets_it(tmp_path):
    svc = _sleepy(tmp_path)
    try:
        a = svc.client()
        sid = svc.session(a, attach=False)
        _wait_state(a, sid, "sleeping")
        mark = len(a.events)
        assert a.request("session.shutdown", {"sessionId": sid}) == {"closing": True}
        a.wait_for(lambda e: e.get("type") == "session.removed" and e.get("session") == sid, after=mark)
        assert _listed(a, sid) is None
        with pytest.raises(ServiceError) as refused:
            a.request("session.attach", {"sessionId": sid})
        assert refused.value.code == protocol.UNKNOWN_SESSION
        a.request("service.shutdown")
        svc.proc.wait(timeout=60)
        svc.start()
        assert all(s["sessionId"] != sid for s in svc.client().hello["sessions"])
    finally:
        svc.stop()


def test_a_live_services_roster_is_not_taken_over(tmp_path, monkeypatch):
    monkeypatch.setenv("DELFIN_RUNTIME_DIR", str(tmp_path / "rt"))
    monkeypatch.setenv("DELFIN_AGENT_SERVICE_DIR", str(tmp_path / "service"))
    (tmp_path / "rt").mkdir()
    entry = {"sessionId": "s1", "workspace": str(tmp_path), "state": "idle", "pid": 4242,
             "settings": {"cwd": str(tmp_path), "fake": True}}
    path = rundir.roster_file()
    path.parent.mkdir(parents=True, exist_ok=True)

    def load(roster: dict) -> dict:
        path.write_text(json.dumps(roster))
        sup = supervisor.Supervisor()
        sup._load_dormant()
        return sup.dormant

    monkeypatch.setattr(supervisor, "_alive_service", lambda pid: pid == 777)
    assert load({"pid": 777, "sessions": [dict(entry)]}) == {}             # alive: its own
    assert load({"pid": 0, "socket": "/elsewhere.sock", "sessions": [dict(entry)]}) == {}
    found = load({"pid": 0, "sessions": [dict(entry)]})
    assert found["s1"]["why"] == "the agent service ended unexpectedly"
    assert found["s1"]["summary"]["state"] == "sleeping" and found["s1"]["oldPid"] == 4242
    old = {"summary": {"sessionId": "s2", "state": "sleeping"}, "settings": {"cwd": str(tmp_path)},
           "since": time.time() - supervisor.DORMANT_KEEP_S - 1, "why": "asleep"}
    assert "s2" not in load({"pid": 0, "dormant": [old]})                   # too old: forgotten


def _cli_env(svc: ServiceProcess) -> dict:
    env = child_env(svc.tmp)
    env.update(DELFIN_RUNTIME_DIR=svc.env["DELFIN_RUNTIME_DIR"], PYTHONPATH=svc.env["PYTHONPATH"])
    return env


@pytest.mark.parametrize("restart", [False, True])
def test_a_stop_on_purpose_is_said_to_the_windows_first(tmp_path, restart):
    svc = ServiceProcess(tmp_path)
    try:
        a, b = svc.client(), svc.client()
        mark = len(b.events)
        a.request("service.shutdown", {"restart": True} if restart else {})
        said, _ = b.wait_for(lambda e: e.get("type") == "service.stopping", after=mark, timeout=30)
        assert said["reason"] == ("restart" if restart else "stop") and said["restart"] is restart
        assert svc.proc.wait(timeout=60) == 0
    finally:
        svc.stop()


def test_service_restart_brings_the_sessions_back_in_the_next_service(tmp_path):
    svc = ServiceProcess(tmp_path)
    try:
        a = svc.client()
        sid = svc.session(a)
        turn(a, sid, "kept across the restart")
        old = svc.proc.pid
        out = subprocess.run([sys.executable, "-m", "delfin.agent.cli", "service", "restart"],
                             env=_cli_env(svc), capture_output=True, text=True, timeout=180)
        assert out.returncode == 0, out.stdout + out.stderr
        assert "runs again" in out.stdout
        svc.proc.wait(timeout=30)
        new = rundir_pid(svc)
        assert new and new != old
        try:
            c = svc.client()
            row = next(s for s in c.hello["sessions"] if s["sessionId"] == sid)
            assert row["state"] == "sleeping"
        finally:
            os.kill(new, signal.SIGTERM)
            _gone(new)
    finally:
        svc.stop()


def test_the_launcher_module_starts_a_service_when_none_runs(tmp_path):
    svc = ServiceProcess(tmp_path, start=False)
    out = subprocess.run([sys.executable, "-m", "delfin.agent.service.launcher"],
                         env=_cli_env(svc), capture_output=True, text=True, timeout=120)
    assert out.returncode == 0, out.stderr
    pid = int(out.stdout.strip())
    try:
        assert svc.client().hello["pid"] == pid
    finally:
        os.kill(pid, signal.SIGTERM)
        assert _gone(pid)
        svc.stop()


def rundir_pid(svc: ServiceProcess) -> int:
    for path in (svc.scratch / ".delfin" / "agent_service" / "run").glob("*.json"):
        try:
            return int(json.loads(path.read_text()).get("pid") or 0)
        except (OSError, ValueError):
            continue
    return 0


def test_a_busy_session_is_not_put_to_sleep(monkeypatch):
    # The worker's answer: what an idle shutdown would cut.
    from types import SimpleNamespace

    from delfin.agent import bash_jobs
    from delfin.agent.service import worker

    class Registry:
        def __init__(self, jobs):
            self.jobs = jobs

        def list_jobs(self, include_finished=True):
            return self.jobs

    idle = SimpleNamespace(busy=False, inbox=[], queued=[])
    no_question = SimpleNamespace(pending=lambda: [])
    monkeypatch.setattr(bash_jobs, "get_registry", lambda: Registry([]))
    assert worker.busy_reason(idle, no_question) == ""
    assert worker.busy_reason(SimpleNamespace(busy=True, inbox=[], queued=[]), no_question) == "a turn runs"
    assert worker.busy_reason(SimpleNamespace(busy=False, inbox=["x"], queued=[]), no_question)
    assert worker.busy_reason(idle, SimpleNamespace(pending=lambda: ["q"])) == "a question is open"
    monkeypatch.setattr(bash_jobs, "get_registry", lambda: Registry([object()]))
    assert worker.busy_reason(idle, no_question) == "1 background shell running"

    # The service: a worker that says it is busy is left awake, and asked
    # again only after another full while.
    import asyncio

    async def run():
        sup = supervisor.Supervisor(sleep_after_s=60)
        asked = []

        class Busy:
            sleeping = True
            sleep_not_before = 0.0
            done = asyncio.get_running_loop().create_future()

            async def request(self, method, params=None, timeout=None):
                asked.append((method, params))
                return {"closing": False, "busy": "1 background shell running"}

        session = Busy()
        before = time.time()
        await sup._put_to_sleep(session)
        return session, asked, before

    session, asked, before = asyncio.run(run())
    assert asked == [("session.shutdown", {"ifIdle": True})]
    assert session.sleeping is False and session.sleep_not_before >= before + 59


def test_an_emergency_stop_is_said_to_the_windows_and_keeps_the_sessions(tmp_path):
    import socket as _socket

    svc = ServiceProcess(tmp_path)
    try:
        a = svc.client()
        sid = svc.session(a)
        turn(a, sid, "before the emergency stop")
        mark = len(a.events)
        # The stop file of this test's own private home -- never the user's.
        stop = svc.home / ".delfin" / "stop_all.json"
        assert str(stop).startswith(str(svc.tmp))
        stop.parent.mkdir(parents=True, exist_ok=True)
        stop.write_text(json.dumps({"at": time.time() + 1, "host": _socket.gethostname().split(".")[0],
                                    "pid": os.getpid(), "reason": "test", "count": 1}))
        said, _ = a.wait_for(lambda e: e.get("type") == "service.stopping", after=mark, timeout=60)
        assert said["reason"] == "emergency stop" and said["restart"] is False
        svc.proc.wait(timeout=90)
        stop.unlink()
        svc.start()
        row = next(s for s in svc.client().hello["sessions"] if s["sessionId"] == sid)
        assert row["state"] == "sleeping"
    finally:
        svc.stop()


def test_an_attach_while_a_busy_session_refuses_to_sleep_is_answered():
    # Asked to sleep, the worker says it is busy: an attach waiting for the
    # session to go to sleep must not wait for an end that will not come.
    import asyncio
    from types import SimpleNamespace

    async def run():
        sup = supervisor.Supervisor()
        loop = asyncio.get_running_loop()
        session = SimpleNamespace(id="s1", sleeping=True, done=loop.create_future())
        sup.sessions["s1"] = session

        async def refuse():
            await asyncio.sleep(0.2)
            session.sleeping = False             # the worker said: busy

        asyncio.ensure_future(refuse())
        return await asyncio.wait_for(sup._awake(SimpleNamespace(may_act=True), "s1"), 5)

    assert asyncio.run(run()).id == "s1"
