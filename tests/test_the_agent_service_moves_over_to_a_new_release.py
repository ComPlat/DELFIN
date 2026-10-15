"""An update: a window of a newer DELFIN finds the agent service of the
release before.

* Idle (no turn running, no question open, and each session, asked,
  runs nothing a restart would cut), it is moved over: it stops, its
  sessions are saved and come back asleep in the new service.
* From before dry runs it cannot be asked (it would end the session it
  is asked about): moved over only when no session of it is live.
* Busy, it is left as it is -- its turns are not cut -- and the window
  says how to move it over later.
* Of another protocol version and idle (by its roster), it is restarted
  the same way.
"""

from __future__ import annotations

import os
import signal
import subprocess
import sys
import time

from conftest import child_env
from tui_service import ServiceProcess, turn

LAUNCH = """
import json
from delfin.agent.service import launcher
try:
    print(json.dumps({"pid": launcher.ensure_service()["pid"]}))
except launcher.ServiceDidNotStart as e:
    print(json.dumps({"error": str(e)}))
"""


def _env(svc: ServiceProcess) -> dict:
    """The test's own home and scratch state (child_env), with the
    service's private runtime directory and this checkout."""
    env = child_env(svc.tmp)
    env.update(DELFIN_RUNTIME_DIR=svc.env["DELFIN_RUNTIME_DIR"], PYTHONPATH=svc.env["PYTHONPATH"])
    return env


#: A service from before dry runs: serviceVersion 1, and a shutdown's
#: dryRun ignored -- the session ends, as it would there.
BEFORE_DRY_RUNS = """
supervisor.SERVICE_VERSION = 1
_forward = supervisor.Supervisor._forward
async def _ignoring_dry_run(self, client, method, params):
    if method == "session.shutdown":
        params = {k: v for k, v in params.items() if k != "dryRun"}
    return await _forward(self, client, method, params)
supervisor.Supervisor._forward = _ignoring_dry_run
"""


def _old_release(svc: ServiceProcess, *, protocol_version: int | None = None,
                 before_dry_runs: bool = False) -> subprocess.Popen:
    """The service of DELFIN 0.0.1 (or of another protocol version)."""
    patch = f"protocol.PROTOCOL_VERSION = {protocol_version}" if protocol_version is not None else ""
    code = "\n".join([
        "import sys",
        "import delfin",
        "delfin.__version__ = '0.0.1'",
        "from delfin.agent.service import protocol",
        patch,
        "import delfin.agent.service.supervisor as supervisor",
        BEFORE_DRY_RUNS if before_dry_runs else "",
        "sys.exit(supervisor.main(['--idle-exit', '120']))",
    ])
    proc = subprocess.Popen([sys.executable, "-c", code], env=_env(svc), start_new_session=True,
                            stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    svc.proc = proc
    deadline = time.monotonic() + 30
    while not svc.socket.exists() and time.monotonic() < deadline:
        time.sleep(0.1)
    return proc


def _launch(svc: ServiceProcess) -> tuple[dict, str]:
    import json

    out = subprocess.run([sys.executable, "-c", LAUNCH], env=_env(svc), capture_output=True,
                         text=True, timeout=180)
    return json.loads(out.stdout.strip().splitlines()[-1]), out.stderr


def _gone(proc: subprocess.Popen, timeout: float = 60) -> bool:
    try:
        proc.wait(timeout=timeout)
        return True
    except subprocess.TimeoutExpired:
        return False


def test_an_idle_service_of_the_release_before_is_moved_over(tmp_path):
    svc = ServiceProcess(tmp_path, start=False)
    old = _old_release(svc)
    try:
        a = svc.client()
        assert a.hello["release"] == "0.0.1"
        sid = svc.session(a)
        turn(a, sid, "made before the update")
        a.request("session.detach", {"sessionId": sid})
        got, said = _launch(svc)
        assert "moving the agent service from DELFIN 0.0.1" in said, said
        assert _gone(old)
        new = got["pid"]
        assert new != old.pid
        try:
            b = svc.client()
            assert b.hello["release"] != "0.0.1"
            row = next(s for s in b.hello["sessions"] if s["sessionId"] == sid)
            assert row["state"] == "sleeping"
            attached = b.request("session.attach", {"sessionId": sid, "turns": 1}, timeout=120)
            assert [e["text"] for e in attached["events"] if e["type"] == "user.message"] == ["made before the update"]
        finally:
            os.kill(new, signal.SIGTERM)
    finally:
        svc.stop()


def test_a_busy_service_of_the_release_before_is_left_as_it_is(tmp_path):
    svc = ServiceProcess(tmp_path, start=False)
    old = _old_release(svc)
    try:
        a = svc.client()
        sid = svc.session(a, permissionMode="default")
        mark = len(a.events)
        a.request("session.prompt", {"sessionId": sid, "text": "please edit the parser"})
        a.wait_for(lambda e: e.get("type") == "permission.asked", after=mark, timeout=30)
        got, said = _launch(svc)
        assert got == {"pid": old.pid}
        assert "`delfin-agent service restart` moves it over" in said, said
        assert old.poll() is None
    finally:
        svc.stop()


def test_a_service_from_before_dry_runs_is_never_asked_and_moves_over_once_nothing_is_live(tmp_path):
    svc = ServiceProcess(tmp_path, start=False)
    old = _old_release(svc, before_dry_runs=True)
    try:
        a = svc.client()
        assert a.hello["serviceVersion"] == 1
        sid = svc.session(a)
        turn(a, sid, "idle, and in view")
        got, said = _launch(svc)
        assert got == {"pid": old.pid}, (got, said)
        assert "cannot say what a restart would cut" in said, said
        # Never asked: asking would have ended it.
        assert next(s for s in a.request("sessions.list")["sessions"] if s["sessionId"] == sid)["state"] == "idle"
        a.request("session.shutdown", {"sessionId": sid})
        deadline = time.monotonic() + 30
        while a.request("sessions.list")["sessions"] and time.monotonic() < deadline:
            time.sleep(0.1)
        got, said = _launch(svc)
        assert "moving the agent service from DELFIN 0.0.1" in said, said
        assert _gone(old) and got["pid"] != old.pid
        os.kill(got["pid"], signal.SIGTERM)
    finally:
        svc.stop()


def test_a_service_of_a_newer_protocol_is_never_stopped_by_an_older_window(tmp_path):
    svc = ServiceProcess(tmp_path, start=False)
    newer = _old_release(svc, protocol_version=99)
    try:
        roster = svc.scratch / ".delfin" / "agent_service" / "run"
        deadline = time.monotonic() + 30
        while not list(roster.glob("*.roster")) and time.monotonic() < deadline:
            time.sleep(0.1)
        got, said = _launch(svc)
        assert "speaks another protocol version" in got.get("error", ""), got
        assert newer.poll() is None
    finally:
        svc.stop()


def test_an_idle_service_of_another_protocol_is_restarted(tmp_path):
    svc = ServiceProcess(tmp_path, start=False)
    old = _old_release(svc, protocol_version=0)
    try:
        roster = svc.scratch / ".delfin" / "agent_service" / "run"
        deadline = time.monotonic() + 30
        while not list(roster.glob("*.roster")) and time.monotonic() < deadline:
            time.sleep(0.1)
        got, said = _launch(svc)
        assert "pid" in got, got
        assert "of another DELFIN version and idle: restarting it" in said, said
        assert _gone(old) and got["pid"] != old.pid
        os.kill(got["pid"], signal.SIGTERM)
    finally:
        svc.stop()


def test_an_update_asks_each_session_what_a_restart_would_cut():
    # "idle" in a summary does not see a background shell or a subagent:
    # the sessions themselves are asked (a dry-run shutdown).
    from delfin.agent.service import launcher

    class Service:
        def __init__(self, answers):
            self.answers = answers

        def request(self, method, params=None, timeout=None):
            assert method == "session.shutdown" and params["dryRun"] is True
            answer = self.answers[params["sessionId"]]
            if isinstance(answer, Exception):
                raise answer
            return answer

    idle = {"sessionId": "a", "state": "idle", "title": "scan"}
    asleep = {"sessionId": "z", "state": "sleeping"}
    assert launcher._busy_sessions(Service({"a": {"closing": False, "busy": ""}}), [idle, asleep]) == []
    assert launcher._busy_sessions(Service({"a": {"closing": False, "busy": "1 background shell running"}}),
                                   [idle]) == ["scan: 1 background shell running"]
    said = launcher._busy_sessions(Service({"a": RuntimeError("no such method")}), [idle])
    assert said and "cannot say" in said[0]                 # an older service: counted busy
