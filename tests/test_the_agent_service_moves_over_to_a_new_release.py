"""An update: a window of a newer DELFIN finds the agent service of the
release before.

* Idle (no turn running, no question open), it is moved over: it stops,
  its sessions are saved and come back asleep in the new service.
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
import textwrap
import time

from tui_service import ServiceProcess, turn

LAUNCH = """
import json
from delfin.agent.service import launcher
try:
    print(json.dumps({"pid": launcher.ensure_service()["pid"]}))
except launcher.ServiceDidNotStart as e:
    print(json.dumps({"error": str(e)}))
"""


def _old_release(svc: ServiceProcess, *, protocol_version: int | None = None) -> subprocess.Popen:
    """The service of DELFIN 0.0.1 (or of another protocol version)."""
    patch = f"protocol.PROTOCOL_VERSION = {protocol_version}" if protocol_version else ""
    code = textwrap.dedent(f"""
        import sys
        import delfin
        delfin.__version__ = "0.0.1"
        from delfin.agent.service import protocol
        {patch}
        import delfin.agent.service.supervisor as supervisor
        sys.exit(supervisor.main(["--idle-exit", "120"]))
    """)
    proc = subprocess.Popen([sys.executable, "-c", code], env=svc.env, start_new_session=True,
                            stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    svc.proc = proc
    deadline = time.monotonic() + 30
    while not svc.socket.exists() and time.monotonic() < deadline:
        time.sleep(0.1)
    return proc


def _launch(svc: ServiceProcess) -> tuple[dict, str]:
    import json

    out = subprocess.run([sys.executable, "-c", LAUNCH], env=svc.env, capture_output=True,
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


def test_an_idle_service_of_another_protocol_is_restarted(tmp_path):
    svc = ServiceProcess(tmp_path, start=False)
    old = _old_release(svc, protocol_version=99)
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
