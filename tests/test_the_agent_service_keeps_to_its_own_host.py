"""Two login nodes, one shared home: each host has its own service, its own
run file and its own roster, and keeps to them.

* A service on another host is named, never stopped or taken over from
  here: its sessions -- running or asleep -- stay its own, also when this
  host's service stops or starts again.
* A window here is told where the other host's sessions are.

The other host is a second service whose host name is patched (one VM:
"loginb"); both share the scratch state (the shared home) and each has
its own runtime directory (/tmp is per node).
"""

from __future__ import annotations

import json
import subprocess
import sys
import tempfile
import textwrap
import time
from pathlib import Path

from conftest import child_env
from delfin.agent.service import rundir
from tui_service import ServiceProcess, turn

OTHER = "loginb"


def _other_host_service(svc: ServiceProcess) -> tuple[subprocess.Popen, Path]:
    """A service of the same user on host OTHER: the same shared home,
    its own /tmp."""
    runtime = Path(tempfile.mkdtemp(prefix="dsvb", dir="/tmp"))
    env = dict(svc.env, DELFIN_RUNTIME_DIR=str(runtime))
    code = textwrap.dedent(f"""
        import socket, sys
        socket.gethostname = lambda: {OTHER!r}
        import delfin.agent.service.supervisor as supervisor
        sys.exit(supervisor.main(["--idle-exit", "120"]))
    """)
    proc = subprocess.Popen([sys.executable, "-c", code], env=env, start_new_session=True,
                            stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    sock = runtime / f"delfin-{__import__('os').getuid()}" / "agent.sock"
    deadline = time.monotonic() + 30
    while not sock.exists() and time.monotonic() < deadline:
        time.sleep(0.1)
    return proc, sock


def _cli(svc: ServiceProcess, *args: str) -> subprocess.CompletedProcess:
    env = child_env(svc.tmp)
    env.update(DELFIN_RUNTIME_DIR=svc.env["DELFIN_RUNTIME_DIR"], PYTHONPATH=svc.env["PYTHONPATH"])
    return subprocess.run([sys.executable, "-m", "delfin.agent.cli", *args], env=env,
                          capture_output=True, text=True, timeout=120)


def test_each_host_keeps_its_own_sessions(tmp_path):
    from delfin.agent.service.client import ServiceClient

    here = ServiceProcess(tmp_path)
    other, other_sock = _other_host_service(here)
    try:
        run = here.scratch / ".delfin" / "agent_service" / "run"
        assert (run / f"{OTHER}.json").exists()
        with ServiceClient(other_sock) as b:
            theirs = b.request("session.create", {"cwd": str(here.workspace()), "fake": True,
                                                  "fakeDelay": 0.003}, timeout=120)["sessionId"]
            turn(b, theirs, "on the other node")
        a = here.client()
        mine = here.session(a)
        # Named, never stopped from here.
        status = _cli(here, "service", "status")
        assert f"one runs on {OTHER}" not in status.stdout       # this host's service answers here
        assert all(s["sessionId"] != theirs for s in a.request("sessions.list")["sessions"])
        # This host's service stops and starts again: it brings back its
        # own session, not the other host's.
        a.request("service.shutdown")
        assert here.proc.wait(timeout=60) == 0
        stopped = _cli(here, "service", "status")
        assert f"one runs on {OTHER}" in stopped.stdout, stopped.stdout
        here.start()
        ids = {s["sessionId"]: s["state"] for s in here.client().hello["sessions"]}
        assert ids == {mine: "sleeping"}
        assert other.poll() is None
        with ServiceClient(other_sock) as b:
            row = next(s for s in b.request("sessions.list")["sessions"] if s["sessionId"] == theirs)
            assert row["state"] == "idle"
        roster = json.loads((run / f"{OTHER}.roster").read_text())
        assert [s["sessionId"] for s in roster["sessions"]] == [theirs]
        # A window here hears where the other host's sessions are.
        notes = subprocess.run(
            [sys.executable, "-c", "import json; from delfin.agent.service import launcher; "
                                   "print(json.dumps(launcher.other_host_notes()))"],
            env=dict(here.env), capture_output=True, text=True, timeout=60)
        said = json.loads(notes.stdout)
        assert said and OTHER in said[0] and "1 session" in said[0], said
    finally:
        try:
            with ServiceClient(other_sock) as b:
                b.request("service.shutdown")
        except Exception:
            other.terminate()
        other.wait(timeout=60)
        here.stop()


def test_a_run_file_names_its_host(tmp_path, monkeypatch):
    monkeypatch.setenv("DELFIN_RUNTIME_DIR", str(tmp_path / "rt"))
    monkeypatch.setenv("DELFIN_AGENT_SERVICE_DIR", str(tmp_path / "service"))
    path = rundir.write_run_file({"pid": 1, "host": rundir.hostname()})
    assert path.name == f"{rundir.hostname()}.json"
    assert rundir.roster_file().name == f"{rundir.hostname()}.roster"
    assert rundir.roster_file(OTHER).name == f"{OTHER}.roster"
