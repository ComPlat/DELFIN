"""The agent service follows DELFIN's cluster jobs and wakes the session a job ended for.

Input: a real agent service (fake model) whose cluster is reached through a
fake ``ssh`` that runs every remote command in the service's home, with
fake ``squeue``/``sacct`` saying what the test writes; a session that went
to sleep after submitting (its job record and watch as cluster_submit
leaves them); the job's folder on the "cluster" holding a DELFIN run's
files. Output: the files in the workspace, the session's state and what its
turn was given.

Semantics pinned here (the milestone's done-when):

* while the job runs nothing comes back and the session sleeps on;
* once Slurm says it ended, the service brings its DELFIN_Data.json and the
  rest back into the job's folder and wakes the sleeping session, whose
  turn is told -- as a system event, not the user -- what ended and where
  its results are;
* also after the PC was off: a job that ended while no service ran is
  brought back, and its session woken, by the next service.
"""

from __future__ import annotations

import json
import os
import sys
import time
from pathlib import Path

import pytest

from delfin.agent import job_monitor
from tui_service import ServiceProcess, turn

_FAKE_SSH = """\
#!{python}
import os, subprocess, sys
WITH_VALUE = set("bcDEeFIiJLlmOoPpQRSWw")
args, i = sys.argv[1:], 0
while i < len(args) and args[i].startswith("-"):
    if args[i] == "-O" and args[i + 1] == "check":
        sys.exit(0 if os.path.exists(os.environ["FAKE_MASTER"]) else 255)
    i += 2 if len(args[i]) == 2 and args[i][1] in WITH_VALUE else 1
sys.exit(subprocess.run(["sh", "-c", " ".join(args[i + 1:])], cwd=os.environ["HOME"]).returncode)
"""
_FAKE_SQUEUE = """#!/bin/sh
[ -f "$FAKE_STATE/squeue" ] && cat "$FAKE_STATE/squeue"
exit 0
"""
_FAKE_SACCT = """#!/bin/sh
case " $* " in
  *" -j "*) [ -f "$FAKE_STATE/sacct_j" ] && cat "$FAKE_STATE/sacct_j" ;;
  *) [ -f "$FAKE_STATE/sacct" ] && cat "$FAKE_STATE/sacct" ;;
esac
exit 0
"""

REMOTE = "delfin/jobs/20261021-160135-89f7-fe"
DONE = {"DELFIN_Data.json": '{"status": "completed", "completed": true}\n',
        "DELFIN.txt": "DELFIN report\n", ".exit_code_4242": "0\n",
        "delfin_4242.out": "DELFIN WRAPPER COMPLETED (exit 0)\n",
        "initial.out": "****ORCA TERMINATED NORMALLY****\n", "initial.gbw": "x" * 2048}


class Cluster:
    """The fake cluster the service reaches, and what its Slurm says."""

    def __init__(self, svc: ServiceProcess, *, agent: dict | None = None):
        self.svc = svc
        self.state = svc.tmp / "fake-slurm"
        self.state.mkdir()
        self.master = svc.tmp / "master"
        self.master.touch()
        bin_dir = svc.tmp / "fake-bin"
        bin_dir.mkdir()
        for name, text in (("ssh", _FAKE_SSH.format(python=sys.executable)),
                           ("squeue", _FAKE_SQUEUE), ("sacct", _FAKE_SACCT)):
            (bin_dir / name).write_text(text, encoding="utf-8")
            (bin_dir / name).chmod(0o755)
        home = svc.home
        ctun = home / ".config" / "cluster-tunnel"
        ctun.mkdir(parents=True, exist_ok=True)
        (ctun / "config.yaml").write_text(
            "clusters:\n  kcist:\n    host: login.example\n    user: me\n"
            "    strict_host_key_checking: 'yes'\n"
            f"    known_hosts_file: {ctun / 'known_hosts'}\n", encoding="utf-8")
        (home / ".delfin_settings.json").write_text(json.dumps({"remote": {"clusters": {
            "kcist": {"default_partition": "small"}}}, **({"agent": agent} if agent else {})}),
            encoding="utf-8")
        svc.env.update(PATH=f"{bin_dir}:/usr/bin:/bin", FAKE_MASTER=str(self.master),
                       FAKE_STATE=str(self.state), XDG_CONFIG_HOME=str(home / ".config"),
                       XDG_CACHE_HOME=str(home / ".cache"),
                       DELFIN_AGENT_CLUSTER_WATCH_S="0.5", DELFIN_AGENT_SLEEP_AFTER_S="1")
        svc.env.pop("CTUN_CONFIG", None)

    def running(self) -> None:
        (self.state / "squeue").write_text("active|4242|RUNNING|0:05|30:00|small|1|N/A|fe\n")
        (self.state / "sacct").write_text("")

    def ended(self) -> None:
        (self.state / "squeue").write_text("")
        end = time.strftime("%Y-%m-%dT%H:%M:%S")
        (self.state / "sacct").write_text(f"{end}|4242|COMPLETED|00:00:30|00:30:00|small|1|{end}|fe\n")
        for rel, text in DONE.items():
            path = self.svc.home / REMOTE / rel
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_text(text, encoding="utf-8")

    def submitted(self, sid: str, ws: Path) -> None:
        """What cluster_submit leaves: the record and the session's watch."""
        (ws / "fe").mkdir(exist_ok=True)
        records = self.svc.scratch / ".delfin" / "remote_jobs.jsonl"
        records.parent.mkdir(parents=True, exist_ok=True)
        with records.open("a", encoding="utf-8") as fh:
            fh.write(json.dumps({
                "cluster": "kcist", "job_id": "4242", "name": "fe", "kind": "delfin",
                "input": "", "remote_dir": f"~/{REMOTE}", "folder": str(ws / "fe"),
                "workspace": str(ws), "session_id": sid, "commit": "c0ffee",
                "submitted_at": time.time()}) + "\n")
        job_monitor.register_agent_job(
            ws, "kcist:4242", description="fe: DELFIN workflow on kcist",
            extra={"kind": "cluster", "cluster": "kcist", "remote_id": "4242",
                   "session_id": sid})


def _row(c, sid: str) -> dict | None:
    return next((s for s in c.request("sessions.list")["sessions"] if s["sessionId"] == sid),
                None)


def _until(cond, timeout: float, what: str):
    deadline = time.monotonic() + timeout
    while time.monotonic() < deadline:
        got = cond()
        if got:
            return got
        time.sleep(0.2)
    raise AssertionError(f"never: {what}")


def _told(c, sid: str) -> list[str]:
    """What the session's turns were given, the user's and the system's."""
    events = c.request("session.replay", {"sessionId": sid})["events"]
    return [str(e.get("text") or "") for e in events if e.get("type") == "user.message"]


@pytest.fixture
def service(tmp_path):
    svc = ServiceProcess(tmp_path, start=False)
    cluster = Cluster(svc)
    svc.start()
    try:
        yield svc, cluster
    finally:
        svc.stop()


def _asleep_after_a_submit(svc, cluster) -> tuple[object, str, Path]:
    c = svc.client()
    ws = svc.workspace()
    sid = svc.session(c, cwd=str(ws))
    turn(c, sid, "submit the batch")
    cluster.submitted(sid, ws)
    c.request("session.detach", {"sessionId": sid})
    _until(lambda: (_row(c, sid) or {}).get("state") == "sleeping", 60, "the session slept")
    return c, sid, ws


def test_an_ended_job_comes_back_and_wakes_its_sleeping_session(service):
    svc, cluster = service
    cluster.running()
    c, sid, ws = _asleep_after_a_submit(svc, cluster)
    time.sleep(2.5)                         # a few ticks of a running job
    assert not (ws / "fe" / "kcist-4242").exists()
    assert (_row(c, sid) or {}).get("state") == "sleeping"
    cluster.ended()
    data = ws / "fe" / "kcist-4242" / "DELFIN_Data.json"
    _until(data.is_file, 60, "DELFIN_Data.json came back")
    assert not (ws / "fe" / "kcist-4242" / "initial.gbw").exists()
    told = _until(lambda: [t for t in _told(c, sid) if "kcist:4242" in t], 90,
                  "the woken session was told")
    assert "[watch — a system event, not the user]" in told[0]
    assert "results in fe/kcist-4242" in told[0]


def test_a_job_that_ended_while_no_service_ran_comes_back_with_the_next(service):
    svc, cluster = service
    cluster.running()
    c, sid, ws = _asleep_after_a_submit(svc, cluster)
    c.request("service.shutdown")           # the PC goes off
    assert svc.proc.wait(timeout=60) == 0
    cluster.ended()
    svc.start()                             # and comes back
    c = svc.client()
    data = ws / "fe" / "kcist-4242" / "DELFIN_Data.json"
    _until(data.is_file, 60, "DELFIN_Data.json came back")
    told = _until(lambda: [t for t in _told(c, sid) if "kcist:4242" in t], 90,
                  "the session was told")
    assert "results in fe/kcist-4242" in told[0]



def test_the_service_waits_for_a_followed_job_and_then_stops(tmp_path):
    """Idle -- no window, no session -- the service stops after idle_exit;
    not while a job it follows is queued or running on a logged-in cluster,
    whose results it is to bring back."""
    svc = ServiceProcess(tmp_path, idle_exit=3.0, start=False)
    cluster = Cluster(svc)
    cluster.running()
    ws = svc.workspace()
    cluster.submitted("no-such-session", ws)
    svc.start()
    try:
        time.sleep(8)
        assert svc.proc.poll() is None, "it stopped while the job ran"
        cluster.ended()
        _until(lambda: svc.proc.poll() is not None, 60, "the service stopped once it was back")
        assert (ws / "fe" / "kcist-4242" / "DELFIN_Data.json").is_file()
    finally:
        svc.stop()


def test_with_wakes_turned_off_a_sleeping_session_sleeps_on(tmp_path):
    """agent.wake_on_job_end: false -- the results still come back; the
    session is not woken for them, and they wait for its next turn."""
    svc = ServiceProcess(tmp_path, start=False)
    cluster = Cluster(svc, agent={"wake_on_job_end": False})
    svc.start()
    try:
        cluster.running()
        c, sid, ws = _asleep_after_a_submit(svc, cluster)
        cluster.ended()
        _until((ws / "fe" / "kcist-4242" / "DELFIN_Data.json").is_file, 60, "it came back")
        time.sleep(3)
        assert (_row(c, sid) or {}).get("state") == "sleeping"
        assert not [t for t in _told(c, sid) if "kcist:4242" in t]
        # Opened again, its worker looks at its watches at once -- and with
        # wakes off starts no turn of its own: the job is left to the next
        # turn the user starts, whose start reads it (job_monitor's cluster
        # watch, tested in test_a_cluster_job_comes_back_to_the_session_...).
        c.request("session.attach", {"sessionId": sid, "turns": 0}, timeout=120)
        time.sleep(3)
        assert not [t for t in _told(c, sid) if "kcist:4242" in t]
        watches = job_monitor.load_watched(job_monitor._agent_watch_path(ws))["jobs"]
        assert "kcist:4242" in watches and not watches["kcist:4242"].get("wake_notified")
    finally:
        svc.stop()


def test_a_login_has_the_watcher_look_at_once(monkeypatch):
    """What waited for a login is followed right after it, not a minute on."""
    import asyncio
    from delfin.agent.service import clusters_feed, supervisor
    monkeypatch.setattr(clusters_feed, "login", lambda name, password, otp: {
        "ok": True, "opened": True, "code": None, "error": None})

    async def run() -> tuple[bool, bool]:
        sup = supervisor.Supervisor()
        sup._watch_now = asyncio.Event()
        await sup._m_cluster_login(None, {"name": "kcist", "password": "pw", "otp": ""})
        looked = sup._watch_now.is_set()
        sup._watch_now.clear()
        monkeypatch.setattr(clusters_feed, "login", lambda name, password, otp: {
            "ok": False, "opened": False, "code": "login_failed", "error": "no"})
        await sup._m_cluster_login(None, {"name": "kcist", "password": "pw", "otp": ""})
        return looked, sup._watch_now.is_set()

    assert asyncio.run(run()) == (True, False)
