"""`delfin-agent service status` lists the sessions, `service end` ends one.

* status: one row per session -- short id, title, workspace (``~`` for the
  home), what it does in words (working, waiting for you, idle, sleeping),
  age, last activity, the windows open on it; ``--json`` the same for
  scripts.
* end: a session ends the way the agent view's ctrl+x twice ends it
  (``session.shutdown``): saved, resumable, the service and the other
  sessions untouched. A unique prefix of the id names it; an ambiguous or
  unknown one ends nothing. ``--all-idle`` ends what is idle and leaves
  what works or waits. Without a terminal to ask on, ``--yes`` is needed.
* Only the user ends a session: from an agent session's own process the
  command is refused before it reaches the service, and the permission gate
  refuses the command under any spelling (as ``service stop``).
"""

from __future__ import annotations

import argparse
import json
import subprocess
import sys
import time

import pytest

from conftest import child_env
from delfin.agent import approval_guard
from delfin.agent.service import cli_service, client
from delfin.agent.service.client import ServiceClient
from tui_service import ServiceProcess, turn


def _env(service: ServiceProcess) -> dict:
    env = child_env(service.tmp)
    env.update(DELFIN_RUNTIME_DIR=service.env["DELFIN_RUNTIME_DIR"],
               PYTHONPATH=service.env["PYTHONPATH"])
    return env


def _cli(service: ServiceProcess, *args: str, env: dict | None = None,
         stdin: str | None = None) -> subprocess.CompletedProcess:
    return subprocess.run([sys.executable, "-m", "delfin.agent.cli", "service", *args],
                          env=env or _env(service), capture_output=True, text=True,
                          timeout=120, input=stdin)


def _listed(c: ServiceClient) -> dict:
    return {s["sessionId"]: s for s in c.request("sessions.list")["sessions"]}


@pytest.fixture
def service(tmp_path):
    svc = ServiceProcess(tmp_path)
    try:
        yield svc
    finally:
        svc.stop()


# -- the words and the rows ---------------------------------------------------------------

@pytest.mark.parametrize("summary, words", [
    ({"state": "running", "pending": 0}, "working"),
    ({"state": "starting", "pending": 0}, "working"),
    ({"state": "needs_input", "pending": 1}, "waiting for you"),
    ({"state": "running", "pending": 2}, "waiting for you"),
    ({"state": "idle", "pending": 0}, "idle"),
    ({"state": "sleeping", "pending": 0}, "sleeping"),
])
def test_a_state_is_said_in_words(summary, words):
    assert cli_service.state_words(summary) == words


def test_a_prefix_names_one_session_or_none():
    rows = [{"id": "abc12345"}, {"id": "abd99999"}, {"id": "ffff0000"}]
    assert cli_service.matching(rows, "abc") == [rows[0]]
    assert cli_service.matching(rows, "ab") == rows[:2]
    assert cli_service.matching(rows, "zz") == []
    assert cli_service.matching(rows, "") == []


def test_a_workspace_is_shortened_from_the_home(monkeypatch):
    monkeypatch.setenv("HOME", "/home/u")
    assert cli_service.short_workspace("/home/u/calc/fe") == "~/calc/fe"
    assert cli_service.short_workspace("/scratch/x") == "/scratch/x"
    long = "/home/u/" + "/".join(["deep"] * 20)
    assert len(cli_service.short_workspace(long, 30)) <= 30
    assert cli_service.short_workspace(long, 30).startswith("~/")


# -- status ---------------------------------------------------------------------------------

def test_status_lists_every_session_with_what_it_does_and_who_looks(service):
    a = service.client()
    home_ws = service.home / "calc" / "fe"
    home_ws.mkdir(parents=True)
    one = service.session(a, cwd=str(home_ws), name="relax Fe complex")
    two = service.session(a, name="scan benzene")
    turn(a, two, "hello")
    a.request("session.detach", {"sessionId": two})
    # A window on the first: the service counts the TUI's connections.
    window = ServiceClient(service.socket, client="delfin-tui")
    service.clients.append(window)
    window.request("session.attach", {"sessionId": one, "turns": 0})
    # The second one asks: it waits for the user.
    start = len(a.events)
    a.request("session.attach", {"sessionId": two, "turns": 0})
    a.request("session.prompt", {"sessionId": two, "text": "ask me"})
    a.wait_for(lambda e: e.get("type") == "permission.asked" and e.get("session") == two,
               after=start, timeout=30)

    out = _cli(service, "status")
    assert out.returncode == 0, out.stderr
    lines = out.stdout.splitlines()
    assert "2 sessions" in lines[0] and "1 window" in lines[0]
    head = next(line for line in lines if "TITLE" in line)
    for column in ("ID", "TITLE", "WORKSPACE", "STATE", "AGE", "ACTIVE", "WINDOWS"):
        assert column in head
    row_one = next(line for line in lines if one[:8] in line)
    row_two = next(line for line in lines if two[:8] in line)
    assert "relax Fe complex" in row_one and "~/calc/fe" in row_one
    assert row_one.rstrip().endswith("1")
    assert "idle" in row_one
    assert "scan benzene" in row_two and "waiting for you" in row_two
    assert row_two.rstrip().endswith("0")

    data = json.loads(_cli(service, "status", "--json").stdout)
    assert data["service"]["pid"] == service.proc.pid
    rows = {r["id"]: r for r in data["sessions"]}
    assert set(rows) == {one, two}
    assert rows[one]["state"] == "idle" and rows[one]["windows"] == 1
    assert rows[one]["workspace"] == str(home_ws)
    assert rows[two]["state"] == "waiting for you" and rows[two]["pending"] == 1
    for row in rows.values():
        assert row["ageS"] >= 0 and row["lastActivityS"] >= 0
        assert row["title"]


def test_status_json_says_when_no_service_runs(tmp_path):
    svc = ServiceProcess(tmp_path, start=False)
    try:
        out = _cli(svc, "status", "--json")
        assert out.returncode == 0, out.stderr
        data = json.loads(out.stdout)
        assert data["service"] is None and data["sessions"] == []
    finally:
        svc.stop()


# -- end ----------------------------------------------------------------------------------------

def test_end_ends_one_session_by_its_prefix_and_keeps_it_resumable(service):
    a = service.client()
    keep = service.session(a, name="keep me")
    end = service.session(a, name="end me", fakeSave=True)
    turn(a, end, "hello")
    mark = len(a.events)
    out = _cli(service, "end", end[:6], "--yes")
    assert out.returncode == 0, out.stdout + out.stderr
    assert "end me" in out.stdout and "saved" in out.stdout
    removed, _ = a.wait_for(lambda e: e.get("type") == "session.removed" and e.get("session") == end,
                            after=mark, timeout=60)
    assert removed["reason"] == "shutdown"
    listed = _listed(a)
    assert end not in listed and keep in listed
    assert service.proc.poll() is None
    # Saved: offered again by the resume picker.
    past = a.request("sessions.past", {"workspace": str(service.workspace())})["sessions"]
    assert end in [r["sessionId"] for r in past]


def test_end_with_an_unknown_or_ambiguous_prefix_ends_nothing(service, monkeypatch):
    a = service.client()
    sid = service.session(a)
    out = _cli(service, "end", "zzzz", "--yes")
    assert out.returncode == 2 and "no session" in out.stderr
    assert sid in _listed(a)
    # Ambiguous: staged in process, two ids with one prefix.
    rows = [{"id": "abc1", "title": "x", "state": "idle"}, {"id": "abc2", "title": "y", "state": "idle"}]
    monkeypatch.setattr(cli_service, "_session_rows", lambda svc: (rows, {}))
    asked: list = []

    class FakeClient:
        def __init__(self, *a, **k): ...
        def __enter__(self): return self
        def __exit__(self, *exc): ...
        def request(self, method, params=None, timeout=30.0):
            asked.append(method)
            return {}

    monkeypatch.setattr(client, "probe", lambda: ({"pid": 1}, ""))
    monkeypatch.setattr(client, "ServiceClient", FakeClient)
    monkeypatch.delenv(approval_guard.ENV_MARKER, raising=False)
    code = cli_service.cmd_service(argparse.Namespace(action="end", target="abc", all_idle=False,
                                                      yes=True, json=False))
    assert code == 2
    assert "session.shutdown" not in asked


def test_end_without_a_terminal_needs_yes(service):
    a = service.client()
    sid = service.session(a)
    out = _cli(service, "end", sid[:8], stdin="y\n")
    assert out.returncode == 2 and "--yes" in out.stderr
    assert sid in _listed(a)


def _on_a_terminal(service: ServiceProcess, answer: str, *args: str) -> tuple[int, str]:
    """`service end` with a pty as its terminal: ``answer`` typed at its
    question; its exit code and what it wrote."""
    import os
    import select

    main, side = os.openpty()
    proc = subprocess.Popen([sys.executable, "-m", "delfin.agent.cli", "service", *args],
                            env=_env(service), stdin=side, stdout=side, stderr=side,
                            close_fds=True)
    os.close(side)
    seen = b""
    deadline = time.monotonic() + 60
    answered = False
    while time.monotonic() < deadline:
        ready, _, _ = select.select([main], [], [], 0.2)
        if ready:
            try:
                chunk = os.read(main, 4096)
            except OSError:
                break
            if not chunk:
                break
            seen += chunk
        if not answered and b"[y/N]" in seen:
            os.write(main, answer.encode() + b"\n")
            answered = True
        if proc.poll() is not None and not ready:
            break
    code = proc.wait(timeout=60)
    os.close(main)
    return code, seen.decode(errors="replace")


def test_end_on_a_terminal_asks_and_no_ends_nothing(service):
    a = service.client()
    sid = service.session(a, name="ask first")
    code, out = _on_a_terminal(service, "n", "end", sid[:8])
    assert "ask first" in out and "[y/N]" in out and "nothing ended" in out, out
    assert code == 1
    assert sid in _listed(a)
    code, out = _on_a_terminal(service, "y", "end", sid[:8])
    assert code == 0 and "saved" in out, out
    assert sid not in _listed(a)


def test_end_all_idle_leaves_what_waits_for_the_user(service):
    a = service.client()
    idle = service.session(a, name="idle one")
    asks = service.session(a, name="asking one")
    start = len(a.events)
    a.request("session.prompt", {"sessionId": asks, "text": "ask me"})
    a.wait_for(lambda e: e.get("type") == "permission.asked" and e.get("session") == asks,
               after=start, timeout=30)
    out = _cli(service, "end", "--all-idle", "--yes")
    assert out.returncode == 0, out.stdout + out.stderr
    deadline = time.monotonic() + 60
    while idle in _listed(a) and time.monotonic() < deadline:
        time.sleep(0.2)
    listed = _listed(a)
    assert idle not in listed
    assert listed[asks]["pending"] == 1


# -- only the user ----------------------------------------------------------------------------

def test_end_is_refused_inside_an_agent_session(service):
    a = service.client()
    sid = service.session(a)
    env = _env(service)
    env[approval_guard.ENV_MARKER] = "4242"          # as in an agent session's shell
    for args in (("end", sid[:8], "--yes"), ("end", "--all-idle", "--yes")):
        out = _cli(service, *args, env=env)
        assert out.returncode == 1, out.stdout + out.stderr
        assert "agent session" in out.stderr and "Traceback" not in out.stderr
    assert _listed(a)[sid]["state"] == "idle"


def test_end_never_falls_back_to_a_signal(monkeypatch, capsys):
    killed: list = []
    monkeypatch.delenv(approval_guard.ENV_MARKER, raising=False)
    monkeypatch.setattr(client, "probe", lambda: (None, "refused"))
    monkeypatch.setattr(cli_service.rundir, "service_pid", lambda: 999_999)
    monkeypatch.setattr(cli_service.os, "kill", lambda pid, sig: killed.append((pid, sig)))
    code = cli_service.cmd_service(argparse.Namespace(action="end", target="abc", all_idle=False,
                                                      yes=True, json=False))
    assert code == 1 and killed == []


@pytest.mark.parametrize("cmd", [
    "delfin-agent service end abc123 --yes",
    "delfin-agent service end --all-idle --yes",
    "python -m delfin.agent.cli service end abc",
    "del'fin'-agent \"service\" e\\nd abc",
    "cd /tmp && delfin-agent service --yes end abc",
])
def test_the_gate_refuses_ending_a_session(cmd):
    assert approval_guard.answer_command(cmd)


@pytest.mark.parametrize("cmd", [
    "delfin-agent service status --json",
    "grep -n endpoint delfin/agent/service/client.py",
    "echo the service endpoint",
])
def test_reading_the_sessions_stays_allowed(cmd):
    assert approval_guard.answer_command(cmd) is None
