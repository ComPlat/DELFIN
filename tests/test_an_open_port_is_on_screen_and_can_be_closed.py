"""An open port is on screen, in the dashboard and the terminal, and can be
closed there; a server is never started in the foreground.

Asked for 2026-10-20: "I want to see in the dashboard that a port is open
and close it there if need be -- in the TUI and the dashboard." A server a
session left running held its port with nothing on screen saying so. And in
session 7cc700dc the agent ran `python3 -m http.server` with bash: the call
waited its full 120 s and was killed before the agent used bash_background.
"""

from __future__ import annotations

import json
import socket
import sys
import time
from types import SimpleNamespace

import pytest

from delfin.agent import api_client as A
from delfin.agent import background_view as BV
from delfin.agent import bash_jobs as BJ
from delfin.agent import repl_commands as RC


@pytest.fixture
def listening_job(tmp_path):
    with socket.socket() as probe:
        probe.bind(("127.0.0.1", 0))
        port = probe.getsockname()[1]
    code = (f"import socket,time; s=socket.socket(); "
            f"s.bind(('127.0.0.1',{port})); s.listen(); time.sleep(60)")
    job = BJ.get_registry().start(
        command=f"python serve.py --port {port}", cwd=str(tmp_path),
        workspace=str(tmp_path), session_id="sess-ports",
        argv=[sys.executable, "-c", code])
    deadline = time.time() + 10
    while time.time() < deadline and port not in BJ.job_listening_ports(job):
        time.sleep(0.1)
    yield job, port
    BJ.get_registry().kill(job.job_id)
    deadline = time.time() + 10
    while job.poll() is None and time.time() < deadline:
        time.sleep(0.05)


def test_the_ports_a_job_listens_on_are_read(listening_job):
    job, port = listening_job
    assert port in BJ.job_listening_ports(job)
    assert any(r["job_id"] == job.job_id and port in r["ports"]
               for r in BJ.open_ports("sess-ports"))
    assert BJ.open_ports("another-session") == [] or all(
        r["job_id"] != job.job_id for r in BJ.open_ports("another-session"))


def test_the_panel_names_the_port_and_its_address(listening_job):
    job, port = listening_job
    extra = BV._ports_and_url(job.job_id, job.command)
    assert extra == {"ports": [port], "url": f"http://127.0.0.1:{port}/"}
    view = {"shells": [{"id": job.job_id, "label": "serve.py",
                        "since": time.time(), **extra}],
            "watches": [], "agents": [], "wakeups": []}
    row = BV.rows(view)[0]
    assert row["kind"] == "Server" and f"port :{port}" in row["detail"]
    assert "close the ports" in row["tip"]
    assert f"http://127.0.0.1:{port}/" in BV.status_line(view)


def test_the_terminal_lists_and_closes_the_port(listening_job):
    job, port = listening_job
    listed = RC._ports(SimpleNamespace(), "").output
    assert f":{port}" in listed and job.job_id in listed
    closed = RC._ports(SimpleNamespace(), f"close {port}").output
    assert f"closed :{port}" in closed
    deadline = time.time() + 10
    while job.poll() is None and time.time() < deadline:
        time.sleep(0.05)
    assert job.poll() is not None
    assert "no background job holds" in RC._ports(
        SimpleNamespace(), f"close {port}").output


@pytest.mark.parametrize("cmd", [
    "python3 -m http.server 8038 --bind 127.0.0.1 --directory site",
    "cd site && python3 -m http.server 8000", "npm run dev", "npx vite",
    "flask run", "uvicorn app:app --port 8000", "jupyter lab",
    "npx vite --port 5173 --host", "vite dev -l info",
])
def test_a_server_in_the_foreground_is_refused_before_anything_runs(
        tmp_path, cmd):
    asked = []
    perms = A.KitToolPermissions(workspace=str(tmp_path), mode="default",
                                 confirm_callback=lambda *a: asked.append(a)
                                 or True)
    out = json.loads(A._doc_executor.execute(
        "bash", {"command": cmd, "description": "d"}, perms))
    assert "bash_background" in out["error"] and "Nothing was run" in out["error"]
    assert asked == []


@pytest.mark.parametrize("cmd", [
    "python3 -m http.server 8000 &", "timeout 5 python3 -m http.server 8000",
    "npm run build", "vite build", "vite build --mode x", "npm test", "curl http://127.0.0.1:8000",
])
def test_a_bounded_or_background_or_ordinary_command_is_not(cmd):
    assert A._foreground_server(cmd) == ""


def test_a_long_option_list_is_judged_in_linear_time():
    """`--?\\S+` read `--a` two ways, so 22 such words took 1.5 s and
    every further word doubled it. The run is refused before it starts, so
    the words are judged by _foreground_server alone."""
    import time
    t = time.monotonic()
    assert A._foreground_server("vite" + " --a" * 5000 + " build") == ""
    assert A._foreground_server("vite" + " --a" * 5000) != ""
    assert time.monotonic() - t < 1.0
