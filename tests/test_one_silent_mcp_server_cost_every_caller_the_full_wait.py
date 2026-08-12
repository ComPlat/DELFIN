"""One MCP server object, one lock, and everybody waited in turn.

The defect, measured before this existed: ``MCPServer._send`` took the
server lock, wrote the request, and then held that lock for the ENTIRE
reply wait — and the HTTP path held it across ``urlopen``. The registry
is a process singleton, so the parent turn and every background
sub-agent share one server object and one lock. A server that accepts a
request and never answers therefore cost each waiter the full 15-second
deadline SERIALLY: four callers measured 1s, 2s, 3s, 4s against a
one-second deadline, so four times the ceiling that was supposed to
bound them. The reader thread had already been moved off the lock for
exactly this reason; the waiter had not.

Two more in the same file: ``stop()`` cleared ``self.proc`` without the
lock while ``_send`` was dereferencing it, and the registry was keyed by
nothing at all — a second session in another repo was silently handed
the first repo's servers, spawned with the first repo's environment,
with its own project config never read and no symptom to see.

What is pinned here: concurrent callers of one server wait concurrently,
each still gets its OWN reply, and the ordinary single-caller path is
unchanged.
"""

from __future__ import annotations

import json
import sys
import threading
import time
from pathlib import Path

import pytest

from delfin.agent import mcp_client


# A server that reads every request and answers none.
_MUTE = "import sys\nfor line in sys.stdin:\n    pass\n"

# A server that answers each request with its own id, after a short
# delay, and answers them OUT OF ORDER — the reply to the second request
# arrives first, so a caller that took whatever line turned up would get
# somebody else's answer.
_ECHO = """
import json, sys, threading, time
def answer(msg, delay):
    time.sleep(delay)
    sys.stdout.write(json.dumps({
        "jsonrpc": "2.0", "id": msg.get("id"),
        "result": {"method": msg.get("method"),
                   "echo": (msg.get("params") or {}).get("echo")}}) + "\\n")
    sys.stdout.flush()
n = 0
for line in sys.stdin:
    line = line.strip()
    if not line:
        continue
    msg = json.loads(line)
    if msg.get("id") is None:
        continue
    n += 1
    threading.Thread(target=answer, args=(msg, 0.30 / n), daemon=True).start()
"""


@pytest.fixture
def server_factory(tmp_path):
    made: list[mcp_client.MCPServer] = []

    def _make(source: str) -> mcp_client.MCPServer:
        script = tmp_path / f"srv{len(made)}.py"
        script.write_text(source, encoding="utf-8")
        srv = mcp_client.MCPServer(
            name=f"srv{len(made)}", command=sys.executable, args=[str(script)])
        srv.start()
        assert srv.proc is not None
        made.append(srv)
        return srv

    yield _make
    for srv in made:
        srv.stop()


# ---------------------------------------------------------------------------
# The measurement
# ---------------------------------------------------------------------------

def test_a_silent_server_costs_each_caller_the_deadline_once(
        monkeypatch, server_factory):
    monkeypatch.setattr(mcp_client, "_RPC_TIMEOUT_S", 1.0)
    srv = server_factory(_MUTE)

    waits: list[float] = []
    started = time.monotonic()

    def caller() -> None:
        t0 = time.monotonic()
        assert "error" in srv._send("tools/list")
        waits.append(time.monotonic() - t0)

    threads = [threading.Thread(target=caller) for _ in range(4)]
    for t in threads:
        t.start()
    for t in threads:
        t.join()

    total = time.monotonic() - started
    # Four callers, one deadline — not four of them end to end.
    assert total < 2.5, waits
    assert max(waits) < 2.0, waits
    assert srv.last_error.startswith("rpc timeout")


def test_each_concurrent_caller_gets_its_own_reply(server_factory):
    """Replies are matched by id in the reader thread. Out of order on
    purpose: whoever answered first must not answer for everybody."""
    srv = server_factory(_ECHO)
    got: dict[str, dict] = {}
    lock = threading.Lock()

    def caller(tag: str) -> None:
        resp = srv._send("tools/call", {"echo": tag})
        with lock:
            got[tag] = resp

    threads = [threading.Thread(target=caller, args=(f"t{i}",))
               for i in range(4)]
    for t in threads:
        t.start()
    for t in threads:
        t.join()

    assert sorted(got) == ["t0", "t1", "t2", "t3"]
    for tag, resp in got.items():
        assert resp["result"]["echo"] == tag, (tag, resp)
    # Distinct ids, one per request.
    assert len({resp["id"] for resp in got.values()}) == 4


def test_a_reply_the_caller_gave_up_on_is_not_handed_to_the_next_one(
        monkeypatch, server_factory):
    monkeypatch.setattr(mcp_client, "_RPC_TIMEOUT_S", 0.05)
    srv = server_factory(_ECHO)
    assert "error" in srv._send("tools/list")     # times out; reply lands late
    time.sleep(0.5)
    monkeypatch.setattr(mcp_client, "_RPC_TIMEOUT_S", 5.0)
    resp = srv._send("tools/call", {"echo": "mine"})
    assert resp["result"]["echo"] == "mine"


# ---------------------------------------------------------------------------
# The lifecycle takes the lock
# ---------------------------------------------------------------------------

def test_stopping_a_server_under_a_waiting_caller_ends_that_caller(
        monkeypatch, server_factory):
    monkeypatch.setattr(mcp_client, "_RPC_TIMEOUT_S", 30.0)
    srv = server_factory(_MUTE)
    out: list[dict] = []

    def caller() -> None:
        out.append(srv._send("tools/list"))

    t = threading.Thread(target=caller)
    t.start()
    time.sleep(0.2)
    srv.stop()
    t.join(timeout=10)
    assert not t.is_alive()
    assert out and "error" in out[0]
    assert srv.proc is None


def test_a_call_after_eof_answers_at_once_rather_than_waiting(
        monkeypatch, server_factory):
    monkeypatch.setattr(mcp_client, "_RPC_TIMEOUT_S", 5.0)
    srv = server_factory("import sys\nsys.exit(0)\n")
    for _ in range(50):
        if srv.proc.poll() is not None:
            break
        time.sleep(0.05)
    started = time.monotonic()
    resp = srv._send("tools/list")
    assert "error" in resp
    assert time.monotonic() - started < 1.0


def test_start_and_stop_hold_the_lock():
    import inspect

    for method in (mcp_client.MCPServer.start, mcp_client.MCPServer.stop):
        assert "self._lock" in inspect.getsource(method), method.__name__


# ---------------------------------------------------------------------------
# The registry is keyed by the workspace it was built for
# ---------------------------------------------------------------------------

def test_a_second_workspace_does_not_inherit_the_first_ones_servers(
        monkeypatch, tmp_path):
    monkeypatch.setattr(Path, "home", staticmethod(lambda: tmp_path))
    repo_a, repo_b = tmp_path / "repo-a", tmp_path / "repo-b"
    for repo, name in ((repo_a, "only-in-a"), (repo_b, "only-in-b")):
        (repo / ".delfin").mkdir(parents=True)
        (repo / ".delfin" / "mcp_servers.json").write_text(json.dumps({
            "servers": {name: {"command": "true", "enabled": True}}}),
            encoding="utf-8")

    mcp_client.reset_registry()
    try:
        reg_a = mcp_client.get_registry(repo_a)
        reg_b = mcp_client.get_registry(repo_b)
        assert reg_a is not reg_b
        assert "only-in-a" in reg_a.servers and "only-in-a" not in reg_b.servers
        assert "only-in-b" in reg_b.servers and "only-in-b" not in reg_a.servers
        # Same workspace, written differently, is the same registry.
        assert mcp_client.get_registry(Path(str(repo_a) + "/")) is reg_a
    finally:
        mcp_client.reset_registry()


def test_resetting_one_workspace_leaves_the_other_alone(monkeypatch, tmp_path):
    monkeypatch.setattr(Path, "home", staticmethod(lambda: tmp_path))
    repo_a, repo_b = tmp_path / "a", tmp_path / "b"
    repo_a.mkdir()
    repo_b.mkdir()
    mcp_client.reset_registry()
    try:
        reg_a = mcp_client.get_registry(repo_a)
        reg_b = mcp_client.get_registry(repo_b)
        mcp_client.reset_registry(repo_a)
        assert mcp_client.get_registry(repo_a) is not reg_a
        assert mcp_client.get_registry(repo_b) is reg_b
    finally:
        mcp_client.reset_registry()
