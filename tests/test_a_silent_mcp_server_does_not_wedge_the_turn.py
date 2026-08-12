"""A stdio MCP server that took a request and never answered used to hang the
whole agent turn, silently and without end.

``MCPServer._send`` was written as if it had a deadline. It computed
``t_end = time.monotonic() + _RPC_TIMEOUT_S``, looped while the clock was
inside that window, and kept an ``rpc timeout`` error ready for the moment it
fell out. But the read inside the loop was ``self.proc.stdout.readline()`` on a
blocking text-mode pipe, and ``readline()`` returns on a newline or on EOF and
on nothing else. The clock was therefore consulted only BETWEEN lines. Against
a server that emitted no lines the loop never came back around to look at it,
so the ceiling the constant's name advertised bounded chatty servers only, and
the ``rpc timeout`` branch was reachable only when a server was busily writing
output that happened not to carry the request id being waited on. Silence --
the cheapest failure a subprocess can produce -- was the one case the deadline
did not cover.

Where that call sits is what made it expensive. ``api_client`` calls
``MCPRegistry.discover_all()`` while it assembles the tool surface for a turn,
before any model output and before the user has seen a single character, and
``_send`` holds ``self._lock`` across the whole read. One misconfigured stdio
server was thus enough to stop every turn dead with no error, no timeout and
no way to tell from the outside what was waiting on what. ``delfin-tools``
ships enabled by default, so the shipped configuration was one bad interpreter
away from it. The module docstring's promise that "servers fail closed: a crash
during discovery just leaves that server's tool set empty" held for a crash and
not for a hang -- and a hang is the easier accident: a server that starts,
accepts ``initialize`` and then blocks on something of its own looks perfectly
healthy right up until the turn stops responding.

These tests hold the deadline to its word. Silence must cost at most the
budget and then report itself, both in the returned error and in the
``last_error`` field that already existed for exactly this purpose. The four
things that were never broken -- a late reply, an unrelated line, an EOF, a
lock that is given back -- have to survive whatever mechanism enforces it.
"""

from __future__ import annotations

import sys
import textwrap
import threading
import time

import pytest

from delfin.agent import mcp_client as M


# ---------------------------------------------------------------------------
# Fake stdio servers. Each is a self-contained script handed to ``python -c``,
# so the tests exercise the real pipe, the real subprocess and the real text
# mode wrapper -- the mechanism under test IS the pipe, and a stub object with
# a sleeping ``readline`` would let a fix pass that only works against stubs.
# ---------------------------------------------------------------------------

# Accepts everything, answers nothing, stays alive: the defect, minimally.
_SILENT = textwrap.dedent('''
    import sys
    sys.stdin.read()
''')

# Answers, but well after the client has started waiting.
_SLOW = textwrap.dedent('''
    import json, sys, time
    while True:
        line = sys.stdin.readline()
        if not line:
            break
        line = line.strip()
        if not line:
            continue
        req = json.loads(line)
        rid = req.get("id")
        if rid is None:
            continue
        time.sleep(0.6)
        sys.stdout.write(json.dumps(
            {"jsonrpc": "2.0", "id": rid, "result": {"slow": True}}) + "\\n")
        sys.stdout.flush()
''')

# Buries the reply under noise: junk, a notification, a stale id, a blank line.
_CHATTY = textwrap.dedent('''
    import json, sys
    while True:
        line = sys.stdin.readline()
        if not line:
            break
        line = line.strip()
        if not line:
            continue
        req = json.loads(line)
        rid = req.get("id")
        if rid is None:
            continue
        noise = [
            "this line is not JSON at all",
            json.dumps({"jsonrpc": "2.0", "method": "notifications/progress",
                        "params": {"pct": 10}}),
            json.dumps({"jsonrpc": "2.0", "id": 99999,
                        "result": {"stale": True}}),
            "",
            json.dumps({"jsonrpc": "2.0", "id": rid,
                        "result": {"chatty": True}}),
        ]
        sys.stdout.write("\\n".join(noise) + "\\n")
        sys.stdout.flush()
''')

# Closes its stdout but keeps running, so the read hits EOF while ``poll()``
# still reports a live process -- otherwise the "server not running" branch
# would answer first and EOF would never be reached. It has to be
# ``os.close(1)``: CPython builds ``sys.stdout`` with ``closefd=False``, so
# closing the Python object leaves the pipe's write end open and the parent
# waits forever instead of seeing the EOF this test is about.
_CLOSES_STDOUT = textwrap.dedent('''
    import os, sys, time
    sys.stdin.readline()
    os.close(1)
    time.sleep(30)
''')


def _server(script: str, name: str = "fake") -> M.MCPServer:
    return M.MCPServer(name=name, command=sys.executable, args=["-c", script])


def _call_within(limit: float, fn):
    """Run ``fn`` on a throwaway thread and fail if it does not come back.

    A regression here is an unbounded wait, which a plain call would turn into
    a suite that hangs until someone notices rather than a test that fails.
    The thread is a daemon so the interpreter can still exit while the old
    ``readline`` is stuck in the kernel; the caller's ``srv.stop()`` releases
    it by killing the process the pipe belongs to.
    """
    box: dict = {}

    def _run() -> None:
        try:
            box["value"] = fn()
        except BaseException as exc:      # re-raised on the calling thread
            box["error"] = exc

    thread = threading.Thread(target=_run, daemon=True)
    thread.start()
    thread.join(limit)
    if thread.is_alive():
        raise AssertionError(
            f"call did not return within {limit}s -- the read is unbounded")
    if "error" in box:
        raise box["error"]
    return box["value"]


@pytest.fixture
def silent_server(monkeypatch):
    """A running silent server with a small RPC budget."""
    monkeypatch.setattr(M, "_RPC_TIMEOUT_S", 0.75)
    srv = _server(_SILENT, name="silent")
    srv.start()
    assert srv.proc is not None
    try:
        yield srv
    finally:
        srv.stop()


def test_a_server_that_never_answers_does_not_block_forever(silent_server):
    started = time.monotonic()
    resp = _call_within(20.0, lambda: silent_server._send("initialize", {}))
    elapsed = time.monotonic() - started
    assert "rpc timeout" in resp["error"]["message"]
    # The budget is spent, not skipped: a wait that returns instantly would
    # break every server that is merely slow.
    assert 0.5 <= elapsed < 10.0, elapsed


def test_a_timed_out_send_gives_the_lock_back(silent_server):
    _call_within(20.0, lambda: silent_server._send("initialize", {}))
    acquired = silent_server._lock.acquire(blocking=False)
    if acquired:
        silent_server._lock.release()
    assert acquired, "the send lock was still held after the timeout"
    # And the next caller gets its own timeout rather than inheriting the
    # first one's wait -- a lock released in name only would show up here.
    second = _call_within(20.0, lambda: silent_server._send("tools/list"))
    assert "rpc timeout" in second["error"]["message"]


def test_a_timeout_is_recorded_in_last_error(silent_server):
    _call_within(20.0, lambda: silent_server._send("tools/list"))
    assert "timeout" in silent_server.last_error.lower()


def test_a_silent_server_leaves_its_tool_set_empty(silent_server):
    """The behaviour the module docstring already promises for a crash."""
    registry = M.MCPRegistry()
    registry.servers["silent"] = silent_server
    tools = _call_within(20.0, registry.discover_all)
    assert tools == []


def test_a_slow_server_that_answers_late_still_gets_its_reply(monkeypatch):
    monkeypatch.setattr(M, "_RPC_TIMEOUT_S", 10.0)
    srv = _server(_SLOW, name="slow")
    srv.start()
    try:
        resp = _call_within(30.0, lambda: srv._send("initialize", {}))
        assert resp.get("result") == {"slow": True}, resp
    finally:
        srv.stop()


def test_a_chatty_server_still_finds_the_matching_reply(monkeypatch):
    monkeypatch.setattr(M, "_RPC_TIMEOUT_S", 10.0)
    srv = _server(_CHATTY, name="chatty")
    srv.start()
    try:
        resp = _call_within(30.0, lambda: srv._send("initialize", {}))
        assert resp.get("result") == {"chatty": True}, resp
    finally:
        srv.stop()


def test_a_server_that_closes_its_output_is_reported_as_eof(monkeypatch):
    monkeypatch.setattr(M, "_RPC_TIMEOUT_S", 10.0)
    srv = _server(_CLOSES_STDOUT, name="eof")
    srv.start()
    try:
        resp = _call_within(30.0, lambda: srv._send("initialize", {}))
        assert "EOF" in resp["error"]["message"], resp
        # EOF is terminal: a second call must not wait out the whole budget
        # hoping for a line from a pipe that is already closed.
        started = time.monotonic()
        again = _call_within(30.0, lambda: srv._send("tools/list"))
        assert "EOF" in again["error"]["message"], again
        assert time.monotonic() - started < 5.0
    finally:
        srv.stop()
