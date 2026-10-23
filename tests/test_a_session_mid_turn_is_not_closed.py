"""A session working through a long turn has not gone away.

``session_presence`` expires a record after 15 minutes without a
refresh, and reaps it after an hour. Measured 2026-09-21, a run of six
sessions: presence was renewed only at the prompt and whenever a
question was asked, so a turn longer than fifteen minutes -- one
model round after another, tools in between -- dropped the session
out of ``open_sessions()``. Another session then got

    no other open session 'runde2c-s4'

although s4 was working, and the operator's overview went quiet about
exactly the sessions that were busiest.

The fix judged here: presence is renewed on every RENDER of the
turn's first event too -- the pump runs the whole turn through one
loop, so the heartbeat travels with whatever the turn is doing, at
the announcement's own throttle (``announce`` itself writes at most
once a minute; the pump's call is cheap between them).

What is judged here, and the instrument:

  the refresh travels with the pump    a scripted RawMode whose
                                       read_ready counts the calls;
                                       three reads inside one turn
                                       refresh presence, while the
                                       turn is still running

  no writes beyond the throttle        announce is the same function
                                       the idle prompt uses; a turn
                                       that renders a hundred items
                                       between two reads writes one
                                       record per heartbeat, not one
                                       per render
"""

from __future__ import annotations

import threading
import time

import pytest

from delfin.agent import repl, repl_keys as rk


class _Theme:
    def dim(self, t): return t
    def bold(self, t): return t
    def red(self, t): return t
    def cyan(self, t): return t


class _Transcript:
    theme = _Theme()
    width = 80

    def chrome(self, line): pass
    def render(self, item): pass
    def finish(self): pass
    def refresh_width(self): pass


class _Engine:
    class _Status:
        def get(self): return {}

    kit_permissions = None
    token_usage = {"input": 0, "output": 0}
    status = _Status()

    def get_status(self): return {}
    def request_stop(self): pass


class _CountingRaw:
    """A RawMode stand-in that never delivers a key and counts reads."""

    active = True

    def __init__(self, reads=5):
        self._reads = reads

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        return False

    def read_ready(self, timeout):
        if self._reads <= 0:
            raise AssertionError("pump still reading after the turn ended")
        self._reads -= 1
        return ""


def _agent(monkeypatch, tmp_path, raw):
    from delfin.agent import session_presence as pres
    monkeypatch.setattr(pres, "_DIR", tmp_path / "presence")
    agent = object.__new__(repl.TerminalAgent)
    agent.engine = _Engine()
    agent.broker = None
    agent.opts = type("O", (), {"cwd": tmp_path, "session_name": "runde2-long",
                                "max_tokens": 0})()
    agent.transcript = _Transcript()
    agent._stdin = object()
    agent._turn_active = threading.Event()
    agent._interrupts = 0
    agent._idle_interrupts = 0
    agent._width_dirty = False
    agent._last_paint = 0.0
    agent._q = repl.queue.Queue()
    agent._wake_last_look = 0.0
    agent._bottom = ""
    agent._repaint_bottom = lambda **_kw: None
    agent._clear_bottom = lambda: None
    agent._show_tasks = False
    agent._bg_view = {"shells": [], "agents": [], "watches": [],
                      "wakeups": [], "errors": []}
    agent._bg_status = ""
    agent._steering_refreshes = 0
    agent._last_result_text = ""
    # Presence is announced before the turn, as run() does.
    agent._announce_presence()
    return agent


def test_a_long_turn_keeps_the_session_listed(tmp_path, monkeypatch):
    """Three quiet reads inside ONE turn -- the shape of a turn that
    spends minutes on one model call -- must refresh the presence
    record while the turn is still running. Before the fix the record
    went stale the moment the turn outlasted _STALE_S and no question
    was asked."""
    from delfin.agent import session_presence as pres
    raw = _CountingRaw(reads=3)
    agent = _agent(monkeypatch, tmp_path, raw)

    monkeypatch.setattr(rk, "RawMode", lambda *a, **k: raw)
    # The worker ends the turn after the reads are spent.
    def _worker():
        while raw._reads > 0:
            time.sleep(0.05)
        agent._q.put(repl.RenderItem("done"))
    t = threading.Thread(target=_worker, daemon=True)
    agent._turn_active.set()
    t.start()
    # Simulate the staleness the fix exists for: nothing else has
    # refreshed the record, and the turn is longer than the idle
    # poll ever ran. The heartbeat window is shortened with it, so
    # the announce the pump makes is one the throttle lets through
    # -- as it is in the real 15-minute case.
    monkeypatch.setattr(pres, "_STALE_S", 0.05)
    monkeypatch.setattr(pres, "_HEARTBEAT_S", 0.05)
    time.sleep(0.1)
    with pytest.raises(AssertionError, match="still reading"):
        agent._pump(t)
    # The pump consumed its reads; presence must still list the
    # session -- refreshed by the pump, not by the prompt.
    assert any(r.get("key") == "runde2-long"
               for r in pres.open_sessions()), (
        "a session mid-turn dropped out of open_sessions; that is the "
        "2026-09-21 finding")


def test_the_refresh_is_throttled_not_per_read(tmp_path, monkeypatch):
    """A turn that renders between reads must not WRITE a presence
    record per read: announce itself throttles unchanged records to
    one per heartbeat, and the pump adds no write of its own. The
    count is of writes, not calls -- the pump may call announce as
    often as it reads."""
    from delfin.agent import state_paths as sp
    from delfin.agent import session_presence as pres
    writes = {"n": 0}
    real_write = sp.write_text_atomic

    def _counting_write(path, text, **kw):
        writes["n"] += 1
        return real_write(path, text, **kw)

    monkeypatch.setattr(sp, "write_text_atomic", _counting_write)
    raw = _CountingRaw(reads=4)
    agent = _agent(monkeypatch, tmp_path, raw)
    monkeypatch.setattr(rk, "RawMode", lambda *a, **k: raw)

    def _worker():
        while raw._reads > 0:
            time.sleep(0.05)
        agent._q.put(repl.RenderItem("done"))
    t = threading.Thread(target=_worker, daemon=True)
    agent._turn_active.set()
    t.start()
    try:
        agent._pump(t)
    except AssertionError:
        pass
    # One write before the turn (from _agent's announce) plus at most
    # one heartbeat per minute from the pump -- never one per read:
    # four reads inside the heartbeat window wrote at most one.
    assert writes["n"] <= 2, writes["n"]
