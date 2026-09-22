"""Operator mail reaches the turn that is running.

The dashboard steers a message into a running turn between rounds
(push_steer, tab_agent._deliver_session_messages); the terminal did
not -- a message that arrived mid-turn waited at the prompt until the
turn ended, and one that arrived during a very long turn was not
seen for exactly as long (the operator's answer to assignment 8,
point b, 2026-09-22: steer between rounds, with sender, never as
user input; prompt delivery stays the fallback).

The delivery point is the pump, which already owns the terminal for
the whole turn: it looks for mail on the same throttle the idle wake
uses, renders it the same way (so it reads as coming from another
session, never from the user), and steers it. A backend without
push_steer leaves the mail in the inbox -- the prompt path delivers
it at the end of the turn, which is what a fallback is for.

What is judged here, and the instrument:

  mail mid-turn is steered       a scripted RawMode inside a running
                                 turn, a real inbox with one
                                 message: engine.steer receives the
                                 rendered text, sender header
                                 included, and the inbox is empty

  not steered as user input      the steered text carries the
                                 "not from the user" header, so no
                                 later reader mistakes the operator
                                 for the person at the keyboard

  no backend, no take            an engine without steer support
                                 leaves the message in the inbox for
                                 the prompt -- taking it would lose
                                 it, which is worse than waiting

  the throttle holds             the pump looks for mail on the
                                 wake's clock, not on every read
"""

from __future__ import annotations

import threading
import time

from delfin.agent import repl, repl_keys as rk


class _Theme:
    def dim(self, t): return t
    def bold(self, t): return t
    def red(self, t): return t
    def cyan(self, t): return t


class _Transcript:
    theme = _Theme()
    width = 80

    def __init__(self):
        self.lines = []

    def chrome(self, line):
        self.lines.append(line)

    def render(self, item): pass
    def finish(self): pass
    def refresh_width(self): pass


class _Engine:
    """Records steers; enough surface for the pump."""

    class _Status:
        def get(self): return {}

    token_usage = {"input": 0, "output": 0}

    def __init__(self, supports_steer=True):
        self.steered = []
        self._supports = supports_steer
        self.kit_permissions = None

    def get_status(self):
        return {}

    def request_stop(self):
        pass

    def steer(self, text):
        if not self._supports:
            return False
        self.steered.append(text)
        return True


class _QuietRaw:
    """A RawMode that never delivers a key and ends after N reads."""

    active = True

    def __init__(self, reads=3):
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


def _agent(engine, raw, tmp_path):
    agent = object.__new__(repl.TerminalAgent)
    agent.engine = engine
    agent.broker = None
    agent.opts = type("O", (), {"cwd": tmp_path, "session_name": "runde2-s8",
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
    agent._bg_status_at = 0.0
    agent._bg_id_at = lambda _row: None
    agent._steering_refreshes = 0
    agent._last_result_text = ""
    return agent


def _run_one_turn(engine, raw, tmp_path, monkeypatch):
    agent = _agent(engine, raw, tmp_path)
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
    return agent


def test_mail_mid_turn_is_steered_with_its_sender(
        tmp_path, monkeypatch):
    from delfin.agent import session_messages as msgs
    monkeypatch.setattr(msgs, "_DIR", tmp_path / "inbox")
    msgs.send("runde2-s8", "leave delfin/api.py alone, session 7 owns it",
              from_title="the operator")
    engine = _Engine()
    agent = _run_one_turn(engine, _QuietRaw(reads=3), tmp_path, monkeypatch)
    assert len(engine.steered) == 1, (
        "mail waiting while a turn ran never reached it; the dashboard "
        "steers it between rounds, the terminal dropped it at the prompt")
    steered = engine.steered[0]
    assert "not from the user" in steered
    assert "session 7 owns it" in steered
    # The message was taken, not copied: the inbox holds nothing for
    # a second delivery.
    assert msgs.take("runde2-s8") == []
    # And the pane said so, because a steer is not echoed either.
    assert any("session 7 owns it" in line for line in
               agent.transcript.lines)


def test_no_steer_support_leaves_the_mail_for_the_prompt(
        tmp_path, monkeypatch):
    """Taking mail a backend cannot deliver would lose it: the
    fallback is the prompt, so the message waits."""
    from delfin.agent import session_messages as msgs
    monkeypatch.setattr(msgs, "_DIR", tmp_path / "inbox")
    msgs.send("runde2-s8", "later", from_title="the operator")
    engine = _Engine(supports_steer=False)
    _run_one_turn(engine, _QuietRaw(reads=3), tmp_path, monkeypatch)
    assert engine.steered == []
    assert msgs.take("runde2-s8"), "the message must still be waiting"


def test_the_look_for_mail_is_throttled(tmp_path, monkeypatch):
    """Three reads inside one throttle window take the inbox at most
    once; a turn may not read the filesystem on every key tick."""
    from delfin.agent import session_messages as msgs
    takes = {"n": 0}
    real_take = msgs.take

    def _counting_take(key):
        takes["n"] += 1
        return real_take(key)

    monkeypatch.setattr(msgs, "_DIR", tmp_path / "inbox")
    monkeypatch.setattr(msgs, "take", _counting_take)
    engine = _Engine()
    _run_one_turn(engine, _QuietRaw(reads=4), tmp_path, monkeypatch)
    # No mail: the take happens once per throttle window, never per
    # read -- and a turn with no mail still costs at most one look.
    assert takes["n"] <= 1, takes["n"]
