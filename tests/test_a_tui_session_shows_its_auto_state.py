"""The window sees auto mode: the commands reach the session, the state
reaches the mode line.

Input: a worker (scripted engine, private HOME) driven over the protocol:
``/auto on``, ``/auto off`` and ``/loop`` typed as prompts. Output: the
notices the commands print and the ``auto`` field of ``session.info`` /
``session.updated`` -- auto mode's status line, "" while it is off -- that
the window draws after the posture.
"""

from __future__ import annotations

from delfin.agent.service import protocol
from tui_worker_client import WorkerClient


def _auto_update(w, after, predicate):
    event, i = w.wait_for(lambda e: e.get("type") == "session.updated"
                          and "auto" in e and predicate(e["auto"]),
                          after=after)
    return event["auto"], i


def test_the_protocol_declares_the_field():
    assert protocol.EVENTS["session.info"]["auto"] == "string?"
    assert protocol.EVENTS["session.updated"]["auto"] == "string?"


def test_auto_on_and_off_reach_the_mode_line(tmp_path):
    with WorkerClient(tmp_path) as w:
        assert w.init_result["session"].get("auto", "") == ""
        mark = len(w.events)
        w.request("session.prompt", {"text": "/auto on"})
        line, at = _auto_update(w, mark, lambda a: a.startswith("auto "))
        assert "0/" in line
        notice, _ = w.wait_for(lambda e: e.get("type") == "notice"
                               and "auto mode on" in e.get("text", ""),
                               after=mark)
        w.request("session.prompt", {"text": "/auto off"})
        line, _ = _auto_update(w, at + 1, lambda a: a == "")
        assert line == ""


def test_loop_is_a_window_command(tmp_path):
    with WorkerClient(tmp_path) as w:
        mark = len(w.events)
        w.request("session.prompt", {"text": "/loop"})
        w.wait_for(lambda e: e.get("type") == "notice"
                   and "no loops defined" in e.get("text", ""), after=mark)
        assert not [e for e in w.events[mark:] if e.get("type") == "turn.started"]
