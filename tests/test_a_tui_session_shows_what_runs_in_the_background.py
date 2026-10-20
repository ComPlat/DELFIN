"""What a session still has out reaches the terminal window as data.

The dashboard draws background_view's rows (shells, watched jobs,
background agents, wake-ups) in a panel; the terminal window had nothing
but `/bash`, which you had to think of asking. The worker now sends the
same rows as ``background`` on ``session.info`` and ``session.updated``:
without ages (the window counts them), so the list goes out when the set
changes and not every second, read at most every few seconds, each label
one plain line.
"""

from __future__ import annotations

import time
import types

from delfin.agent import background_view
from delfin.agent.service import protocol
from delfin.agent.service.worker import Worker


def _view(now: float) -> dict:
    return {"shells": [{"id": "235d02ab", "label": "pytest \x1b[31mtests\n -q",
                        "since": now - 252}],
            "watches": [{"id": "w1", "kind": "ci", "label": "CI on main",
                         "state": "RUNNING", "since": now - 60}],
            "agents": [], "errors": [],
            "wakeups": [{"id": "k1", "kind": "once", "label": "look again",
                         "at": now + 600}]}


def test_the_snapshot_is_the_panels_rows_without_the_clock(monkeypatch):
    now = time.time()
    view = _view(now)
    snap = background_view.snapshot(view)
    assert [r["id"] for r in snap] == [r["id"] for r in background_view.rows(view)]
    assert [r["kind"] for r in snap] == ["Shell", "Watch · CI", "Wake-up"]
    assert snap[0]["since"] == now - 252 and snap[2]["at"] == now + 600
    assert snap[1]["state"] == "running"
    # Nothing in it moves with the clock: two reads a minute apart agree,
    # where the panel's rows do not.
    assert all("detail" not in r for r in snap)
    monkeypatch.setattr(background_view.time, "time", lambda: now + 60)
    assert background_view.snapshot(view) == snap
    assert background_view.rows(view)[0]["detail"] != background_view.rows(view, now=now)[0]["detail"]


def _worker(tmp_path, _unused):
    w = types.SimpleNamespace(
        workspace=str(tmp_path), engine=types.SimpleNamespace(session_id="s1"),
        _background=[], _background_at=-1e9,
        _BACKGROUND_EVERY_S=Worker._BACKGROUND_EVERY_S,
        _BACKGROUND_MAX=Worker._BACKGROUND_MAX)
    return w


def test_the_worker_sends_this_sessions_rows_as_plain_lines(tmp_path, monkeypatch):
    now = time.time()
    asked = []

    def collect(ws, *, session_id=None, now=None):
        asked.append((ws, session_id))
        return _view(time.time())

    monkeypatch.setattr(background_view, "collect", collect)
    w = _worker(tmp_path, None)
    rows = Worker._background_rows(w)
    assert asked == [(str(tmp_path), "s1")]           # scoped to the session
    assert [r["id"] for r in rows] == ["235d02ab", "w1", "k1"]
    label = rows[0]["label"]
    assert "\x1b" not in label and "\n" not in label
    assert len(label) <= protocol.TASK_LABEL_MAX
    assert abs(rows[0]["since"] - (now - 252)) < 5
    # Read again only after the interval: the tick runs every 0.25 s.
    Worker._background_rows(w)
    assert len(asked) == 1
    w._background_at -= Worker._BACKGROUND_EVERY_S
    Worker._background_rows(w)
    assert len(asked) == 2


def test_a_list_that_cannot_be_read_keeps_the_last_one(tmp_path, monkeypatch):
    monkeypatch.setattr(background_view, "collect",
                        lambda ws, *, session_id=None, now=None: _view(time.time()))
    w = _worker(tmp_path, None)
    first = Worker._background_rows(w)
    assert first

    def broken(ws, *, session_id=None, now=None):
        raise OSError("registry unreadable")

    monkeypatch.setattr(background_view, "collect", broken)
    w._background_at -= Worker._BACKGROUND_EVERY_S
    assert Worker._background_rows(w) == first      # not "nothing runs"


def test_the_protocol_declares_the_field():
    assert protocol.EVENTS["session.info"]["background"] == "[object]?"
    assert protocol.EVENTS["session.updated"]["background"] == "[object]?"


def test_the_session_info_a_window_opens_with_carries_the_list(tmp_path):
    """session.info lists its fields one by one (Worker.info); a field
    added to the facts only reached session.updated. Found by a probe
    against a real worker, not by the tests above."""
    from tui_worker_client import WorkerClient

    with WorkerClient(tmp_path) as w:
        assert w.init_result["session"]["background"] == []
        header = w.request("session.replay", {"turns": 1})["events"][0]
        assert header["type"] == "session.info"
        assert header["background"] == []
