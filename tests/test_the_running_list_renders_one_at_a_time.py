"""The running list is drawn by one thread at a time.

RunningPanel.render runs from the event timer (a summary arrived) and from
the session list's 3 s tick (the ages). Two renders at once each found a
new session's row missing and built one: the list then showed one row's
widgets while the panel's buttons and Stops belonged to the other, and a
click on the shown row did nothing.

Driven with the panel's own render; the row builder is slowed so the two
renders overlap the way a tick and an event do. Control: one render alone
is consistent.
"""

from __future__ import annotations

import threading
import time

import pytest

pytest.importorskip("ipywidgets")


class _Client:
    def subscribe(self, fn):
        return lambda: None

    def sessions(self):
        return {"sessions": [{"sessionId": "s" * 16, "title": "one", "workspace": "/w/one",
                              "model": "m", "state": "idle", "pending": 0,
                              "createdAt": time.time(), "lastActivity": time.time()}]}


def _consistent(panel) -> bool:
    shown = list(panel.list_box.children)
    rows = panel._rows
    return (len(shown) == len(rows) == 1
            and all(rows[sid]["box"] in shown for sid in rows)
            and all(panel.buttons[sid] is rows[sid]["title"] for sid in rows)
            and all(panel.stops[sid] in rows[sid]["box"].children for sid in rows))


def _panel(monkeypatch, slow: bool):
    from delfin.dashboard import agent_sessions as AS

    panel = AS.RunningPanel(None, on_open=lambda sid: None, is_open=lambda sid: False)
    panel.built = []
    if slow:
        # A widget build inside the row builder, slowed: between finding
        # the row missing and keeping the new one.
        real_text = AS.widgets.Text

        def slow_text(*a, **kw):
            panel.built.append(threading.current_thread().name)
            time.sleep(0.2)
            return real_text(*a, **kw)
        monkeypatch.setattr(AS.widgets, "Text", slow_text)
    panel.running.load(_Client().sessions())
    return panel


def test_one_render_alone_is_consistent(monkeypatch):
    panel = _panel(monkeypatch, slow=True)
    panel.render()
    assert _consistent(panel) and len(panel.built) == 1


def test_two_renders_at_once_leave_one_consistent_list(monkeypatch):
    panel = _panel(monkeypatch, slow=True)
    threads = [threading.Thread(target=panel.render) for _ in range(2)]
    for t in threads:
        t.start()
    for t in threads:
        t.join(timeout=10)
    # One row for one session: a second one built beside it is a set of
    # widgets whose buttons nothing shows (or the shown ones answer nothing).
    assert len(panel.built) == 1, panel.built
    assert _consistent(panel)
