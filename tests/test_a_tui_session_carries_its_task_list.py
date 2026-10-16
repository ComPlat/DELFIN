"""The session's task list reaches the window as data.

The engine keeps a session's tasks in the workspace's task store
(agent_tasks.TaskStore); the dashboard and the REPL draw them through
task_ticker.rows. The worker sends the same rows -- num, status, label,
in display order -- as ``tasks`` on ``session.info`` and, whenever they
change, on ``session.updated``, so a window draws the list without asking
and a window attaching later gets the current list in the replay header.

Labels are written by the model: they arrive as one plain line, escapes
and control characters made spaces, cut to TASK_LABEL_MAX.
"""

from __future__ import annotations

from delfin.agent.agent_tasks import TaskStore
from delfin.agent.service import protocol
from tui_worker_client import WorkerClient


def _tasks_update(w: WorkerClient, after: int, predicate=lambda rows: True):
    event, i = w.wait_for(lambda e: e.get("type") == "session.updated" and "tasks" in e
                          and predicate(e["tasks"]), after=after)
    return event["tasks"], i


def test_the_list_starts_empty_and_follows_the_store(tmp_path):
    with WorkerClient(tmp_path) as w:
        sid = w.init_result["session"]["sessionId"]
        assert w.init_result["session"]["tasks"] == []
        store = TaskStore(w.workspace)
        mark = len(w.events)
        a = store.create("read the input files", session_id=sid)
        b = store.create("write the parser", active_form="writing the parser", session_id=sid)
        store.create("someone else's work", session_id="another-session")
        rows, at = _tasks_update(w, mark, lambda r: len(r) == 2)
        assert rows == [{"num": 1, "status": "pending", "label": "read the input files"},
                        {"num": 2, "status": "pending", "label": "write the parser"}]
        store.update(a["id"], status="completed")
        store.update(b["id"], status="in_progress")
        rows, at = _tasks_update(w, at + 1, lambda r: r and r[0]["status"] == "in_progress"
                                 and any(x["status"] == "completed" for x in r))
        # The one in hand first, under its active form; the done ones last.
        assert rows == [{"num": 2, "status": "in_progress", "label": "writing the parser"},
                        {"num": 1, "status": "completed", "label": "read the input files"}]
        # Nothing changed, nothing sent: one update per change, not per tick.
        n = len(w.events)
        w.request("session.status")
        w.request("session.status")
        assert not [e for e in w.events[n:] if e.get("type") == "session.updated" and "tasks" in e]


def test_a_label_is_one_plain_line(tmp_path):
    with WorkerClient(tmp_path) as w:
        sid = w.init_result["session"]["sessionId"]
        mark = len(w.events)
        TaskStore(w.workspace).create("fix \x1b[31mred\x1b[0m\nnext line​ " + "x" * 400,
                                      session_id=sid)
        rows, _ = _tasks_update(w, mark, lambda r: len(r) == 1)
        label = rows[0]["label"]
        assert "\x1b" not in label and "\n" not in label and "​" not in label
        assert label.startswith("fix red next line")
        assert len(label) <= protocol.TASK_LABEL_MAX


def test_a_replay_header_carries_the_current_list(tmp_path):
    with WorkerClient(tmp_path) as w:
        sid = w.init_result["session"]["sessionId"]
        mark = len(w.events)
        TaskStore(w.workspace).create("one task", session_id=sid)
        _tasks_update(w, mark, lambda r: len(r) == 1)
        header = w.request("session.replay", {"turns": 1})["events"][0]
        assert header["type"] == "session.info"
        assert header["tasks"] == [{"num": 1, "status": "pending", "label": "one task"}]


def test_the_protocol_declares_the_field():
    assert protocol.EVENTS["session.info"]["tasks"] == "[object]?"
    assert protocol.EVENTS["session.updated"]["tasks"] == "[object]?"


def test_the_list_moves_while_the_turn_runs(tmp_path):
    """What the window draws during a run: each step of the list arrives
    as it happens, before the turn ends -- not one update at its end."""
    with WorkerClient(tmp_path, delay=0.02) as w:
        start, end = w.turn("work through the tasks")
        updates = [(i, e["tasks"]) for i, e in enumerate(w.events[start:end], start)
                   if e.get("type") == "session.updated" and e.get("tasks")]
        in_hand = [rows[0]["label"] for _, rows in updates if rows[0]["status"] == "in_progress"]
        # Each task was seen in hand, in turn, under its active form.
        assert [l for i, l in enumerate(in_hand) if l not in in_hand[:i]] == [
            "reading the input files", "fixing the energy parser", "running the tests"]
        ended = next(i for i in range(start, end + 1) if w.events[i].get("type") == "turn.ended")
        assert updates and all(i < ended for i, _ in updates[:-1])
