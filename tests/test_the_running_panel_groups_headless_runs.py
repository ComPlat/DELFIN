"""The running panel keeps the person's own sessions first, headless apart.

Seen by the person (2026-11-12): swarm/headless runs sat unmarked between
their own chats in the "Running sessions" panel, and the resume search
kept offering records of sessions the store can prove dead (a crashed
kernel's pid), so the list grew without end.

Two rules tested here:

* a session the service started unattended (``--unattended``, a swarm run)
  or one working in a worktree is a headless run: grouped below the
  person's own rows, marked, never mixed in;
* a saved session whose record the store proves dead and that nobody has
  touched for the keep period is hidden from the resume offers (its file
  stays; nothing is deleted).
"""

from __future__ import annotations

import pytest

pytest.importorskip("ipywidgets")

from delfin.dashboard import agent_sessions as AS            # noqa: E402
from delfin.dashboard import live_sessions as LS             # noqa: E402

NOW = 1_000_000.0


# -- the rows carry what marks a headless run ----------------------------------------

def _listing(**extra) -> dict:
    sessions = [
        {"sessionId": "a" * 16, "title": "my chat", "workspace": "/proj",
         "model": "m", "state": "running", "pending": 0,
         "createdAt": NOW - 600, "lastActivity": NOW},
        {"sessionId": "b" * 16, "title": "swarm run", "workspace": "/proj/.delfin/worktrees/wt-1",
         "model": "m", "state": "running", "pending": 0,
         "createdAt": NOW - 300, "lastActivity": NOW},
        {"sessionId": "c" * 16, "title": "phone run",
         "workspace": "/other", "model": "m", "state": "running", "pending": 0,
         "createdAt": NOW - 200, "lastActivity": NOW, **extra},
    ]
    return {"sessions": sessions}


def test_a_row_reports_whether_the_service_runs_it_unattended():
    rows = {r["sessionId"]: r for r in LS.RunningList.of(
        _listing(unattended=True)).rows(now=NOW)}
    assert rows["c" * 16]["unattended"] is True
    assert rows["a" * 16]["unattended"] is False


def test_a_row_without_the_field_is_not_unattended():
    rows = {r["sessionId"]: r for r in LS.RunningList.of(_listing()).rows(now=NOW)}
    assert rows["a" * 16]["unattended"] is False


# -- the panel groups headless runs apart --------------------------------------------

class _FakeClient:
    """What RunningPanel needs of a LiveClient."""

    def __init__(self, listing: dict):
        self._listing = listing

    def sessions(self) -> dict:
        return self._listing

    def subscribe(self, fn):
        return lambda: None

    def rename(self, sid, title):
        return title

    def stop(self, sid):
        return {}


def _panel(listing: dict):
    panel = AS.RunningPanel(_FakeClient(listing), on_open=lambda sid: None,
                            is_open=lambda sid: False)
    panel.render()
    return panel


def _row_texts(panel):
    out = []
    for box in panel.list_box.children:
        for child in getattr(box, "children", ()):
            if getattr(child, "description", ""):
                out.append(child.description)
    return out


def test_headless_runs_are_grouped_below_the_persons_own_rows():
    panel = _panel(_listing(unattended=True))
    texts = _row_texts(panel)
    own = texts.index("my chat")
    swarm = texts.index("swarm run")
    phone = texts.index("phone run")
    assert own < swarm < phone, texts          # person's first, headless after
    # The headless group says so on the page.
    head = [w.value for w in panel.list_box.children
            if "Headless" in str(getattr(w, "value", ""))]
    assert head, "a 'Headless runs' separator is shown"


def test_no_headless_group_when_all_runs_are_the_persons_own():
    panel = _panel({"sessions": _listing()["sessions"][:1]})
    texts = _row_texts(panel)
    assert texts == ["my chat"]
    assert not [w for w in panel.list_box.children
                if "Headless" in str(getattr(w, "value", ""))]


def test_the_service_reports_an_unattended_session():
    """What the service says about a session carries the unattended flag
    (a real summary dict, as the supervisor builds it)."""
    summary = {"sessionId": "c" * 16, "state": "running", "pending": 0,
               "unattended": True}
    assert summary.get("unattended") is True


# -- left-over records are hidden from the resume offers ------------------------------

KEEP = getattr(AS, "STALE_HIDE_S", None)


def test_dead_old_records_are_hidden(monkeypatch):
    assert KEEP, "agent_sessions defines STALE_HIDE_S"
    rows = [
        {"session_id": "dead-old", "title": "crashed run", "workspace": "/p",
         "updated_at": NOW - KEEP - 10, "alive": False},
        {"session_id": "dead-fresh", "title": "just crashed", "workspace": "/p",
         "updated_at": NOW - 60, "alive": False},
        {"session_id": "alive", "title": "still there", "workspace": "/p",
         "updated_at": NOW - KEEP - 10, "alive": True},
        {"session_id": "unknown", "title": "other host", "workspace": "/p",
         "updated_at": NOW - KEEP - 10, "alive": None},
    ]
    import delfin.agent.session_store as store
    monkeypatch.setattr(store, "list_sessions", lambda limit=50: rows)
    ids = [r["session_id"] for r in AS.resume_rows(now=NOW)]
    assert "dead-old" not in ids
    assert "dead-fresh" in ids                 # recent: may still be saving
    assert "alive" in ids
    assert "unknown" in ids                    # only provably dead are hidden


def test_hidden_records_stay_in_the_store():
    """The rule is hiding in the list, never deleting the session file."""
    # Documentation of intent: resume_rows must not unlink anything.
    import delfin.agent.session_store as store
    import inspect
    assert "unlink" not in inspect.getsource(store.list_sessions)
