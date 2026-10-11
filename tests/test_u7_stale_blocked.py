"""U7 / phase 2 — a blocked reason older than the last message reads as stale.

Finding (wave 13, package U7): the task ticker renders a blocked task's
``waiting on <reason>`` identically whether the block was set a turn ago or
a day ago. A user returning to the panel cannot tell an up-to-date reason
("waiting on the cluster queue") from a stale one that no longer applies
("waiting on the approval from Tuesday"). The reading: a blocked reason set
before the session's last message is stale and is shown so.

The reference is passed in (``last_message_at``, epoch seconds) rather than
read inside the ticker, so the render decision is a pure function of the two
timestamps and the task store — testable without building a session record.
Wiring the real session reference at the call sites (repl / tab_agent) ships
as a `.gate` patch alongside the mark_read callers.

Default: ``last_message_at=None`` marks nothing — every existing caller keeps
today's exact output.
"""

from __future__ import annotations

import pytest

from delfin.agent import task_ticker as TT
from delfin.agent.agent_tasks import get_store


@pytest.fixture
def store(tmp_path):
    from delfin.agent import agent_tasks
    agent_tasks._STORES.clear()
    return get_store(tmp_path)


def _blocked(store, subject, reason="waiting on approval"):
    t = store.create(subject, session_id="s")
    return store.update(t["id"], status="blocked", blocked_reason=reason)


def test_a_blocked_task_older_than_the_last_message_is_stale(store, tmp_path):
    _blocked(store, "push the branch")
    row = TT.rows(tmp_path, session_id="s", last_message_at=10**13)
    assert " (stale)" in row[0]["label"]


def test_a_blocked_task_newer_than_the_last_message_is_not_stale(store, tmp_path):
    _blocked(store, "push the branch")
    row = TT.rows(tmp_path, session_id="s", last_message_at=0.0)
    assert " (stale)" not in row[0]["label"]
    assert "waiting on approval" in row[0]["label"]


def test_render_text_marks_the_stale_reason(store, tmp_path):
    _blocked(store, "run the suite", reason="waiting on the queue")
    text = TT.render_text(tmp_path, session_id="s", last_message_at=10**13)
    assert "waiting on the queue" in text and "(stale)" in text


def test_non_blocked_tasks_never_carry_the_stale_marker(store, tmp_path):
    t = store.create("plan the migration", session_id="s")
    store.update(t["id"], status="in_progress")
    store.update(t["id"], status="completed")
    rows = TT.rows(tmp_path, session_id="s", last_message_at=10**13)
    for r in rows:
        assert " (stale)" not in r["label"]
    assert rows


def test_default_reference_marks_nothing_so_callers_stay_unchanged(store, tmp_path):
    _blocked(store, "push the branch")
    rows = TT.rows(tmp_path, session_id="s")
    assert not rows or all(" (stale)" not in r["label"] for r in rows)


def test_a_fresh_block_equals_the_reference_is_not_stale(store, tmp_path):
    # A block set at the same moment as the last message is current, not
    # stale: only a block that strictly precedes it reads as outdated.
    _blocked(store, "push the branch")
    updated = store.get(1)["updated_at"]
    import datetime as _dt
    equal = _dt.datetime.fromisoformat(updated.rstrip("Z")).replace(
        tzinfo=_dt.timezone.utc).timestamp()
    row = TT.rows(tmp_path, session_id="s", last_message_at=equal)
    assert " (stale)" not in row[0]["label"]


def test_a_non_numeric_last_message_at_never_takes_the_render_down(store, tmp_path):
    """Reviewer finding (nacht-s12): ``_stale_suffix`` ran ``float()`` on the
    reference unguarded, so a non-numeric ``last_message_at`` (a string
    arriving over a wire) raised out of rows()/render_text()/render_html() and
    killed the whole panel — breaking ``_block_epoch``'s "never raises: a
    ticker is a courtesy" contract. The reference must be treated as "nothing
    stale" when it cannot parse, not as a crash."""
    _blocked(store, "push the branch")
    row = TT.rows(tmp_path, session_id="s", last_message_at="not-a-number")
    assert " (stale)" not in row[0]["label"]
    text = TT.render_text(tmp_path, session_id="s", last_message_at="oops")
    assert " (stale)" not in text
    html = TT.render_html(tmp_path, session_id="s", last_message_at=object())
    assert "(stale)" not in html
