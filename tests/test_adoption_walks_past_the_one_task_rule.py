"""Two tasks in progress in one session, and neither tool refused it.

The store enforces one ``in_progress`` per session, and it enforced it on
the STATUS TRANSITION only::

    if new_status == "in_progress":
        ... refuse if this session already has one

``task_adopt`` changes no status. It writes ``session_id`` and nothing
else, to take over a task a previous session left open. A task that was
already ``in_progress`` under its old owner therefore arrived in a
session that had one of its own, and the invariant every reader depends
on -- the ticker, the per-turn reminder, the completion window -- was
false with no tool having refused anything.

The rule was written as a condition on one event instead of as a
property of the store, which is the shape that keeps recurring here.
"""

from __future__ import annotations

import pytest

from delfin.agent import agent_tasks


@pytest.fixture
def store(tmp_path):
    agent_tasks._STORES.clear()
    return agent_tasks.TaskStore(tmp_path)


def _running(store, *, session, subject):
    t = store.create(subject=subject, session_id=session)
    store.update(t["id"], status="in_progress")
    return t["id"]


# ---------------------------------------------------------------------------
# The hole
# ---------------------------------------------------------------------------

def test_adopting_a_running_task_into_a_busy_session_is_refused(store):
    mine = _running(store, session="s2", subject="what I am on")
    theirs = _running(store, session="s1", subject="what they were on")

    with pytest.raises(ValueError) as caught:
        store.update(theirs, session_id="s2")

    assert str(mine) in str(caught.value)
    assert "in_progress" in str(caught.value)


def test_the_refusal_leaves_the_task_where_it_was(store):
    _running(store, session="s2", subject="mine")
    theirs = _running(store, session="s1", subject="theirs")
    with pytest.raises(ValueError):
        store.update(theirs, session_id="s2")
    assert store.get(theirs)["session_id"] == "s1"


def test_the_session_still_has_exactly_one_running_task(store):
    _running(store, session="s2", subject="mine")
    theirs = _running(store, session="s1", subject="theirs")
    with pytest.raises(ValueError):
        store.update(theirs, session_id="s2")
    running = [t for t in store.list(session_id="s2")
               if t["status"] == "in_progress"]
    assert len(running) == 1


# ---------------------------------------------------------------------------
# What adoption is FOR still works
# ---------------------------------------------------------------------------

def test_adopting_a_running_task_into_an_idle_session_is_allowed(store):
    """The case the tool exists for: a previous session died mid-task."""
    orphan = _running(store, session="s1", subject="left open")
    store.update(orphan, session_id="s2")
    assert store.get(orphan)["session_id"] == "s2"
    assert store.get(orphan)["status"] == "in_progress"


def test_adopting_a_pending_task_is_always_allowed(store):
    _running(store, session="s2", subject="mine")
    parked = store.create(subject="not started", session_id="s1")["id"]
    store.update(parked, session_id="s2")
    assert store.get(parked)["session_id"] == "s2"


def test_adopting_a_blocked_task_is_always_allowed(store):
    _running(store, session="s2", subject="mine")
    stuck = _running(store, session="s1", subject="waiting on a key")
    store.update(stuck, status="blocked", blocked_reason="needs a credential")
    store.update(stuck, session_id="s2")
    assert store.get(stuck)["session_id"] == "s2"


def test_a_completed_task_can_be_adopted(store):
    _running(store, session="s2", subject="mine")
    done = _running(store, session="s1", subject="finished")
    store.update(done, status="completed")
    store.update(done, session_id="s2")
    assert store.get(done)["session_id"] == "s2"


# ---------------------------------------------------------------------------
# Not a new refusal for anything that was already fine
# ---------------------------------------------------------------------------

def test_rewriting_a_task_to_its_own_session_is_not_an_adoption(store):
    """An update that mentions the session it already has must not be
    read as a takeover and refused against itself."""
    mine = _running(store, session="s2", subject="mine")
    store.update(mine, session_id="s2", subject="renamed")
    assert store.get(mine)["subject"] == "renamed"


def test_an_unscoped_store_is_unaffected(store):
    """With no session ids there is nothing to own and nothing to refuse."""
    a = _running(store, session="", subject="a")
    b = store.create(subject="b", session_id="")["id"]
    store.update(b, session_id="")
    assert store.get(a)["status"] == "in_progress"


def test_the_ordinary_second_start_is_still_refused(store):
    """The rule this widens, unchanged."""
    _running(store, session="s1", subject="first")
    second = store.create(subject="second", session_id="s1")["id"]
    with pytest.raises(ValueError):
        store.update(second, status="in_progress")
