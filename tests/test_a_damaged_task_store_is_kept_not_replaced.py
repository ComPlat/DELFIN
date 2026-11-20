"""A damaged task store is quarantined, never emptied by a write.

``TaskStore`` (``delfin.agent.agent_tasks``) persists the session's task
list. ``create`` and ``update`` are the writers: each read the store,
changed one task, and wrote the whole file back. The shared read answered
"no tasks" for both an absent file and one that cannot be read, so a
store truncated the way a full disk or a kill mid-write truncates it
became, on the next task tool call, a store holding only the one task
just touched -- every other task in the plan gone, with no copy anywhere
(the old single ``.bak`` slot did not even survive a second damage).

The fix splits the two reads. Readers (``get``/``list``) keep the lenient
answer -- a damaged store must not crash the tool. Writers go through
``_load_for_write``, which on damage moves the file to the dated
``.unreadable-`` quarantine of ``agent.state_file`` (bytes kept) and
raises ``StateFileUnreadable``: a store that cannot be read cannot be
safely replaced, so the write is refused and the reason carried to the
caller (the task tools turn it into a visible tool error). The store is
gone, so the next create writes a fresh one safely.
"""

from __future__ import annotations

import pytest

from delfin.agent.agent_tasks import TaskStore
from delfin.agent.state_file import StateFileUnreadable


def _store(tmp_path) -> TaskStore:
    return TaskStore(base_dir=tmp_path)


def _path(store: TaskStore):
    return store.path


def _aside(store: TaskStore):
    return sorted(q for q in _path(store).parent.iterdir()
                  if "unreadable" in q.name)


def test_create_refuses_to_replace_a_damaged_store_and_keeps_its_bytes(tmp_path):
    store = _store(tmp_path)
    # Two real tasks, then the file truncated the way a full disk or a
    # kill mid-write truncates it.
    store.create("first task")
    store.create("second task")
    damaged = _path(store).read_text(encoding="utf-8")[:30]
    assert not damaged.endswith("}")   # genuinely torn, not valid JSON
    _path(store).write_text(damaged, encoding="utf-8")

    # Readers stay lenient while the damaged file is still present: a broken
    # store does not crash the tool, it answers "nothing is tracked" and
    # reports the read failure rather than hiding it.
    assert store.list() == []
    assert store.last_load_error, "the read failure is reported, not hidden"

    with pytest.raises(StateFileUnreadable):
        store.create("third task")

    # The write is refused: the damaged bytes are not rewritten, truncated
    # further, or silently replaced by a store that drops the other tasks.
    moved = _aside(store)
    assert len(moved) == 1, "the unreadable store must be kept, not deleted"
    assert moved[0].read_text(encoding="utf-8") == damaged, (
        "the quarantined file must still hold exactly what was in it")
    assert not _path(store).exists(), (
        "nothing may be written over a store that could not be read")

    # The store is gone, so the next create writes a fresh one safely.
    store.create("fresh task")
    assert len(store.list()) == 1
    assert len(_aside(store)) == 1, "the second create quarantines nothing"


def test_update_refuses_to_write_over_a_damaged_store(tmp_path):
    store = _store(tmp_path)
    store.create("a task")
    damaged = _path(store).read_text(encoding="utf-8")[:30]
    _path(store).write_text(damaged, encoding="utf-8")

    with pytest.raises(StateFileUnreadable):
        store.update(1, status="completed")

    assert len(_aside(store)) == 1
    assert _aside(store)[0].read_text(encoding="utf-8") == damaged
    assert not _path(store).exists()


def test_an_intact_store_is_untouched_by_any_of_this(tmp_path):
    """The control: nothing above may change the normal path."""
    store = _store(tmp_path)
    store.create("a task")
    store.create("another task")
    store.update(1, status="completed")
    tasks = store.list()
    assert len(tasks) == 2
    assert not _aside(store), "nothing is quarantined when the store parses"
    assert not store.last_load_error
