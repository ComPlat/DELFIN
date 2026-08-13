"""A cross-process lock that gives up must not give up quietly.

Two independent implementations of the same lock met in one merge, with
opposite failure behaviour: one blocked for ever on ``LOCK_EX``, the other
waited five seconds and then carried on unlocked. Blocking for ever turns
a bookkeeping file into a way to stop the agent, so the bounded wait is
the one that survived -- but carrying on unlocked is a lost update, and
the deadline expires under exactly the contention that makes one likely.

What a lost update costs here is not abstract: the registry write that
loses is an acknowledged flag or an exit code, so a live job becomes
unaddressable (``bash_status`` answering "unknown job_id" about a process
running right now) or a twelve-hour calculation announces its completion
twice. Neither is something to find out about by accident.
"""

from __future__ import annotations

import fcntl
import types

import pytest

from delfin.agent import bash_jobs as bj


@pytest.fixture(autouse=True)
def _clear_notes():
    bj.take_lock_timeouts()
    yield
    bj.take_lock_timeouts()


def test_the_lock_is_taken_when_nobody_holds_it(tmp_path):
    """The ordinary path records nothing -- a note on every call would be
    noise, and noise is how the real one gets skipped."""
    state = tmp_path / "bash_jobs.json"
    with bj.cross_process_lock(state):
        pass
    assert bj.take_lock_timeouts() == []


def test_giving_up_on_the_lock_is_recorded(tmp_path, monkeypatch):
    """Held by someone else past the deadline: the body still runs -- that
    is the deliberate part -- and the fact that it ran unlocked survives."""
    monkeypatch.setattr(bj, "_LOCK_TIMEOUT_S", 0.05)
    state = tmp_path / "bash_jobs.json"
    sidecar = tmp_path / ("bash_jobs.json" + bj._LOCK_SUFFIX)
    sidecar.parent.mkdir(parents=True, exist_ok=True)

    # A second open file description conflicts with the first even inside
    # one process: flock is per-description, not per-process.
    holder = open(sidecar, "a+")
    fcntl.flock(holder.fileno(), fcntl.LOCK_EX)
    try:
        ran = False
        with bj.cross_process_lock(state):
            ran = True
        assert ran, "the cycle must proceed; a wait that cannot end is worse"
    finally:
        fcntl.flock(holder.fileno(), fcntl.LOCK_UN)
        holder.close()

    noted = bj.take_lock_timeouts()
    assert noted, "the lock was not taken and nothing said so"
    assert str(sidecar) in noted[0]
    # Taken means taken: the next reader must not see it a second time and
    # report a degradation that already reached a turn.
    assert bj.take_lock_timeouts() == []


def test_the_note_is_bounded(tmp_path, monkeypatch):
    """A wedged holder would otherwise grow this list for as long as the
    process lives."""
    monkeypatch.setattr(bj, "_LOCK_TIMEOUT_S", 0.0)
    state = tmp_path / "bash_jobs.json"
    sidecar = tmp_path / ("bash_jobs.json" + bj._LOCK_SUFFIX)
    holder = open(sidecar, "a+")
    fcntl.flock(holder.fileno(), fcntl.LOCK_EX)
    try:
        for _ in range(bj._LOCK_TIMEOUT_CAP + 10):
            with bj.cross_process_lock(state):
                pass
    finally:
        fcntl.flock(holder.fileno(), fcntl.LOCK_UN)
        holder.close()
    assert len(bj.take_lock_timeouts()) == bj._LOCK_TIMEOUT_CAP


def test_the_degradation_reaches_the_turn(tmp_path, monkeypatch):
    """Recording it in a module list is only half a mechanism. The engine
    block that carries finished jobs into the next turn is where somebody
    reads, so that is where it has to appear -- and it has to appear even
    when no job finished, because an unlocked write is a fact about the
    registry, not about any one job."""
    from delfin.agent.engine import AgentEngine

    monkeypatch.setattr(bj, "_LOCK_TIMEOUT_S", 0.05)
    state = tmp_path / "bash_jobs.json"
    sidecar = tmp_path / ("bash_jobs.json" + bj._LOCK_SUFFIX)
    holder = open(sidecar, "a+")
    fcntl.flock(holder.fileno(), fcntl.LOCK_EX)
    try:
        with bj.cross_process_lock(state):
            pass
    finally:
        fcntl.flock(holder.fileno(), fcntl.LOCK_UN)
        holder.close()

    engine = AgentEngine.__new__(AgentEngine)
    # kit_permissions is a read-only view onto the client's permissions.
    engine.client = types.SimpleNamespace(
        _permissions=types.SimpleNamespace(workspace=str(tmp_path)))
    block = engine._build_finished_jobs_block()
    assert "WITHOUT the cross-process lock" in block
    assert "unknown job id" in block
