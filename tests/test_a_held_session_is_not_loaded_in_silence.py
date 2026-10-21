"""A conversation held by another writer opened anyway, and saved nothing.

Three defects on one path.

The dashboard never released the writer lock. ``release_session_lock``
had exactly one caller, the CLI's ``atexit`` hook, so every dashboard
session ever run left its lock file behind. Until the lock went stale an
hour later it refused its own resume -- and on another login node, where
a holder pid names nothing checkable, by age alone.

The guard against opening a held conversation consulted only the windows
of the dashboard asking. A holder in a different process, or on another
node sharing the directory, was invisible to it, so the conversation
opened.

And then every save was refused into ``except Exception: pass``. The
session ran for the lock's full hour with each turn dropped, and the loss
only became visible on the next reload.

``session_lock_holder`` is the read-only form of the question
``acquire_session_lock`` already answers, in the same module, so there is
one answer and not a second copy that can drift.
"""

from __future__ import annotations

import inspect
import json
import os
import time

import pytest

from delfin.agent import session_store as SS


@pytest.fixture(autouse=True)
def _store(monkeypatch, tmp_path):
    monkeypatch.setattr(SS, "_SESSIONS_DIR", tmp_path / "sessions")
    (tmp_path / "sessions").mkdir()


def _write_lock(sid, *, pid, host=None, age_s=0.0):
    host = SS._this_host() if host is None else host
    (SS._SESSIONS_DIR / f"{sid}.lock").write_text(json.dumps(
        {"pid": pid, "ts": time.time() - age_s, "host": host}))


def _a_dead_pid() -> int:
    """A pid that has certainly exited, without assuming a number."""
    import subprocess
    import sys
    p = subprocess.Popen([sys.executable, "-c", "pass"])
    p.wait()
    return p.pid


# ---------------------------------------------------------------------------
# The read-only question
# ---------------------------------------------------------------------------

def test_no_lock_means_no_holder():
    assert SS.session_lock_holder("s1") is None


def test_our_own_lock_is_not_a_holder():
    """Refreshing our own lock is not a conflict, so it must not read as one."""
    _write_lock("s1", pid=os.getpid())
    assert SS.session_lock_holder("s1") is None


def test_a_live_foreign_pid_on_this_host_is_a_holder():
    # pid 1 is alive on every POSIX host and is not this process.
    _write_lock("s1", pid=1)
    held = SS.session_lock_holder("s1")
    assert held is not None
    assert held["pid"] == 1
    assert held["elsewhere"] is False


def test_a_dead_pid_is_not_a_holder():
    _write_lock("s1", pid=_a_dead_pid())
    assert SS.session_lock_holder("s1") is None


def test_an_expired_lock_is_not_a_holder():
    _write_lock("s1", pid=1, age_s=SS._LOCK_MAX_AGE_S + 60)
    assert SS.session_lock_holder("s1") is None


def test_a_fresh_lock_from_another_node_is_a_holder_whatever_its_pid():
    """A pid from another login node names nothing here, so age is the only
    judge available -- and a dead-looking number must not break the lock."""
    _write_lock("s1", pid=_a_dead_pid(), host="some-other-node")
    held = SS.session_lock_holder("s1")
    assert held is not None
    assert held["elsewhere"] is True
    assert held["host"] == "some-other-node"


def test_an_expired_lock_from_another_node_is_released():
    _write_lock("s1", pid=4242, host="some-other-node",
                age_s=SS._LOCK_MAX_AGE_S + 60)
    assert SS.session_lock_holder("s1") is None


def test_a_corrupt_lock_file_is_not_a_holder_and_does_not_raise():
    (SS._SESSIONS_DIR / "s1.lock").write_text("{not json")
    assert SS.session_lock_holder("s1") is None


def test_asking_does_not_take_the_lock():
    """The difference from acquire_session_lock: this one writes nothing."""
    _write_lock("s1", pid=1)
    before = (SS._SESSIONS_DIR / "s1.lock").read_text()
    SS.session_lock_holder("s1")
    assert (SS._SESSIONS_DIR / "s1.lock").read_text() == before


def test_it_agrees_with_acquire_on_the_same_lock():
    """One answer, not two: where the holder reads as held, acquiring must
    raise, and where it reads as free, acquiring must succeed."""
    _write_lock("held", pid=1)
    assert SS.session_lock_holder("held") is not None
    with pytest.raises(SS.SessionLockedError):
        SS.acquire_session_lock("held")

    _write_lock("free", pid=_a_dead_pid())
    assert SS.session_lock_holder("free") is None
    assert SS.acquire_session_lock("free").exists()


# ---------------------------------------------------------------------------
# Releasing
# ---------------------------------------------------------------------------

def test_a_release_drops_our_own_lock():
    SS.acquire_session_lock("s1")
    assert (SS._SESSIONS_DIR / "s1.lock").exists()
    SS.release_session_lock("s1")
    assert not (SS._SESSIONS_DIR / "s1.lock").exists()


def test_a_release_leaves_another_writers_lock_alone():
    _write_lock("s1", pid=1)
    SS.release_session_lock("s1")
    assert (SS._SESSIONS_DIR / "s1.lock").exists()


# ---------------------------------------------------------------------------
# The dashboard wiring
# ---------------------------------------------------------------------------

def _tab_source() -> str:
    from delfin.dashboard import tab_agent as T
    return inspect.getsource(T)


def test_the_dashboard_releases_the_lock_when_it_shuts_down():
    """Registered with atexit, so it covers clean close and kernel exit."""
    src = _tab_source()
    shutdown = src[src.index("def _shutdown_tab"):]
    shutdown = shutdown[:shutdown.index("_atexit.register(_shutdown_tab")]
    assert "release_session_lock" in shutdown
    assert "active_session_id" in shutdown
    assert "_atexit.register(_shutdown_tab" in src


def test_a_refused_save_is_reported_and_not_swallowed():
    src = _tab_source()
    body = src[src.index("def _auto_save_session"):]
    body = body[:body.index("def _load_saved_session")]
    i = body.index("except SessionLockedError")
    assert i < body.index("except Exception:\n            pass"), (
        "the catch-all must come after the one case that means data loss")
    arm = body[i:body.index("except Exception:\n            pass")]
    assert "_append_system_message" in arm, "a refused save must be visible"
    assert "_save_lock_warned" in arm, "once per session, not once per turn"


def test_opening_a_conversation_consults_the_lock_on_disk():
    src = _tab_source()
    body = src[src.index("def _load_saved_session"):]
    body = body[:body.index("from delfin.agent.session_store import load_session")]
    assert "session_lock_holder" in body, (
        "the in-memory window list cannot see another process or node")
    assert "_append_system_message" in body
    assert "return" in body
