"""A kernel is not ended while a turn is running in it.

The culler reads windows, not activity, and for good reason: a
dashboard kernel is idle to Jupyter while its agent works in a thread.
But a window is also gone when a tab is closed by accident or a
WebSocket ping times out on a slow link, and the grace then ended a
kernel with a run still in flight -- seen on 2026-09-17, where a kernel
was ended 99 seconds after its last window while the agents were still
working.

So the kernel leaves a record while a turn runs and the culler reads
it. What is pinned here:

  a running turn holds it      the record is there -> not ended
  the grace is served after    the record goes -> the next check ends it
  the clock starts again       the grace runs after the turn, not during
  a dead kernel cannot hold    the writing process is gone -> ignored
  a stale record cannot hold   older than the freshness window -> ignored
  another machine's record     honoured on age alone; its pid is not ours
  it crosses processes         the server reads what a kernel wrote
  the tab marks both ends      a turn starts and ends the record
"""

from __future__ import annotations

import asyncio
import json
import logging
import os
import subprocess
import sys
import time
from pathlib import Path

import pytest

from delfin.dashboard import resume_server as R
from delfin.dashboard import session as S
from delfin.dashboard import turn_record as T


class _Base:
    """The parts of a MappingKernelManager the wrapper touches."""

    def __init__(self):
        self._kernel_connections = {}
        self._kernels = {}
        self.log = logging.getLogger("test")
        self.ended = []

    def __contains__(self, kid):
        return kid in self._kernels

    async def start_kernel(self, **kwargs):
        kid = kwargs.get("kernel_id") or f"k{len(self._kernels) + 1}"
        self._kernels[kid] = object()
        self._kernel_connections[kid] = 0
        return kid

    def notify_connect(self, kid):
        self._kernel_connections[kid] += 1

    def notify_disconnect(self, kid):
        self._kernel_connections[kid] -= 1

    async def shutdown_kernel(self, kid, *a, **k):
        self.ended.append(kid)
        self._kernels.pop(kid, None)
        self._kernel_connections.pop(kid, None)

    async def shutdown_all(self, *a, **k):
        for kid in list(self._kernels):
            await self.shutdown_kernel(kid)


def _run(coro):
    return asyncio.run(coro)


@pytest.fixture
def turns(tmp_path, monkeypatch):
    """The turn records, in a directory of this test's own."""
    root = tmp_path / "turns"
    monkeypatch.setattr(T, "RECORD_DIR", str(root))
    monkeypatch.delenv("DELFIN_TURN_RECORD_DIR", raising=False)
    yield root
    T._reset_for_tests()


@pytest.fixture
def manager(tmp_path, monkeypatch, turns):
    monkeypatch.setattr(S, "RECORD_DIR", str(tmp_path / "kept"))
    monkeypatch.delenv(R.GRACE_ENV, raising=False)
    cls = R.resume_kernel_manager_class(_Base)
    return cls()


def _leave(manager, kid):
    """The last window goes."""
    manager.notify_connect(kid)
    manager.notify_disconnect(kid)


def _age(manager, kid, seconds):
    """Move a kernel's unwatched-since back in time.

    The window announced its closing as it went, the way a page being
    unloaded does: that is what makes the short grace apply at all. A
    connection that merely dropped is the subject of
    test_a_closed_window_is_not_a_dropped_connection.py.
    """
    since, had = manager._delfin_unwatched()[kid]
    manager._delfin_unwatched()[kid] = (since - seconds, had)
    if had:
        manager.delfin_window_closed(kid)
        manager._delfin_closed()[kid] = since - seconds


def _unwatched(manager, kid, seconds):
    """Take the window away and put it that long ago."""
    _leave(manager, kid)
    _age(manager, kid, seconds)


def _write_record(root: Path, kid: str, **fields) -> Path:
    root.mkdir(parents=True, exist_ok=True)
    payload = {
        "kernel_id": kid,
        "pid": os.getpid(),
        "host": T._hostname(),
        "updated_at": time.time(),
    }
    payload.update(fields)
    path = root / f"{kid}.json"
    path.write_text(json.dumps(payload), encoding="utf-8")
    return path


# -- the culler -------------------------------------------------------------

def test_a_kernel_with_a_turn_running_is_not_ended(manager, turns):
    kid = _run(manager.start_kernel())
    _unwatched(manager, kid, R.GRACE_SECONDS + 30)
    T.mark(True, kid=kid)
    _run(manager.cull_kernel_if_idle(kid))
    assert manager.ended == []


def test_without_the_record_the_same_kernel_is_ended(manager, turns):
    """The control: everything the same, no turn running."""
    kid = _run(manager.start_kernel())
    _unwatched(manager, kid, R.GRACE_SECONDS + 30)
    _run(manager.cull_kernel_if_idle(kid))
    assert manager.ended == [kid]


def test_the_grace_is_served_after_the_turn_not_during_it(manager, turns):
    kid = _run(manager.start_kernel())
    _unwatched(manager, kid, R.GRACE_SECONDS + 30)
    T.mark(True, kid=kid)
    _run(manager.cull_kernel_if_idle(kid))
    assert manager.ended == []

    # The turn ends. The kernel is not ended on the spot: the clock the
    # check above restarted has to run out first, so somebody who comes
    # back right after a long run still finds the page.
    T.mark(False, kid=kid)
    _run(manager.cull_kernel_if_idle(kid))
    assert manager.ended == []

    seconds, _had = manager._delfin_seconds_unwatched(kid)
    assert seconds < R.GRACE_SECONDS

    _age(manager, kid, R.GRACE_SECONDS + 1)
    _run(manager.cull_kernel_if_idle(kid))
    assert manager.ended == [kid]


def test_the_kernel_is_kept_for_as_long_as_the_turn_lasts(manager, turns):
    """Poll after poll, a long run is not ended."""
    kid = _run(manager.start_kernel())
    T.mark(True, kid=kid)
    _leave(manager, kid)
    for _ in range(5):
        _age(manager, kid, R.GRACE_SECONDS * 10)
        _run(manager.cull_kernel_if_idle(kid))
    assert manager.ended == []


def test_a_window_that_comes_back_reports_the_next_turn_again(manager, turns,
                                                              caplog):
    """The report is made once an episode, not at every poll."""
    kid = _run(manager.start_kernel())
    T.mark(True, kid=kid)
    with caplog.at_level(logging.WARNING, logger="test"):
        _unwatched(manager, kid, R.GRACE_SECONDS + 1)
        _run(manager.cull_kernel_if_idle(kid))
        _age(manager, kid, R.GRACE_SECONDS + 1)
        _run(manager.cull_kernel_if_idle(kid))
        kept = [r for r in caplog.records if "keeping kernel" in r.getMessage()]
        assert len(kept) == 1

        # Somebody looks again, and leaves again: a new episode, and
        # worth saying once more.
        manager.notify_connect(kid)
        caplog.clear()
        manager.notify_disconnect(kid)
        _age(manager, kid, R.GRACE_SECONDS + 1)
        _run(manager.cull_kernel_if_idle(kid))
        again = [r for r in caplog.records if "keeping kernel" in r.getMessage()]
        assert len(again) == 1


def test_only_the_kernel_that_runs_a_turn_is_kept(manager, turns):
    working = _run(manager.start_kernel())
    other = _run(manager.start_kernel())
    T.mark(True, kid=working)
    for kid in (working, other):
        _unwatched(manager, kid, R.GRACE_SECONDS + 5)
        _run(manager.cull_kernel_if_idle(kid))
    assert manager.ended == [other]


# -- what a record may not do -----------------------------------------------

def test_a_record_from_a_dead_kernel_holds_nothing(manager, turns):
    gone = subprocess.Popen([sys.executable, "-c", "pass"])
    gone.wait()
    kid = _run(manager.start_kernel())
    _write_record(turns, kid, pid=gone.pid)
    assert T.running_kernel_ids() == set()

    _unwatched(manager, kid, R.GRACE_SECONDS + 5)
    _run(manager.cull_kernel_if_idle(kid))
    assert manager.ended == [kid]


def test_a_record_nobody_refreshed_holds_nothing(manager, turns):
    kid = _run(manager.start_kernel())
    _write_record(turns, kid, updated_at=time.time() - T.FRESH_SECONDS - 1)
    assert T.running_kernel_ids() == set()

    _unwatched(manager, kid, R.GRACE_SECONDS + 5)
    _run(manager.cull_kernel_if_idle(kid))
    assert manager.ended == [kid]


def test_a_record_from_another_machine_is_read_by_its_age(turns):
    """A pid means nothing across machines, so age is the whole test --
    and a fresh record from elsewhere still counts."""
    _write_record(turns, "far-away", host="some-other-node", pid=999999)
    assert "far-away" in T.running_kernel_ids()

    _write_record(turns, "long-ago", host="some-other-node", pid=999999,
                  updated_at=time.time() - T.FRESH_SECONDS - 1)
    assert "long-ago" not in T.running_kernel_ids()


def test_junk_in_the_directory_is_not_a_running_turn(turns):
    turns.mkdir(parents=True, exist_ok=True)
    (turns / "not-json.json").write_text("{", encoding="utf-8")
    (turns / "a-list.json").write_text("[1, 2]", encoding="utf-8")
    (turns / "no-id.json").write_text(json.dumps({"pid": os.getpid()}),
                                      encoding="utf-8")
    (turns / "notes.txt").write_text("hello", encoding="utf-8")
    assert T.running_kernel_ids() == set()


def test_a_missing_directory_is_simply_no_turn(tmp_path, monkeypatch):
    monkeypatch.setattr(T, "RECORD_DIR", str(tmp_path / "never-made"))
    assert T.running_kernel_ids() == set()


def test_outside_a_kernel_there_is_nothing_to_announce(turns, monkeypatch):
    monkeypatch.setattr(S, "kernel_id", lambda: "")
    assert T.mark(True) == ""
    assert T.running_kernel_ids() == set()


# -- the two ends of it -----------------------------------------------------

def test_the_server_reads_what_another_process_wrote(turns):
    """The kernel and the server are different processes; the record has
    to cross that gap or it protects nothing."""
    code = (
        "import os, sys, time\n"
        "os.environ['DELFIN_TURN_RECORD_DIR'] = sys.argv[1]\n"
        "sys.path.insert(0, sys.argv[3])\n"
        "from delfin.dashboard import turn_record as T\n"
        "T.RECORD_DIR = sys.argv[1]\n"
        "print(T.mark(True, kid=sys.argv[2]), flush=True)\n"
        "time.sleep(30)\n"
    )
    repo = str(Path(__file__).resolve().parents[1])
    kernel = subprocess.Popen(
        [sys.executable, "-c", code, str(turns), "from-a-kernel", repo],
        stdout=subprocess.PIPE, text=True)
    try:
        assert kernel.stdout.readline().strip()
        assert "from-a-kernel" in T.running_kernel_ids()
    finally:
        kernel.kill()
        kernel.wait()
    # The kernel is gone; its record no longer speaks for it.
    assert T.running_kernel_ids() == set()


def test_the_agent_tab_marks_the_start_and_the_end_of_a_turn():
    """Both ends, or the record is either useless or never released."""
    source = (Path(__file__).resolve().parents[1]
              / "delfin" / "dashboard" / "tab_agent.py").read_text(encoding="utf-8")
    assert source.count("_turns.mark(True)") == 1
    # The worker's finally, and the stop button, which bumps the
    # generation so that finally leaves the record alone.
    assert source.count("_turns.mark(False)") == 2
