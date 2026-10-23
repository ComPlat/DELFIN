"""The cross-process lock holds where flock does not reach.

Login nodes share the home directory, and flock stays on one node on NFS
mounted nolock, Lustre with localflock and BeeGFS by default (audit
2026-09-16). Two nodes then both load, change and write the job registry,
the memory index or a session inbox, and one node's update is gone -- in
the registry, a live job nobody can address. The lock now adds a lease
file created with O_EXCL. These tests switch flock off entirely, the way
such a filesystem behaves across nodes, and count lost updates.
"""
import json
import multiprocessing
import os
import time

import pytest

from delfin.agent import bash_jobs as BJ


def _worker(path, rounds, use_lease, flock_works):
    import fcntl as _f
    if not flock_works:
        BJ.fcntl.flock = lambda *a, **k: None
    if not use_lease:
        BJ._take_lease = lambda p, d: None
        BJ._note_lock_timeout = lambda p: None
    for _ in range(rounds):
        with BJ.cross_process_lock(path):
            data = json.loads(path.read_text())
            n = data["n"]
            time.sleep(0.0005)
            data["n"] = n + 1
            path.write_text(json.dumps(data))


def _count(tmp_path, *, use_lease, procs=6, rounds=40):
    path = tmp_path / "state.json"
    path.write_text(json.dumps({"n": 0}))
    ctx = multiprocessing.get_context("fork")
    ps = [ctx.Process(target=_worker, args=(path, rounds, use_lease, False)) for _ in range(procs)]
    for p in ps:
        p.start()
    for p in ps:
        p.join(120)
    return json.loads(path.read_text())["n"], procs * rounds


def test_without_a_working_flock_the_lease_loses_no_update(tmp_path):
    got, want = _count(tmp_path, use_lease=True)
    assert got == want


def _meet(path, barrier, use_lease):
    BJ.fcntl.flock = lambda *a, **k: None
    if not use_lease:
        BJ._take_lease = lambda p, d: None
        BJ._note_lock_timeout = lambda p: None
    with BJ.cross_process_lock(path):
        n = json.loads(path.read_text())["n"]
        try:
            barrier.wait(timeout=3)     # both inside at once, if the lock lets them
        except Exception:
            pass
        path.write_text(json.dumps({"n": n + 1}))


def _meeting(tmp_path, use_lease):
    path = tmp_path / "state.json"
    path.write_text(json.dumps({"n": 0}))
    ctx = multiprocessing.get_context("fork")
    barrier = ctx.Barrier(2)
    ps = [ctx.Process(target=_meet, args=(path, barrier, use_lease)) for _ in range(2)]
    for p in ps:
        p.start()
    for p in ps:
        p.join(60)
    return json.loads(path.read_text())["n"]


def test_the_control_without_the_lease_loses_an_update(tmp_path):
    """The old lock with flock gone: two writers meet inside and one update
    is lost -- every time, not by the luck of the scheduler."""
    assert _meeting(tmp_path, use_lease=False) == 1


def test_with_the_lease_the_second_writer_waits(tmp_path):
    assert _meeting(tmp_path, use_lease=True) == 2


def test_a_lease_left_by_a_dead_process_here_is_broken_at_once(tmp_path):
    target = tmp_path / "registry.json"
    lease = tmp_path / "registry.json.lease"
    here = BJ.this_machine()
    lease.write_text(json.dumps({"host": here["host"], "boot": here["boot_id"],
                                 "pid": 999999999, "ts": time.time(), "nonce": "old"}))
    t0 = time.monotonic()
    with BJ.cross_process_lock(target):
        assert json.loads(lease.read_text())["pid"] == os.getpid()
    assert time.monotonic() - t0 < 2
    assert not lease.exists()


def test_a_fresh_lease_of_another_node_is_waited_for(tmp_path, monkeypatch):
    monkeypatch.setattr(BJ, "_LOCK_TIMEOUT_S", 0.3)
    notes = []
    monkeypatch.setattr(BJ, "_note_lock_timeout", lambda p: notes.append(p))
    target = tmp_path / "registry.json"
    lease = tmp_path / "registry.json.lease"
    lease.write_text(json.dumps({"host": "another-login-node", "boot": "b",
                                 "pid": os.getpid(), "ts": time.time(), "nonce": "theirs"}))
    t0 = time.monotonic()
    with BJ.cross_process_lock(target):
        pass
    assert time.monotonic() - t0 >= 0.3 and notes
    assert json.loads(lease.read_text())["nonce"] == "theirs"      # not ours to remove


def test_an_expired_lease_of_another_node_is_broken(tmp_path):
    target = tmp_path / "registry.json"
    lease = tmp_path / "registry.json.lease"
    lease.write_text(json.dumps({"host": "another-login-node", "boot": "b", "pid": 1,
                                 "ts": time.time() - 600, "nonce": "old"}))
    past = time.time() - 600
    os.utime(lease, (past, past))
    with BJ.cross_process_lock(target):
        assert json.loads(lease.read_text())["pid"] == os.getpid()


def test_the_session_inbox_takes_the_lock(tmp_path, monkeypatch):
    from delfin.agent import session_messages as SM
    monkeypatch.setattr(SM, "_DIR", tmp_path)
    seen = []
    real = BJ.cross_process_lock

    import contextlib

    @contextlib.contextmanager
    def spy(path):
        seen.append(str(path))
        with real(path):
            yield

    monkeypatch.setattr(BJ, "cross_process_lock", spy)
    SM.send("abc", "hello", from_key="me")
    assert SM.take("abc")[0]["text"] == "hello"
    assert len(seen) == 2
