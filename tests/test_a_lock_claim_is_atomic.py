"""Two writers could both believe they held the single-writer lock.

``acquire_session_lock`` read the lock file, decided, and then wrote. When
the lock did not exist yet, every concurrent caller read "no holder",
every one of them wrote, and every one of them returned holding it. The
whole point of the lock is that a second dashboard's saves must fail
loudly instead of overwriting the first's turns.

Measured on this installation, 24 and 48 processes released by one
barrier and every one kept ALIVE while the others decided (so the holder
is demonstrably live and cannot be mistaken for stale):

    main      N=24   3 of 10 runs had more than one holder
    main      N=48   3 of 8 runs, up to FOUR simultaneous holders
    patched   N=24   0 of 10
    patched   N=48   0 of 8

Two earlier instruments said otherwise and both were wrong: the first
used ``Pool.map``, which chunks the iterable so several claims ran in one
worker (same pid -- a refresh, not a race), and the second let a claimer
exit as soon as it had the lock, so the next reader found a dead pid and
broke the lock, which is correct succession rather than a defect. The
numbers here come from the third.

The claim is now a temp file ``os.link``-ed into place: atomic, and the
content is complete before the name exists. ``O_CREAT|O_EXCL`` alone
publishes an empty file and fills it after, and a reader finding a lock
with no pid treats it as stale -- which hands the lock to two callers
again by another route.

Universal: real subprocesses through ``sys.executable``, rendezvous
through files rather than a shared Barrier object, and no assumption
about the start method, the host or the home directory.
"""

from __future__ import annotations

import json
import os
import subprocess
import sys
import textwrap
import time
from pathlib import Path

import pytest

from delfin.agent import session_store as SS


@pytest.fixture(autouse=True)
def _store(monkeypatch, tmp_path):
    monkeypatch.setattr(SS, "_SESSIONS_DIR", tmp_path / "sessions")
    (tmp_path / "sessions").mkdir()


def _lock(sid="s1") -> Path:
    return SS._SESSIONS_DIR / f"{sid}.lock"


def _a_dead_pid() -> int:
    p = subprocess.Popen([sys.executable, "-c", "pass"])
    p.wait()
    return p.pid


# ---------------------------------------------------------------------------
# The primitive
# ---------------------------------------------------------------------------

def test_a_claim_succeeds_once_and_then_reports_the_file(tmp_path):
    p = tmp_path / "x.lock"
    assert SS._claim_lock_file(p, '{"pid": 1}') is True
    assert SS._claim_lock_file(p, '{"pid": 2}') is False
    assert json.loads(p.read_text())["pid"] == 1, "a claim overwrote a holder"


def test_a_claimed_file_is_never_seen_empty(tmp_path):
    """The reason it is a link and not an exclusive create: a lock with no
    pid in it reads as stale, and a stale lock gets broken."""
    p = tmp_path / "x.lock"
    SS._claim_lock_file(p, '{"pid": 4242, "ts": 1.0, "host": "h"}')
    assert SS._read_lock(p) == (4242, 1.0, "h")


def test_a_failed_claim_leaves_no_litter(tmp_path):
    """The temp file a refused claim wrote must not stay in the sessions
    directory -- it is scanned for locks."""
    p = tmp_path / "x.lock"
    SS._claim_lock_file(p, '{"pid": 1}')
    SS._claim_lock_file(p, '{"pid": 2}')
    leftovers = [q.name for q in tmp_path.iterdir()
                 if q.name.endswith(".claim") or ".claim." in q.name]
    assert leftovers == [], leftovers
    assert p.is_file()


def test_threads_claiming_one_path_produce_one_winner(tmp_path):
    """Thread-level check of the primitive: the function must not hand True
    to two callers for the same path."""
    import threading

    p = tmp_path / "x.lock"
    gate = threading.Barrier(16)
    wins: list[bool] = []
    lock = threading.Lock()

    def _go():
        gate.wait()
        got = SS._claim_lock_file(p, '{"pid": 1}')
        with lock:
            wins.append(got)

    threads = [threading.Thread(target=_go) for _ in range(16)]
    for t in threads:
        t.start()
    for t in threads:
        t.join(30)
    assert sum(1 for w in wins if w) == 1, wins


# ---------------------------------------------------------------------------
# Absent and unreadable are different answers
# ---------------------------------------------------------------------------

def test_an_absent_lock_reads_as_none(tmp_path):
    assert SS._read_lock(tmp_path / "nothing.lock") is None


def test_a_corrupt_lock_reads_as_present_but_unowned(tmp_path):
    p = tmp_path / "x.lock"
    p.write_text("{not json")
    assert SS._read_lock(p) == (0, 0.0, "")


def test_an_empty_lock_reads_as_present_but_unowned(tmp_path):
    p = tmp_path / "x.lock"
    p.write_text("")
    assert SS._read_lock(p) == (0, 0.0, "")


def test_a_lock_that_vanishes_is_not_broken(monkeypatch):
    """The residual race in the first version of the fix.

    A loser whose read came back empty treated the lock as stale and
    unlinked it -- so it took a lock a live holder was already holding.
    "I could not read it" is evidence of another claimer, not of
    staleness, and must cost a retry rather than an unlink.
    """
    _lock().write_text(json.dumps(
        {"pid": 1, "ts": time.time(), "host": SS._this_host()}))
    unlinked: list[str] = []
    real_unlink = os.unlink
    monkeypatch.setattr(
        os, "unlink",
        lambda p, *a, **k: (unlinked.append(str(p)), real_unlink(p))[1]
        if str(p).endswith(".lock") else real_unlink(p, *a, **k))
    # The file is there for the claim and gone for the read, once.
    calls = {"n": 0}
    real_read = SS._read_lock

    def flaky(path):
        calls["n"] += 1
        return None if calls["n"] == 1 else real_read(path)

    monkeypatch.setattr(SS, "_read_lock", flaky)
    with pytest.raises(SS.SessionLockedError):
        SS.acquire_session_lock("s1")
    assert not unlinked, "a lock was broken because a read came back empty"
    assert _lock().is_file()


# ---------------------------------------------------------------------------
# The decisions the lock already made, unchanged
# ---------------------------------------------------------------------------

def test_a_free_id_is_acquired():
    assert SS.acquire_session_lock("s1").is_file()
    assert SS._read_lock(_lock())[0] == os.getpid()


def test_our_own_lock_is_refreshed_not_refused():
    SS.acquire_session_lock("s1")
    before = SS._read_lock(_lock())[1]
    time.sleep(0.01)
    SS.acquire_session_lock("s1")
    assert SS._read_lock(_lock())[1] > before


def test_a_live_foreign_holder_is_refused():
    _lock().write_text(json.dumps(
        {"pid": 1, "ts": time.time(), "host": SS._this_host()}))
    with pytest.raises(SS.SessionLockedError):
        SS.acquire_session_lock("s1")


def test_a_dead_holder_is_broken():
    _lock().write_text(json.dumps(
        {"pid": _a_dead_pid(), "ts": time.time(), "host": SS._this_host()}))
    assert SS.acquire_session_lock("s1").is_file()
    assert SS._read_lock(_lock())[0] == os.getpid()


def test_an_expired_lock_is_broken():
    _lock().write_text(json.dumps(
        {"pid": 1, "ts": time.time() - SS._LOCK_MAX_AGE_S - 60,
         "host": SS._this_host()}))
    assert SS.acquire_session_lock("s1").is_file()


def test_a_fresh_lock_from_another_node_holds_whatever_its_pid():
    _lock().write_text(json.dumps(
        {"pid": _a_dead_pid(), "ts": time.time(), "host": "some-other-node"}))
    with pytest.raises(SS.SessionLockedError):
        SS.acquire_session_lock("s1")


def test_a_corrupt_lock_is_broken():
    _lock().write_text("{not json")
    assert SS.acquire_session_lock("s1").is_file()


def test_the_read_only_holder_question_uses_the_same_reader():
    """One answer: session_lock_holder must agree with acquire."""
    _lock().write_text(json.dumps(
        {"pid": 1, "ts": time.time(), "host": SS._this_host()}))
    assert SS.session_lock_holder("s1") is not None
    with pytest.raises(SS.SessionLockedError):
        SS.acquire_session_lock("s1")


# ---------------------------------------------------------------------------
# The whole function, under real concurrency
# ---------------------------------------------------------------------------

_CHILD = textwrap.dedent('''
    import os, sys, time, pathlib
    repo, sessions, sid, n = sys.argv[1], sys.argv[2], sys.argv[3], int(sys.argv[4])
    sys.path.insert(0, repo)
    from delfin.agent import session_store as SS
    SS._SESSIONS_DIR = pathlib.Path(sessions)
    base = pathlib.Path(sessions).parent
    ready, out = base / "ready", base / "out"
    (ready / str(os.getpid())).write_text("1")
    deadline = time.time() + 60
    while len(list(ready.iterdir())) < n and time.time() < deadline:
        pass
    try:
        SS.acquire_session_lock(sid)
        verdict = "held"
    except SS.SessionLockedError:
        verdict = "refused"
    except Exception as exc:
        verdict = "error:" + type(exc).__name__
    (out / str(os.getpid())).write_text(verdict)
    # Stay alive until every sibling has decided, so a holder is never
    # mistaken for a dead one.
    while len(list(out.iterdir())) < n and time.time() < deadline:
        pass
''')


def test_only_one_of_many_processes_holds_a_fresh_lock(tmp_path):
    """The defect itself: separate processes, one fresh id, one winner.

    Separate processes because the holder is identified by pid -- two
    threads of one process are both "our own pid", which is a refresh and
    not a race. Each child stays alive until all have decided, so the
    winner cannot be read as a dead holder.
    """
    n = 8
    repo = str(Path(SS.__file__).resolve().parents[2])
    script = tmp_path / "claim.py"
    script.write_text(_CHILD)
    (tmp_path / "ready").mkdir()
    (tmp_path / "out").mkdir()
    procs = [subprocess.Popen(
        [sys.executable, str(script), repo, str(SS._SESSIONS_DIR),
         "race-1", str(n)]) for _ in range(n)]
    for p in procs:
        p.wait(timeout=120)
    verdicts = [q.read_text() for q in (tmp_path / "out").iterdir()]
    assert len(verdicts) == n, verdicts
    assert not [v for v in verdicts if v.startswith("error")], verdicts
    assert verdicts.count("held") == 1, (
        f"{verdicts.count('held')} processes believe they hold the lock: "
        f"{verdicts}")
