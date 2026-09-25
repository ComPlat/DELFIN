"""No suite git call may take a lock in the real checkout's index.

Measured twice on 2026-09-22 (full suite runs, 16 parallel pytest
processes, two worktrees at once): an EMPTY ``.git/worktrees/<name>/
index.lock`` was left in BOTH worktrees with the same timestamp to the
second, no git process alive, and every later ``git add`` failed with
"Unable to create ... index.lock: File exists". In the same runs the
OOM killer ended several test processes.

The one suite git call that targets the real checkout is the guard's
own ``git -C <checkout> status --porcelain`` (conftest,
``_unexpected_under_a_generated_root``): ``git status`` opportunistically
refreshes the index, which takes exactly that lock; a process killed
mid-refresh leaves it behind. Killed here means the OOM killer, which
does not give git the chance to clean up.

The tests below pin three things: the mechanism (a plain ``git status``
rewrites a stat-dirty index; ``--no-optional-locks`` does not), the
incident's effect (a stale lock blocks every later write), and the
suite's own call (it must carry the flag), plus a scan that keeps every
future git call in tests/ either inside a tmp repository or lock-free.
"""

from __future__ import annotations

import inspect
import os
import subprocess
import sys
from pathlib import Path

import pytest


def _repo(tmp_path: Path) -> Path:
    """A committed repository with one tracked file."""
    repo = tmp_path / "repo"
    repo.mkdir()
    def git(*args):
        subprocess.run(["git", "-C", str(repo), *args], check=True,
                       capture_output=True, text=True,
                       env={**os.environ, "GIT_OPTIONAL_LOCKS": "0"})
    git("init", "-q", "-b", "main")
    git("config", "user.email", "t@example.invalid")
    git("config", "user.name", "t")
    (repo / "f.txt").write_text("one\n")
    git("add", "f.txt")
    git("commit", "-qm", "initial")
    return repo


def _index(repo: Path) -> Path:
    out = subprocess.run(
        ["git", "-C", str(repo), "rev-parse", "--git-path", "index"],
        check=True, capture_output=True, text=True).stdout.strip()
    # --git-path answers relative to the repo when the repo's git dir is
    # relative, which a plain `git init` makes it.
    p = Path(out)
    return p if p.is_absolute() else repo / p


def _make_stat_dirty(repo: Path) -> None:
    """Same content, new mtime: the index wants a refresh, not a change."""
    f = repo / "f.txt"
    f.write_text("one\n")
    os.utime(f, None)


class TestTheMechanism:
    def test_plain_git_status_rewrites_a_stat_dirty_index(self, tmp_path):
        """Why the flag matters: an ordinary ``git status`` writes the
        index (under the lock) to cache the refresh. A process killed
        during that write is the incident."""
        repo = _repo(tmp_path)
        _make_stat_dirty(repo)
        before = _index(repo).stat().st_mtime_ns
        subprocess.run(["git", "-C", str(repo), "status", "--porcelain"],
                       check=True, capture_output=True)
        assert _index(repo).stat().st_mtime_ns != before

    def test_no_optional_locks_leaves_the_index_alone(self, tmp_path):
        """The fix's mechanism, measured: with the flag the index file is
        not rewritten at all, so there is no lock to leave behind."""
        repo = _repo(tmp_path)
        _make_stat_dirty(repo)
        before = _index(repo).stat().st_mtime_ns
        subprocess.run(
            ["git", "-C", str(repo), "--no-optional-locks", "status",
             "--porcelain"],
            check=True, capture_output=True)
        assert _index(repo).stat().st_mtime_ns == before


class TestTheIncidentsEffect:
    def test_a_stale_empty_lock_blocks_every_later_write(self, tmp_path):
        """The empty lock file the incident left, reproduced: whatever
        wrote it died before filling it, and git refuses to proceed."""
        repo = _repo(tmp_path)
        lock = _index(repo).with_name("index.lock")
        lock.write_bytes(b"")               # empty, as measured
        (repo / "g.txt").write_text("two\n")
        done = subprocess.run(["git", "-C", str(repo), "add", "g.txt"],
                              capture_output=True, text=True)
        assert done.returncode != 0
        assert "index.lock" in done.stderr and "exists" in done.stderr


class TestTheSuitesOwnCall:
    def test_the_checkout_scan_takes_no_lock(self):
        """The control. The one git call the suite aims at the real
        checkout must be lock-free, or a killed suite run leaves the
        checkout unusable for every later ``git add``."""
        import conftest
        src = inspect.getsource(conftest._unexpected_under_a_generated_root)
        assert "--no-optional-locks" in src, (
            "the suite runs `git status` against the real checkout "
            "without --no-optional-locks; a process killed mid-refresh "
            "leaves an index.lock behind (measured 2026-09-22, both "
            "worktrees, same second)")
