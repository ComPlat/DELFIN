"""Closing a preview always ends in a final state, within a bounded time.

Input: a preview record (``branch_preview.Preview``) for a worktree under
``<parent>/delfin-preview/`` of a checkout made here, in four conditions
seen on a cluster home on a parallel filesystem: the server is already
gone (a stale record), ``worktree_use.in_use`` raises, ``git worktree
remove`` hangs (a ``git`` on PATH that sleeps past the step timeout), and
the record names another host. Output of :func:`branch_preview.close`: a
:class:`CloseResult` whose ``state`` is final (``closed``, ``kept``,
``ask`` or ``failed``), whose ``text`` is plain (no port, token, pid or
absolute path), and whose ``details`` carry the commands that finish the
job by hand when a step failed. ``progress`` receives the step lines.

Rejected: a close that waits for whatever a step does (the field report:
"Stopping the preview ..." never replaced), and killing a process group
named in a record written on another host (the pid means nothing here).
"""

from __future__ import annotations

import os
import socket
import subprocess
import sys
import time
from pathlib import Path

import pytest

from delfin.dashboard import branch_preview as BP

_ENV = {**os.environ, "GIT_CONFIG_GLOBAL": os.devnull,
        "GIT_CONFIG_SYSTEM": os.devnull}


def _git(cwd, *args) -> str:
    return subprocess.run(
        ["git", "-c", "user.email=t@t", "-c", "user.name=t", *args],
        cwd=str(cwd), env=_ENV, check=True, capture_output=True,
        text=True).stdout.strip()


@pytest.fixture
def rec(tmp_path, monkeypatch):
    monkeypatch.setenv("HOME", str(tmp_path / "home"))
    (tmp_path / "home").mkdir()
    co = tmp_path / "repo" / "DELFIN"
    co.mkdir(parents=True)
    _git(co, "init", "-q", "-b", "main")
    (co / "a.txt").write_text("a\n")
    _git(co, "add", "a.txt")
    _git(co, "commit", "-q", "-m", "a")
    _git(co, "branch", "feat")
    wt = BP.worktree_path(co, "feat")
    _git(co, "worktree", "add", "-q", str(wt), "feat")
    dead = subprocess.Popen([sys.executable, "-c", "pass"])
    dead.wait()
    p = BP.Preview(
        slug="feat", ref="feat", sha=_git(co, "rev-parse", "feat"),
        branch="feat", repo=str(co), worktree=str(wt), port=45999,
        token="tok-secret-1", url_path="/voila/render/x.ipynb", pid=dead.pid,
        proc_start="gone", ticks=None, host=socket.gethostname(),
        started_at=time.time(), delfin_file="", log="")
    BP._save(p)
    return co, p


def _plain(text: str, p) -> None:
    assert str(p.port) not in text and p.token not in text
    assert str(p.pid) not in text and p.worktree not in text


def test_a_stale_record_closes_at_once(rec):
    co, p = rec
    lines = []
    t0 = time.monotonic()
    r = BP.close(p.slug, progress=lines.append)
    assert time.monotonic() - t0 < 5
    assert r.ok and r.state == "closed", r
    assert not BP._record_path(p.slug).exists()
    assert not Path(p.worktree).exists()
    assert lines and all(isinstance(x, str) for x in lines)
    _plain(r.text, p)


def test_in_use_raising_ends_kept_with_the_manual_command(rec, monkeypatch):
    co, p = rec
    from delfin.dashboard import worktree_use

    def boom(*_a, **_k):
        raise OSError("stale file handle")
    monkeypatch.setattr(worktree_use, "in_use", boom)
    r = BP.close(p.slug)
    assert r.state in ("kept", "failed") and r.state != "", r
    assert "stale file handle" in r.text + r.details
    assert "worktree remove" in r.details
    assert not BP._record_path(p.slug).exists()     # the server part is done
    assert Path(p.worktree).exists()                # nothing removed blindly
    _plain(r.text, p)


def test_a_hanging_git_worktree_remove_is_bounded(rec, monkeypatch, tmp_path):
    co, p = rec
    real = subprocess.run(["sh", "-c", "command -v git"], capture_output=True,
                          text=True).stdout.strip()
    bindir = tmp_path / "bin"
    bindir.mkdir()
    fake = bindir / "git"
    fake.write_text("#!/bin/sh\n"
                    'case " $* " in *" worktree remove "*) sleep 30;; esac\n'
                    f'exec {real} "$@"\n')
    fake.chmod(0o755)
    monkeypatch.setenv("PATH", f"{bindir}{os.pathsep}{os.environ['PATH']}")
    t0 = time.monotonic()
    r = BP.close(p.slug, step_timeout=2)
    assert time.monotonic() - t0 < 15
    assert r.state == "failed" and not r.ok, r
    assert "worktree remove" in r.details
    _plain(r.text, p)


def test_a_record_from_another_host_is_not_killed(rec, monkeypatch):
    co, p = rec
    import dataclasses
    other = dataclasses.replace(p, host="elsewhere-node-7", pid=os.getpid())
    BP._save(other)
    killed = []
    monkeypatch.setattr(BP, "_kill_group", lambda *a, **k: killed.append(a))
    r = BP.close(p.slug)
    assert killed == []
    assert r.state == "failed" and "elsewhere-node-7" in r.text
    assert BP._record_path(p.slug).exists()


def test_stop_keeps_its_tuple_and_never_raises(rec, monkeypatch):
    co, p = rec
    monkeypatch.setattr(BP, "changes", lambda _p: (_ for _ in ()).throw(
        RuntimeError("git status timed out")))
    ok, msg, ask = BP.stop(p.slug)
    assert isinstance(ok, bool) and isinstance(msg, str) and ask is False
    assert "git status timed out" in msg


def test_a_worktree_still_in_use_right_after_the_stop_is_waited_for(rec, monkeypatch):
    """The server's kernels end a moment after the server (field: a close
    from the preview's own tab left the folder behind, and nobody was told)."""
    co, p = rec
    from delfin.dashboard import worktree_use
    calls = []

    def busy_twice(paths, **_k):
        calls.append(1)
        return {list(paths)[0]: "a kernel"} if len(calls) <= 2 else {}
    monkeypatch.setattr(worktree_use, "in_use", busy_twice)
    monkeypatch.setattr(BP, "_alive", lambda _p: True)
    monkeypatch.setattr(BP, "_kill_group", lambda *_a, **_k: None)
    monkeypatch.setattr(BP, "_group_alive", lambda _pid: False)
    r = BP.close(p.slug, settle_s=5)
    assert r.state == "closed", r
    assert len(calls) == 3
    assert not Path(p.worktree).exists()


def test_a_worktree_that_stays_in_use_is_left_as_a_question(rec, monkeypatch):
    co, p = rec
    from delfin.dashboard import worktree_use
    monkeypatch.setattr(worktree_use, "in_use",
                        lambda paths, **_k: {list(paths)[0]: "a kernel"})
    r = BP.close(p.slug, settle_s=1)
    assert r.state == "kept" and Path(p.worktree).exists()
    assert [x.slug for x in BP.pending_closes(co)] == [p.slug]
