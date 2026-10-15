"""A worktree that any process is using is never offered or removed.

Input: worktrees, and processes of this user. Output of
``worktree_use.in_use(paths)``: ``{path: why}`` for every path a process
uses -- its working directory, an open file, a mapped file or its
executable lies under it -- read from ``/proc``; on a host without
``/proc``, from ``lsof``; without either, a file modified within the last
hour counts as use. DELFIN's own job registry is not consulted: a benchmark
started from a terminal registers nowhere.

Clean up keeps such a worktree (with the process named), re-checks at
deletion time, and Close preview keeps the preview's worktree while
something else still uses it.

Universal: the process is a child this test starts itself; the open-file
and mapped-file cases use a /proc tree written by the test; lsof is a
stand-in script.
"""

from __future__ import annotations

import os
import stat
import subprocess
import sys
import time
from pathlib import Path

import pytest

from conftest import child_env
from delfin.dashboard import branch_cleanup as BC
from delfin.dashboard import worktree_use as WU

_ENV = {**os.environ, "GIT_CONFIG_GLOBAL": os.devnull,
        "GIT_CONFIG_SYSTEM": os.devnull}


def _git(cwd, *args) -> str:
    return subprocess.run(
        ["git", "-c", "user.email=t@t", "-c", "user.name=t", *args],
        cwd=str(cwd), env=_ENV, check=True, capture_output=True,
        text=True).stdout.strip()


@pytest.fixture
def repo(tmp_path, monkeypatch):
    monkeypatch.setenv("HOME", str(tmp_path / "home"))
    (tmp_path / "home").mkdir()
    origin = tmp_path / "origin.git"
    _git(tmp_path, "init", "-q", "--bare", "-b", "main", str(origin))
    co = tmp_path / "checkout"
    _git(tmp_path, "clone", "-q", str(origin), str(co))
    _git(co, "checkout", "-q", "-b", "main")
    (co / "a.txt").write_text("a\n")
    _git(co, "add", "a.txt")
    _git(co, "commit", "-q", "-m", "a")
    _git(co, "push", "-q", "-u", "origin", "main")
    return co


@pytest.fixture
def busy(repo, tmp_path):
    """A clean worktree detached at main (stale by every git rule) with a
    child process of this test working in it."""
    wt = tmp_path / "wt-bench"
    _git(repo, "worktree", "add", "-q", "--detach", str(wt), "main")
    child = subprocess.Popen([sys.executable, "-c", "import time; time.sleep(120)"],
                             cwd=str(wt), env=child_env(tmp_path))
    time.sleep(0.3)
    yield wt, child
    child.kill()
    child.wait()


def test_a_process_working_in_a_worktree_marks_it_in_use(busy):
    wt, child = busy
    use = WU.in_use([str(wt)])
    assert str(wt) in use
    assert str(child.pid) in use[str(wt)]


def test_clean_up_keeps_it_and_names_the_process(repo, busy):
    wt, child = busy
    cands, kept = BC.candidates(repo, gh=None)
    assert str(wt) not in {c.path for c in cands}
    why = {k.name: k.why for k in kept}[str(wt)]
    assert str(child.pid) in why


def test_the_rule_is_checked_again_at_deletion(repo, tmp_path):
    wt = tmp_path / "wt-later"
    _git(repo, "worktree", "add", "-q", "--detach", str(wt), "main")
    cands, _ = BC.candidates(repo, gh=None)
    mine = [c for c in cands if c.path == str(wt)]
    assert mine                                   # offered while idle
    child = subprocess.Popen([sys.executable, "-c", "import time; time.sleep(120)"],
                             cwd=str(wt), env=child_env(tmp_path))
    try:
        time.sleep(0.3)
        report = BC.delete(repo, mine)
        assert report.deleted == []
        assert "in use" in report.skipped[0][1]
        assert wt.is_dir()
    finally:
        child.kill()
        child.wait()


def _fake_proc(tmp_path, wt):
    proc = tmp_path / "proc"
    for pid, how in (("4242", "fd"), ("4343", "maps"), ("4444", "exe"),
                     ("4545", "none")):
        d = proc / pid
        (d / "fd").mkdir(parents=True)
        (d / "comm").write_text(f"worker{pid}\n")
        os.symlink("/", d / "cwd")
        os.symlink(sys.executable if how != "exe" else str(wt / "bin" / "tool"),
                   d / "exe")
        if how == "fd":
            os.symlink(str(wt / "data" / "results.db"), d / "fd" / "7")
        maps = "7f00-7f01 r-xp 0 08:01 1 /usr/lib/libc.so.6\n"
        if how == "maps":
            maps += f"7f02-7f03 r-xp 0 08:01 2 {wt}/build/libx.so\n"
        (d / "maps").write_text(maps)
    (proc / "self").mkdir()
    return proc


def test_an_open_file_a_mapped_file_and_an_executable_count(tmp_path):
    wt = tmp_path / "wt"
    wt.mkdir()
    proc = _fake_proc(tmp_path, wt)
    other = tmp_path / "other"
    other.mkdir()
    use = WU.in_use([str(wt), str(other)], proc_root=str(proc))
    assert set(use) == {str(wt)}
    why = use[str(wt)]
    assert "4242" in why and "open file" in why
    for pid in ("4242", "4343", "4444"):
        assert WU.in_use([str(wt)], proc_root=str(proc), pids=[pid]), pid
    assert not WU.in_use([str(wt)], proc_root=str(proc), pids=["4545"])


def test_a_sibling_with_a_common_prefix_is_not_inside(tmp_path):
    wt = tmp_path / "wt"
    wt.mkdir()
    proc = _fake_proc(tmp_path, tmp_path / "wt-other")
    assert WU.in_use([str(wt)], proc_root=str(proc)) == {}


def test_without_proc_lsof_answers(tmp_path):
    wt = tmp_path / "wt"
    wt.mkdir()
    lsof = tmp_path / "lsof"
    lsof.write_text(f"#!{sys.executable}\nimport sys\nprint('p777')\nprint('cbench')\n")
    lsof.chmod(lsof.stat().st_mode | stat.S_IXUSR)
    use = WU.in_use([str(wt)], proc_root=str(tmp_path / "no-proc"), lsof=str(lsof))
    assert "777" in use[str(wt)]


def test_without_proc_or_lsof_a_recent_change_counts(tmp_path):
    wt = tmp_path / "wt"
    (wt / "sub").mkdir(parents=True)
    f = wt / "sub" / "out.log"
    f.write_text("x")
    old = time.time() - 3 * 3600
    os.utime(f, (old, old))
    os.utime(wt / "sub", (old, old))
    os.utime(wt, (old, old))
    kw = dict(proc_root=str(tmp_path / "no-proc"), lsof="")
    assert WU.in_use([str(wt)], **kw) == {}
    f.write_text("y")                                   # touched now
    assert "modified" in WU.in_use([str(wt)], **kw)[str(wt)]


def test_close_preview_keeps_a_worktree_still_in_use(repo, tmp_path, monkeypatch):
    from delfin.dashboard import branch_preview as BP
    target, _ = BP.resolve_target(repo, "main")
    ok, _m, wt, _a = BP.ensure_worktree(repo, target)
    assert ok
    child = subprocess.Popen([sys.executable, "-c", "import time; time.sleep(120)"],
                             cwd=str(wt), env=child_env(tmp_path))
    try:
        time.sleep(0.3)
        p = BP.Preview(slug="main", ref="main", sha=target.sha, branch="",
                       repo=str(repo), worktree=str(wt), port=1, token="t",
                       url_path="/", pid=999999, proc_start="", ticks=None,
                       host="", started_at=time.time(), delfin_file="", log="")
        ok, msg, ask = BP._release_worktree(p, "Stopped.", "discard")
        assert Path(wt).is_dir()
        assert "in use" in msg and str(child.pid) in msg
    finally:
        child.kill()
        child.wait()
