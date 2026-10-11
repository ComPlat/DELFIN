"""U8 review: git/patch write targets the first version did not name.

Each case is a command whose write the reader returned no target for (or a
target in the wrong tree), or a reader input that blocked or refused
ordinary work. Pure: ``patch_targets.targets`` on command text and patch
files under tmp_path.
"""

import os
import threading

import pytest

from delfin.agent import patch_targets as pt


@pytest.fixture
def ws(tmp_path):
    (tmp_path / "sub").mkdir()
    (tmp_path / "f.txt").write_text("a\n")
    (tmp_path / "a.diff").write_text(
        "diff --git a/a.py b/a.py\n--- a/a.py\n+++ b/a.py\n"
        "@@ -1 +1 @@\n-x\n+y\n")
    return tmp_path


# --- patch contents ----------------------------------------------------

def test_a_traditional_delete_names_the_old_file(ws):
    (ws / "d.diff").write_text(
        "--- a/victim.txt\n+++ /dev/null\n@@ -1 +0,0 @@\n-a\n")
    assert str(ws / "victim.txt") in pt.targets("git apply d.diff", ws)
    assert str(ws / "victim.txt") in pt.targets("patch -p1 -i d.diff", ws)


def test_a_context_diff_names_its_files(ws):
    (ws / "c.diff").write_text(
        "*** victim.txt\n--- victim.txt\n***************\n"
        "*** 1 ****\n! a\n--- 1 ----\n! b\n")
    assert pt.targets("patch -p0 -i c.diff", ws) == [str(ws / "victim.txt")]


def test_gnu_patch_names_the_old_name_too(ws):
    (ws / "o.diff").write_text(
        "--- victim.txt\n+++ other.txt\n@@ -1 +1 @@\n-a\n+b\n")
    got = pt.targets("patch -p0 -i o.diff", ws)
    assert str(ws / "victim.txt") in got and str(ws / "other.txt") in got


def test_a_rename_names_the_rename_line(ws):
    (ws / "r.diff").write_text(
        "diff --git a/x b/y\nsimilarity index 100%\n"
        "rename from x\nrename to elsewhere/z\n")
    assert str(ws / "elsewhere/z") in pt.targets("git apply r.diff", ws)


def test_a_fifo_patch_does_not_block_the_reader(ws):
    os.mkfifo(ws / "p.diff")
    box = []
    t = threading.Thread(
        target=lambda: box.append(pt.targets("git apply p.diff", ws)),
        daemon=True)
    t.start()
    t.join(5)
    assert box and box[0] is pt.UNKNOWABLE


def test_an_oversized_patch_is_unknowable(ws, monkeypatch):
    monkeypatch.setattr(pt, "_MAX_PATCH_BYTES", 10)
    assert pt.targets("git apply a.diff", ws) is pt.UNKNOWABLE


# --- command line ------------------------------------------------------

@pytest.mark.parametrize("cmd", [
    "env git apply a.diff", "command git apply a.diff",
    "nice -n 5 git apply a.diff", "timeout 5 git apply a.diff",
    "nohup git apply a.diff", "env -i FOO=1 git apply a.diff",
    "sudo git apply a.diff", "time git apply a.diff",
    "exec git apply a.diff", "stdbuf -oL git apply a.diff",
])
def test_a_command_prefix_does_not_hide_git(ws, cmd):
    assert pt.targets(cmd, ws) == [str(ws / "a.py")], cmd


def test_xargs_feeding_a_patch_writer_is_unknowable(ws):
    assert pt.targets("xargs git apply < list", ws) is pt.UNKNOWABLE
    assert pt.targets("ls | xargs git checkout --", ws) is pt.UNKNOWABLE


@pytest.mark.parametrize("cmd", [
    "git apply a.diff 2>&1", "git apply a.diff 2>/dev/null",
    "git apply a.diff > log.txt", "git apply < a.diff",
    "git apply - < a.diff",
])
def test_a_redirection_is_not_a_patch_file(ws, cmd):
    assert pt.targets(cmd, ws) == [str(ws / "a.py")], cmd


def test_patch_reads_its_patch_from_a_stdin_redirect(ws):
    assert pt.targets("patch -p1 < a.diff", ws) == [str(ws / "a.py")]
    assert pt.targets("patch -p1 -i a.diff 2>&1", ws) == [str(ws / "a.py")]


@pytest.mark.parametrize("cmd,want", [
    ("git checkout -f /x/nx", "/x/nx"),
    ("git checkout --theirs /x/nx", "/x/nx"),
    ("git checkout -p /x/nx", "/x/nx"),
    ("git checkout -q /x/nx", "/x/nx"),
    ("git checkout HEAD~1 newfile", "newfile"),
    ("git restore --worktree /x/nx", "/x/nx"),
    ("git restore -q /x/nx", "/x/nx"),
    ("git restore -p /x/nx", "/x/nx"),
    ("git restore --staged --worktree /x/nx", "/x/nx"),
])
def test_a_checkout_or_restore_flag_does_not_swallow_the_path(ws, cmd, want):
    assert str(ws / want) in pt.targets(cmd, ws), cmd


@pytest.mark.parametrize("cmd", [
    "git restore --staged a b", "git restore -S a",
    "git restore --source=HEAD --staged a",
])
def test_a_staged_only_restore_writes_no_file(ws, cmd):
    assert pt.targets(cmd, ws) == [], cmd


@pytest.mark.parametrize("cmd", [
    "GIT_WORK_TREE=/x git checkout -- f",
    "git -c core.worktree=/x checkout -- f",
    "git --work-tree /x checkout -- f",
    "git -C sub -C /x checkout -- f",
])
def test_the_tree_git_is_pointed_at_is_where_paths_land(ws, cmd):
    assert pt.targets(cmd, ws) == ["/x/f"], cmd


def test_repeated_dash_c_composes(ws):
    assert pt.targets("git -C sub -C deeper checkout -- f", ws) == [
        str(ws / "sub/deeper/f")]


@pytest.mark.parametrize("cmd,want", [
    ("patch --strip=1 --input=a.diff", "a.py"),
    ("patch -p1 -i a.diff -r /x/r.rej", "/x/r.rej"),
    ("patch -p1 -i a.diff --reject-file=/x/r.rej", "/x/r.rej"),
    ("patch -p1 -i a.diff -B /x/bk/", "/x/bk/a.py"),
    ("patch -p1 -i a.diff --prefix=/x/bk/", "/x/bk/a.py"),
])
def test_gnu_patch_long_options_and_side_files(ws, cmd, want):
    assert str(ws / want) in pt.targets(cmd, ws), cmd


def test_gnu_patch_directory_long_form_reads_there(ws):
    (ws / "sub" / "s.diff").write_text(
        "--- a/s.py\n+++ b/s.py\n@@ -1 +1 @@\n-x\n+y\n")
    assert pt.targets("patch --directory=sub -p1 -i s.diff", ws) == [
        str(ws / "sub/s.py")]


@pytest.mark.parametrize("cmd", [
    "GIT_WORK_TREE=/x git reset --hard",
    "git --work-tree=/x stash pop",
    "git -c core.worktree=/x merge feature",
    "git --work-tree /x clean -fd",
])
def test_a_tree_moved_away_is_named_for_a_tree_wide_write(ws, cmd):
    assert pt.targets(cmd, ws) == ["/x"], cmd


@pytest.mark.parametrize("cmd", [
    "GIT_WORK_TREE=/x git status", "git --work-tree=/x log",
    "git reset --hard", "git stash pop",
])
def test_a_read_or_an_unmoved_tree_names_nothing(ws, cmd):
    assert pt.targets(cmd, ws) == [], cmd
