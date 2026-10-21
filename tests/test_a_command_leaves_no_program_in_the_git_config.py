"""A shell command leaves no program in the repository's git config.

Where the sandbox cannot hold .git/config read-only (Landlock only, no
layer, a layer that failed at run time), the config is compared before and
after each shell command and background job. A program-naming entry the
command added is undone; any other config change stays. Nothing here runs
the program: the entries are only written and compared.
"""
from __future__ import annotations

import os
import subprocess

import pytest

from delfin.agent import api_client as A
from delfin.agent import git_config_watch as W


@pytest.fixture
def repo(tmp_path, monkeypatch):
    monkeypatch.setattr(A, "_record_security_event", lambda *a, **k: None)
    A._GIT_ROOTS_CACHE.clear()
    r = tmp_path / "repo"
    r.mkdir()
    env = {**os.environ, "GIT_CONFIG_GLOBAL": "/dev/null"}
    subprocess.run(["git", "init", "-q"], cwd=r, env=env, check=True)
    return r


def _append(repo, text):
    cfg = repo / ".git" / "config"
    cfg.write_text(cfg.read_text() + text)
    return cfg


def test_a_program_entry_the_command_added_is_undone(repo):
    snap = W.take(repo)
    before = (repo / ".git" / "config").read_bytes()
    _append(repo, "[core]\n\tpager = cat\n")
    notes = W.check(snap)
    assert (repo / ".git" / "config").read_bytes() == before
    assert notes and "core.pager" in notes[0]


@pytest.mark.parametrize("entry", [
    "[filter \"x\"]\n\tclean = cat\n",
    "[include]\n\tpath = elsewhere.cfg\n",
    "[alias]\n\tst = !cat\n",
    "[diff \"pic\"]\n\ttextconv = cat\n",
    "[core]\n\tsshCommand = ssh\n",
])
def test_every_program_naming_key_counts(repo, entry):
    snap = W.take(repo)
    before = (repo / ".git" / "config").read_bytes()
    _append(repo, entry)
    assert W.check(snap)
    assert (repo / ".git" / "config").read_bytes() == before


def test_an_ordinary_config_change_stays(repo):
    snap = W.take(repo)
    cfg = _append(repo, "[user]\n\tname = x\n[remote \"origin\"]\n\turl = /tmp/x\n")
    after = cfg.read_bytes()
    assert W.check(snap) == []
    assert cfg.read_bytes() == after


def test_a_value_change_of_an_existing_program_key_counts(repo):
    _append(repo, "[core]\n\tpager = cat\n")
    snap = W.take(repo)
    before = (repo / ".git" / "config").read_bytes()
    cfg = repo / ".git" / "config"
    cfg.write_text(cfg.read_text().replace("pager = cat", "pager = less"))
    assert W.check(snap)
    assert cfg.read_bytes() == before


def test_the_restore_keeps_the_file_mode(repo):
    cfg = repo / ".git" / "config"
    os.chmod(cfg, 0o600)
    snap = W.take(repo)
    _append(repo, "[core]\n\tpager = cat\n")
    W.check(snap)
    assert (cfg.stat().st_mode & 0o777) == 0o600


def test_the_shell_journal_runs_the_watch(repo):
    """The executor's snapshot and journal carry the watch, so foreground
    commands, timeouts and background jobs (on_finish) all get it."""
    perms = A.KitToolPermissions(workspace=str(repo))
    ex = A._doc_executor
    snap = ex._shell_snapshot(perms)
    assert getattr(snap, "git_watch", None) is not None
    before = (repo / ".git" / "config").read_bytes()
    _append(repo, "[core]\n\tfsmonitor = cat\n")
    notes = ex._journal_shell_changes("bash", snap, perms)
    assert any("core.fsmonitor" in n for n in notes), notes
    assert (repo / ".git" / "config").read_bytes() == before


def test_a_snapshot_with_only_the_watch_still_checks(repo, monkeypatch):
    perms = A.KitToolPermissions(workspace=str(repo))
    ex = A._doc_executor
    monkeypatch.setattr(type(ex), "_shell_snapshot_files", lambda self, p: None)
    snap = ex._shell_snapshot(perms)
    assert isinstance(snap, A._GitWatchOnly)
    _append(repo, "[core]\n\tpager = cat\n")
    assert ex._journal_shell_changes("bash", snap, perms)
