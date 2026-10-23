"""A test result belongs to the state it ran against -- freshness controls.

Night run 2026-09-26 (assignment Y): a session reported "2 passed, 1
skipped" for a file that had long had 6 tests, and a commit message said
"live runs report nothing" while the last real run had a finding. Same
pattern both times: a test result is quoted after the code beneath it
changed. These tests pin the countermeasure -- a state fingerprint on
every evidence entry, and ``is_stale`` judging whether a quoted result
still describes the tree it ran on.

All scenarios use a throwaway git repo in tmp_path: run, then change.
"""

from __future__ import annotations

import os
import subprocess

import pytest

from delfin.agent.evidence_freshness import (
    fingerprint,
    is_stale,
    note,
    stamp,
)


def _git(repo, *args):
    # Fully self-contained git environment: no global/system config is
    # read (the gate cage has a private /tmp and may lack ~/.gitconfig),
    # dubious-ownership checks are off, and identity is fixed.
    env = {
        "PATH": os.environ.get("PATH", "/usr/bin:/bin"),
        "GIT_CONFIG_GLOBAL": "/dev/null",
        "GIT_CONFIG_SYSTEM": "/dev/null",
        "GIT_CONFIG_NOSYSTEM": "1",
        "GIT_CONFIG_COUNT": "1",
        "GIT_CONFIG_KEY_0": "safe.directory",
        "GIT_CONFIG_VALUE_0": "*",
        "GIT_AUTHOR_NAME": "t", "GIT_AUTHOR_EMAIL": "t@t",
        "GIT_COMMITTER_NAME": "t", "GIT_COMMITTER_EMAIL": "t@t",
    }
    subprocess.run(
        ["git", "-C", str(repo), *args], check=True,
        stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True,
        env=env)


def _write(repo, rel, text):
    path = repo / rel
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text)
    return rel


def _commit_all(repo, msg="c"):
    _git(repo, "add", "-A")
    _git(repo, "commit", "--allow-empty", "-m", msg)


@pytest.fixture
def repo(tmp_path):
    """A repo with a module, its test, and an unrelated stranger file."""
    _write(tmp_path, "pkg/module.py", "X = 1\n")
    _write(tmp_path, "tests/test_module.py", "def test_x():\n    assert True\n")
    _write(tmp_path, "other/stranger.py", "Y = 2\n")
    _git(tmp_path, "init", "-q")
    _commit_all(repo := tmp_path, "initial")
    return repo


def _run(repo, target="tests/test_module.py"):
    """A finished test run observed against the CURRENT tree state."""
    return stamp({"tool": "run_tests", "command": target, "exit_code": 0,
                  "status": "ok", "passed": 3, "failed": 0, "ts": 100.0},
                 repo)


def _touch(repo, rel, text=None, mtime=None):
    path = repo / rel
    if text is not None:
        path.write_text(text)
    if mtime is not None:
        os.utime(path, (mtime, mtime))


# --- fingerprint ---------------------------------------------------------

def test_fingerprint_names_commit_and_dirty_state(repo):
    fp = fingerprint(repo)
    assert fp["commit"]
    assert fp["dirty"] is not None
    _touch(repo, "pkg/module.py", "X = 2\n", mtime=200)
    assert fingerprint(repo)["dirty"] != fp["dirty"]


def test_stamp_attaches_the_fingerprint(repo):
    ev = _run(repo)
    assert ev["fingerprint"] == fingerprint(repo)
    # the observed entry itself is untouched
    assert ev["command"] == "tests/test_module.py"


# --- is_stale: the four contract scenarios -------------------------------

def test_change_of_the_tested_file_makes_it_stale(repo):
    ev = _run(repo)
    _touch(repo, "pkg/module.py", "X = 2\n", mtime=200)
    reason = is_stale(ev, fingerprint(repo))
    assert reason is not None


def test_change_of_a_stranger_file_keeps_it_fresh(repo):
    ev = _run(repo)
    _touch(repo, "other/stranger.py", "Y = 3\n", mtime=200)
    assert is_stale(ev, fingerprint(repo)) is None


def test_a_commit_that_changes_the_tested_module_makes_it_stale(repo):
    ev = _run(repo)
    _touch(repo, "pkg/module.py", "X = 9\n", mtime=300)
    _commit_all(repo, "second")
    assert is_stale(ev, fingerprint(repo)) is not None


def test_committing_what_was_tested_keeps_the_run_fresh(repo):
    # The usual flow -- change, test green, commit, mark done. The
    # session's first version judged every commit stale, which would
    # have turned exactly this into "unmet". The old decision was "any
    # commit move invalidates"; content decides now.
    _touch(repo, "pkg/module.py", "X = 5\n", mtime=50)
    ev = _run(repo)
    _commit_all(repo, "the tested change")
    assert is_stale(ev, fingerprint(repo)) is None
    _commit_all(repo, "an empty commit")
    assert is_stale(ev, fingerprint(repo)) is None


def test_a_commit_to_a_stranger_does_not_accuse(repo):
    ev = _run(repo)
    _touch(repo, "other/stranger.py", "Y = 3\n", mtime=300)
    _commit_all(repo, "elsewhere")
    assert is_stale(ev, fingerprint(repo)) is None


def test_unchanged_state_is_fresh(repo):
    ev = _run(repo)
    assert is_stale(ev, fingerprint(repo)) is None


# --- is_stale: reach of the staleness ------------------------------------

def test_a_change_before_the_run_does_not_accuse(repo):
    # dirty file at RUN time: fingerprint recorded it, so a matching
    # dirty state later is the SAME state, not staleness.
    _touch(repo, "pkg/module.py", "X = 5\n", mtime=50)
    ev = _run(repo)
    assert is_stale(ev, fingerprint(repo)) is None


def test_a_run_from_an_older_commit_is_stale_even_if_files_match(repo):
    ev = _run(repo)
    _commit_all(repo, "second")
    _touch(repo, "pkg/module.py", "X = 5\n", mtime=50)   # dirt differs too
    assert is_stale(ev, fingerprint(repo)) is not None


def test_evidence_without_fingerprint_is_judged_stale(repo):
    ev = {"tool": "run_tests", "command": "tests/test_module.py",
          "status": "ok", "passed": 3, "failed": 0, "ts": 100.0}
    # an old observation from before the fingerprint existed must not be
    # trusted as fresh -- the state it ran on is unknown
    assert is_stale(ev, fingerprint(repo)) is not None


# --- note ----------------------------------------------------------------

def test_note_names_the_state_and_the_command(repo):
    ev = _run(repo)
    _touch(repo, "pkg/module.py", "X = 2\n", mtime=200)
    reason = is_stale(ev, fingerprint(repo))
    text = note(reason)
    assert "tests/test_module.py" in text     # re-run the command
    assert ev["fingerprint"]["commit"][:7] in text


def test_note_on_a_stale_commit_names_the_new_commit(repo):
    ev = _run(repo)
    _touch(repo, "tests/test_module.py", "def test_x():\n    assert 1\n",
           mtime=300)
    _commit_all(repo, "second")
    text = note(is_stale(ev, fingerprint(repo)))
    assert "re-run" in text

