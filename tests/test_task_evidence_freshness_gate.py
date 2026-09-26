"""The completion check demands a green run for the CURRENT tree state.

Night run 2026-09-26 (assignment Y): a green count was quoted after the
code beneath it had changed. ``check_completion_claim`` now takes a
``current_fingerprint``; a green run whose own fingerprint is stale does
not verify a test task -- and the note names the state it came from and
the command to re-run. Without a fingerprint argument the check behaves
exactly as before (the characterization tests pin that).

Uses the same throwaway-repo pattern as test_evidence_freshness.py.
"""

from __future__ import annotations

import os
import subprocess

import pytest

from delfin.agent.evidence_freshness import fingerprint, stamp
from delfin.agent.task_evidence import check_completion_claim


def _git(repo, *args):
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


@pytest.fixture
def repo(tmp_path):
    (tmp_path / "pkg").mkdir()
    (tmp_path / "tests").mkdir()
    (tmp_path / "pkg" / "module.py").write_text("X = 1\n")
    (tmp_path / "tests" / "test_module.py").write_text(
        "def test_x():\n    assert True\n")
    _git(tmp_path, "init", "-q")
    _git(tmp_path, "add", "-A")
    _git(tmp_path, "commit", "--allow-empty", "-m", "initial")
    return tmp_path


def _green_run(repo, target="tests/test_module.py"):
    return stamp({"tool": "run_tests", "command": target, "exit_code": 0,
                  "status": "ok", "passed": 3, "failed": 0, "ts": 100.0},
                 repo)


SUBJECT = "Run the test suite and verify everything passes"
# "tests ... verify" matches the test-task vocabulary; deliberately no
# path token so the claim reaches the tests branch, not the path branch.


def test_green_run_for_the_current_state_verifies(repo):
    res = check_completion_claim(
        SUBJECT, tests=[_green_run(repo)], window_start=0.0,
        current_fingerprint=fingerprint(repo))
    assert res["verdict"] == "verified"
    assert res["kind"] == "tests"


def test_green_run_from_an_older_state_does_not_verify(repo):
    run = _green_run(repo)
    (repo / "pkg" / "module.py").write_text("X = 2\n")   # tested file moved
    res = check_completion_claim(
        SUBJECT, tests=[run], window_start=0.0,
        current_fingerprint=fingerprint(repo))
    assert res["verdict"] == "unmet"
    assert res["kind"] == "tests_stale_state"
    # the note names the state and asks for the re-run
    assert "re-run tests/test_module.py" in res["note"]


def test_stale_run_of_a_stranger_file_still_verifies(repo):
    run = _green_run(repo)
    stranger = repo / "other"
    stranger.mkdir()
    (stranger / "unrelated.py").write_text("Y = 2\n")
    res = check_completion_claim(
        SUBJECT, tests=[run], window_start=0.0,
        current_fingerprint=fingerprint(repo))
    assert res["verdict"] == "verified"


def test_without_a_fingerprint_argument_nothing_changes(repo):
    # Old callers pass no fingerprint: a green run verifies as before,
    # even when the tree has moved on since the run.
    run = _green_run(repo)
    (repo / "pkg" / "module.py").write_text("X = 2\n")
    res = check_completion_claim(
        SUBJECT, tests=[run], window_start=0.0)
    assert res["verdict"] == "verified"
    assert res["kind"] == "tests"


def test_unstamped_ledger_entries_are_not_trusted_as_fresh(repo):
    # Ledgers recorded before stamping existed carry no fingerprint:
    # unknown state is never quoted as fresh (safe direction).
    old = {"tool": "run_tests", "command": "tests/test_module.py",
           "exit_code": 0, "status": "ok", "passed": 3, "failed": 0,
           "ts": 100.0}
    res = check_completion_claim(
        SUBJECT, tests=[old], window_start=0.0,
        current_fingerprint=fingerprint(repo))
    assert res["verdict"] == "unmet"
    assert res["kind"] == "tests_stale_state"
