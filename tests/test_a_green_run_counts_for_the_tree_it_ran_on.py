"""A green test run counts for the tree it ran on -- and not after that tree moved.

The wiring of evidence_freshness: every test run the agent is seen doing
is stamped with the work tree's state, and the completion check judges
the stamped runs against the state at the moment a test task is marked
done. Committing what was tested does not move the state; editing the
tested module afterwards does.
"""

from __future__ import annotations

import os
import subprocess

from delfin.agent import api_client
from delfin.agent.task_evidence import check_completion_claim


def _git(repo, *args):
    subprocess.run(["git", "-C", str(repo), *args], check=True,
                   capture_output=True,
                   env={**os.environ, "GIT_AUTHOR_NAME": "t",
                        "GIT_AUTHOR_EMAIL": "t@t", "GIT_COMMITTER_NAME": "t",
                        "GIT_COMMITTER_EMAIL": "t@t"})


def _repo(tmp_path):
    (tmp_path / "pkg").mkdir()
    (tmp_path / "tests").mkdir()
    (tmp_path / "pkg" / "module.py").write_text("X = 1\n")
    (tmp_path / "tests" / "test_module.py").write_text("def test_x(): pass\n")
    _git(tmp_path, "init", "-q")
    _git(tmp_path, "add", "-A")
    _git(tmp_path, "commit", "-qm", "init")
    return tmp_path


def _green(target="tests/test_module.py"):
    return {"tool": "run_tests", "command": target, "exit_code": 0,
            "status": "ok", "passed": 1, "failed": 0, "ts": 100.0}


def _claim(tests, workspace):
    return check_completion_claim(
        "Run the tests", "", changes=[], observed=None,
        tests=tests, window_start=0.0,
        current_fingerprint=api_client._current_fingerprint(tests, workspace))


def test_only_the_new_runs_are_stamped(tmp_path):
    repo = _repo(tmp_path)
    ledger = [_green()]
    ledger.append(_green())
    api_client._stamp_new_evidence(ledger, 1, repo)
    assert "fingerprint" not in ledger[0]
    assert ledger[1]["fingerprint"]["commit"]


def test_outside_git_nothing_is_stamped_or_judged(tmp_path):
    ledger = [_green()]
    api_client._stamp_new_evidence(ledger, 0, tmp_path)
    assert "fingerprint" not in ledger[0]
    assert api_client._current_fingerprint(ledger, tmp_path) is None


def test_commit_then_done_is_verified_and_a_later_edit_is_not(tmp_path):
    repo = _repo(tmp_path)
    (repo / "pkg" / "module.py").write_text("X = 2\n")
    ledger = [_green()]
    api_client._stamp_new_evidence(ledger, 0, repo)
    _git(repo, "commit", "-qam", "the tested change")
    assert _claim(ledger, repo)["verdict"] == "verified"
    (repo / "pkg" / "module.py").write_text("X = 3  # after the run\n")
    verdict = _claim(ledger, repo)
    assert verdict["verdict"] == "unmet"
    assert verdict["kind"] == "tests_stale_state"
