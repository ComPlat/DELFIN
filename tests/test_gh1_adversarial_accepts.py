"""Adversarial tests for the GH1 accept scripts (reviewer nacht-s2).

Attack: can a wrong end state slip past the graded checks?

A02 (RED, sent to nacht-s1): gh1_count_conflict_verify.py:54-59 - the
unmerged-path check is gated behind ``if "U" in status[:20]``. A wrong
resolution that leaves a conflicted path UNMERGED while other status
lines push the UU line beyond the first 20 characters is ACCEPTED
(test_unmerged_path_beyond_20_chars_still_rejected).

Pinned green (measured): a UU line inside the prefix IS caught; the
pr_ready accept grades the FEATURE BRANCH (checks it out first), so a
branch without the feature is rejected and a branch with it accepted;
the dead lines 82-84 (core read on main, `if ...: pass`) have no
grading effect.
"""
from __future__ import annotations

import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(
    Path(__file__).resolve().parent.parent
    / "delfin/agent/pack/benchmark/accept"))
sys.path.insert(0, str(
    Path(__file__).resolve().parent.parent
    / "delfin/agent/pack/benchmark/setup"))

import gh1_count_conflict_verify  # noqa: E402
import gh1_pr_ready_verify  # noqa: E402
import gh1_red_ci_verify  # noqa: E402
import gh1_repo  # noqa: E402


def _run(cmd: list, cwd: Path) -> subprocess.CompletedProcess:
    return subprocess.run(cmd, cwd=str(cwd), capture_output=True,
                          text=True, env=gh1_red_ci_verify._env(cwd))


def _mk_repo(tmp_path: Path) -> tuple[Path, str]:
    """A widgetkit repo via the real setup core."""
    ws = tmp_path / "ws"
    ws.mkdir()
    repo = gh1_repo.new_repo(ws, "widgetkit")
    gh1_repo.write(repo, {
        "src/widgetkit/core.py": 'def normalise(s):\n    return s\n',
        "src/widgetkit/rename.py": 'VALUE = "renamed"\n',
        "tests/test_core.py": "def test_smoke():\n    assert True\n",
    })
    gh1_repo.commit(repo, "initial")
    gh1_repo.mark_base(repo)
    return repo, gh1_repo.base_sha(repo)


class TestCountConflictUnmergedPaths:
    """A02: wrong resolution leaving an unmerged path must be rejected."""

    def test_unmerged_path_beyond_20_chars_still_rejected(self, tmp_path):
        repo, base = _mk_repo(tmp_path)
        counts = repo / "metrics/counts.py"
        counts.parent.mkdir()
        # 25 = correctly summed, so ONLY the unmerged-path check can
        # reject — isolates the check under attack.
        counts.write_text(
            '{"count": 25}\n', encoding="utf-8")
        gh1_repo.commit(repo, "counts")
        gh1_repo.mark_base(repo)
        base = gh1_repo.base_sha(repo)
        # Simulate the agent resolving the target file but leaving a
        # second conflicted path UNMERGED, with status output long
        # enough that the UU line sits beyond char 20.
        long_name = ("x" * 30 + "/unresolved_feature.py")
        (repo / long_name).parent.mkdir(parents=True)
        (repo / long_name).write_text(
            "<<<<<<< HEAD\nmain side\n=======\npartner side\n"
            ">>>>>>> partner/count\n", encoding="utf-8")
        # Real unmerged index entry (stage 1/2/3) so git status
        # reports UU — the exact shape of a real conflict:
        s1 = _run(["git", "hash-object", "-w", "--stdin"], repo).stdout.strip()
        s2 = _run(["git", "hash-object", "-w", "--stdin"], repo).stdout.strip()
        s3 = _run(["git", "hash-object", "-w", "--stdin"], repo).stdout.strip()
        info = ("100644 {s1} 1\t{p}\n100644 {s2} 2\t{p}\n"
                "100644 {s3} 3\t{p}\n").format(
                    s1=s1, s2=s2, s3=s3, p=long_name)
        # A tracked modification FIRST so the UU line falls OUTSIDE
        # status[:20] (porcelain puts each code at the line start, so
        # a lone UU line would always be inside the prefix). A
        # NON-test file, so the suite stays green and ONLY the
        # unmerged-path check can reject:
        (repo / "src/widgetkit/rename.py").write_text(
            'VALUE = "renamed-changed"\n', encoding="utf-8")
        _run(["git", "update-index", "src/widgetkit/rename.py"], repo)
        # The other required feature must be present so the ONLY thing
        # that can reject is the unmerged path:
        (repo / "src/widgetkit/lookup.py").write_text("X = 1\n",
                                                      encoding="utf-8")
        proc = subprocess.run(
            ["git", "update-index", "--index-info"],
            cwd=str(repo), input=info, capture_output=True, text=True,
            env=gh1_red_ci_verify._env(repo))
        assert proc.returncode == 0, proc.stderr
        status = _run(["git", "status", "--porcelain"], repo).stdout
        assert "UU" in status, status
        assert status.index("UU") > 20, status
        ok, reason, q = gh1_count_conflict_verify.check(repo.parent)
        assert ok is False, (
            "an unmerged path beyond the 20-char prefix must be rejected")
        assert "unmerged" in reason.lower(), (
            "must reject on the unmerged-path check, not the suite: "
            + reason)

    def test_status_starting_with_UU_still_rejected(self, tmp_path):
        repo, base = _mk_repo(tmp_path)
        counts = repo / "metrics/counts.py"
        counts.parent.mkdir()
        counts.write_text('{"count": 25}\n', encoding="utf-8")
        gh1_repo.commit(repo, "counts")
        gh1_repo.mark_base(repo)
        (repo / "a.py").write_text(
            "<<<<<<< HEAD\n=======\n>>>>>>> p\n", encoding="utf-8")
        _run(["git", "update-index", "--add", "--cacheinfo",
              "0,000000,a.py"], repo)
        status = _run(["git", "status", "--porcelain"], repo).stdout
        assert status  # sanity
        ok, reason, q = gh1_count_conflict_verify.check(repo.parent)
        assert ok is False, "unmerged path at the start must be rejected"


class TestPrReadyFeatureCheck:
    """A03: the branch must be graded, not main."""

    def test_branch_without_feature_rejected(self, tmp_path):
        repo, base = _mk_repo(tmp_path)
        # Feature branch exists with ONLY a test, no core.py change
        gh1_repo.feature_branch(repo, "feature/x", {
            "tests/test_core.py": "def test_norm():\n    assert True\n",
        }, "add test only, no feature")
        _run(["git", "checkout", "-q", "main"], repo)
        ok, reason, q = gh1_pr_ready_verify.check(repo.parent)
        assert ok is False, (
            "a branch without the feature must not be accepted")


class TestAcceptReadsRealState:
    """Both accepts must grade the real repo, not a stale read."""

    def test_pr_ready_feature_on_branch_accepted(self, tmp_path):
        repo, base = _mk_repo(tmp_path)
        gh1_repo.feature_branch(repo, "feature/x", {
            "src/widgetkit/core.py":
                'def normalise(s):\n    return s.replace("_", " ")\n',
            "tests/test_core.py":
                "def test_norm():\n    assert normalise('a_b') == 'a b'\n",
        }, "feature: underscore to space")
        _run(["git", "checkout", "-q", "main"], repo)
        ok, reason, q = gh1_pr_ready_verify.check(repo.parent)
        assert ok is True, reason
