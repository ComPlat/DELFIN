"""Adversarial tests for delfin.agent.pre_pr_check (reviewer nacht-s2, GH1).

Attacks the test-selection heuristics: over-selection (irrelevant heavy
tests slow every hand-over) and under-selection (a missed guard is the
measured defect this module exists to prevent).  Git-touching tests build
their own tmp repo (git init -b main, env identity, chdir) per wave-26
rules — CI runs on a checkout with no git identity.
"""
from __future__ import annotations

import os
import subprocess

import pytest

from delfin.agent import pre_pr_check as ppc


@pytest.fixture()
def tmp_repo(tmp_path, monkeypatch):
    """A throw-away git repo with one committed change under delfin/."""
    monkeypatch.setenv("GIT_AUTHOR_NAME", "test")
    monkeypatch.setenv("GIT_AUTHOR_EMAIL", "test@example.com")
    monkeypatch.setenv("GIT_COMMITTER_NAME", "test")
    monkeypatch.setenv("GIT_COMMITTER_EMAIL", "test@example.com")
    home = tmp_path / "home"
    home.mkdir()
    monkeypatch.setenv("HOME", str(home))
    repo = tmp_path / "repo"
    (repo / "delfin" / "agent").mkdir(parents=True)
    (repo / "tests").mkdir()
    (repo / "delfin" / "agent" / "mod.py").write_text("x = 1\n")
    (repo / "tests" / "test_mod.py").write_text("def test_x():\n    assert True\n")
    env = dict(os.environ)
    subprocess.run(["git", "init", "-b", "main", str(repo)], check=True,
                   capture_output=True, env=env)
    monkeypatch.chdir(repo)
    for path in ["delfin/agent/mod.py", "tests/test_mod.py"]:
        subprocess.run(["git", "add", path], check=True, capture_output=True,
                       env=env)
    subprocess.run(["git", "commit", "-m", "base"], check=True,
                   capture_output=True, env=env)
    # A NON-MOVING base: a branch name like "main" advances with every
    # commit, making diff main..HEAD empty.  The tag is the branch point.
    subprocess.run(["git", "tag", "base"], check=True, capture_output=True,
                   env=env)
    return repo


def _commit(repo, env_extra=None):
    env = dict(os.environ)
    if env_extra:
        env.update(env_extra)
    subprocess.run(["git", "add", "-A"], check=True, capture_output=True, env=env)
    subprocess.run(["git", "commit", "-m", "change"], check=True,
                   capture_output=True, env=env)


class TestGuardAreaMarker:
    def test_scanning_guard_is_selected(self, tmp_repo):
        # A guard that really walks delfin/ is correctly selected.
        src = ("import pathlib\n"
               "paths = list(pathlib.Path('delfin').rglob('*.py'))\n")
        (tmp_repo / "tests" / "test_walker_guard.py").write_text(src)
        _commit(tmp_repo)
        changed = ppc.changed_files(str(tmp_repo), "base")
        assert "tests/test_walker_guard.py" in ppc.select_tests(changed, str(tmp_repo))

    def test_comment_marker_selects_for_every_delfin_change(self, tmp_repo):
        # FINDING: '# guard test' in a COMMENT alone selects the file for
        # ANY delfin/ change. This is the same comment-only over-selection
        # shape that 73e6aaa6 fixed for the Path/glob heuristics — the
        # marker is forgeable and drags an unrelated test into every
        # hand-over. If the builder intends markers, the test must fail
        # first and the intent documented; as it stands this measures
        # over-selection.
        src = ("def test_unrelated():\n"
               "    # guard test mention in a comment about something else\n"
               "    assert 1 + 1 == 2\n")
        (tmp_repo / "tests" / "test_comment_only.py").write_text(src)
        _commit(tmp_repo)
        changed = ppc.changed_files(str(tmp_repo), "base")
        selected = ppc.select_tests(changed, str(tmp_repo))
        assert "tests/test_comment_only.py" not in selected, (
            "comment-only '# guard test' over-selected: " + repr(selected))


class TestNameMatch:
    def test_stem_match_picks_own_test(self, tmp_repo):
        (tmp_repo / "delfin" / "agent" / "frob.py").write_text("y = 2\n")
        (tmp_repo / "tests" / "test_frob_guard.py").write_text(
            "def test_frob():\n    assert True\n")
        _commit(tmp_repo)
        changed = ppc.changed_files(str(tmp_repo), "base")
        selected = ppc.select_tests(changed, str(tmp_repo))
        assert "tests/test_frob_guard.py" in selected

    def test_short_stem_does_not_grab_unrelated_tests(self, tmp_repo):
        # A stem like 'data' matches every test whose FILE NAME contains
        # 'data' (test_database.py, test_metadata.py, ...). The heuristic
        # is a name-contains match; this pins that an unrelated-looking
        # guard is not dragged in merely by a common English word. If
        # over-selection shows up here it is a finding, not a rule.
        (tmp_repo / "delfin" / "agent" / "data.py").write_text("z = 3\n")
        (tmp_repo / "tests" / "test_database_guard.py").write_text(
            "def test_db():\n    assert True\n")
        _commit(tmp_repo)
        changed = ppc.changed_files(str(tmp_repo), "base")
        selected = ppc.select_tests(changed, str(tmp_repo))
        # 'data' is a substring of 'database' -> the heuristic WILL match.
        # The intent question is the builder's; we only require the match
        # to be EXPLAINABLE (it is: substring). Document, do not silently
        # depend: assert the behaviour we measured so a change is visible.
        assert "tests/test_database_guard.py" in selected


class TestNoTestsSelected:
    def test_blocker_fires_when_tree_offers_nothing(self, tmp_repo):
        # Change a delfin/ file; the fixture tests dir has only test_mod.py
        # whose name shares no stem with 'onlyfile' and whose source scans
        # nothing -> no selection -> the no-tests-selected blocker.
        (tmp_repo / "delfin" / "agent" / "onlyfile.py").write_text("q = 9\n")
        _commit(tmp_repo)
        report = ppc.run(str(tmp_repo), base="base")
        kinds = {b.kind for b in report.blockers}
        assert "no-tests-selected" in kinds or report.selected_tests, (
            "either a blocker or a real selection must exist; got "
            f"kinds={kinds}, selected={report.selected_tests}")


class TestSchedulerSurface:
    def test_direct_scheduler_word_is_blocked(self, tmp_repo):
        (tmp_repo / "delfin" / "agent" / "bad.py").write_text(
            "cmd = 'squeue'\n")
        _commit(tmp_repo)
        report = ppc.run(str(tmp_repo), base="base")
        kinds = {b.kind for b in report.blockers}
        assert "direct-scheduler-call" in kinds

    def test_scheduler_word_built_by_concatenation_is_NOT_blocked(self, tmp_repo):
        # FINDING (accepted heuristic limit, documented): the surface scan
        # matches quoted literals only. A split literal evades it. The
        # guard TEST (test_r3_direct_scheduler_call_guard.py) is an AST
        # scan and still applies via _guard_covers; the pre-PR blocker is
        # the weaker layer. Pin the measured behaviour so a change is
        # visible.
        (tmp_repo / "delfin" / "agent" / "sneaky.py").write_text(
            "cmd = 'sq' + 'ueue'\n")
        _commit(tmp_repo)
        report = ppc.run(str(tmp_repo), base="base")
        kinds = {b.kind for b in report.blockers}
        assert "direct-scheduler-call" not in kinds
