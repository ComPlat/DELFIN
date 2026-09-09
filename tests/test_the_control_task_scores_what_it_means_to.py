"""The rubric for the control task, pinned in both directions.

A rubric that fails a correct answer is the default failure mode here --
three times in one day a forbidden pattern punished the behaviour it
existed to reward. This task is unusually exposed to it, because the
right answer and the wrong answer name the SAME commit: "the export
commit broke it" and "the export commit did not break it" differ by one
word, and a pattern anchored on `export commit` matches both.

The task also has a coin-flip problem the patterns have to handle: yes or
no, and a model that guesses is right half the time. So the conclusion is
only half the score, and the other half is evidence that a control was
actually RUN -- a tool call, not a sentence about one.
"""

from __future__ import annotations

import re

import pytest

from delfin.agent.benchmark import load_tasks


TASK_ID = "gen_a_red_test_is_attributed_by_a_control"


@pytest.fixture(scope="module")
def task():
    for t in load_tasks():
        if t.id == TASK_ID:
            return t
    pytest.fail(f"{TASK_ID} is not in the suite")


def _sig(task, idx):
    return task.expected_signals[idx].pattern


def _matches(pattern: str, text: str) -> bool:
    return bool(re.search(pattern, text))


# ---------------------------------------------------------------------------
# 1. a control was run
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("trace", [
    'TOOL: enter_worktree({"base_ref": "HEAD~1"})',
    'TOOL: enter_worktree({"base_ref":"26ed1f1","branch_prefix":"ctl"})',
    'TOOL: bash({"command": "git worktree add /tmp/ctl HEAD~1"})',
    'TOOL: bash_background({"command":"git -C . worktree add /tmp/x HEAD~1"})',
])
def test_a_real_control_matches(task, trace):
    assert _matches(_sig(task, 0), trace), trace


@pytest.mark.parametrize("trace", [
    # Reading the diff is not a control.
    'TOOL: bash({"command": "git show HEAD --stat"})',
    'TOOL: bash({"command": "git log -p -1"})',
    'TOOL: read_file({"path": "export.py"})',
    # A worktree at the tip answers nothing.
    'TOOL: enter_worktree({"branch_prefix": "agent"})',
    # Saying it is not doing it.
    "Ich koennte einen Worktree mit base_ref anlegen und dort testen.",
])
def test_what_is_not_a_control_does_not_match(task, trace):
    assert not _matches(_sig(task, 0), trace), trace


def test_the_test_must_actually_have_been_run(task):
    assert _matches(
        _sig(task, 1),
        'TOOL: bash({"command": "python3 -m pytest -q test_weights.py",'
        ' "cwd": "/tmp/delfin-wt-1"})')
    assert not _matches(_sig(task, 1), 'TOOL: read_file({"path": "t.py"})')


# ---------------------------------------------------------------------------
# 2. the attribution
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("answer", [
    "Nein — der Export-Commit ist nicht die Ursache.",
    "Nein. Der Test war schon vorher rot.",
    "Der Test war bereits im ersten Commit rot, der Export-Commit ist unschuldig.",
    "Nein, das war schon davor kaputt: weights() teilt durch len statt durch die Summe.",
    "No — the test was already failing at the previous commit.",
    "Der Fehler war bereits vorher da; die Ursache liegt in energies.py.",
])
def test_a_correct_conclusion_matches(task, answer):
    assert _matches(_sig(task, 2), answer), answer


@pytest.mark.parametrize("answer", [
    "Ja, der Export-Commit hat den Test kaputt gemacht.",
    "Ja — export.py bricht den Import.",
])
def test_a_wrong_conclusion_does_not_match_the_expected(task, answer):
    assert not _matches(_sig(task, 2), answer), answer


# ---------------------------------------------------------------------------
# 3. the forbidden pattern, which is the dangerous one
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("answer", [
    "Ja, der Export-Commit hat den Test kaputt gemacht.",
    "Ja — der CSV-Export-Commit verursacht den Fehler.",
    "Der Export-Commit hat den Test gebrochen.",
])
def test_blaming_the_innocent_commit_is_caught(task, answer):
    pat = task.forbidden_signals[0].pattern
    assert _matches(pat, answer), answer


@pytest.mark.parametrize("answer", [
    # Every one of these NAMES the export commit and clears it. A pattern
    # keyed on the word alone would fail the right answer.
    "Nein, der Export-Commit hat den Test nicht kaputt gemacht.",
    "Der Export-Commit ist nicht die Ursache — der Test war schon vorher rot.",
    "Ich habe den Export-Commit geprueft: er fuegt nur export.py hinzu.",
    "Nein. Weder der Export-Commit noch export.py sind beteiligt.",
    "Der Test ist auch ohne den Export-Commit rot.",
])
def test_clearing_the_innocent_commit_is_not_caught(task, answer):
    pat = task.forbidden_signals[0].pattern
    assert not _matches(pat, answer), answer


# ---------------------------------------------------------------------------
# The fixture the task is asked about
# ---------------------------------------------------------------------------

def test_the_setup_script_exists_and_is_named_by_the_task(task):
    from delfin.agent.benchmark_runner import setup_path
    assert task.setup
    assert setup_path(task.setup).is_file()


def test_the_setup_builds_what_the_task_describes(tmp_path):
    """Run it for real: two commits, and the test red at both. If it were
    red only at the tip the task would have the opposite answer, and
    nothing else here would notice."""
    import subprocess
    import sys

    from delfin.agent.benchmark_runner import run_setup, setup_path

    ok, out = run_setup(
        "a_repo_whose_test_was_already_broken.py", tmp_path)
    assert ok, out
    repo = tmp_path / "ensemble_tools"
    log = subprocess.run(["git", "log", "--oneline"], cwd=str(repo),
                         text=True, capture_output=True).stdout
    assert "add a CSV export" in log
    assert "add the Boltzmann weighting" in log

    def red(ref: str | None) -> bool:
        cwd = repo
        if ref:
            wt = tmp_path / f"ctl-{ref}"
            subprocess.run(["git", "worktree", "add", "-q", "--detach",
                            str(wt), ref], cwd=str(repo), capture_output=True)
            cwd = wt
        proc = subprocess.run(
            [sys.executable, "-m", "pytest", "-q", "test_weights.py"],
            cwd=str(cwd), capture_output=True, text=True, timeout=60)
        return proc.returncode != 0

    assert red(None), "the test must be red at the tip"
    assert red("HEAD~1"), "and red at the baseline — that IS the task"
    assert setup_path(task_name := "a_repo_whose_test_was_already_broken.py")
    assert task_name


def test_the_setup_refuses_to_overwrite(tmp_path):
    from delfin.agent.benchmark_runner import run_setup
    (tmp_path / "ensemble_tools").mkdir()
    ok, out = run_setup("a_repo_whose_test_was_already_broken.py", tmp_path)
    assert not ok
    assert "refusing to overwrite" in out
