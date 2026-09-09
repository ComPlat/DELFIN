"""A task about a working program is graded on the program.

The coding side of this suite was four tasks, each a ~30-line script over
one small JSON file, and every one graded on what the model WROTE ABOUT
its work. For "build something that works" that is the wrong question,
and no amount of pattern care fixes it: a rubric that reads prose about a
program can be satisfied by prose. The same session that added this had
just finished repairing ten rubrics that measured spelling.

So `Task.verify` names a script under `pack/benchmark/accept/`. The
runner hands it the workspace after the turn and BEFORE the fixture guard
puts it back — the artifacts exist only there — and its exit status is
the verdict, its output the reason. The script is ours, never the
model's, and it is resolved inside that directory so a task file, which
is data, cannot name a path outside it.

Three states, not two: passed, failed, and NOT RUN. A missing or timing
out script is not a model failure, the same distinction `denials` makes
between unobserved and zero and `is_unmeasured` makes about a turn that
never reached the model.
"""

from __future__ import annotations

import shutil
import subprocess
import sys
from pathlib import Path

import pytest

from delfin.agent.benchmark import Task, Trajectory, load_tasks, score_outcome
from delfin.agent.benchmark_runner import (
    _BEHAVIOR_WS_RELS, acceptance_path, run_acceptance, workspace_for,
)

_ROOT = Path(__file__).resolve().parents[1]
_FIXTURE = _ROOT / "tests" / "fixtures" / "science_workspace"
_ACCEPT = _ROOT / "delfin" / "agent" / "pack" / "benchmark" / "accept"

_TASK_ID = "science_summarises_only_the_runs_that_converged"

# A correct answer and the everyday wrong one, kept here rather than in
# the fixture: they are what this file asserts ABOUT the acceptance
# script, and a solution sitting in the workspace would be handed to the
# model it is meant to measure.
_CORRECT = '''import csv, re
from pathlib import Path
E = re.compile(r"TOTAL ENERGY\\s+(-?\\d+\\.\\d+)\\s+Eh")
G = re.compile(r"HOMO-LUMO GAP\\s+(\\d+\\.\\d+)\\s+eV")
rows = []
for p in sorted(Path(".").glob("*.out")):
    t = p.read_text()
    if "abnormal termination" in t or "normal termination" not in t:
        continue
    rows.append([p.stem, E.search(t).group(1), G.search(t).group(1)])
with open("summary.csv", "w", newline="") as fh:
    w = csv.writer(fh)
    w.writerow(["run", "energy_hartree", "gap_ev"])
    w.writerows(rows)
'''

_NAIVE = _CORRECT.replace(
    '    if "abnormal termination" in t or "normal termination" not in t:\n'
    '        continue\n', "")


@pytest.fixture
def workspace(tmp_path):
    if not _FIXTURE.is_dir():
        pytest.skip("the science fixture is not in this checkout")
    ws = tmp_path / "science_workspace"
    shutil.copytree(_FIXTURE, ws)
    return ws


# ---------------------------------------------------------------------------
# The acceptance runner
# ---------------------------------------------------------------------------

def test_a_script_that_accepts_says_so():
    ok, out = run_acceptance("_selftest_ok.py", Path.cwd(), root=_ROOT)
    assert ok is True, out


def test_a_script_that_rejects_gives_its_reason():
    ok, out = run_acceptance("_selftest_fail.py", Path.cwd(), root=_ROOT)
    assert ok is False
    assert "reason" in out, out


def test_a_missing_script_is_not_run_rather_than_failed():
    """The distinction the score depends on: our harness being broken is
    not the model getting it wrong."""
    ok, out = run_acceptance("_no_such_script.py", Path.cwd(), root=_ROOT)
    assert ok is None
    assert "not found" in out


@pytest.mark.parametrize("escape", [
    "../../../../etc/passwd", "/etc/passwd", "../../conftest.py",
])
def test_a_task_cannot_name_a_script_outside_the_pack(escape):
    """A task file is data. It must not be able to point the runner at an
    arbitrary path and have it executed."""
    with pytest.raises(ValueError):
        acceptance_path(escape, _ROOT)


def test_the_scripts_live_where_the_resolver_looks():
    assert _ACCEPT.is_dir(), f"{_ACCEPT} is missing"
    assert acceptance_path("_selftest_ok.py", _ROOT).is_file()


# ---------------------------------------------------------------------------
# The verdict reaches the score
# ---------------------------------------------------------------------------

def _scored(verify_ok, output=""):
    task = Task(id="probe", task_class="science_analysis", mode="solo",
                prompt="p", verify="_selftest_ok.py")
    return score_outcome(task, Trajectory(
        text="fertig", verify_ok=verify_ok, verify_output=output))


def test_an_accepted_artifact_passes_the_task():
    assert _scored(True).success


def test_a_rejected_artifact_fails_it_and_keeps_the_reason():
    res = _scored(False, "summary.csv carries the energy of run_d")
    assert not res.success
    assert any("verify:failed" in m for m in res.missing_signals)
    assert "run_d" in " ".join(res.signal_evidence.values())


def test_a_script_that_never_ran_is_not_scored_as_a_pass():
    res = _scored(None)
    assert not res.success
    assert any("verify:not_run" in m for m in res.missing_signals)


def test_a_task_without_a_verify_is_untouched():
    task = Task(id="p", task_class="generic_project", mode="solo", prompt="p")
    res = score_outcome(task, Trajectory(text="ok"))
    assert res.success and res.verify_ok is None


# ---------------------------------------------------------------------------
# The first task, and what its acceptance script actually decides
# ---------------------------------------------------------------------------

def _task():
    return next(t for t in load_tasks() if t.id == _TASK_ID)


def test_the_task_names_a_script_that_exists():
    task = _task()
    assert task.verify
    assert acceptance_path(task.verify, _ROOT).is_file()


def test_the_fixture_holds_a_run_that_did_not_converge():
    """The premise. Without it the task measures nothing."""
    texts = {p.name: p.read_text(encoding="utf-8")
             for p in _FIXTURE.glob("*.out")}
    assert len(texts) >= 3
    bad = [n for n, t in texts.items() if "abnormal termination" in t]
    assert len(bad) == 1, f"expected exactly one failed run, got {bad}"


def _accept(workspace) -> subprocess.CompletedProcess:
    return subprocess.run(
        [sys.executable, str(acceptance_path(_task().verify, _ROOT)),
         str(workspace)],
        capture_output=True, text=True, timeout=120)


def test_the_correct_answer_is_accepted(workspace):
    (workspace / "summary.py").write_text(_CORRECT)
    subprocess.run([sys.executable, "summary.py"], cwd=str(workspace),
                   check=True, capture_output=True)
    proc = _accept(workspace)
    assert proc.returncode == 0, proc.stdout + proc.stderr


def test_summarising_the_run_that_failed_is_rejected(workspace):
    """The everyday mistake: the numbers are there, so they get used."""
    (workspace / "summary.py").write_text(_NAIVE)
    subprocess.run([sys.executable, "summary.py"], cwd=str(workspace),
                   check=True, capture_output=True)
    proc = _accept(workspace)
    assert proc.returncode == 1
    assert "abnormal termination" in proc.stdout


def test_a_summary_that_was_never_produced_is_rejected(workspace):
    proc = _accept(workspace)
    assert proc.returncode == 1
    assert "summary.py" in proc.stdout


def test_normal_termination_is_not_read_inside_abnormal():
    """The acceptance script's own trap, found by running it against the
    untouched fixture: "normal termination" is a substring of "abnormal
    termination", so the naive test reads every run as converged --
    including the only one this task is about."""
    source = acceptance_path(_task().verify, _ROOT).read_text(encoding="utf-8")
    assert "abnormal termination" in source, (
        "the converged check no longer excludes the abnormal case")


# ---------------------------------------------------------------------------
# The workspace it writes into
# ---------------------------------------------------------------------------

def test_the_science_class_gets_its_own_fixture():
    ws = workspace_for(_ROOT, task_class="science_analysis")
    assert ws == _FIXTURE


def test_the_science_workspace_is_restored_between_attempts():
    """It was not. A task graded on what it BUILDS always leaves
    something behind, and adding the workspace to workspace_for without
    adding it here left the built script in the checkout after the run."""
    assert any(rel.name == "science_workspace" for rel in _BEHAVIOR_WS_RELS), (
        "the science fixture is not under the guard, so its artifacts "
        "survive the run")
