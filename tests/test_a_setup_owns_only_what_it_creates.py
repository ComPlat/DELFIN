"""A setup script must not touch a fixture it does not own.

`Task.setup` prepares a workspace before an attempt, and for a
science_analysis task that workspace is
`tests/fixtures/science_workspace/` — shared by every other science
task, and already holding `run_a.out`..`run_d.out`, `pipeline.py`,
`ensemble.csv` and a `README.md` describing them.

The first version of `two_timings_that_overlap.py` wrote its own
`README.md` there. The guard around an attempt snapshots and restores
the fixture directories, so nothing would have survived to the next
task — which is not a licence to write over somebody else's fixture. A
guard that has to undo a write is one failure away from not undoing it,
and the `_PristineWorkspace.failed` flag exists because a snapshot that
did not happen used to be swallowed.

This is asserted for EVERY setup script the suite declares, not just
that one, so the next script added inherits the rule.
"""

from __future__ import annotations

import shutil
import subprocess
import sys
from pathlib import Path

import pytest

from delfin.agent.benchmark import load_tasks
from delfin.agent.benchmark_runner import setup_path

_ROOT = Path(__file__).resolve().parents[1]
_WS_FOR_CLASS = {
    "science_analysis": "science_workspace",
    "generic_project": "user_project_workspace",
}


def _setups():
    seen = {}
    for t in load_tasks():
        if t.setup:
            seen.setdefault(t.setup, t.task_class)
    return sorted(seen.items())


def test_the_suite_still_has_setup_scripts():
    """The premise: a parametrised test over an empty list passes."""
    assert _setups(), "no task declares a setup script any more"


@pytest.mark.parametrize("name,task_class", _setups())
def test_a_setup_adds_files_and_changes_none(name, task_class, tmp_path):
    rel = _WS_FOR_CLASS.get(task_class)
    if rel is None:
        pytest.skip(f"no fixture workspace mapped for {task_class}")
    source = _ROOT / "tests" / "fixtures" / rel
    if not source.is_dir():
        pytest.skip(f"{source} is not in this checkout")

    ws = tmp_path / rel
    shutil.copytree(source, ws)
    before = {p.relative_to(ws): p.read_bytes()
              for p in ws.rglob("*") if p.is_file()}

    proc = subprocess.run([sys.executable, str(setup_path(name)), str(ws)],
                          capture_output=True, text=True, timeout=300)
    assert proc.returncode == 0, proc.stderr[:400]

    after = {p.relative_to(ws): p.read_bytes()
             for p in ws.rglob("*") if p.is_file()}
    changed = sorted(str(k) for k, v in before.items() if after.get(k) != v)
    assert not changed, f"{name} rewrote fixture files it does not own: {changed}"
    assert set(after) - set(before), f"{name} created nothing"


@pytest.mark.parametrize("name,task_class", _setups())
def test_a_setup_refuses_to_run_twice(name, task_class, tmp_path):
    """Non-zero on a second run, so a workspace that was not cleaned is
    reported as unmeasured rather than silently half-rebuilt."""
    rel = _WS_FOR_CLASS.get(task_class)
    if rel is None:
        pytest.skip(f"no fixture workspace mapped for {task_class}")
    ws = tmp_path / rel
    ws.mkdir(parents=True)
    first = subprocess.run([sys.executable, str(setup_path(name)), str(ws)],
                           capture_output=True, text=True, timeout=300)
    assert first.returncode == 0, first.stderr[:400]
    second = subprocess.run([sys.executable, str(setup_path(name)), str(ws)],
                            capture_output=True, text=True, timeout=300)
    assert second.returncode != 0, (
        f"{name} rebuilt over its own output instead of refusing")
    assert "refusing to overwrite" in second.stderr


@pytest.mark.parametrize("name,task_class", _setups())
def test_a_setup_says_what_it_wants_when_given_nothing(name, task_class):
    proc = subprocess.run([sys.executable, str(setup_path(name))],
                          capture_output=True, text=True, timeout=60)
    assert proc.returncode != 0
    assert "usage" in (proc.stdout + proc.stderr).lower()
