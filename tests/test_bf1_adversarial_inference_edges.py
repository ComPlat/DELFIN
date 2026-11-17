"""BF1 adversarial edge cases (reviewer nacht-s14) against the fixed
1-argument task inference: wrong-family deliverables, ambiguity, unknown
folders, unicode, and the no-deliverable workspace. All through the REAL
run_acceptance. On the fixed code these must behave as graded failures
(exit 1 / ok=False), never usage errors (exit 2 / "usage")."""
from __future__ import annotations

import json
from pathlib import Path

from delfin.agent.benchmark_runner import run_acceptance


def _ws(tmp_path: Path) -> Path:
    ws = tmp_path / "fixture_workspace"
    ws.mkdir()
    return ws


def _spec(dir_: Path) -> None:
    dir_.mkdir(parents=True)
    (dir_ / "spec.json").write_text("{}", encoding="utf-8")


def _result(dir_: Path) -> None:
    dir_.mkdir(parents=True)
    (dir_ / "result.json").write_text(json.dumps({}), encoding="utf-8")


def test_w21_build_wrong_family_folder_is_not_graded_as_a_task(tmp_path):
    """A deliverable in a folder outside the known task ids must NOT be
    adopted as the task (that would let an agent pass by writing to a
    folder the task never names). It grades as nothing -> exit 1."""
    ws = _ws(tmp_path)
    _spec(ws / "i_am_not_a_task")
    ok, out = run_acceptance("w21_build_verify.py", ws, root=Path.cwd())
    assert ok is False, (ok, out)
    assert "usage" not in out, out


def test_w21_run_two_task_folders_is_a_failure_not_a_crash(tmp_path):
    """Ambiguity (several task folders with deliverables) grades as a
    failure (exit 1) -- the grade must be deterministic, not lucky."""
    ws = _ws(tmp_path)
    for name in ("w21_run_methanol", "w21_run_water_ir"):
        _spec(ws / name)
    ok, out = run_acceptance("w21_run_verify.py", ws, root=Path.cwd())
    assert ok is False, (ok, out)
    assert "ambiguous" in out, out


def test_w20_workflow_deliverable_only_in_unknown_folder_fails(tmp_path):
    ws = _ws(tmp_path)
    _result(ws / "mystery_folder")
    ok, out = run_acceptance("w20_workflow_verify.py", ws, root=Path.cwd())
    assert ok is False, (ok, out)
    assert "usage" not in out, out


def test_w20_workflow_unicode_folder_name_does_not_crash(tmp_path):
    """A unicode-named deliverable folder that is NOT a task id must be
    handled: not adopted, not a crash, a graded failure."""
    ws = _ws(tmp_path)
    _result(ws / "wörkfløw")
    ok, out = run_acceptance("w20_workflow_verify.py", ws, root=Path.cwd())
    assert ok is False, (ok, out)


def test_w21_build_empty_workspace_fails_not_usage(tmp_path):
    ws = _ws(tmp_path)
    ok, out = run_acceptance("w21_build_verify.py", ws, root=Path.cwd())
    assert ok is False, (ok, out)
    assert "usage" not in out, out


def test_w21_build_missing_workspace_arg_is_usage(tmp_path):
    """The usage error is reserved for a malformed invocation: with NO
    argument at all the usage line is correct behaviour."""
    import subprocess
    import sys
    script = (Path("delfin") / "agent" / "pack" / "benchmark" / "accept"
              / "w21_build_verify.py").resolve()
    proc = subprocess.run(
        [sys.executable, str(script)], capture_output=True, text=True,
        cwd=str(Path.cwd()))
    assert proc.returncode == 2, (proc.returncode, proc.stderr)
    assert "usage" in proc.stderr
