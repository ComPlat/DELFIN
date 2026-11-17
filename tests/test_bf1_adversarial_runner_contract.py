"""BF1 adversarial (reviewer nacht-s14): every verify: script named in the
benchmark yaml must grade through the REAL caller, benchmark_runner.
run_acceptance, which hands each script exactly ONE argument -- the
workspace -- and runs it with cwd = workspace.

Red on main (foundings, package BF1):
- cas1_verify.py, w20_workflow_verify.py, w21_build_verify.py,
  w21_run_verify.py require len(argv) == 3 and exit 2 ("usage") under
  the real caller, so every CAS1 / W20-workflow / W21 task fails
  regardless of what the agent did.
- w17_longrun_wait_verify.py takes an optional second argument (a
  longrun_inputs JSON path) the runner never passes; with one argument
  the waiting rubric reads "unmeasured" and the grade is False.
"""
from __future__ import annotations

import json
from pathlib import Path

import pytest

from delfin.agent.benchmark_runner import run_acceptance


def _ws(tmp_path: Path) -> Path:
    ws = tmp_path / "fixture_workspace"
    ws.mkdir()
    return ws


# --- cas1_verify.py -------------------------------------------------------

def _cas1_good_answer(ws: Path, task_id: str, trap: bool = False) -> None:
    d = ws / task_id
    d.mkdir(parents=True)
    if trap:
        (d / "answer.txt").write_text(
            "Nactel is too large: NEL exceeds feasible space. "
            "Use RAS or ICE or DMRG with nel and norb set.",
            encoding="utf-8")
    else:
        (d / "answer.txt").write_text(
            "proposal validated", encoding="utf-8")
        (d / "proposal.json").write_text(json.dumps(
            {"n_electrons": 6, "n_orbitals": 5,
             "justification": "pi system",
             "named_orbital_classes":
                 ["2 pi bonding + 3 pi* antibonding electrons in the "
                  "pi/pi* orbitals"]}),
            encoding="utf-8")


def test_cas1_trap_does_not_exit_with_usage(tmp_path):
    """The trap answer is deliberately WRONG content-wise; the point is
    the FAILURE MESSAGE: on main the script dies with its usage line
    before it ever grades anything."""
    ws = _ws(tmp_path)
    _cas1_good_answer(ws, "cas1_trap_oversized", trap=True)
    ok, out = run_acceptance("cas1_verify.py", ws, root=Path.cwd())
    assert ok is not None, f"not run at all: {out}"
    assert "usage" not in out, out


def test_cas1_proposal_grades_ok_through_the_runner(tmp_path):
    ws = _ws(tmp_path)
    _cas1_good_answer(ws, "cas1_propose_ozone")
    ok, out = run_acceptance("cas1_verify.py", ws, root=Path.cwd())
    assert ok is True, out


# --- w20_workflow_verify.py ----------------------------------------------

def test_w20_workflow_does_not_exit_with_usage(tmp_path):
    ws = _ws(tmp_path)
    d = ws / "w20_opt_1"
    d.mkdir()
    ok, out = run_acceptance("w20_workflow_verify.py", ws, root=Path.cwd())
    assert ok is not None, f"not run at all: {out}"
    assert "usage" not in out, out


# --- w21_build_verify.py / w21_run_verify.py ------------------------------

def test_w21_build_does_not_exit_with_usage(tmp_path):
    ws = _ws(tmp_path)
    ok, out = run_acceptance("w21_build_verify.py", ws, root=Path.cwd())
    assert ok is not None, f"not run at all: {out}"
    assert "usage" not in out, out


def test_w21_run_does_not_exit_with_usage(tmp_path):
    ws = _ws(tmp_path)
    ok, out = run_acceptance("w21_run_verify.py", ws, root=Path.cwd())
    assert ok is not None, f"not run at all: {out}"
    assert "usage" not in out, out


# --- w17_longrun_wait_verify.py ------------------------------------------

def test_w17_longrun_wait_grades_fail_not_usage(tmp_path):
    """The waiting rubric is unmeasurable through the real caller (QS
    ruling, 2026-11): the script must keep grading HONESTLY unmeasured —
    a grade failure (exit 1 with 'unmeasured'), never a usage error
    (exit 2). Measurement itself is a protected-hook request to the
    operator, not an inference inside this script. Supersedes the
    earlier case that demanded 'unmeasured' not appear at all."""
    ws = _ws(tmp_path)
    d = ws / "w17_longrun" / "lr_wait_for_job"
    d.mkdir(parents=True)
    (d / "result.json").write_text(json.dumps(
        {"n": 6, "mean_kj_per_mol": 42.19}), encoding="utf-8")
    ok, out = run_acceptance(
        "w17_longrun_wait_verify.py", ws, root=Path.cwd())
    assert ok is not None, f"not run at all: {out}"
    assert ok is False, out
    assert "unmeasured" in out, out
    assert "usage" not in out, out
