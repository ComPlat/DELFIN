"""BF1 phase 2: each fixed script grades right from wrong THROUGH THE
REAL run_acceptance -- one argument, the workspace -- and the old
two-argument form still works.
"""

from __future__ import annotations

import json
from pathlib import Path

from delfin.agent.benchmark_runner import run_acceptance

ROOT = Path(__file__).resolve().parents[1]
ACCEPT = ROOT / "delfin" / "agent" / "pack" / "benchmark" / "accept"

_GOOD_PROPOSAL = {
    "n_electrons": 18,
    "n_orbitals": 12,
    "justification": (
        "ozone: the delocalised pi system and the open-shell resonance "
        "forms need both pi bonding and antibonding orbitals active"),
    "named_orbital_classes": ["pi bonding", "pi antibonding"],
}

_TRAP_ANSWER = (
    "The keywords are NOT real %casscf keywords (ORCA manual 3.13.3: "
    "the real keys are nel and norb), and 24 orbitals exceed what "
    "CASSCF can handle. Use RAS or ICE instead, or DMRG for larger "
    "spaces.")


def test_cas1_grades_a_proposal_through_the_runner(tmp_path):
    ws = tmp_path / "ws"
    d = ws / "cas1_propose_ozone"
    d.mkdir(parents=True)
    (d / "proposal.json").write_text(json.dumps(_GOOD_PROPOSAL),
                                     encoding="utf-8")
    (d / "answer.txt").write_text("the space validates through the "
                                  "casscf workflow.", encoding="utf-8")
    ok, out = run_acceptance("cas1_verify.py", ws, root=ROOT)
    assert ok is True, out


def test_cas1_grades_a_refused_trap_through_the_runner(tmp_path):
    ws = tmp_path / "ws"
    (ws / "cas1_trap_oversized").mkdir(parents=True)
    (ws / "cas1_trap_oversized" / "answer.txt").write_text(
        _TRAP_ANSWER, encoding="utf-8")
    ok, out = run_acceptance("cas1_verify.py", ws, root=ROOT)
    assert ok is True, out


def test_cas1_refuses_a_passed_through_trap_through_the_runner(tmp_path):
    ws = tmp_path / "ws"
    (ws / "cas1_trap_oversized").mkdir(parents=True)
    (ws / "cas1_trap_oversized" / "answer.txt").write_text(
        "Nactel 24 and Nactorb 24 are fine keywords.", encoding="utf-8")
    ok, out = run_acceptance("cas1_verify.py", ws, root=ROOT)
    assert ok is False
    assert "usage" not in out.lower()


def test_cas1_fails_on_an_empty_workspace_without_a_usage_error(tmp_path):
    ws = tmp_path / "ws"
    ws.mkdir()
    ok, out = run_acceptance("cas1_verify.py", ws, root=ROOT)
    assert ok is False
    assert "usage" not in out.lower()


def test_cas1_old_two_argument_form_still_works(tmp_path):
    """The shipped contract: `<workspace> <task_id>` keeps grading the
    named folder (extend, never lose)."""
    ws = tmp_path / "ws"
    d = ws / "cas1_trap_oversized"
    d.mkdir(parents=True)
    (d / "answer.txt").write_text(_TRAP_ANSWER, encoding="utf-8")
    script = ACCEPT / "cas1_verify.py"
    import subprocess
    import sys
    proc = subprocess.run(
        [sys.executable, str(script), str(ws), "cas1_trap_oversized"],
        capture_output=True, text=True, cwd=str(ws))
    assert proc.returncode == 0, proc.stderr


def test_w20_workflow_grades_right_and_wrong_through_the_runner(tmp_path):
    good = {
        "energy_hartree": -154.09, "n_imag": 0,
        "thermal_corr_298k": 0.098,
    }
    bad = {"energy_hartree": -154.09, "n_imag": 1,
           "thermal_corr_298k": 0.098}
    for case, data, expect in (("good", good, True), ("bad", bad, False)):
        ws = tmp_path / f"ws_{case}"
        d = ws / "w20_classic_opt_freq_1"
        d.mkdir(parents=True)
        (d / "result.json").write_text(json.dumps(data), encoding="utf-8")
        (d / "answer.txt").write_text(
            "the optimised structure is a saddle point"
            if expect is False and data["n_imag"] else
            "no negative curvature: not a saddle point", encoding="utf-8")
        ok, out = run_acceptance("w20_workflow_verify.py", ws, root=ROOT)
        assert ok is expect, out


def test_w20_workflow_two_argument_form_still_works(tmp_path):
    ws = tmp_path / "ws"
    d = ws / "w20_redox_1"
    d.mkdir(parents=True)
    (d / "result.json").write_text(
        json.dumps({"energy_oxidized_hartree": -154.0,
                    "energy_reduced_hartree": -154.3,
                    "potential_v": 0.45}), encoding="utf-8")
    script = ACCEPT / "w20_workflow_verify.py"
    import subprocess
    import sys
    proc = subprocess.run(
        [sys.executable, str(script), str(ws), "w20_redox_1"],
        capture_output=True, text=True, cwd=str(ws))
    assert proc.returncode == 0, proc.stderr


def test_w21_build_grades_right_and_wrong_through_the_runner(tmp_path):
    good = {
        "name": "ethanol_conformer_thermo",
        "steps": [
            {"step": "crest_conformers", "kwargs": {"charge": 0}},
            {"step": "orca_opt_freq",
             "kwargs": {"charge": 0, "mult": 1, "method": "PBE0",
                        "basis": "def2-TZVP", "freq_temp": 298.15,
                        "freq_pressure": 1.0}},
        ],
    }
    bad = {
        "name": "flexible_noconf",
        "steps": [{"step": "orca_opt_freq",
                   "kwargs": {"charge": 0, "mult": 1, "method": "PBE0",
                              "basis": "def2-TZVP"}}],
    }
    for case, spec, expect in (("good", good, True), ("bad", bad, False)):
        ws = tmp_path / f"ws_{case}"
        d = ws / "w21_build_conformers"
        d.mkdir(parents=True)
        (d / "brief.txt").write_text(
            "ethanol is flexible. I need its conformer landscape mapped "
            "cheaply, then the thermochemistry of the lowest conformer "
            "at one consistent level of theory.", encoding="utf-8")
        (d / "spec.json").write_text(json.dumps(spec), encoding="utf-8")
        (d / "answer.txt").write_text(
            "no deviations from the defaults needed.", encoding="utf-8")
        ok, out = run_acceptance("w21_build_verify.py", ws, root=ROOT)
        assert ok is expect, out


def test_w21_run_requires_run_evidence_through_the_runner(tmp_path):
    """The negative grade: a spec without real run evidence fails -- the
    run half exists in the graded contract, not only the spec half."""
    ws = tmp_path / "ws"
    d = ws / "w21_run_methanol"
    d.mkdir(parents=True)
    (d / "brief.txt").write_text(
        "Optimise methanol in the gas phase and confirm the geometry "
        "is a true minimum.", encoding="utf-8")
    (d / "spec.json").write_text(json.dumps({
        "name": "methanol_opt_freq",
        "steps": [{"step": "orca_opt_freq",
                   "kwargs": {"charge": 0, "mult": 1, "method": "PBE0",
                              "basis": "def2-TZVP"}}],
    }), encoding="utf-8")
    (d / "answer.txt").write_text("no deviations.", encoding="utf-8")
    ok, out = run_acceptance("w21_run_verify.py", ws, root=ROOT)
    assert ok is False
    assert "usage" not in out.lower()
    assert "run" in out.lower()
