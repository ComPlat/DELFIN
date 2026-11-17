"""Adversarial tests for package A4 — `check_goat_e2e` (reviewer s13).

We drive the REAL checker `check_goat_e2e` and the REAL reader
`read_mlips_goat` / `read_goat_ensemble` through synthetic run directories
built in the test (tempfile), never ORCA-generated output in the repo
(licence).  These tests pin the checker's contract on broken runs: the
checker is the thing that must NOT crash on the exact cases it exists to
report.

The two red findings on the current head (effa0e4c5):
- F1: a frame comment with NO energy in the final ensemble makes
  `check_goat_e2e` RAISE ValueError instead of returning ok=False with a
  finding.  `read_goat_ensemble` (orca_goat_opi.py:314-318) raises on a
  no-energy frame; `check_goat_e2e` (goat_e2e_check.py:64) calls it with
  no try/except.
- F2: mixing the two comment forms (orca_goat_opi.py:319-323) raises the
  same way.
"""

from __future__ import annotations

from pathlib import Path
import math

import pytest

from delfin.mlp_tools.goat_e2e_check import check_goat_e2e

_TERM = "****ORCA TERMINATED NORMALLY****"


class _RunDir:
    """Build a synthetic finished-run directory shaped like a GOAT run."""

    def __init__(self, tmp_path: Path, base: str = "extopt"):
        self.dir = tmp_path / "run"
        self.dir.mkdir()
        self.base = base

    def write_out(self, text: str) -> None:
        (self.dir / f"{self.base}.out").write_text(text)

    def write_ensemble(self, frames: list[tuple[int, str, list[str]]]) -> None:
        """frames: (natoms, comment, [atomline,...]) in file order."""
        lines: list[str] = []
        for natoms, comment, atoms in frames:
            lines.append(str(natoms))
            lines.append(comment)
            lines.extend(atoms)
        (self.dir / f"{self.base}.finalensemble.xyz").write_text("\n".join(lines) + "\n")

    def write_global_minimum(self) -> None:
        (self.dir / f"{self.base}.globalminimum.xyz").write_text(
            "5\n<E> converged=true\nC 0 0 0 0\n" * 0 + "\n"
        )


#: two plausible conformers, ascending energy (lowest first) as GOAT writes.
_A1 = "C   0.0  0.0  0.0"
_A2 = "C   0.0  0.0  1.0"
_A3 = "C   0.0  1.0  1.0"
#: the eh comment form the reader accepts: <float> converged=true.
def _E(v):
    return f"{v} converged=true"


@pytest.fixture
def healthy_run(tmp_path):
    r = _RunDir(tmp_path)
    r.write_out(f"Total Job Time (sec):     42.0\n{_TERM}\n")
    r.write_ensemble([
        (3, _E(-50.0), [_A1, _A2, _A3]),
        (3, _E(-49.5), [_A1, _A2, _A3]),
    ])
    r.write_global_minimum()
    return r


def test_successful_run_is_ok(healthy_run):
    """A clean, well-formed run is ok=True (the positive control)."""
    g = check_goat_e2e(healthy_run.dir)
    assert g["ok"] is True
    assert g["terminated_normally"] is True
    assert g["n_conformers"] == 2
    assert g["energy_source"] == "mlip"
    assert g["final"] is False


# ---------------------------------------------------------------- findings

def test_missing_out_is_error(tmp_path):
    """No .out at all -> an error finding, no crash."""
    r = _RunDir(tmp_path)
    r.write_ensemble([(3, _E(-50.0), [_A1, _A2, _A3])])
    g = check_goat_e2e(r.dir)
    assert g["ok"] is False
    assert any(f["level"] == "error" and ".out" in f["message"] for f in g["findings"])


def test_missing_termination_marker_is_error(healthy_run):
    """No TERMINATED marker in the .out -> error, even with no failure line."""
    (healthy_run.dir / f"{healthy_run.base}.out").write_text("Total Job Time (sec): 1.0\n")
    g = check_goat_e2e(healthy_run.dir)
    assert g["ok"] is False
    assert any("TERMINATED" in f["message"] for f in g["findings"])


def test_scf_not_converged_is_error(healthy_run):
    """A _failure_lines match in the .out -> error even if marker present."""
    (healthy_run.dir / f"{healthy_run.base}.out").write_text(
        f"SCF NOT CONVERGED AFTER 50 ITERATIONS\n{_TERM}\n"
    )
    g = check_goat_e2e(healthy_run.dir)
    assert g["ok"] is False
    assert any("SCF NOT CONVERGED" in f["message"].upper() for f in g["findings"])


def test_missing_ensemble_is_error(tmp_path):
    """No finalensemble.xyz -> error 'is missing'."""
    r = _RunDir(tmp_path)
    r.write_out(f"{_TERM}\n")
    r.write_global_minimum()
    g = check_goat_e2e(r.dir)
    assert g["ok"] is False
    assert any("is missing" in f["message"] for f in g["findings"])


def test_empty_ensemble_is_error(tmp_path):
    """A present but frame-less ensemble -> error 'no conformer frame'."""
    r = _RunDir(tmp_path)
    r.write_out(f"{_TERM}\n")
    (r.dir / f"{r.base}.finalensemble.xyz").write_text("")
    r.write_global_minimum()
    g = check_goat_e2e(r.dir)
    assert g["ok"] is False
    assert any("no conformer frame" in f["message"] for f in g["findings"])


def test_missing_global_minimum_is_error(tmp_path):
    r = _RunDir(tmp_path)
    r.write_out(f"{_TERM}\n")
    r.write_ensemble([(3, _E(-50.0), [_A1, _A2, _A3])])
    g = check_goat_e2e(r.dir)
    assert g["ok"] is False
    assert any("global minimum" in f["message"].lower() for f in g["findings"])


# ---------------------------------------------------------- RED findings F1/F2

def test_frame_comment_with_no_energy_returns_finding_not_raise(tmp_path):
    """A frame whose comment carries no energy is a broken ensemble.

    F1: on the current head `check_goat_e2e` RAISES ValueError here
    (read_goat_ensemble -> ValueError, no-raise path absent in the
    checker).  The E2E checker must report it as ok=False with a finding,
    not crash -- a run it is supposed to flag.  RED before the fix.
    """
    r = _RunDir(tmp_path)
    r.write_out(f"{_TERM}\n")
    # comment line is blank -> 'tokens' empty -> no energy form -> ValueError
    r.write_ensemble([(3, "", [_A1, _A2, _A3])])
    r.write_global_minimum()
    g = check_goat_e2e(r.dir)          # must not raise
    assert g["ok"] is False
    assert any(f["level"] == "error" for f in g["findings"])


def test_mixed_energy_forms_return_finding_not_raise(tmp_path):
    """An ensemble mixing '<E> converged=true' and 'Energy <x>' forms.

    F2: the real reader raises ValueError on mixed forms; the checker must
    surface it as ok=False, not propagate the exception.  RED before the fix.
    """
    r = _RunDir(tmp_path)
    r.write_out(f"{_TERM}\n")
    r.write_ensemble([
        (3, _E(-50.0), [_A1, _A2, _A3]),
        (3, "Energy -0.5", [_A1, _A2, _A3]),
    ])
    r.write_global_minimum()
    g = check_goat_e2e(r.dir)          # must not raise
    assert g["ok"] is False
    assert any(f["level"] == "error" for f in g["findings"])


# ------------------------------------------------------ energy / GPU / timing

def test_nan_energy_is_error(tmp_path):
    """A NaN ensemble energy -> ok=False (finite check is a hard error)."""
    r = _RunDir(tmp_path)
    r.write_out(f"{_TERM}\n")
    r.write_ensemble([
        (3, _E(float("nan")), [_A1, _A2, _A3]),
        (3, _E(-49.5), [_A1, _A2, _A3]),
    ])
    r.write_global_minimum()
    g = check_goat_e2e(r.dir)
    assert g["ok"] is False
    assert any("finite" in f["message"] for f in g["findings"])


def test_unsorted_energies_produce_a_signal(tmp_path):
    """Out-of-order energies (not lowest-first) are reported.

    Advisory note: the task lists 'energies finite and sorted' as checks;
    the builder reports both finite (hard error) and sorted (warn).  This
    pins that unsorted is at least SIGNALLED -- a finding exists -- while
    ok stays True (MLIP is a screen).  If the builder argues sorted must be
    a hard gate too, that is a design decision for the QS lead, not for us.
    """
    r = _RunDir(tmp_path)
    r.write_out(f"{_TERM}\n")
    r.write_ensemble([
        (3, _E(-49.5), [_A1, _A2, _A3]),
        (3, _E(-50.0), [_A1, _A2, _A3]),
    ])
    r.write_global_minimum()
    g = check_goat_e2e(r.dir)
    assert g["energies_sorted"] is False
    assert any("sorted" in f["message"] for f in g["findings"])


def test_walltime_recorded_from_out(tmp_path):
    """Wall time is read from the ORCA .out timing line."""
    r = _RunDir(tmp_path)
    r.write_out(f"Total Job Time (sec):     123.75\n{_TERM}\n")
    r.write_ensemble([(3, _E(-50.0), [_A1, _A2, _A3])])
    r.write_global_minimum()
    g = check_goat_e2e(r.dir)
    assert g["walltime_s"] == 123.75


def test_no_gpu_evidence_is_a_warning_not_error(tmp_path):
    """No CUDA line in the logs is advisory only; ok is unaffected."""
    r = _RunDir(tmp_path)
    r.write_out(f"{_TERM}\n")
    r.write_ensemble([(3, _E(-50.0), [_A1, _A2, _A3])])
    r.write_global_minimum()
    g = check_goat_e2e(r.dir)
    assert g["ok"] is True               # GPU absence never fails the run
    assert g["gpu_evidence"] is None
    assert any("GPU" in f["message"] for f in g["findings"])


def test_stale_pid_marker_is_a_warning(tmp_path):
    """A leftover .pid marker is flagged, not a hard failure."""
    r = _RunDir(tmp_path)
    r.write_out(f"{_TERM}\n")
    r.write_ensemble([(3, _E(-50.0), [_A1, _A2, _A3])])
    r.write_global_minimum()
    (r.dir / "extopt.pid").write_text("12345\n")
    g = check_goat_e2e(r.dir)
    assert g["ok"] is True
    assert any("pid" in f["message"].lower() for f in g["findings"])
