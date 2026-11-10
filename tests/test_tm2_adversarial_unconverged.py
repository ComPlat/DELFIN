"""Adversarial reviewer tests for package TM2 (TURBOMOLE engine).

Target-3 defect under test ("an unconverged run is never reported as a
result / success criteria for the run are CHECKED after it finishes"):

    delfin/tools/adapters/turbomole.py:335 -- TurbomoleOptAdapter.execute
        success = proc.returncode == 0

`jobex` returns 0 when its driver script exits cleanly even if the SCF or
the geometry optimisation never converged (no `energy` file, jobex stopped
at --maxcycles without satisfying the convergence criteria). The adapter
reports that run as SUCCESS with empty data, so a user receives a "result"
for a calculation that produced none.

The physical-real contract (wave-20 priorities #1, #3; license rules) is
that an unconverged run is FAILED, never SUCCESS.

These tests are synthetic by design (TURBOMOLE licence rule: no actual
TURBOMOLE output in the repo): a fake `jobex` binary writes the artefacts
the parser reads. No real TURBOMOLE binary or output is involved.
"""

from __future__ import annotations

import os
import stat
from pathlib import Path

import pytest

from delfin.tools.adapters.turbomole import TurbomoleOptAdapter
from delfin.tools._types import ErrorKind, StepStatus

# Minimal valid XYZ that the adapter's manual coord fallback parses.
XYZ = """3

O 0.0 0.0 0.1173
H 0.0 0.7572 -0.4692
H 0.0 -0.7572 -0.4692
"""


def _make_fake_jobex(tmp_path: Path, energy: bool) -> Path:
    """Write a fake `jobex` on PATH that exits 0.

    energy=False -> the fake jobex produces NO `energy` file, exactly what
    a real jobex of an unconverged SCF leaves behind; it still exits 0.
    energy=True -> it writes a real two-column TURBOMOLE `energy` file
    (``<iteration>  <energy>``, the format delfin.tools.adapters.turbomole
    _extract_energy_tm parses -- `len(parts) >= 2`, value in column 1).
    """
    bin_dir = tmp_path / "bin"
    bin_dir.mkdir()
    script = bin_dir / "jobex"
    body = "#!/bin/sh\n"
    if energy:
        body += 'echo "1   -2.891329348" > energy\n'
    body += "exit 0\n"
    script.write_text(body)
    script.chmod(script.stat().st_mode | stat.S_IEXEC)
    return bin_dir


def _run_opt(monkeypatch: pytest.MonkeyPatch, tmp_path: Path,
             use_energy: bool) -> object:
    """Drive the real public path: TurbomoleOptAdapter.execute()."""
    work_dir = tmp_path / "work"
    work_dir.mkdir()
    geom = work_dir / "mol.xyz"
    geom.write_text(XYZ)

    fake_bin = _make_fake_jobex(tmp_path, energy=use_energy)
    old_path = os.environ.get("PATH", "")
    monkeypatch.setenv("PATH", f"{fake_bin}:{old_path}")

    adapter = TurbomoleOptAdapter()
    return adapter.execute(work_dir, geometry=geom, charge=0, mult=1)


def test_unconverged_jobex_is_failed_not_success(monkeypatch, tmp_path):
    """A jobex that exits 0 but produced no energy (unconverged) is FAILED.

    Red on the current code: the adapter's only success signal is
    `proc.returncode == 0`, so today it reports SUCCESS for this run.
    """
    result = _run_opt(monkeypatch, tmp_path, use_energy=False)
    assert result.status is StepStatus.FAILED
    assert not result.ok


def test_converged_jobex_is_success(monkeypatch, tmp_path):
    """Guard: a jobex that did produce an energy file is still SUCCESS.

    Ensures the convergence check does not reject a genuinely finished
    run (no false negative introduced by the fix).
    """
    result = _run_opt(monkeypatch, tmp_path, use_energy=True)
    assert result.status is StepStatus.SUCCESS
    assert result.ok


def test_unconverged_run_reports_convergence_error_kind(monkeypatch, tmp_path):
    """An unconverged run is diagnosable: FAILED with error_kind CONVERGENCE.

    Red on the current code: the adapter reports error_kind NONE (and
    SUCCESS) for a run whose SCF never converged, so a downstream caller
    or recovery policy cannot tell "converged" from "no result".
    """
    result = _run_opt(monkeypatch, tmp_path, use_energy=False)
    assert result.status is StepStatus.FAILED
    assert result.error_kind is ErrorKind.CONVERGENCE
