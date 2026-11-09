"""Convergence guard for the Turbomole adapters (package TM2, target 3).

Defect TM2-1 (reviewer nacht-s19): ``TurbomoleOptAdapter.execute`` reported
SUCCESS for an UNCONVERGED jobex run — ``success = proc.returncode == 0``
was the only signal, so a jobex that exits 0 but produced no converged
energy (converged energy file) came back as StepResult(status=SUCCESS,
data={}, error=None, error_kind=NONE).

That violates target 3: "success criteria are checked after the run; an
unconverged run is never reported as a result." The same hole existed in
turbomole_sp and turbomole_freq.

This test drives the real registered adapters (the public path) with a
synthetic fake ``jobex``/``ridft`` on PATH — no real TURBOMOLE engine and no
TURBOMOLE-generated output in the repo (licence rule honored).
"""

from __future__ import annotations

import os
import stat
from pathlib import Path

import pytest

from delfin.tools._registry import get
from delfin.tools._types import ErrorKind


def _make_bin(bin_dir: Path, name: str, body: str) -> Path:
    p = bin_dir / name
    p.write_text(body)
    p.chmod(p.stat().st_mode | stat.S_IEXEC)
    return p


def _xyz(work_dir: Path) -> Path:
    p = work_dir / "mol.xyz"
    p.write_text("3\nwater\nO 0.0 0.0 0.0\nH 0.9 0.0 0.0\nH -0.3 0.8 0.0\n")
    return p


def _run(name: str, work_dir: Path, **kwargs):
    adapter = get(name)
    assert adapter is not None, f"adapter {name} not registered"
    return adapter.execute(work_dir, geometry=_xyz(work_dir), cores=1, **kwargs)


class TestOptConvergence:
    """jobex exiting 0 without a converged energy file is NOT success."""

    def _make_env(self, tmp_path, monkeypatch, jobex_body: str) -> Path:
        work_dir = tmp_path / "work"
        work_dir.mkdir()
        bin_dir = tmp_path / "fakebin"
        bin_dir.mkdir()
        # Shadow PATH so only our fakes are found (deterministic).
        _make_bin(bin_dir, "jobex", jobex_body)
        _make_bin(bin_dir, "ridft", "#!/bin/sh\nexit 0\n")
        _make_bin(bin_dir, "aoforce", "#!/bin/sh\nexit 0\n")
        monkeypatch.setenv("PATH", str(bin_dir) + os.pathsep + os.environ.get("PATH", ""))
        return work_dir

    def test_unconverged_jobex_is_failed_not_success(self, tmp_path, monkeypatch):
        # jobex exits 0 and deliberately writes NO energy file.
        work_dir = self._make_env(
            tmp_path, monkeypatch, "#!/bin/sh\nexit 0\n",
        )
        result = _run("turbomole_opt", work_dir, charge=0)
        assert result.status.name == "FAILED", (
            f"unconverged jobex must FAIL, got {result.status.name} "
            f"error={result.error!r} error_kind={result.error_kind}"
        )
        assert result.ok is False

    def test_unconverged_run_reports_convergence_error_kind(self, tmp_path, monkeypatch):
        work_dir = self._make_env(
            tmp_path, monkeypatch, "#!/bin/sh\nexit 0\n",
        )
        result = _run("turbomole_opt", work_dir, charge=0)
        assert result.error_kind == ErrorKind.CONVERGENCE, (
            f"expected error_kind CONVERGENCE, got {result.error_kind}"
        )
        # The message must say WHY, not just fail silently.
        assert result.error and "not converge" in result.error.lower()

    def test_converged_jobex_is_success(self, tmp_path, monkeypatch):
        # jobex exits 0 AND writes the converged energy file.
        work_dir = self._make_env(
            tmp_path, monkeypatch,
            f"#!/bin/sh\nprintf '1  -76.4\\n$end\\n' > {str(tmp_path / 'work' / 'energy')}\n",
        )
        result = _run("turbomole_opt", work_dir, charge=0)
        assert result.status.name == "SUCCESS", result.error
        assert result.data.get("energy_Eh") is not None


class TestSpConvergence:
    """ridft exiting 0 without energy is an unconverged SCF, not SUCCESS."""

    def test_unconverged_sp_is_failed_convergence(self, tmp_path, monkeypatch):
        work_dir = tmp_path / "work"
        work_dir.mkdir()
        bin_dir = tmp_path / "fakebin"
        bin_dir.mkdir()
        _make_bin(bin_dir, "ridft", "#!/bin/sh\nexit 0\n")  # no energy file
        _make_bin(bin_dir, "dscf", "#!/bin/sh\nexit 0\n")   # no energy file
        monkeypatch.setenv("PATH", str(bin_dir) + os.pathsep + os.environ.get("PATH", ""))
        result = _run("turbomole_sp", work_dir, charge=0, ri=False)  # dscf path
        assert result.status.name == "FAILED"
        assert result.error_kind == ErrorKind.CONVERGENCE
        # Neither path (ridft nor dscf) reports an empty SUCCESS.
        result2 = _run("turbomole_sp", work_dir, charge=0, ri=True)
        assert result2.status.name == "FAILED"
        assert result2.error_kind == ErrorKind.CONVERGENCE


def test_energy_file_is_written_by_ridft(monkeypatch, tmp_path):
    """Authoritative check that the converged fixture path does write energy."""
    work_dir = tmp_path / "work"
    work_dir.mkdir()
    bin_dir = tmp_path / "fakebin"
    bin_dir.mkdir()
    _make_bin(bin_dir, "ridft", f"#!/bin/sh\nprintf '1  -1.0\\n$end\\n' > {work_dir}/energy\n")
    monkeypatch.setenv("PATH", str(bin_dir) + os.pathsep + os.environ.get("PATH", ""))
    result = _run("turbomole_sp", work_dir, charge=0)
    assert result.status.name == "SUCCESS", result.error
    assert result.data["energy_Eh"] == -1.0


def test_multicore_run_sets_para_arch(monkeypatch, tmp_path):
    """A cores>1 run through the adapter path must set PARA_ARCH, not just
    PARNODES.  The w21 workflow calls the adapters directly (bypassing the
    dashboard submit path that exports TM_PARA_ARCH), so without PARA_ARCH
    here a multicore jobex/ridft never tells Turbomole's runtime to
    parallelise — the T2 target's \"PARA_ARCH, PARNODES\" is met only if the
    adapter itself provides the parallel-architecture variable.
    """
    work_dir = tmp_path / "work"
    work_dir.mkdir()
    bin_dir = tmp_path / "fakebin"
    bin_dir.mkdir()
    env_dump = tmp_path / "env.dump"
    # The fake records the parallel env the child actually inherits (double
    # quotes so $VAR expands to the child's real value at runtime).
    _make_bin(bin_dir, "ridft",
              f"#!/bin/sh\nprintf 'PARA_ARCH=%s\\nPARNODES=%s\\n' \"$PARA_ARCH\" \"$PARNODES\" > {env_dump}\n"
              "printf '1  -1.0\\n$end\\n' > energy\n")
    monkeypatch.setenv("PATH", str(bin_dir) + os.pathsep + os.environ.get("PATH", ""))
    adapter = get("turbomole_sp")
    assert adapter is not None
    result = adapter.execute(work_dir, geometry=_xyz(work_dir), cores=2)
    assert result.status.name == "SUCCESS", result.error
    dumped = env_dump.read_text()
    # The first line is PARA_ARCH=$PARA_ARCH; it must not be the empty default.
    para_arch_line = dumped.splitlines()[0] if dumped.splitlines() else ""
    assert para_arch_line == "PARA_ARCH=SMP", (
        f"PARA_ARCH must be set (got {para_arch_line!r}); env was:\n{dumped}")
    assert "PARNODES=2" in dumped, dumped


# ---------------------------------------------------------------------------
# Round-4 finding (operator real-run, 09.11.): the freq adapter reported
# imaginary_frequencies=3 for a true minimum (vibspectrum: 6 x 0.00,
# 1637.76, 3789.09, 3884.13 cm-1). The old counter read the WHOLE aoforce
# text and counted any line whose 2nd token parsed as a negative float --
# that hits the coordinate block ('0.00000000 -1.43078088 1.26069347 h')
# and the normal-mode eigenvector rows ('z -0.07021 ...'). The clean source
# is the frequency table: the $vibrational spectrum block (the vibspectrum
# file) whose column 3 is the wave number; TURBOMOLE prints imaginary
# modes there as negative numbers. Fixtures are hand-written in the
# parser's format (licence rule: no TURBOMOLE-generated output in the
# repo).
# ---------------------------------------------------------------------------

# Hand-written aoforce.out shape: coordinate rows, eigenvector rows with
# negatives, and a vibspectrum ($vibrational spectrum) table whose
# wave-number column is all positive -- a true minimum. The old counter
# reported 3 imaginary modes from the coordinate/eigenvector negatives.
_AOFORCE_MINIMUM_TEXT = """\
$hessian
symmetry    c1

     1       0.00000000      0.00000000      0.00000000  o
     2       0.00000000      1.43078088      1.26069347  h
     3       0.00000000     -1.43078088      1.26069347  h

 vibration    1    frequency   1637.76
  x         0.00000      y        -0.43249      z         0.00000
  x        -0.07021      y         0.43249      z         0.00000

 vibration    2    frequency   3789.09
  x         0.00000      y         0.11234      z         0.00000

$vibrational spectrum
       1        -1     1637.76
       2        -1     3789.09
       3        -1     3884.13
$end
"""

# Same file, but the vibspectrum table carries one negative wave number --
# one imaginary mode (TURBOMOLE prints imaginary modes there as negatives).
_AOFORCE_ONE_IMAGINARY_TEXT = """\
$vibrational spectrum
       1        -1      -512.34
       2        -1      1637.76
       3        -1      3789.09
$end
"""


class TestImaginaryCountReadsVibspectrumTable:
    def test_minimum_with_negative_eigenvectors_counts_zero(self, tmp_path):
        from delfin.tools.adapters.turbomole import _count_imaginary_tm

        assert _count_imaginary_tm(_AOFORCE_MINIMUM_TEXT) == 0, (
            "coordinate/eigenvector negatives must NOT be counted; only the "
            "vibspectrum wave-number column decides")

    def test_one_negative_wave_number_counts_one(self, tmp_path):
        from delfin.tools.adapters.turbomole import _count_imaginary_tm

        assert _count_imaginary_tm(_AOFORCE_ONE_IMAGINARY_TEXT) == 1

