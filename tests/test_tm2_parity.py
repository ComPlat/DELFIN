"""Engine-parity defect for the Turbomole freq adapter (package TM2, target 3).

Defect TM2-T3-A: ``TurbomoleFreqAdapter.execute`` runs ``aoforce`` through
``_run_tm_command`` with ``capture_output`` and never persists the output to
``work_dir/aoforce.out``. The workflow's IMAG logic (``_count_imaginary`` in
w21_turbomole_workflow.py) reads ``work_dir / "aoforce.out"`` — a file the
adapter never writes — so imaginary modes are never surfaced in the real
path. A non-minimum therefore comes back as a clean SUCCESS with no
``imaginary_frequencies`` data, violating target 3 ("aoforce imaginary modes
via DELFIN's IMAG logic").

This drives the real registered adapter (the public path) with a synthetic
fake ``aoforce``/``ridft`` on PATH — no real TURBOMOLE engine and no
TURBOMOLE-generated output in the repo (licence rule honored).
"""

from __future__ import annotations

import os
import stat
from pathlib import Path

from delfin.tools._registry import get
from delfin.tools._types import ErrorKind


def _make_bin(bin_dir: Path, name: str, body: str) -> Path:
    p = bin_dir / name
    p.write_text(body)
    p.chmod(p.stat().st_mode | stat.S_IEXEC)
    return p


def _env(tmp_path, monkeypatch):
    work_dir = tmp_path / "work"
    work_dir.mkdir()
    bin_dir = tmp_path / "fakebin"
    bin_dir.mkdir()
    aoforce = _make_bin(
        bin_dir, "aoforce",
        "#!/bin/sh\n"
        "printf '%s\\n' '$title aoforce' \\\n"
        "  '      1   -123.45i     0.3' \\\n"
        "  '      2    456.78      ostic' \\\n"
        "  '      3   1234.56      not(ating)' \\\n"
        "  '$end'\n"
        "printf '1  -1.0\\n$end\\n' > energy\n"
        "printf 'synthetic hessian\\n' > hessian\n",
    )
    _make_bin(bin_dir, "ridft",
              "#!/bin/sh\nprintf '1  -1.0\\n$end\\n' > energy\n")
    monkeypatch.setenv(
        "PATH", str(bin_dir) + os.pathsep + os.environ.get("PATH", ""))
    return work_dir, aoforce


def test_freq_persists_aoforce_output(tmp_path, monkeypatch):
    # Target 3: DELFIN's IMAG logic reads the aoforce output. The adapter
    # must leave it on disk, not swallow it into a CompletedProcess.
    work_dir, _ = _env(tmp_path, monkeypatch)
    result = _run("turbomole_freq", work_dir, charge=0)
    assert result.status.name == "SUCCESS", result.error
    out = work_dir / "aoforce.out"
    assert out.is_file(), "aoforce output must be persisted for IMAG logic"
    assert "-123.45i" in out.read_text()


def test_freq_surfaces_imaginary_frequencies_in_data(tmp_path, monkeypatch):
    # A synthetic aoforce table with one imaginary mode must be counted and
    # surfaced in the adapter data (ORCA parity: extract_imaginary_frequencies),
    # never hidden behind a clean SUCCESS.
    work_dir, _ = _env(tmp_path, monkeypatch)
    result = _run("turbomole_freq", work_dir, charge=0)
    assert result.status.name == "SUCCESS", result.error
    assert result.data.get("imaginary_frequencies") == 1, (
        f"one imaginary mode expected in data, got "
        f"{result.data.get('imaginary_frequencies')!r} (data={result.data})")


def test_freq_zero_imaginary_is_a_clean_minimum(tmp_path, monkeypatch):
    # A run with no imaginary modes reports 0 and stays SUCCESS (a converged
    # minimum, not an error) -- parity with ORCA, where imag-freq is a datum.
    work_dir, aoforce = _env(tmp_path, monkeypatch)
    aoforce.write_text(
        "#!/bin/sh\n"
        "printf '%s\\n' '      1  120.00      0.1' \\\n"
        "  '      2  456.78      os' \\\n"
        "  '$end'\n"
        "printf '1  -1.0\\n$end\\n' > energy\n"
        "printf 'synthetic hessian\\n' > hessian\n")
    aoforce.chmod(aoforce.stat().st_mode | stat.S_IEXEC)
    result = _run("turbomole_freq", work_dir, charge=0)
    assert result.status.name == "SUCCESS", result.error
    assert result.data.get("imaginary_frequencies") == 0


def test_freq_aoforce_failure_is_not_a_result(tmp_path, monkeypatch):
    # aoforce exits 0 but writes no hessian: the declared product is absent --
    # fail honestly (reviewer finding F4), never SUCCESS with partial data.
    # The fake writes energy + aoforce.out but drops the hessian file.
    work_dir, aoforce = _env(tmp_path, monkeypatch)
    aoforce.write_text(
        "#!/bin/sh\n"
        "printf '%s\\n' '$title aoforce' \\\n"
        "  '      1   -123.45i     0.3' \\\n"
        "  '$end'\n"
        "printf '1  -1.0\\n$end\\n' > energy\n"
        # deliberately NO hessian file
    )
    aoforce.chmod(aoforce.stat().st_mode | stat.S_IEXEC)
    result = _run("turbomole_freq", work_dir, charge=0)
    assert result.status.name == "FAILED"
    assert "hessian" in result.error.lower(), result.error
    assert result.error_kind == ErrorKind.PARSE


def _run(name: str, work_dir: Path, **kwargs):
    adapter = get(name)
    assert adapter is not None, f"adapter {name} not registered"
    return adapter.execute(work_dir, geometry=_xyz(work_dir), **kwargs)


def _xyz(work_dir: Path) -> Path:
    p = work_dir / "mol.xyz"
    p.write_text(
        "3\nwater\nO 0.0 0.0 0.0\nH 0.9 0.0 0.0\nH -0.3 0.8 0.0\n")
    return p
