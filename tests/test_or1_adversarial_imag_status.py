"""Adversarial OR1: IMAG adapter must not report SUCCESS while modes remain.

delfin/tools/adapters/imag.py:82-87 builds its StepResult with
StepStatus.SUCCESS unconditionally, even when eliminate_imaginary_modes
returned ImagResult.remaining != [] (the "still a saddle after
IMAG_max_rounds" case, delfin/imag.py:487-489).  A pipeline consumer that
reads only the step status is then told the structure is a clean minimum
when it is not.  The data dict carries n_imaginary, but nothing forces a
consumer to read it; status is the contract.

Red test: with modes remaining, the adapter must NOT report SUCCESS.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from delfin.tools.adapters.imag import ImagFixAdapter


class _FakeImagResult:
    """Shape-compatible stand-in for delfin.imag.ImagResult."""

    def __init__(self, remaining):
        self.label = "imag_fix"
        self.rounds = 2
        self.remaining = remaining
        self.reason = "still a saddle after IMAG_max_rounds=2"


def _fake_hess_calc(tmp_path: Path) -> Path:
    """Write a hess file with sibling .inp/.out, as the ORCA adapter does."""
    (tmp_path / "calc.inp").write_text("! B3LYP def2-SVP Freq\n* xyz 0 1\n")
    (tmp_path / "calc.out").write_text("fake output\n")
    (tmp_path / "calc.hess").write_text("$hessian\n")
    return tmp_path / "calc.hess"


def test_adapter_reports_failure_when_imag_modes_remain(tmp_path, monkeypatch):
    """n_imaginary > 0 after the IMAG budget must not be StepStatus.SUCCESS."""
    hess = _fake_hess_calc(tmp_path)

    monkeypatch.setattr(
        "delfin.imag.eliminate_imaginary_modes",
        lambda **kwargs: _FakeImagResult(remaining=[(0, -123.4), (1, -45.6)]),
    )

    adapter = ImagFixAdapter()
    work_dir = tmp_path / "imag_work"
    result = adapter.execute(
        work_dir, geometry=None, cores=1,
        hess_file=str(hess), config={}, step_name="imag_fix",
    )

    assert result.status != __import__(
        "delfin.tools._types", fromlist=["StepStatus"]
    ).StepStatus.SUCCESS, (
        "IMAG adapter reported SUCCESS while 2 imaginary modes remain; "
        "a status-only consumer is told the structure is a clean minimum"
    )
    assert "n_imaginary" in (result.data or {})
    assert result.data["n_imaginary"] == 2


def test_adapter_reports_success_when_clean(tmp_path, monkeypatch):
    """Ordinary action still works: zero remaining modes stays SUCCESS."""
    hess = _fake_hess_calc(tmp_path)

    monkeypatch.setattr(
        "delfin.imag.eliminate_imaginary_modes",
        lambda **kwargs: _FakeImagResult(remaining=[]),
    )

    adapter = ImagFixAdapter()
    work_dir = tmp_path / "imag_work_clean"
    result = adapter.execute(
        work_dir, geometry=None, cores=1,
        hess_file=str(hess), config={}, step_name="imag_fix",
    )

    assert result.status == __import__(
        "delfin.tools._types", fromlist=["StepStatus"]
    ).StepStatus.SUCCESS
