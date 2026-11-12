"""SV1 adversarial tests: validation must hold on the REAL call path.

Finding R1 (2026-10-25, reviewer nacht-s19): the phase-1 realism limits
live only in ``validate_params``, but the only real call site in the
product, ``delfin/agent/campaign.py:373``, calls ``adapter.execute(...)``
directly. Neither adapter's ``execute`` runs the validation first, so an
invalid method/mult/charge/n_solvent still reaches the ORCA input
writer (queue time spent on a guaranteed-failing input) or the packmol
driver (n_solvent=0 passed straight through).

These tests drive the adapters the way the product does — through
``execute()`` with the downstream driver (run_orca / packmol_wrapper.
solvate) monkeypatched to a recorder — so they measure the real path.
"""

from pathlib import Path
from typing import Any, Optional

import pytest

from delfin.tools.adapters.packmol import PackmolSolvateAdapter
from delfin.tools.adapters.xtb_solvator import XtbSolvatorAdapter


class _Recorder:
    """Stands in for run_orca / solvate and records that it was called."""

    def __init__(self) -> None:
        self.calls: list[dict[str, Any]] = []

    def __call__(self, *args: Any, **kwargs: Any) -> Any:
        self.calls.append({"args": args, "kwargs": kwargs})
        return None


@pytest.fixture()
def orca_recorder(monkeypatch: pytest.MonkeyPatch) -> _Recorder:
    rec = _Recorder()
    import delfin.orca as orca_mod

    monkeypatch.setattr(orca_mod, "run_orca", rec)
    return rec


@pytest.fixture()
def solvate_recorder(monkeypatch: pytest.MonkeyPatch) -> _Recorder:
    rec = _Recorder()
    import delfin.analysis_tools.packmol_wrapper as pw

    monkeypatch.setattr(pw, "solvate", rec)
    return rec


# ---------------------------------------------------------------------------
# xtb_solvator: execute() must validate before touching ORCA
# ---------------------------------------------------------------------------

def _xyz(tmp_path: Path) -> Path:
    p = tmp_path / "mol.xyz"
    p.write_text("2\nwater\nO 0.0 0.0 0.0\nH 0.9 0.0 0.0\n", encoding="utf-8")
    return p


def test_xtb_solvator_execute_rejects_dft_method_before_orca(
    tmp_path: Path, orca_recorder: _Recorder
) -> None:
    """A DFT method name must be refused in execute(), not written into
    the ORCA input (campaign.py:373 calls execute() directly)."""
    result = XtbSolvatorAdapter().execute(
        tmp_path, geometry=_xyz(tmp_path), solvent="water", method="B3LYP"
    )
    assert result.status.name == "FAILED"
    assert "method" in (result.error or "")
    assert orca_recorder.calls == [], "run_orca must not be called for an invalid method"


def test_xtb_solvator_execute_rejects_multiplicity_zero_before_orca(
    tmp_path: Path, orca_recorder: _Recorder
) -> None:
    result = XtbSolvatorAdapter().execute(
        tmp_path, geometry=_xyz(tmp_path), solvent="water", mult=0
    )
    assert result.status.name == "FAILED"
    assert "mult" in (result.error or "")
    assert orca_recorder.calls == []


def test_xtb_solvator_execute_rejects_nonpositive_n_solvent_before_orca(
    tmp_path: Path, orca_recorder: _Recorder
) -> None:
    result = XtbSolvatorAdapter().execute(
        tmp_path, geometry=_xyz(tmp_path), solvent="water", n_solvent=0
    )
    assert result.status.name == "FAILED"
    assert "n_solvent" in (result.error or "")
    assert orca_recorder.calls == []


def test_xtb_solvator_execute_rejects_whitespace_solvent(
    tmp_path: Path, orca_recorder: _Recorder
) -> None:
    result = XtbSolvatorAdapter().execute(
        tmp_path, geometry=_xyz(tmp_path), solvent="   "
    )
    assert result.status.name == "FAILED"
    assert orca_recorder.calls == []


# ---------------------------------------------------------------------------
# packmol_solvate: execute() must validate before touching packmol
# ---------------------------------------------------------------------------

def test_packmol_execute_rejects_zero_n_solvent_before_packmol(
    tmp_path: Path, solvate_recorder: _Recorder
) -> None:
    solute = tmp_path / "solute.xyz"
    solute.write_text("2\nsolute\nO 0.0 0.0 0.0\nH 0.9 0.0 0.0\n", encoding="utf-8")
    solvent = tmp_path / "water.xyz"
    solvent.write_text("3\nwater\nO 0.0 0.0 0.0\nH 0.9 0.0 0.0\nH -0.9 0.0 0.0\n", encoding="utf-8")
    result = PackmolSolvateAdapter().execute(
        tmp_path, geometry=solute, solvent_xyz=str(solvent), n_solvent=0
    )
    assert result.status.name == "FAILED"
    assert "n_solvent" in (result.error or "")
    assert solvate_recorder.calls == [], "solvate() must not be called with n_solvent=0"


def test_packmol_execute_rejects_nonpositive_box_and_tolerance_before_packmol(
    tmp_path: Path, solvate_recorder: _Recorder
) -> None:
    solute = tmp_path / "solute.xyz"
    solute.write_text("2\nsolute\nO 0.0 0.0 0.0\nH 0.9 0.0 0.0\n", encoding="utf-8")
    solvent = tmp_path / "water.xyz"
    solvent.write_text("3\nwater\nO 0.0 0.0 0.0\nH 0.9 0.0 0.0\nH -0.9 0.0 0.0\n", encoding="utf-8")
    for bad in ({"box_side": 0.0}, {"tolerance": -1.0}):
        result = PackmolSolvateAdapter().execute(
            tmp_path, geometry=solute, solvent_xyz=str(solvent), **bad
        )
        assert result.status.name == "FAILED", bad
        assert solvate_recorder.calls == [], bad


# ---------------------------------------------------------------------------
# Does-not-fire direction: valid params must still reach the driver
# ---------------------------------------------------------------------------

def test_xtb_solvator_execute_accepts_valid_params(
    tmp_path: Path, orca_recorder: _Recorder
) -> None:
    result = XtbSolvatorAdapter().execute(
        tmp_path, geometry=_xyz(tmp_path), solvent="water",
        method="GFN2-xTB", charge=0, mult=1, n_solvent=12,
    )
    # run_orca is monkeypatched, so the run "fails" to produce the
    # solvated xyz — but the ORCA input must have been handed over and
    # must contain exactly the validated parameters.
    assert len(orca_recorder.calls) == 1
    assert result.status.name == "FAILED"  # recorder produces no output file
    inp = tmp_path / "XTB_SOLVATOR.inp"
    text = inp.read_text(encoding="utf-8")
    assert "!GFN2-xTB ALPB(water)" in text
    assert "%SOLVATOR NSOLV 12 END" in text
    assert "*xyz 0 1" in text
