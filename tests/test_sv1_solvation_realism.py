"""SV1: solvation adapters declare and enforce physical realism limits.

Red tests for package SV1 phase 1. The two solvation adapters
(packmol_solvate, xtb_solvator) accepted any parameter combination
silently before the fix: a DFT method name in SOLVATOR, a multiplicity
of 0, an empty solvent string, a box that cannot hold the solute. Every
check here exists so a physically wrong input fails BEFORE compute is
spent.

Sources:
- Packmol: Martínez, Andrade, Birgin, Martínez, J. Comput. Chem. 2009,
  30, 2157-2164, DOI 10.1002/jcc.21224 (packing tolerance is a
  hard lower bound on inter-molecular distance; the box must hold the
  solute plus at least one solvent molecule).
- ORCA SOLVATOR: ORCA 6 manual, %SOLVATOR block — NSOLV counts explicit
  solvent molecules; the method must be an xTB Hamiltonian (GFN0/GFN1/
  GFN2-xTB) because SOLVATOR is an xTB feature.
- Spin multiplicity 2S+1 >= 1 (delfin/tools/_science.py:87-94 states
  the same limit).
"""

import textwrap
from pathlib import Path

import pytest

from delfin.tools.adapters.packmol import PackmolSolvateAdapter
from delfin.tools.adapters.xtb_solvator import XtbSolvatorAdapter


# ---------------------------------------------------------------------------
# xtb_solvator: declared limits
# ---------------------------------------------------------------------------

def test_xtb_solvator_rejects_dft_method_name():
    """SOLVATOR is an xTB feature — a DFT method string must be refused
    before an ORCA input is written (it would fail in ORCA, after
    minutes of queue time)."""
    adapter = XtbSolvatorAdapter()
    with pytest.raises(ValueError, match="method"):
        adapter.validate_params(solvent="water", method="B3LYP")


def test_xtb_solvator_accepts_declared_xtb_methods():
    adapter = XtbSolvatorAdapter()
    for method in ("GFN0-xTB", "GFN1-xTB", "GFN2-xTB"):
        adapter.validate_params(solvent="water", method=method)


def test_xtb_solvator_rejects_multiplicity_below_one():
    adapter = XtbSolvatorAdapter()
    with pytest.raises(ValueError, match="mult"):
        adapter.validate_params(solvent="water", mult=0)


def test_xtb_solvator_rejects_negative_charge_type_and_n_solvent():
    adapter = XtbSolvatorAdapter()
    with pytest.raises(ValueError, match="n_solvent"):
        adapter.validate_params(solvent="water", n_solvent=0)
    with pytest.raises(ValueError, match="charge"):
        adapter.validate_params(solvent="water", charge="zero")


def test_xtb_solvator_rejects_empty_or_nonstring_solvent():
    adapter = XtbSolvatorAdapter()
    with pytest.raises(ValueError, match="solvent"):
        adapter.validate_params(solvent="")
    with pytest.raises(ValueError, match="solvent"):
        adapter.validate_params(solvent=123)


# ---------------------------------------------------------------------------
# packmol_solvate: declared limits
# ---------------------------------------------------------------------------

def _write_solute(path: Path, extent_a: float) -> Path:
    """A 3-atom collinear solute of about `extent_a` Å length."""
    coords = textwrap.dedent(f"""\
        3
        line solute
        O 0.000000 0.000000 0.000000
        H {extent_a / 2:.6f} 0.000000 0.000000
        H {extent_a:.6f} 0.000000 0.000000
        """)
    path.write_text(coords, encoding="utf-8")
    return path


def test_packmol_rejects_nonpositive_n_solvent_and_box():
    adapter = PackmolSolvateAdapter()
    with pytest.raises(ValueError, match="n_solvent"):
        adapter.validate_params(solvent_xyz="x.xyz", n_solvent=0)
    with pytest.raises(ValueError, match="box_side"):
        adapter.validate_params(solvent_xyz="x.xyz", box_side=0.0)
    with pytest.raises(ValueError, match="tolerance"):
        adapter.validate_params(solvent_xyz="x.xyz", tolerance=0.0)


def test_packmol_rejects_box_smaller_than_solute(tmp_path):
    """A box that cannot hold the solute is a physically meaningless
    packing problem; Packmol would churn to its 300 s timeout and the
    run would fail with 'Packmol solvation failed' — the adapter must
    refuse it up front with the numbers in the message."""
    adapter = PackmolSolvateAdapter()
    solute = _write_solute(tmp_path / "solute.xyz", extent_a=15.0)
    result = adapter.execute(
        tmp_path,
        geometry=solute,
        solvent_xyz=str(tmp_path / "water.xyz"),
        n_solvent=50,
        box_side=10.0,
    )
    assert result.status.name == "FAILED"
    assert "box_side" in (result.error or "")


# ---------------------------------------------------------------------------
# R1: execute() must enforce the limits itself — the product call site
# (delfin/agent/campaign.py:373) calls execute() directly, validate_params
# is not in that path. Red tests below drive execute() with the SAME bad
# parameters the validate_params tests use.
# ---------------------------------------------------------------------------

def test_xtb_solvator_execute_rejects_dft_method(tmp_path):
    """A DFT method must not reach the ORCA input writer: the input file
    must never be written with !B3LYP for an xTB feature."""
    adapter = XtbSolvatorAdapter()
    solute = tmp_path / "solute.xyz"
    solute.write_text("1\nx\nO 0.0 0.0 0.0\n", encoding="utf-8")
    result = adapter.execute(
        tmp_path, geometry=solute, method="B3LYP", solvent="water")
    assert result.status.name == "FAILED"
    assert "xTB Hamiltonian" in (result.error or "")


def test_xtb_solvator_execute_rejects_bad_mult_and_n_solvent(tmp_path):
    adapter = XtbSolvatorAdapter()
    solute = tmp_path / "solute.xyz"
    solute.write_text("1\nx\nO 0.0 0.0 0.0\n", encoding="utf-8")
    for kwargs in ({"mult": 0}, {"n_solvent": 0}, {"solvent": "   "}):
        result = adapter.execute(tmp_path, geometry=solute, **kwargs)
        assert result.status.name == "FAILED", kwargs
        # An ORCA input must not be written for a refused combination
        assert not (tmp_path / "XTB_SOLVATOR.inp").exists(), kwargs


def test_xtb_solvator_execute_accepts_valid_params(tmp_path, monkeypatch):
    """Does-not-fire direction: with valid params the adapter proceeds —
    the ORCA input is written with the declared values (run_orca is
    stubbed; this tests the adapter, not ORCA)."""
    calls = {}

    def fake_run_orca(inp, out, working_dir=None):
        calls["inp"] = Path(inp).read_text(encoding="utf-8")

    monkeypatch.setattr(
        "delfin.orca.run_orca", fake_run_orca)
    adapter = XtbSolvatorAdapter()
    solute = tmp_path / "solute.xyz"
    solute.write_text("1\nx\nO 0.0 0.0 0.0\n", encoding="utf-8")
    result = adapter.execute(
        tmp_path, geometry=solute, solvent="water", charge=0, mult=1,
        n_solvent=12, method="GFN2-xTB")
    # With a stubbed run_orca no solvator xyz appears -> FAILED is the
    # honest outcome; what matters is the input carries the right physics.
    assert "GFN2-xTB" in calls["inp"]
    assert "%SOLVATOR NSOLV 12" in calls["inp"]
    assert "*xyz 0 1" in calls["inp"]


def test_packmol_execute_rejects_bad_params(tmp_path):
    """A refused combination must not reach the packmol wrapper."""
    adapter = PackmolSolvateAdapter()
    solute = tmp_path / "solute.xyz"
    solute.write_text("1\nx\nO 0.0 0.0 0.0\n", encoding="utf-8")
    for kwargs in ({"n_solvent": 0}, {"tolerance": -1.0}, {"box_side": 0}):
        result = adapter.execute(
            tmp_path, geometry=solute,
            solvent_xyz=str(solute), **kwargs)
        assert result.status.name == "FAILED", kwargs
        assert "invalid" in (result.error or ""), kwargs
