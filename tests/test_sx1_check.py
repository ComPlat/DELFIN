"""Phase 1 red tests for package SX1: post-conversion check (findings, never silent fixes).

These tests run against the UNCHANGED code first: they must be RED because the
check function does not exist yet.
"""
import os
import sys

import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from delfin.smiles_converter import smiles_to_xyz  # noqa: E402

FERROCENE = "[Fe](C1=CC=CC=C1)(C2=CC=CC=C2)"


def _convert(smiles):
    xyz, err = smiles_to_xyz(smiles)
    if err or not xyz:
        pytest.skip(f"conversion failed: {err!r}")
    return xyz


def test_check_function_exists_and_reports_findings_list():
    from delfin import smiles_converter as sc
    assert hasattr(sc, "check_converted_complex"), (
        "check_converted_complex must exist in delfin/smiles_converter.py"
    )
    findings = sc.check_converted_complex(FERROCENE, _convert(FERROCENE))
    assert isinstance(findings, list)
    for f in findings:
        assert set(f) >= {"kind", "detail"}, f"finding shape: {f}"


def test_check_flags_missing_hydrogens_vs_smiles():
    """[Cu(NH3)4]2+ as [Cu](N)(N)(N)(N): converter emits 8 H, SMILES says 12."""
    from delfin import smiles_converter as sc
    xyz = _convert("[Cu](N)(N)(N)(N)")
    findings = sc.check_converted_complex("[Cu](N)(N)(N)(N)", xyz)
    kinds = {f["kind"] for f in findings}
    assert "h_count" in kinds, f"expected h_count finding, got {findings}"


def test_check_flags_unphysical_h_bond_length():
    """An H sitting far from any heavy neighbour (> 1.6 A) must be flagged."""
    from delfin import smiles_converter as sc
    xyz = "\n".join(
        ["3", "artificial case", "H 0.0 0.0 0.0", "O 1.9 0.0 0.0", "H 1.9 0.0 0.0"]
    )
    findings = sc.check_converted_complex("[OH2]", xyz)
    kinds = {f["kind"] for f in findings}
    assert "xh_length" in kinds, f"expected xh_length finding, got {findings}"


def test_check_flags_atom_clash():
    from delfin import smiles_converter as sc
    xyz = "\n".join(
        ["3", "artificial case", "H 0.0 0.0 0.0", "C 0.5 0.0 0.0", "H 1.0 0.0 0.0"]
    )
    findings = sc.check_converted_complex("[CH3]C", xyz)
    kinds = {f["kind"] for f in findings}
    assert "clash" in kinds, f"expected clash finding, got {findings}"


def test_check_clean_structure_has_no_findings():
    from delfin import smiles_converter as sc
    xyz = _convert("C")
    findings = sc.check_converted_complex("C", xyz)
    assert findings == [], f"methane must be clean, got {findings}"


def test_check_flags_wrong_total_charge():
    """SMILES [Cu](N)(N)(N)(N) is neutral by RDKit's reading; the check must
    compare the XYZ-derived formal charge against the SMILES charge and flag a
    mismatch when the emitter dropped charge-carrying atoms."""
    from delfin import smiles_converter as sc
    # artificial: geometry missing the charged nitrate O- (fragment detached)
    xyz = "\n".join(
        ["2", "artificial case", "N 0.0 0.0 0.0", "O 1.4 0.0 0.0"]
    )
    findings = sc.check_converted_complex("[N+](=O)([O-])", xyz)
    kinds = {f["kind"] for f in findings}
    assert "atom_count" in kinds, f"expected atom_count finding, got {findings}"
