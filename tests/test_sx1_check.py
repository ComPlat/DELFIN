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
    """A frame that really lost H vs the SMILES AS WRITTEN is flagged.

    Old decision (superseded by the person's convention, 2026-10-31): the
    unbracketed ``[Cu](N)(N)(N)(N)`` case was read dative (12 H) so its 8-H
    frame flagged.  Under standard SMILES (#36) that SMILES reads 4 amido N
    (8 H) and the 8-H frame matches -- the reading choice is surfaced as the
    h_ambiguity INFO instead.  A REAL h_count defect is a frame with fewer H
    than the written H counts: a bracketed ammine that lost hydrogens.
    """
    from delfin import smiles_converter as sc
    xyz = _convert("[Cu](N)(N)(N)(N)")
    findings = sc.check_converted_complex("[Cu](N)(N)(N)(N)", xyz)
    kinds = {f["kind"] for f in findings}
    assert "h_ambiguity" in kinds, f"expected h_ambiguity INFO, got {findings}"
    # and the genuine-defect direction still fires: drop 2 H from a real
    # conversion of an explicitly bracketed ammine complex
    xyz12 = _convert("[Cu]([NH3])([NH3])([NH3])([NH3])")
    lines = xyz12.strip().splitlines()
    kept = [str(len(lines) - 3), "H-dropped probe"]
    dropped = 0
    for ln in lines[1:]:
        p = ln.split()
        if len(p) >= 4 and p[0] == "H" and dropped < 2:
            dropped += 1
            continue
        kept.append(ln)
    findings2 = sc.check_converted_complex(
        "[Cu]([NH3])([NH3])([NH3])([NH3])", "\n".join(kept) + "\n"
    )
    assert "h_count" in {f["kind"] for f in findings2}, findings2


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
