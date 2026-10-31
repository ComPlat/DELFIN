"""Ferrocene ships FLAGGED (QS s11 ruling, 2026-10-31): the eta5 Cp embed
does not place the rings correctly today, and the fix is a feature, not a
hunk (measured: 12 C from a 10-C SMILES, Fe-C spread 1.78-4.60 A with four
beyond 4 A, UFF "Unrecognized atom type: Fe5+2" log evidence — probe
retired to .gate/).  This test pins the STANDING STATE so any future embed
fix turns it green; the ml_length tolerance is NOT loosened here.
Parsable spelling; the dotted aromatic spelling fails RDKit parse.
"""
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from delfin.smiles_converter import (  # noqa: E402
    check_converted_complex,
    smiles_to_xyz,
)

SMI = "[Fe](C1=CC=CC=C1)(C2=CC=CC=C2)"


def _run():
    xyz, err = smiles_to_xyz(SMI)
    assert xyz is not None, f"ferrocene conversion failed: {err!r}"
    return check_converted_complex(SMI, xyz)


def test_ferrocene_is_flagged_with_ml_length_findings():
    findings = _run()
    ml = [f for f in findings if f["kind"] == "ml_length"]
    assert ml, f"ferrocene must stay flagged until the eta5 embed fix: {findings}"
    for f in ml:
        assert "detail" in f and "exceeds" in f["detail"], f
