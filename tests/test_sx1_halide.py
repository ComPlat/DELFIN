"""F1 (halide phantom-H): the checker's H expectation must not give a halide
bound to a metal a phantom hydrogen.

Semantics (QS arbitration, nacht-s11): the dative-donor reading (default
valence minus non-metal bonds) applies ONLY to unbracketed metal neighbours
carrying implicit H, mirroring the pin's own guard
(delfin/manta/single_structure.py:161-163, ``if n_h <= 0: continue``).  A
halide (n_h == 0) is read by its own H count, so no phantom H.

For cisplatin ``[Pt](Cl)(Cl)(N)N`` the expectation is therefore 6 H
(2 x NH3 dative donors), while the converter emits 4 -- a REAL h_count
finding, and the detail must say 6, never the old phantom-8 reading.
"""
from delfin.smiles_converter import check_converted_complex, smiles_to_xyz_quick


def test_cisplatin_halide_no_phantom_h_in_detail():
    """The detail says 6 implied H (2 x NH3), not 8 (with phantom halide H)."""
    xyz, err = smiles_to_xyz_quick("[Pt](Cl)(Cl)(N)N")
    assert xyz is not None, err
    findings = check_converted_complex("[Pt](Cl)(Cl)(N)N", xyz)
    h = [f for f in findings if f["kind"] == "h_count"]
    assert len(h) == 1, f"expected exactly one h_count finding, got {h}"
    detail = h[0]["detail"]
    assert "implies 6 H" in detail, f"phantom-H reading in detail: {detail}"


def test_missing_h_still_flagged():
    """A genuinely lost hydrogen is still flagged (no over-correction)."""
    xyz, err = smiles_to_xyz_quick("[Pt](Cl)(Cl)(N)N")
    assert xyz is not None, err
    lines = xyz.strip().splitlines()
    kept = [str(len(lines) - 1), "XYZ probe"]
    dropped = False
    for ln in lines:
        p = ln.split()
        if len(p) < 4 or not p[0].isalpha():
            continue
        if not dropped and p[0] == "H":
            dropped = True
            continue
        kept.append(ln)
    findings = check_converted_complex("[Pt](Cl)(Cl)(N)N", "\n".join(kept) + "\n")
    assert any(f["kind"] == "h_count" for f in findings), findings


def test_ammonia_bracketed_cisplatin_no_finding():
    """Explicit [NH3] donors convert with all H and raise no h_count finding."""
    xyz, err = smiles_to_xyz_quick("[Pt](Cl)(Cl)([NH3])[NH3]")
    assert xyz is not None, err
    findings = check_converted_complex("[Pt](Cl)(Cl)([NH3])[NH3]", xyz)
    assert [f for f in findings if f["kind"] == "h_count"] == [], findings
