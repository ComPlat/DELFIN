"""F1 (halide phantom-H): the checker's H expectation must not give a halide
bound to a metal a phantom hydrogen.

Semantics (person's convention decision, 2026-10-31, superseding the phase-2a
dative reading): the SMILES reads as WRITTEN under standard OpenSMILES (#36) --
an unbracketed N/O on a metal is AMIDO/HYDROXO (the M-L bond counts against
the donor's valence, i.e. H4 for ``[Pt](Cl)(Cl)(N)N``); a dative arrow that
ORIGINATES at the donor (``N->[Pt]``) shifts to the donor reading (NH3).  The
default-valence donor rule fires only on such a donor-originated arrow; a
halide (no implicit H) is read by its own H count, so no phantom H.

Consequence of the convention: cisplatin ``[Pt](Cl)(Cl)(N)N`` reads H4 and the
converter emits 4 -- NO h_count finding any more.  The real remaining defect
(Converter loses 2 H for the ammine-donor intent, [Cu](N)... = 8 vs 12) is
reported by the h_ambiguity INFO, which names both readings; the h_count
finding fires when the XYZ disagrees with the SMILES as written.
"""
from delfin.smiles_converter import check_converted_complex, smiles_to_xyz_quick


def test_cisplatin_halide_no_phantom_h_in_detail():
    """Amido reading (#36): H4 SMILES vs 4 H emitted -> no h_count finding.

    Old decision (superseded): the dative reading expected 6 and misflagged a
    correct H4 frame.  The ambiguity between the two readings is now carried
    by the h_ambiguity INFO finding instead of a wrong h_count.
    """
    xyz, err = smiles_to_xyz_quick("[Pt](Cl)(Cl)(N)N")
    assert xyz is not None, err
    findings = check_converted_complex("[Pt](Cl)(Cl)(N)N", xyz)
    assert [f for f in findings if f["kind"] == "h_count"] == [], findings
    # and the ambiguity is named for the user (amido NH2 vs ammine [NH3])
    amb = [f for f in findings if f["kind"] == "h_ambiguity"]
    assert len(amb) == 1 and "NH2" in amb[0]["detail"] and "NH3" in amb[0]["detail"]


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
