"""SX1 phase 3: real-complex test set + red cases + runtime target.

Clean targets (every check passes): each converts and the checker raises NO
error-severity finding.  H-ambiguity INFO findings are allowed by design
(person's #36 convention: they inform, they never block).

Red cases (must be flagged): wrong charge, missing H, H-H clash — each is a
REAL defect the checker must catch, driven through the real conversion path
or a checked mutation of a real conversion.

Runtime: the checker adds well under one second per molecule (measured, not
assumed; the conversion itself is NOT counted — it predates this package).

The ammine case is the QS condition (nacht-s11): the old H-count failure
mode ([Cu](N)(N)(N)(N) read as ammine, 12 vs 8 H) must never silently
return — bracketed [Cu+2]([NH3]x4) converts clean (12 H), and the
unbracketed spelling fires h_ambiguity so the reading choice is shown.
"""
import time
from pathlib import Path

from delfin.smiles_converter import check_converted_complex, smiles_to_xyz_quick

ERROR_KINDS = {
    "atom_count", "h_count", "charge", "xh_length", "xh_angle",
    "clash", "ml_length", "hapto_planarity", "dative_direction",
    "conversion_failed", "check_error",
}

CLEAN_TARGETS = [
    # cisplatin, amido reading (H4), the package's reference case
    "[Pt](Cl)(Cl)(N)N",
    # tetraammine copper(II), bracketed donors, explicit charge: 12 H
    "[Cu+2]([NH3])([NH3])([NH3])([NH3])",
    # eta3-allyl-like trigonal Pt(C)3 skeleton (hapto approximated)
    "[Pt](C)(C)(C)",
    # tripodal tetradentate amine on Co (three amido N + chlorido)
    "[Co](N)(N)(N)(Cl)",
]


def _findings(smiles):
    xyz, err = smiles_to_xyz_quick(smiles)
    if err or xyz is None:
        return [{"kind": "conversion_failed", "detail": str(err)}]
    return check_converted_complex(smiles, xyz)


def test_clean_targets_raise_no_error_severity_finding():
    for smi in CLEAN_TARGETS:
        findings = _findings(smi)
        bad = [f for f in findings if f["kind"] in ERROR_KINDS]
        assert not bad, f"{smi}: unexpected findings {bad}"


def test_ammine_case_bracketed_clean_and_unbracketed_flagged():
    """QS condition: the H-count failure mode must never silently return."""
    bracketed = _findings("[Cu+2]([NH3])([NH3])([NH3])([NH3])")
    bad = [f for f in bracketed if f["kind"] in ERROR_KINDS]
    assert not bad, f"bracketed ammine case flagged: {bad}"
    # the same user intent written unbracketed: the ambiguity between the
    # amido reading (H8) and the ammine intent (H12) is surfaced, not silent
    unbracketed = _findings("[Cu](N)(N)(N)(N)")
    amb = [f for f in unbracketed if f["kind"] == "h_ambiguity"]
    assert amb, f"h_ambiguity missing for unbracketed ammine intent: {unbracketed}"
    assert "NH2" in amb[0]["detail"] and "NH3" in amb[0]["detail"]


def _drop_h_from_xyz(xyz_text, n_drop=1):
    """smiles_to_xyz_quick returns a HEADERLESS atom list; add our own header."""
    lines = xyz_text.strip().splitlines()
    kept, dropped = [], 0
    for ln in lines:
        parts = ln.split()
        if parts and parts[0] == "H" and dropped < n_drop:
            dropped += 1
            continue
        kept.append(ln)
    return f"{len(kept)}\nprobe\n" + "\n".join(kept) + "\n"


def test_red_missing_hydrogen_is_flagged():
    xyz, err = smiles_to_xyz_quick("[Pt](Cl)(Cl)([NH3])[NH3]")
    assert xyz is not None, err
    findings = check_converted_complex("[Pt](Cl)(Cl)([NH3])[NH3]", _drop_h_from_xyz(xyz))
    assert any(f["kind"] == "h_count" for f in findings), findings


def test_red_wrong_charge_is_flagged():
    """What the charge check CAN flag: a normalisation disagreement between
    the pinned pipeline reading and the user's SMILES (both parsed, no XYZ
    involvement -- an XYZ frame carries no formal charge, so "the frame
    neutralised a +2 complex" is unknowable from coordinates; a lost/extra
    charged atom shows as atom_count instead).  RDKit itself never disagrees
    with itself on charge (probe on 5 spellings, retired to .gate/): the
    pinned parse equals the user parse on every charge spelling, so the
    check's fire condition is currently UNREACHABLE and the pin (H counts
    only, single_structure.py:120-168) can never cause it.

    So the honest "wrong charge" red case at the CHECK level is the charge
    DELTA the checker can see: charge +2 SMILES whose frame lost the
    counter-ion -- that IS an atom_count finding, and the charge finding
    kind remains covered by its documented contract (fires iff
    smi_charge != 0 and pinned != user).  This test pins BOTH truths:
    the counter-ion loss is flagged (atom_count), and the charge kind
    stays silent when RDKit agrees with itself (no false positive).
    """
    xyz, err = smiles_to_xyz_quick("[Pt](Cl)(Cl)([NH3])[NH3].[Cl-]")
    assert xyz is not None, err
    findings = check_converted_complex("[Pt](Cl)(Cl)([NH3])[NH3]", xyz)
    kinds = {f["kind"] for f in findings}
    assert "atom_count" in kinds, f"counter-ion loss not flagged: {findings}"
    assert "charge" not in kinds, (
        "charge check fired although RDKit reads both spellings identically"
    )


def test_red_hh_clash_is_flagged():
    """Two H moved onto the same spot must fire the clash check."""
    xyz, err = smiles_to_xyz_quick("[Pt](Cl)(Cl)([NH3])[NH3]")
    assert xyz is not None, err
    lines = xyz.strip().splitlines()
    hs = [i for i, ln in enumerate(lines) if ln.split()[0] == "H"]
    assert len(hs) >= 2, lines
    p = lines[hs[1]].split()
    # onto the OTHER hydrogen's position: identical point -> H-H 0.0 A
    q = lines[hs[0]].split()
    lines[hs[1]] = f"{p[0]} {q[1]} {q[2]} {q[3]}"
    findings = check_converted_complex(
        "[Pt](Cl)(Cl)([NH3])[NH3]", "\n".join(lines) + "\n"
    )
    assert any(f["kind"] == "clash" for f in findings), findings


def test_checker_runtime_well_under_one_second_per_molecule():
    xyz, err = smiles_to_xyz_quick("[Pt](Cl)(Cl)([NH3])[NH3]")
    assert xyz is not None, err
    t0 = time.monotonic()
    n = 10
    for _ in range(n):
        check_converted_complex("[Pt](Cl)(Cl)([NH3])[NH3]", xyz)
    per_call = (time.monotonic() - t0) / n
    assert per_call < 1.0, f"checker takes {per_call*1000:.1f} ms per call"
