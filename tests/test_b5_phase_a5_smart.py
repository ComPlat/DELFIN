"""Tests for Patch P-A5-SMART: Jahn-Teller / multiple-bond / cyclometal /
mer-tridentate aware Phase A.5 projection in the B5 v2 post-optimizer.

Validates that, when ``DELFIN_B5_PHASE_A5_SMART=1``:

* Cu(II) hexaaqua / Cu(II) Cl4 → classified as ``"jahn_teller_cu2"`` and
  Phase A.5 polyhedron projection is skipped, preserving the 4+2 elongation.
* Vanadyl V=O (geometric short bond) → classified as ``"metal_oxo"``.
* Cyclometal Ir(ppy)-like 5-ring → classified as ``"cyclometal_5ring"``.
* mer-N3 tridentate (collinear i-M-j pair with k perpendicular) →
  classified as ``"mer_tridentate"``.
* Regular Cr(NH3)6 octahedral → NOT classified (returns ``None``), so
  standard projection runs.
* ``DELFIN_B5_PHASE_A5_SMART=0`` (default) → bit-exact behaviour to the
  legacy P-A5-ENV gate (Phase A.5 disabled).
* ``DELFIN_B5_PHASE_A5_SMART=1`` + Cu(II) → distortion preserved (output
  geometry does NOT collapse to a symmetric polyhedron).
* Malformed / None mol → safe fallback (no crash, returns None).

All tests skip cleanly when RDKit or the B5 module are missing.
"""
from __future__ import annotations

from typing import List, Sequence, Tuple

import numpy as np
import pytest

_post_optimizer = pytest.importorskip(
    "delfin.manta._post_optimizer",
    reason="Baustein 5 post-optimizer module not yet implemented",
)
Chem = pytest.importorskip("rdkit.Chem", reason="RDKit required for tests")
_rdGeometry = pytest.importorskip("rdkit.Geometry", reason="RDKit required")

_classify_distortion_mandatory = getattr(
    _post_optimizer, "_classify_distortion_mandatory", None
)
_project_donors_to_polyhedron = getattr(
    _post_optimizer, "_project_donors_to_polyhedron", None
)
post_optimize_geometry = getattr(_post_optimizer, "post_optimize_geometry", None)
if (
    _classify_distortion_mandatory is None
    or _project_donors_to_polyhedron is None
    or post_optimize_geometry is None
):
    pytest.skip(
        "Smart Phase A.5 helpers not wired", allow_module_level=True
    )


Atom = Tuple[str, float, float, float]
Bond = Tuple[int, int]


# ---------------------------------------------------------------------------
# Fixture helpers
# ---------------------------------------------------------------------------


def _xyz(atoms: Sequence[Atom]) -> str:
    lines = [str(len(atoms)), "test"]
    for sym, x, y, z in atoms:
        lines.append(f"{sym:<3s} {x:14.8f} {y:14.8f} {z:14.8f}")
    return "\n".join(lines) + "\n"


def _build_mol(
    atoms: Sequence[Atom],
    bonds: Sequence[Bond],
    formal_charges: Sequence[Tuple[int, int]] = (),
    double_bonds: Sequence[Bond] = (),
):
    """Build an RDKit Mol with explicit topology + 3D positions.

    formal_charges: list of (atom_idx, charge) pairs.
    double_bonds:   list of (i, j) pairs that should be DOUBLE-bonded.
    """
    mol = Chem.RWMol()
    for sym, _, _, _ in atoms:
        a = Chem.Atom(sym)
        a.SetNoImplicit(True)
        mol.AddAtom(a)
    double_set = {tuple(sorted(b)) for b in double_bonds}
    for i, j in bonds:
        bt = (
            Chem.BondType.DOUBLE
            if tuple(sorted((i, j))) in double_set
            else Chem.BondType.SINGLE
        )
        mol.AddBond(i, j, bt)
    for idx, ch in formal_charges:
        mol.GetAtomWithIdx(int(idx)).SetFormalCharge(int(ch))
    conf = Chem.Conformer(len(atoms))
    for i, (_, x, y, z) in enumerate(atoms):
        conf.SetAtomPosition(
            i, _rdGeometry.Point3D(float(x), float(y), float(z))
        )
    mol.AddConformer(conf, assignId=True)
    return mol


def _coords_of(atoms: Sequence[Atom]) -> np.ndarray:
    return np.array([[float(x), float(y), float(z)] for _, x, y, z in atoms],
                    dtype=float)


# ---------------------------------------------------------------------------
# Fixture builders — chemically meaningful but minimal
# ---------------------------------------------------------------------------


def _cu2_hexaaqua_4plus2() -> Tuple[List[Atom], List[Bond], List[Tuple[int, int]]]:
    """Cu(II) hexaaqua with 4+2 Jahn-Teller elongation.

    Equatorial Cu-O at 2.00 Å along ±x and ±y; axial Cu-O elongated to
    2.50 Å along ±z. We model only the O donors (heavy atoms) — H of the
    aqua ligands is omitted, which is sufficient for the classifier.
    """
    eq = 2.00
    ax = 2.50
    atoms: List[Atom] = [
        ("Cu", 0.0, 0.0, 0.0),
        ("O",  eq,  0.0, 0.0),
        ("O", -eq,  0.0, 0.0),
        ("O",  0.0,  eq, 0.0),
        ("O",  0.0, -eq, 0.0),
        ("O",  0.0,  0.0,  ax),
        ("O",  0.0,  0.0, -ax),
    ]
    bonds: List[Bond] = [(0, i) for i in range(1, 7)]
    charges: List[Tuple[int, int]] = [(0, 2)]   # Cu2+
    return atoms, bonds, charges


def _cu2_cl4_square_planar() -> Tuple[List[Atom], List[Bond], List[Tuple[int, int]]]:
    """Cu(II) Cl4 distorted square planar (also JT-prone)."""
    atoms: List[Atom] = [
        ("Cu", 0.0, 0.0, 0.0),
        ("Cl",  2.40, 0.0, 0.0),
        ("Cl", -2.40, 0.0, 0.0),
        ("Cl",  0.0,  2.20, 0.0),
        ("Cl",  0.0, -2.20, 0.0),
        # 5th and 6th donors at long axial distance to make CN>=5
        ("Cl",  0.0,  0.0,  2.80),
        ("Cl",  0.0,  0.0, -2.80),
    ]
    bonds: List[Bond] = [(0, i) for i in range(1, 7)]
    charges: List[Tuple[int, int]] = [(0, 2)]
    return atoms, bonds, charges


def _vanadyl_vo() -> Tuple[List[Atom], List[Bond], List[Tuple[int, int]], List[Bond]]:
    """V=O vanadyl with one short M=O bond (1.60 Å) + 4 equatorial donors.

    The short V=O distance is the geometric trigger for the metal_oxo
    classification (< 1.85 Å threshold). Equatorial donors are normal
    V-N bonds at ~2.10 Å.
    """
    atoms: List[Atom] = [
        ("V",  0.0, 0.0, 0.0),
        ("O",  0.0, 0.0, 1.60),    # apical V=O — short
        ("N",  2.10, 0.0, 0.0),
        ("N", -2.10, 0.0, 0.0),
        ("N",  0.0,  2.10, 0.0),
        ("N",  0.0, -2.10, 0.0),
    ]
    bonds: List[Bond] = [(0, 1), (0, 2), (0, 3), (0, 4), (0, 5)]
    double_bonds: List[Bond] = [(0, 1)]   # explicit V=O
    charges: List[Tuple[int, int]] = [(0, 4)]   # V(IV) for VO^2+
    return atoms, bonds, charges, double_bonds


def _ir_cyclometal_5ring() -> Tuple[List[Atom], List[Bond], List[Tuple[int, int]]]:
    """Minimal Ir-(C^N) cyclometalation: 5-membered ring Ir-C-C-C-N closure.

    Layout (in xy-plane around Ir at origin):
        Ir — C_aryl — C — C — N — back to Ir.
    The C_aryl is the cyclometal site (sp2 aromatic carbon). N is the
    pyridyl donor. Add 4 more sigma donors to give Ir an octahedral CN=6.
    """
    atoms: List[Atom] = [
        ("Ir", 0.0,  0.0, 0.0),
        ("C",  2.00, 0.0, 0.0),    # aryl C (cyclometal site)
        ("C",  2.50, 1.20, 0.0),   # ring-internal C
        ("C",  1.50, 2.20, 0.0),   # ring-internal C
        ("N",  0.30, 1.80, 0.0),   # pyridyl N (closes the 5-ring with Ir)
        # 4 spectator donors
        ("N",  0.0, -2.10, 0.0),
        ("N",  0.0,  0.0,  2.10),
        ("N",  0.0,  0.0, -2.10),
        ("Cl", -2.30, 0.0, 0.0),
    ]
    bonds: List[Bond] = [
        (0, 1),     # Ir-C
        (1, 2),     # C-C (aryl)
        (2, 3),     # C-C
        (3, 4),     # C-N
        (0, 4),     # Ir-N — closes the 5-ring
        (0, 5), (0, 6), (0, 7), (0, 8),   # spectator donors
    ]
    charges: List[Tuple[int, int]] = [(0, 3)]   # Ir(III)
    # Mark the C-aryl atom aromatic so the classifier picks it up.
    return atoms, bonds, charges


def _ir_cyclometal_with_aromatic_flag(
    atoms: List[Atom], bonds: List[Bond],
    charges: List[Tuple[int, int]],
):
    """Build mol and mark the cyclometal-C aromatic + sp2."""
    mol = _build_mol(atoms, bonds, charges)
    # mark atom 1 (aryl C) as aromatic so the classifier's aromatic branch fires
    aryl_c = mol.GetAtomWithIdx(1)
    aryl_c.SetIsAromatic(True)
    return mol


def _mer_tridentate_n3() -> Tuple[List[Atom], List[Bond], List[Tuple[int, int]]]:
    """Co(III)-N3 mer-tridentate: planar tridentate ligand N-C-N-C-N backbone
    coordinating Co with N1-Co-N3 ~180° and middle N2 perpendicular. The
    three N donors are chemically linked through a chelate chain (mimicking
    terpyridine / pincer geometry). Pad to CN=6 with 3 spectator Cl donors
    off the mer plane.
    """
    atoms: List[Atom] = [
        ("Co", 0.0, 0.0, 0.0),
        ("N",  0.0, 0.0,  1.95),   # mer terminal N1 (along +z)
        ("N",  0.0, 0.0, -1.95),   # mer terminal N3 (along -z, trans to N1)
        ("N",  1.95, 0.0, 0.0),    # mer middle N2 (along +x)
        # Backbone carbons linking the three N donors (off-axis to keep
        # chemistry plausible but the backbone connects N1-Cb1-N2-Cb2-N3).
        ("C",  1.30, 0.0, 1.20),   # Cb1: bonded to N1 and N2
        ("C",  1.30, 0.0, -1.20),  # Cb2: bonded to N2 and N3
        ("Cl", -2.30, 0.0, 0.0),
        ("Cl", 0.0,  2.30, 0.0),
        ("Cl", 0.0, -2.30, 0.0),
    ]
    bonds: List[Bond] = [
        # Co-donor bonds (3 N + 3 Cl = CN 6)
        (0, 1), (0, 2), (0, 3),
        (0, 6), (0, 7), (0, 8),
        # Backbone chelate links (mer-tridentate ligand)
        (1, 4), (3, 4),    # N1-Cb1-N2
        (3, 5), (2, 5),    # N2-Cb2-N3
    ]
    charges: List[Tuple[int, int]] = [(0, 3)]
    return atoms, bonds, charges


def _cr_nh3_6() -> Tuple[List[Atom], List[Bond], List[Tuple[int, int]]]:
    """Regular Cr(III) hexaammine octahedral — NO mandatory distortion.

    Cr3+ is d^3, not Jahn-Teller; all N donors are spectator sp3; no rings;
    no near-trans triple. Classifier must return None.
    """
    eq = 2.10
    atoms: List[Atom] = [
        ("Cr", 0.0, 0.0, 0.0),
        ("N",  eq, 0.0, 0.0),
        ("N", -eq, 0.0, 0.0),
        ("N",  0.0,  eq, 0.0),
        ("N",  0.0, -eq, 0.0),
        ("N",  0.0,  0.0,  eq * 0.95),   # slight asymmetry to break mer triple
        ("N",  0.05, 0.0, -eq),
    ]
    bonds: List[Bond] = [(0, i) for i in range(1, 7)]
    charges: List[Tuple[int, int]] = [(0, 3)]
    return atoms, bonds, charges


# ============================================================================
# 1. Direct classifier unit tests
# ============================================================================


def test_classify_cu2_jahn_teller():
    """Cu(II) hexaaqua → ``jahn_teller_cu2``."""
    atoms, bonds, charges = _cu2_hexaaqua_4plus2()
    mol = _build_mol(atoms, bonds, charges)
    coords = _coords_of(atoms)
    donor_idxs = list(range(1, 7))
    label = _classify_distortion_mandatory(mol, coords, 0, donor_idxs)
    assert label == "jahn_teller_cu2", (
        f"Cu(II) with CN=6 must be classified as Jahn-Teller, got {label!r}"
    )


def test_classify_metal_oxo_short_bond():
    """V=O with d(V-O)=1.60 Å (geometric criterion) → ``metal_oxo``."""
    atoms, bonds, charges, double_bonds = _vanadyl_vo()
    mol = _build_mol(atoms, bonds, charges, double_bonds)
    coords = _coords_of(atoms)
    donor_idxs = [1, 2, 3, 4, 5]
    label = _classify_distortion_mandatory(mol, coords, 0, donor_idxs)
    assert label == "metal_oxo", (
        f"V=O short bond must be classified metal_oxo, got {label!r}"
    )


def test_classify_cyclometal_5ring():
    """Ir cyclometalation with 5-membered Ir-C-C-C-N ring → ``cyclometal_5ring``.

    The fixture has an explicit Ir-C aryl + Ir-N pyridyl forming a 5-ring;
    we mark the C as aromatic so the classifier's sp2/aromatic branch fires.
    """
    atoms, bonds, charges = _ir_cyclometal_5ring()
    mol = _ir_cyclometal_with_aromatic_flag(atoms, bonds, charges)
    coords = _coords_of(atoms)
    donor_idxs = [1, 4, 5, 6, 7, 8]   # all heavy donors bonded to Ir
    label = _classify_distortion_mandatory(mol, coords, 0, donor_idxs)
    assert label == "cyclometal_5ring", (
        f"Ir 5-ring cyclometal must be classified cyclometal_5ring, got {label!r}"
    )


def test_classify_mer_tridentate():
    """mer-N3 (180° + 90° + 90°) with chelate backbone → ``mer_tridentate``.

    Six donors total (3 mer-N + 3 spectator Cl); the three N's are linked
    by a backbone (N1-Cb1-N2-Cb2-N3) so they belong to the same chelate
    ligand fragment, satisfying the graph-connectivity requirement of the
    CN >= 4 mer detector.
    """
    atoms, bonds, charges = _mer_tridentate_n3()
    mol = _build_mol(atoms, bonds, charges)
    coords = _coords_of(atoms)
    # Donors are the 3 N (idx 1, 2, 3) + 3 Cl (idx 6, 7, 8). Atoms 4, 5
    # are backbone carbons, NOT bonded to the metal.
    donor_idxs = [1, 2, 3, 6, 7, 8]
    label = _classify_distortion_mandatory(mol, coords, 0, donor_idxs)
    assert label == "mer_tridentate", (
        f"mer-N3 must be classified mer_tridentate, got {label!r}"
    )


def test_classify_regular_octahedral_is_none():
    """Cr(NH3)6 — NO mandatory distortion → ``None``."""
    atoms, bonds, charges = _cr_nh3_6()
    mol = _build_mol(atoms, bonds, charges)
    coords = _coords_of(atoms)
    donor_idxs = list(range(1, 7))
    label = _classify_distortion_mandatory(mol, coords, 0, donor_idxs)
    assert label is None, (
        f"Cr(NH3)6 must NOT be classified as distortion-mandatory, got {label!r}"
    )


def test_classify_malformed_mol_safe():
    """``mol=None`` or empty donor list → ``None`` (no crash)."""
    coords = np.zeros((1, 3))
    assert _classify_distortion_mandatory(None, coords, 0, [1]) is None
    # empty donors
    mol = _build_mol(
        [("Cu", 0.0, 0.0, 0.0)],
        [],
        [(0, 2)],
    )
    assert _classify_distortion_mandatory(mol, np.zeros((1, 3)), 0, []) is None


# ============================================================================
# 2. Env-flag wiring + bit-exact default-OFF contract
# ============================================================================


def test_smart_env_default_off_bit_exact_to_legacy(monkeypatch):
    """``DELFIN_B5_PHASE_A5_SMART`` unset → bit-exact to legacy P-A5-ENV gate.

    With legacy env=0 (default), Phase A.5 is disabled regardless of param;
    smart env=0 must reproduce the same disabled output for the same input.
    """
    monkeypatch.delenv("DELFIN_B5_PHASE_A5_ENABLE", raising=False)
    monkeypatch.delenv("DELFIN_B5_PHASE_A5_SMART", raising=False)
    atoms, bonds, charges = _cu2_hexaaqua_4plus2()
    xyz = _xyz(atoms)
    mol_a = _build_mol(atoms, bonds, charges)
    mol_b = _build_mol(atoms, bonds, charges)

    # Smart-mode disabled: enable_symmetry=True must STILL be no-op (legacy gate).
    out_smart_off, _ = post_optimize_geometry(
        xyz, mol_a, class_label="sigma", max_iter=3,
        step_size=0.3, bond_tol=0.15, enable_symmetry=True,
    )
    out_legacy_off, _ = post_optimize_geometry(
        xyz, mol_b, class_label="sigma", max_iter=3,
        step_size=0.3, bond_tol=0.15, enable_symmetry=False,
    )
    assert out_smart_off == out_legacy_off, (
        "Smart env unset must be bit-exact to legacy disabled state"
    )


def test_smart_env_on_cu2_skips_projection(monkeypatch):
    """``DELFIN_B5_PHASE_A5_SMART=1`` + Cu(II) → projection skipped via
    classifier; report records the skip.

    We instrument :func:`_best_polyhedron_for_metal` to count invocations:
    in smart-mode, the Cu(II) classifier fires BEFORE the polyhedron picker,
    so _best_polyhedron_for_metal must NOT be called.
    """
    monkeypatch.setenv("DELFIN_B5_PHASE_A5_SMART", "1")
    monkeypatch.delenv("DELFIN_B5_PHASE_A5_ENABLE", raising=False)

    pick_log: List[bool] = []
    real_pick = _post_optimizer._best_polyhedron_for_metal

    def _spy(coords, m_idx, donor_idxs):
        pick_log.append(True)
        return real_pick(coords, m_idx, donor_idxs)

    monkeypatch.setattr(
        _post_optimizer, "_best_polyhedron_for_metal", _spy
    )

    atoms, bonds, charges = _cu2_hexaaqua_4plus2()
    xyz = _xyz(atoms)
    mol = _build_mol(atoms, bonds, charges)
    _, report = post_optimize_geometry(
        xyz, mol, class_label="sigma", max_iter=2,
        step_size=0.3, bond_tol=0.15, enable_symmetry=True,
    )
    # Polyhedron picker must NOT have been called for the JT-Cu metal.
    assert len(pick_log) == 0, (
        f"Cu(II) Jahn-Teller must skip _best_polyhedron_for_metal entirely; "
        f"got {len(pick_log)} invocations"
    )
    # Skip metadata recorded.
    skipped = report.get("smart_a5_skipped_metals", [])
    assert skipped, "smart_a5_skipped_metals must be populated for Cu(II) JT"
    assert any(lbl == "jahn_teller_cu2" for _, lbl in skipped), (
        f"Expected jahn_teller_cu2 label in skipped list, got {skipped!r}"
    )


def test_smart_env_on_regular_oct_runs_projection(monkeypatch):
    """``DELFIN_B5_PHASE_A5_SMART=1`` + Cr(NH3)6 (no distortion) → standard
    projection runs (no skip).
    """
    monkeypatch.setenv("DELFIN_B5_PHASE_A5_SMART", "1")
    monkeypatch.delenv("DELFIN_B5_PHASE_A5_ENABLE", raising=False)

    pick_log: List[bool] = []
    real_pick = _post_optimizer._best_polyhedron_for_metal

    def _spy(coords, m_idx, donor_idxs):
        pick_log.append(True)
        return real_pick(coords, m_idx, donor_idxs)

    monkeypatch.setattr(
        _post_optimizer, "_best_polyhedron_for_metal", _spy
    )

    atoms, bonds, charges = _cr_nh3_6()
    xyz = _xyz(atoms)
    mol = _build_mol(atoms, bonds, charges)
    _, report = post_optimize_geometry(
        xyz, mol, class_label="sigma", max_iter=2,
        step_size=0.3, bond_tol=0.15, enable_symmetry=True,
    )
    # Polyhedron picker MUST have been called for Cr.
    assert len(pick_log) >= 1, (
        f"Cr(NH3)6 should run standard projection (>=1 picker call); "
        f"got {len(pick_log)}"
    )
    # No skip metadata.
    skipped = report.get("smart_a5_skipped_metals", [])
    assert not skipped, (
        f"Cr(NH3)6 must not produce skip entries, got {skipped!r}"
    )


def test_smart_preserves_cu2_4plus2_geometry(monkeypatch):
    """With smart-mode ON, Cu(II) 4+2 axial elongation must be preserved:
    after post_optimize_geometry, the two axial Cu-O distances must remain
    distinct from the equatorial Cu-O distances (no projection to a
    symmetric polyhedron).
    """
    monkeypatch.setenv("DELFIN_B5_PHASE_A5_SMART", "1")
    monkeypatch.delenv("DELFIN_B5_PHASE_A5_ENABLE", raising=False)
    atoms, bonds, charges = _cu2_hexaaqua_4plus2()
    xyz = _xyz(atoms)
    mol = _build_mol(atoms, bonds, charges)
    out_xyz, _ = post_optimize_geometry(
        xyz, mol, class_label="sigma", max_iter=5,
        step_size=0.3, bond_tol=0.15, enable_symmetry=True,
    )
    # Parse output XYZ.
    lines = out_xyz.strip().splitlines()
    coords = np.array(
        [[float(x) for x in line.split()[1:4]] for line in lines[2:2 + 7]],
        dtype=float,
    )
    m = coords[0]
    eq_dists = [float(np.linalg.norm(coords[i] - m)) for i in (1, 2, 3, 4)]
    ax_dists = [float(np.linalg.norm(coords[i] - m)) for i in (5, 6)]
    # Axial must remain significantly longer than equatorial (Δ > 0.30 Å).
    assert min(ax_dists) - max(eq_dists) > 0.30, (
        f"Cu(II) Jahn-Teller 4+2 must be preserved by smart-mode: "
        f"axial {ax_dists}, equatorial {eq_dists}"
    )


def test_smart_env_invalid_value_safe(monkeypatch):
    """Garbage env value → treated as 0 (off, fail-safe)."""
    monkeypatch.setenv("DELFIN_B5_PHASE_A5_SMART", "not-an-int")
    monkeypatch.delenv("DELFIN_B5_PHASE_A5_ENABLE", raising=False)
    atoms, bonds, charges = _cu2_hexaaqua_4plus2()
    xyz = _xyz(atoms)
    mol_a = _build_mol(atoms, bonds, charges)
    mol_b = _build_mol(atoms, bonds, charges)
    out_invalid, _ = post_optimize_geometry(
        xyz, mol_a, class_label="sigma", max_iter=3,
        step_size=0.3, bond_tol=0.15, enable_symmetry=True,
    )
    out_explicit_off, _ = post_optimize_geometry(
        xyz, mol_b, class_label="sigma", max_iter=3,
        step_size=0.3, bond_tol=0.15, enable_symmetry=False,
    )
    assert out_invalid == out_explicit_off, (
        "Invalid smart-env value must fall back to OFF — no Phase A.5"
    )


def test_smart_a5_skipped_metals_not_set_when_no_skip(monkeypatch):
    """Standard projection path with no distortion → no skip key in report."""
    monkeypatch.setenv("DELFIN_B5_PHASE_A5_SMART", "1")
    monkeypatch.delenv("DELFIN_B5_PHASE_A5_ENABLE", raising=False)
    atoms, bonds, charges = _cr_nh3_6()
    xyz = _xyz(atoms)
    mol = _build_mol(atoms, bonds, charges)
    _, report = post_optimize_geometry(
        xyz, mol, class_label="sigma", max_iter=2,
        step_size=0.3, bond_tol=0.15, enable_symmetry=True,
    )
    assert "smart_a5_skipped_metals" not in report, (
        f"smart_a5_skipped_metals must NOT appear when no metal is skipped, "
        f"got {report.get('smart_a5_skipped_metals')!r}"
    )
