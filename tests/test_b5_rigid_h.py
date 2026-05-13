"""Tests for Patch P-H-TRACK: rigid-H tracking in B5 v2 (Baustein 5).

Validates that bonded H atoms are dragged with their heavy parent during
single-atom moves in :func:`_stage_bonds` (M, D shifts) and
:func:`_stage_clashes` Strategy 2 (clash push), preserving C-H / N-H / O-H
bond lengths that the heavy-atom-only topology gate does not check.

All tests skip cleanly if the module or RDKit is missing.
"""
from __future__ import annotations

from typing import List, Sequence, Tuple

import pytest

_post_optimizer = pytest.importorskip(
    "delfin._post_optimizer",
    reason="Baustein 5 post-optimizer module not yet implemented",
)
Chem = pytest.importorskip("rdkit.Chem", reason="RDKit required for tests")
_rdGeometry = pytest.importorskip("rdkit.Geometry", reason="RDKit required")

_compute_h_neighbors = getattr(_post_optimizer, "_compute_h_neighbors", None)
post_optimize_geometry = getattr(_post_optimizer, "post_optimize_geometry", None)
if _compute_h_neighbors is None or post_optimize_geometry is None:
    pytest.skip("rigid-H tracking not wired", allow_module_level=True)


Atom = Tuple[str, float, float, float]
Bond = Tuple[int, int]


def _xyz(atoms: Sequence[Atom]) -> str:
    lines = [str(len(atoms)), "test"]
    for sym, x, y, z in atoms:
        lines.append(f"{sym:<3s} {x:14.8f} {y:14.8f} {z:14.8f}")
    return "\n".join(lines) + "\n"


def _build_mol(atoms: Sequence[Atom], bonds: Sequence[Bond]):
    mol = Chem.RWMol()
    for sym, _, _, _ in atoms:
        a = Chem.Atom(sym)
        a.SetNoImplicit(True)
        mol.AddAtom(a)
    for i, j in bonds:
        mol.AddBond(i, j, Chem.BondType.SINGLE)
    conf = Chem.Conformer(len(atoms))
    for i, (_, x, y, z) in enumerate(atoms):
        conf.SetAtomPosition(i, _rdGeometry.Point3D(float(x), float(y), float(z)))
    mol.AddConformer(conf, assignId=True)
    return mol


def _parse(xyz: str) -> List[Atom]:
    lines = xyz.strip().splitlines()
    n = int(lines[0])
    out: List[Atom] = []
    for line in lines[2:2 + n]:
        parts = line.split()
        out.append((parts[0], float(parts[1]), float(parts[2]), float(parts[3])))
    return out


def _dist(a: Atom, b: Atom) -> float:
    return ((a[1] - b[1]) ** 2 + (a[2] - b[2]) ** 2 + (a[3] - b[3]) ** 2) ** 0.5


# ----------------------------------------------------------------------
# 1. _compute_h_neighbors basic correctness
# ----------------------------------------------------------------------
def test_h_neighbors_methane():
    """C atom in CH4 must list all 4 bonded H atoms."""
    atoms = [
        ("C", 0.0, 0.0, 0.0),
        ("H", 1.09, 0.0, 0.0),
        ("H", -0.36, 1.03, 0.0),
        ("H", -0.36, -0.51, 0.89),
        ("H", -0.36, -0.51, -0.89),
    ]
    bonds = [(0, 1), (0, 2), (0, 3), (0, 4)]
    mol = _build_mol(atoms, bonds)
    h_nbrs = _compute_h_neighbors(mol)
    assert sorted(h_nbrs[0]) == [1, 2, 3, 4]
    # H atoms themselves get empty lists.
    for h in (1, 2, 3, 4):
        assert h_nbrs[h] == []


def test_h_neighbors_no_hydrogens():
    """Molecule without H atoms gets all-empty lists (safe to use unconditionally)."""
    atoms = [("Pt", 0.0, 0.0, 0.0), ("Cl", 2.3, 0.0, 0.0), ("Cl", -2.3, 0.0, 0.0)]
    bonds = [(0, 1), (0, 2)]
    mol = _build_mol(atoms, bonds)
    h_nbrs = _compute_h_neighbors(mol)
    assert len(h_nbrs) == 3
    assert all(v == [] for v in h_nbrs)


def test_h_neighbors_mixed_donors():
    """N, O, C donors each report their H children correctly."""
    atoms = [
        ("Pt", 0.0, 0.0, 0.0),
        ("N", 2.05, 0.0, 0.0),
        ("H", 2.55, 0.85, 0.0),
        ("H", 2.55, -0.85, 0.0),
        ("O", -2.05, 0.0, 0.0),
        ("H", -2.55, 0.85, 0.0),
        ("C", 0.0, 2.05, 0.0),
        ("H", 0.85, 2.55, 0.0),
    ]
    bonds = [(0, 1), (1, 2), (1, 3), (0, 4), (4, 5), (0, 6), (6, 7)]
    mol = _build_mol(atoms, bonds)
    h_nbrs = _compute_h_neighbors(mol)
    assert sorted(h_nbrs[1]) == [2, 3]   # N has 2 H
    assert sorted(h_nbrs[4]) == [5]      # O has 1 H
    assert sorted(h_nbrs[6]) == [7]      # C has 1 H
    assert h_nbrs[0] == []               # Pt has no H


# ----------------------------------------------------------------------
# 2. Stage 1 bond shift preserves N-H when rigid_h=True (vs broken when False)
# ----------------------------------------------------------------------
def _pt_nh3_distorted() -> Tuple[List[Atom], List[Bond]]:
    """Pt-NH3 with Pt-N stretched from 2.05 to 2.15 Å (inside [0.93, 1.07] gate
    so Stage 1 actually fires). N has 3 H children at ~1.01 Å each."""
    atoms = [
        ("Pt", 0.0, 0.0, 0.0),
        ("N",  2.15, 0.0, 0.0),
        ("H",  2.45, 0.95, 0.0),
        ("H",  2.45, -0.475, 0.823),
        ("H",  2.45, -0.475, -0.823),
    ]
    bonds = [(0, 1), (1, 2), (1, 3), (1, 4)]
    return atoms, bonds


def _two_donor_clash() -> Tuple[List[Atom], List[Bond]]:
    """Bimetallic Pt2 with two CH3 donors whose C atoms clash mildly.

    Two independent Pt-CH3 fragments, C atoms 2.7 Å apart — inside the
    VdW clash threshold (~2.89 Å for C-C) but loose enough that Strategy 2
    (mass-weighted pair push) succeeds without breaking the Pt-C topology
    gate. Strategy 1 is skipped because the two atoms belong to fragments
    under different metals, and Strategy 2 fires before Strategy 3.

    Without rigid_h, the heavy C atoms move but their bonded H stay put →
    C-H drift detectable. With rigid_h=True, H follow rigidly.
    """
    atoms = [
        # Left fragment: Pt0 - C1 - {H2, H3, H4}
        ("Pt", -2.05,  0.0,    0.0),
        ("C",   0.0,   0.0,    0.0),
        ("H",   0.30,  1.03,   0.0),
        ("H",   0.30, -0.515,  0.892),
        ("H",   0.30, -0.515, -0.892),
        # Right fragment: C6 - Pt5  (Pt at +4.75 so Pt-C ideal 2.05)
        ("Pt",  4.75,  0.0,    0.0),
        ("C",   2.70,  0.0,    0.0),
        ("H",   2.40,  1.03,   0.0),
        ("H",   2.40, -0.515,  0.892),
        ("H",   2.40, -0.515, -0.892),
    ]
    bonds = [(0, 1), (1, 2), (1, 3), (1, 4),
             (5, 6), (6, 7), (6, 8), (6, 9)]
    return atoms, bonds


def _ch_dists(atoms_now, c_idx: int, h_idxs):
    return [_dist(atoms_now[c_idx], atoms_now[h]) for h in h_idxs]


def test_rigid_h_off_breaks_ch_bond():
    """Without rigid_h: Stage 3 Strategy 2 pushes the two C atoms apart but
    leaves their H children behind → C-H bonds stretch. Pins the pre-patch bug.
    """
    atoms, bonds = _two_donor_clash()
    initial_ch1 = _ch_dists(atoms, 1, (2, 3, 4))
    initial_ch2 = _ch_dists(atoms, 6, (7, 8, 9))
    xyz = _xyz(atoms)
    mol = _build_mol(atoms, bonds)
    out_xyz, report = post_optimize_geometry(
        xyz, mol, class_label="sigma", max_iter=5, step_size=0.5,
        bond_tol=0.05, md_drift_max=0.5, rigid_h=False,
    )
    out_atoms = _parse(out_xyz)
    c_moved = abs(_dist(out_atoms[1], out_atoms[6])
                  - _dist(atoms[1], atoms[6])) > 1e-3
    assert c_moved, "Test setup must trigger Stage 3 clash push"
    final_ch1 = _ch_dists(out_atoms, 1, (2, 3, 4))
    final_ch2 = _ch_dists(out_atoms, 6, (7, 8, 9))
    max_drift = max(
        max(abs(f - i) for f, i in zip(final_ch1, initial_ch1)),
        max(abs(f - i) for f, i in zip(final_ch2, initial_ch2)),
    )
    assert max_drift > 1e-3, (
        f"Without rigid_h, C moves but H stays → C-H drift expected, got {max_drift:.6f}"
    )


def test_rigid_h_on_preserves_ch_bonds():
    """With rigid_h=True: Stage 3 Strategy 2 drags bonded H of each C →
    C-H bond lengths are preserved exactly. Core acceptance test for Patch P-H-TRACK.
    """
    atoms, bonds = _two_donor_clash()
    initial_ch1 = _ch_dists(atoms, 1, (2, 3, 4))
    initial_ch2 = _ch_dists(atoms, 6, (7, 8, 9))
    xyz = _xyz(atoms)
    mol = _build_mol(atoms, bonds)
    out_xyz, report = post_optimize_geometry(
        xyz, mol, class_label="sigma", max_iter=5, step_size=0.5,
        bond_tol=0.05, md_drift_max=0.5, rigid_h=True,
    )
    out_atoms = _parse(out_xyz)
    final_ch1 = _ch_dists(out_atoms, 1, (2, 3, 4))
    final_ch2 = _ch_dists(out_atoms, 6, (7, 8, 9))
    for f, i in zip(final_ch1, initial_ch1):
        assert abs(f - i) < 1e-6, (
            f"C1-H drifted under rigid_h=True: {i:.6f} → {f:.6f}"
        )
    for f, i in zip(final_ch2, initial_ch2):
        assert abs(f - i) < 1e-6, (
            f"C2-H drifted under rigid_h=True: {i:.6f} → {f:.6f}"
        )


# ----------------------------------------------------------------------
# 3. Default-OFF guarantee
# ----------------------------------------------------------------------
def test_rigid_h_default_off_bit_exact():
    """rigid_h omitted == rigid_h=False (bit-exact pre-patch default)."""
    atoms, bonds = _pt_nh3_distorted()
    xyz = _xyz(atoms)
    mol_a = _build_mol(atoms, bonds)
    mol_b = _build_mol(atoms, bonds)
    out_a, _ = post_optimize_geometry(
        xyz, mol_a, class_label="sigma", max_iter=5,
        step_size=0.5, bond_tol=0.05,
    )
    out_b, _ = post_optimize_geometry(
        xyz, mol_b, class_label="sigma", max_iter=5,
        step_size=0.5, bond_tol=0.05, rigid_h=False,
    )
    assert out_a == out_b, "default and rigid_h=False must produce identical XYZ"


# ----------------------------------------------------------------------
# 4. No-H molecule: rigid_h flag has zero effect (no spurious behaviour change)
# ----------------------------------------------------------------------
def test_rigid_h_no_op_when_no_hydrogens():
    """If the molecule has no H atoms, rigid_h=True must equal rigid_h=False."""
    atoms = [
        ("Pt", 0.0, 0.0, 0.0),
        ("Cl", 2.6, 0.0, 0.0),
        ("Cl", -2.6, 0.0, 0.0),
        ("Cl", 0.0, 2.6, 0.0),
        ("Cl", 0.0, -2.6, 0.0),
    ]
    bonds = [(0, 1), (0, 2), (0, 3), (0, 4)]
    xyz = _xyz(atoms)
    mol_a = _build_mol(atoms, bonds)
    mol_b = _build_mol(atoms, bonds)
    out_a, _ = post_optimize_geometry(
        xyz, mol_a, class_label="sigma", max_iter=5,
        step_size=0.3, bond_tol=0.05, rigid_h=False,
    )
    out_b, _ = post_optimize_geometry(
        xyz, mol_b, class_label="sigma", max_iter=5,
        step_size=0.3, bond_tol=0.05, rigid_h=True,
    )
    assert out_a == out_b, "rigid_h must be no-op when molecule has no H"


# ----------------------------------------------------------------------
# 5. Topology gate still respected when rigid_h=True
# ----------------------------------------------------------------------
def test_rigid_h_preserves_topology_gate():
    """rigid_h=True must not break topology — output Pt-N must stay in [0.85, 1.10] × ideal
    (the heavy-atom-only topology gate enforced by _passes_topology) and the
    optimizer must report topology_preserved=True. Convergence into the tighter
    [0.93, 1.07] v2 quality band is not required from a single 0.30 Å distortion."""
    atoms, bonds = _pt_nh3_distorted()
    xyz = _xyz(atoms)
    mol = _build_mol(atoms, bonds)
    out_xyz, report = post_optimize_geometry(
        xyz, mol, class_label="sigma", max_iter=10, step_size=0.5,
        bond_tol=0.02, md_drift_max=0.5, rigid_h=True,
    )
    out_atoms = _parse(out_xyz)
    pt_n = _dist(out_atoms[0], out_atoms[1])
    # v2 topology gate is [0.93, 1.07] × ideal (Pt-N ideal ~2.0-2.10 Å fallback)
    assert pt_n <= 1.07 * 2.10, f"Pt-N exceeded topology gate: {pt_n:.3f}"
    assert pt_n >= 0.93 * 1.95, f"Pt-N below topology gate: {pt_n:.3f}"
    assert report.get("topology_preserved", False), "topology_preserved flag must be True"


# ----------------------------------------------------------------------
# 6. Multi-iter correctness: H tracks heavy parent across many sweeps
# ----------------------------------------------------------------------
def test_rigid_h_stable_over_multiple_iters():
    """N-H lengths stay constant across many Gauss-Seidel sweeps with rigid_h=True."""
    atoms, bonds = _pt_nh3_distorted()
    initial_nh = [_dist(atoms[1], atoms[i]) for i in (2, 3, 4)]
    xyz = _xyz(atoms)
    mol = _build_mol(atoms, bonds)
    out_xyz, _ = post_optimize_geometry(
        xyz, mol, class_label="sigma", max_iter=20, step_size=0.2,
        bond_tol=0.02, md_drift_max=0.5, rigid_h=True,
    )
    out_atoms = _parse(out_xyz)
    final_nh = [_dist(out_atoms[1], out_atoms[i]) for i in (2, 3, 4)]
    for f, i in zip(final_nh, initial_nh):
        assert abs(f - i) < 1e-6, (
            f"N-H drift after 20 sweeps with rigid_h=True: {i:.6f} → {f:.6f}"
        )
