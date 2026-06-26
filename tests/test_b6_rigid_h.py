"""Tests for Patch P-H-TRACK (B6 port): rigid-H tracking in the B6
variational refiner (``delfin._variational_refiner``).

Mirrors :mod:`tests.test_b5_rigid_h`.  L-BFGS-B minimises an 8-term
energy functional where each H atom carries its own 3 DoFs; under
angle / clash / topology pressure the heavy parent migrates and the H
child lags behind, stretching the X-H bond beyond the heavy-atom-only
topology gate's notice.  Patch P-H-TRACK ties terminal H atoms rigidly
to their unique heavy parent via coordinate substitution, eliminating
the lag.

All tests skip cleanly if the module, RDKit, or SciPy is missing.
"""
from __future__ import annotations

import os
from typing import List, Sequence, Tuple

import pytest

_variational_refiner = pytest.importorskip(
    "delfin._variational_refiner",
    reason="Baustein 6 variational refiner module not yet available",
)
Chem = pytest.importorskip("rdkit.Chem", reason="RDKit required for tests")
_rdGeometry = pytest.importorskip("rdkit.Geometry", reason="RDKit required")
pytest.importorskip("scipy", reason="SciPy required for L-BFGS-B")

variational_refine = getattr(_variational_refiner, "variational_refine", None)
_compute_h_neighbors = getattr(_variational_refiner, "_compute_h_neighbors", None)
_compute_rigid_h_map = getattr(_variational_refiner, "_compute_rigid_h_map", None)
if (variational_refine is None
        or _compute_h_neighbors is None
        or _compute_rigid_h_map is None):
    pytest.skip("rigid-H tracking not wired in B6", allow_module_level=True)


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
# 1. _compute_h_neighbors basic correctness (mirrors B5 contract)
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
    for h in (1, 2, 3, 4):
        assert h_nbrs[h] == []


def test_h_neighbors_no_hydrogens():
    """Molecule without H atoms gets all-empty lists."""
    atoms = [("Pt", 0.0, 0.0, 0.0), ("Cl", 2.3, 0.0, 0.0), ("Cl", -2.3, 0.0, 0.0)]
    bonds = [(0, 1), (0, 2)]
    mol = _build_mol(atoms, bonds)
    h_nbrs = _compute_h_neighbors(mol)
    assert len(h_nbrs) == 3
    assert all(v == [] for v in h_nbrs)


# ----------------------------------------------------------------------
# 2. _compute_rigid_h_map identifies eligible terminal H
# ----------------------------------------------------------------------
def test_rigid_h_map_terminal_only():
    """Terminal C-H, N-H, O-H are rigid; the heavy parent index is recorded."""
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
    rigid, parent = _compute_rigid_h_map(mol)
    assert rigid == [2, 3, 5, 7]
    assert parent[2] == 1
    assert parent[3] == 1
    assert parent[5] == 4
    assert parent[7] == 6
    # heavy atoms keep -1
    assert parent[0] == -1
    assert parent[1] == -1
    assert parent[4] == -1
    assert parent[6] == -1


# ----------------------------------------------------------------------
# 3. rigid_h=True preserves X-H bond lengths bit-exactly
# ----------------------------------------------------------------------
def _methylamine_distorted() -> Tuple[List[Atom], List[Bond]]:
    """Methylamine (CH3-NH2) with a slightly squashed N-H and a bent C-N-H.

    L-BFGS-B will pull the heavy C and N back to ideal geometry; bonded H
    must be dragged rigidly with their parents when rigid_h=True.
    """
    atoms = [
        ("C",  0.000, 0.000, 0.000),
        ("N",  1.470, 0.000, 0.000),
        ("H", -0.515, 0.892, 0.355),
        ("H", -0.515, -0.892, 0.355),
        ("H", -0.515, 0.000, -0.890),
        ("H",  1.870, 0.820, 0.300),
        ("H",  1.870, -0.820, 0.300),
    ]
    bonds = [(0, 1), (0, 2), (0, 3), (0, 4), (1, 5), (1, 6)]
    return atoms, bonds


def test_rigid_h_on_preserves_ch_nh_exact():
    """With rigid_h=True every C-H and N-H length is preserved bit-exactly
    regardless of how far heavy atoms drift during L-BFGS-B.
    """
    atoms, bonds = _methylamine_distorted()
    initial_ch = [_dist(atoms[0], atoms[i]) for i in (2, 3, 4)]
    initial_nh = [_dist(atoms[1], atoms[i]) for i in (5, 6)]
    xyz = _xyz(atoms)
    mol = _build_mol(atoms, bonds)
    out_xyz, report = variational_refine(
        xyz, mol, class_label="no_metal",
        max_iter=80, ftol=1e-7, enable_global_pg=False,
        rigid_h=True,
    )
    assert not report.get("fallback_used", True), report
    assert report.get("rigid_h", False) is True
    assert int(report.get("rigid_h_count", 0)) == 5
    out_atoms = _parse(out_xyz)
    final_ch = [_dist(out_atoms[0], out_atoms[i]) for i in (2, 3, 4)]
    final_nh = [_dist(out_atoms[1], out_atoms[i]) for i in (5, 6)]
    for f, i in zip(final_ch, initial_ch):
        assert abs(f - i) < 1.0e-6, (
            f"C-H drifted under rigid_h=True: {i:.6f} -> {f:.6f}"
        )
    for f, i in zip(final_nh, initial_nh):
        assert abs(f - i) < 1.0e-6, (
            f"N-H drifted under rigid_h=True: {i:.6f} -> {f:.6f}"
        )


# ----------------------------------------------------------------------
# 4. rigid_h=False == default OFF == pre-patch bit-exact
# ----------------------------------------------------------------------
def test_rigid_h_default_off_bit_exact_when_env_unset():
    """rigid_h omitted equals rigid_h=False when DELFIN_B6_RIGID_H is unset."""
    atoms, bonds = _methylamine_distorted()
    xyz = _xyz(atoms)
    mol_a = _build_mol(atoms, bonds)
    mol_b = _build_mol(atoms, bonds)
    env_backup = os.environ.pop("DELFIN_B6_RIGID_H", None)
    try:
        out_a, rep_a = variational_refine(
            xyz, mol_a, class_label="no_metal",
            max_iter=50, ftol=1e-6, enable_global_pg=False,
        )
        out_b, rep_b = variational_refine(
            xyz, mol_b, class_label="no_metal",
            max_iter=50, ftol=1e-6, enable_global_pg=False,
            rigid_h=False,
        )
    finally:
        if env_backup is not None:
            os.environ["DELFIN_B6_RIGID_H"] = env_backup
    assert rep_a.get("rigid_h", None) is False
    assert rep_b.get("rigid_h", None) is False
    assert out_a == out_b, "default (env unset) must equal explicit rigid_h=False"


# ----------------------------------------------------------------------
# 5. DELFIN_B6_RIGID_H environment flag picks up when arg is None
# ----------------------------------------------------------------------
def test_rigid_h_env_flag_enables():
    """Setting DELFIN_B6_RIGID_H=1 with rigid_h=None turns the patch on."""
    atoms, bonds = _methylamine_distorted()
    xyz = _xyz(atoms)
    mol = _build_mol(atoms, bonds)
    env_backup = os.environ.get("DELFIN_B6_RIGID_H")
    os.environ["DELFIN_B6_RIGID_H"] = "1"
    try:
        out_xyz, report = variational_refine(
            xyz, mol, class_label="no_metal",
            max_iter=40, ftol=1e-6, enable_global_pg=False,
            # rigid_h argument omitted -> env wins
        )
    finally:
        if env_backup is None:
            os.environ.pop("DELFIN_B6_RIGID_H", None)
        else:
            os.environ["DELFIN_B6_RIGID_H"] = env_backup
    assert report.get("rigid_h", False) is True
    assert int(report.get("rigid_h_count", 0)) == 5


# ----------------------------------------------------------------------
# 6. No-H molecule: rigid_h flag is a no-op
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
    out_a, rep_a = variational_refine(
        xyz, mol_a, class_label="sigma",
        max_iter=30, ftol=1e-6, enable_global_pg=False,
        rigid_h=False,
    )
    out_b, rep_b = variational_refine(
        xyz, mol_b, class_label="sigma",
        max_iter=30, ftol=1e-6, enable_global_pg=False,
        rigid_h=True,
    )
    assert int(rep_b.get("rigid_h_count", 0)) == 0
    assert out_a == out_b, "rigid_h must be no-op when molecule has no H"


# ----------------------------------------------------------------------
# 7. Topology hard-gate respected with rigid_h=True
# ----------------------------------------------------------------------
def test_rigid_h_preserves_topology_gate():
    """With rigid_h=True, the refiner must still report
    topology_preserved=True and the X-H bonds must not have drifted.
    """
    atoms, bonds = _methylamine_distorted()
    xyz = _xyz(atoms)
    mol = _build_mol(atoms, bonds)
    out_xyz, report = variational_refine(
        xyz, mol, class_label="no_metal",
        max_iter=60, ftol=1e-7, enable_global_pg=False,
        rigid_h=True,
    )
    assert report.get("topology_preserved", False)
    assert not report.get("fallback_used", True)
    # X-H bonds preserved exactly (already tested in #3, repeated here as
    # part of the topology-gate acceptance bundle).
    out_atoms = _parse(out_xyz)
    for h_idx, parent in [(2, 0), (3, 0), (4, 0), (5, 1), (6, 1)]:
        d_old = _dist(atoms[parent], atoms[h_idx])
        d_new = _dist(out_atoms[parent], out_atoms[h_idx])
        assert abs(d_old - d_new) < 1.0e-6, (
            f"X-H drift on atom {h_idx} despite rigid_h=True: "
            f"{d_old:.6f} -> {d_new:.6f}"
        )


# ----------------------------------------------------------------------
# 8. Energy non-increase: rigid_h must still lower (or keep) U
# ----------------------------------------------------------------------
def test_rigid_h_lowers_energy_like_full_run():
    """With rigid_h=True the optimizer still drives U downward — the patch
    does not strand the refiner at a higher-energy minimum.  We require
    only U_final <= U_initial + 1e-3 (tolerating round-off near zero).
    """
    atoms, bonds = _methylamine_distorted()
    xyz = _xyz(atoms)
    mol = _build_mol(atoms, bonds)
    _, report = variational_refine(
        xyz, mol, class_label="no_metal",
        max_iter=60, ftol=1e-7, enable_global_pg=False,
        rigid_h=True,
    )
    assert report["energy_final"] <= report["energy_initial"] + 1.0e-3
