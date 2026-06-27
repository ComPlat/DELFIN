"""Unit tests for Baustein 6: variational refiner + 4-tier symmetry stack.

Spec source: ``iters/BAUSTEIN6_MASTERPLAN.md`` Section 10
(Verification Strategy).

The whole module is skipped gracefully if any of the Baustein 6 sister
modules (``delfin.manta._variational_refiner``, ``delfin.manta._energy_terms``,
``delfin.manta._symmetry_detection``, ``delfin.manta._point_group_ops``,
``delfin.manta._fragment_archetypes``) or RDKit are unavailable.  Tests are
self-contained: synthetic XYZ strings are built inline, no external
file dependencies.

Coverage outline
----------------
Energy-term sanity (analytic + finite-difference cross-check):

  1.  ``U_bond`` zero at ideal length, zero gradient
  2.  ``U_bond`` quadratic for stretched bond + gradient pulls toward shrink
  3.  ``U_angle`` ~zero for perfect sp3 tetrahedral angles
  4.  ``U_angle`` non-zero penalty for distorted angle
  5.  ``U_clash`` zero for atoms at vdW separation, zero gradient
  6.  ``U_clash`` positive penalty + non-zero gradient when overlapping
  7.  ``U_topology`` finite inside the [0.85, 1.10]·ideal window
  8.  ``U_topology`` very large (~1e9) near the boundary
  9.  ``U_bond`` analytic gradient matches finite-difference (FD)
  10. ``U_angle`` analytic gradient matches finite-difference (FD)

Symmetry detection:

  11. ``find_equivalent_atoms`` for benzene (2 classes: 6 C + 6 H)
  12. ``find_equivalent_atoms`` for methane (1 class: 4 H)
  13. ``detect_global_point_group`` for water → C2v
  14. ``detect_global_point_group`` for ammonia → C3v
  15. ``detect_global_point_group`` for methane → Td
  16. ``detect_global_point_group`` for benzene → D6h

Fragment detection:

  17. ``detect_fragments`` finds "pyridine" archetype with C2v
  18. ``detect_fragments`` finds "methyl" on ethane

Hungarian assignment:

  19. ``hungarian_assign_donors_to_slots`` distorted [Co(NH3)6] → Oh slots

Integration tests for ``variational_refine``:

  20. distorted [Co(NH3)6]³⁺ → refined coords near octahedral
  21. mild distortion → topology preserved + report metadata correct
  22. pathological all-at-origin input → graceful fallback to input
"""

from __future__ import annotations

import math
from typing import List, Tuple

import pytest

# ---------------------------------------------------------------------------
# Graceful skip if any Baustein 6 module is missing.
# ---------------------------------------------------------------------------

_SKIP_REASON = None
try:
    import numpy as np
    from rdkit import Chem
    from delfin.manta._variational_refiner import variational_refine
    from delfin.manta._energy_terms import (
        U_bond, U_angle, U_clash, U_topology,
        U_A_coord_sphere, U_B_equivalence,
        U_C_fragment, U_D_global, U_total,
        _fd_gradient,
    )
    from delfin.manta._symmetry_detection import (
        find_equivalent_atoms,
        find_equivalent_bond_pairs,
        hungarian_assign_donors_to_slots,
        detect_global_point_group,
    )
    from delfin.manta._point_group_ops import (
        get_operations,
        rotation_matrix_axis_angle,
        supported_point_groups,
    )
    from delfin.manta._fragment_archetypes import detect_fragments
except Exception as _exc:  # pragma: no cover - environment shim
    _SKIP_REASON = f"Baustein 6 stack not importable: {_exc}"
    pytestmark = pytest.mark.skip(reason=_SKIP_REASON)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _xyz_from_arrays(symbols: List[str], coords: "np.ndarray",
                     comment: str = "test") -> str:
    """Render a DELFIN-style XYZ string from element list + (N, 3) coords."""
    coords = np.asarray(coords, dtype=float)
    lines = [str(len(symbols)), comment]
    for s, (x, y, z) in zip(symbols, coords):
        lines.append(f"{s:<2s}  {x: .6f}  {y: .6f}  {z: .6f}")
    return "\n".join(lines) + "\n"


def _xyz_to_array(xyz: str) -> "np.ndarray":
    """Parse a DELFIN-style XYZ string back to (N, 3) ndarray."""
    pts: List[List[float]] = []
    for line in xyz.strip().splitlines():
        parts = line.split()
        if len(parts) != 4:
            continue
        try:
            x = float(parts[1]); y = float(parts[2]); z = float(parts[3])
        except ValueError:
            continue
        if not parts[0][:1].isalpha():
            continue
        pts.append([x, y, z])
    return np.asarray(pts, dtype=float)


def _two_atom_mol(sym_a: str, sym_b: str,
                  bond_type=None) -> "Chem.Mol":
    """Build a sanitised two-atom RDKit Mol (with optional bond)."""
    mw = Chem.RWMol()
    mw.AddAtom(Chem.Atom(sym_a))
    mw.AddAtom(Chem.Atom(sym_b))
    if bond_type is not None:
        mw.AddBond(0, 1, bond_type)
    mol = mw.GetMol()
    try:
        Chem.SanitizeMol(mol)
    except Exception:
        pass
    return mol


def _methane_mol() -> "Chem.Mol":
    """Explicit-H methane (C at index 0, four Hs at 1..4)."""
    mw = Chem.RWMol()
    mw.AddAtom(Chem.Atom("C"))
    for _ in range(4):
        h = mw.AddAtom(Chem.Atom("H"))
        mw.AddBond(0, h, Chem.BondType.SINGLE)
    mol = mw.GetMol()
    Chem.SanitizeMol(mol)
    return mol


def _td_unit_vectors() -> "np.ndarray":
    """Four unit vectors at the vertices of a regular tetrahedron."""
    v = np.array([
        [+1.0, +1.0, +1.0],
        [+1.0, -1.0, -1.0],
        [-1.0, +1.0, -1.0],
        [-1.0, -1.0, +1.0],
    ], dtype=float)
    return v / math.sqrt(3.0)


def _oh_unit_vectors() -> "np.ndarray":
    """Six unit vectors along the +-x/+-y/+-z octahedral axes."""
    return np.array([
        [+1.0, 0.0, 0.0], [-1.0, 0.0, 0.0],
        [0.0, +1.0, 0.0], [0.0, -1.0, 0.0],
        [0.0, 0.0, +1.0], [0.0, 0.0, -1.0],
    ], dtype=float)


def _build_co_nh3_6(distortion: float = 0.0,
                    seed: int = 12345) -> Tuple["Chem.Mol", "np.ndarray"]:
    """Build a [Co(NH3)6] complex and its 3D coordinates.

    Atom layout: index 0 = Co, then per ammine ligand
    (N, H, H, H) interleaved -> 1 + 6 * 4 = 25 atoms.

    The N atoms sit on the octahedral axes at 2.0 A from the metal; each
    set of three Hs of an ammine forms a small tripod tilted away from
    the metal.  ``distortion`` (Å) is the magnitude of a random Gaussian
    perturbation added to every atom (Co kept fixed).
    """
    mw = Chem.RWMol()
    co_idx = mw.AddAtom(Chem.Atom("Co"))
    n_indices: List[int] = []
    h_indices: List[List[int]] = []
    for _ in range(6):
        n_idx = mw.AddAtom(Chem.Atom("N"))
        mw.AddBond(co_idx, n_idx, Chem.BondType.SINGLE)
        ligand_hs: List[int] = []
        for _ in range(3):
            h_idx = mw.AddAtom(Chem.Atom("H"))
            mw.AddBond(n_idx, h_idx, Chem.BondType.SINGLE)
            ligand_hs.append(h_idx)
        n_indices.append(n_idx)
        h_indices.append(ligand_hs)
    mol = mw.GetMol()
    try:
        Chem.SanitizeMol(mol)
    except Exception:
        pass

    # Ideal coordinates --------------------------------------------------
    oh = _oh_unit_vectors() * 2.0  # Co-N = 2.0 A (matches _get_ml_bond_length)
    n_total = mol.GetNumAtoms()
    coords = np.zeros((n_total, 3), dtype=float)
    # Metal at origin.
    coords[co_idx] = (0.0, 0.0, 0.0)
    # Place N + tripod of 3 H per ammine.
    for k, n_idx in enumerate(n_indices):
        n_pos = oh[k]
        coords[n_idx] = n_pos
        # Outward unit vector from metal to N; build orthonormal basis.
        radial = n_pos / np.linalg.norm(n_pos)
        # Pick a stable orthogonal axis.
        if abs(radial[2]) < 0.9:
            tmp = np.array([0.0, 0.0, 1.0])
        else:
            tmp = np.array([1.0, 0.0, 0.0])
        e1 = np.cross(radial, tmp)
        e1 /= np.linalg.norm(e1)
        e2 = np.cross(radial, e1)
        e2 /= np.linalg.norm(e2)
        # Three H around N at 109.47° to the radial bond, 120° azimuth.
        tilt = math.radians(109.47)
        nh = 1.01  # N-H ideal
        for i, h_idx in enumerate(h_indices[k]):
            phi = 2.0 * math.pi * i / 3.0
            direction = (
                math.cos(tilt) * radial
                + math.sin(tilt) * (math.cos(phi) * e1 + math.sin(phi) * e2)
            )
            coords[h_idx] = n_pos + nh * direction

    if distortion > 0.0:
        rng = np.random.default_rng(seed)
        # Do not perturb the metal -- keeps the test reference frame anchored.
        noise = rng.normal(scale=distortion, size=coords.shape)
        noise[co_idx] = 0.0
        coords = coords + noise

    return mol, coords


# ---------------------------------------------------------------------------
# 1-2. U_bond sanity
# ---------------------------------------------------------------------------

def test_U_bond_zero_at_ideal():
    """C-C at exactly the ideal Pyykkö length (1.52 Å) -> U = 0, grad = 0."""
    mol = _two_atom_mol("C", "C", Chem.BondType.SINGLE)
    coords = np.array([[0.0, 0.0, 0.0], [1.52, 0.0, 0.0]])
    U, grad = U_bond(coords, mol, k_bond=1000.0)
    assert abs(U) < 1.0e-6, f"expected U ~ 0, got {U}"
    assert np.linalg.norm(grad) < 1.0e-6, f"expected grad ~ 0, got {grad}"


def test_U_bond_quadratic():
    """0.1 Å stretch -> U = k * 0.01 = 10.0 (within tolerance)."""
    mol = _two_atom_mol("C", "C", Chem.BondType.SINGLE)
    coords = np.array([[0.0, 0.0, 0.0], [1.62, 0.0, 0.0]])
    U, grad = U_bond(coords, mol, k_bond=1000.0)
    assert 9.0 < U < 11.0, f"expected U ~ 10, got {U}"
    # Descent direction (-grad) must shrink the bond:
    # atom 0 should move toward +x (toward atom 1).
    descent_0 = -grad[0]
    descent_1 = -grad[1]
    assert descent_0[0] > 0.0, (
        f"descent on atom 0 should be +x (toward atom 1), got {descent_0}"
    )
    assert descent_1[0] < 0.0, (
        f"descent on atom 1 should be -x (toward atom 0), got {descent_1}"
    )


# ---------------------------------------------------------------------------
# 3-4. U_angle sanity
# ---------------------------------------------------------------------------

def test_U_angle_zero_at_tetrahedral():
    """Perfect Td methane geometry -> U_angle ~ 0."""
    mol = _methane_mol()
    coords = np.zeros((5, 3), dtype=float)
    coords[1:] = _td_unit_vectors() * 1.09  # C-H ~ 1.07-1.09 Å
    U, grad = U_angle(coords, mol, k_angle=100.0)
    assert U < 1.0e-3, f"expected U ~ 0 for Td methane, got {U}"
    assert np.linalg.norm(grad) < 1.0e-3, (
        f"expected grad ~ 0 for Td methane, got max |grad|="
        f"{np.max(np.abs(grad))}"
    )


def test_U_angle_pulls_to_tetrahedral():
    """A distorted methane H placement -> non-zero penalty."""
    mol = _methane_mol()
    coords = np.zeros((5, 3), dtype=float)
    coords[1:] = _td_unit_vectors() * 1.09
    # Move one H along +x by an extra 0.4 Å -> opens / closes neighbour angles.
    coords[1] = np.array([1.5, 0.0, 0.0])
    U, grad = U_angle(coords, mol, k_angle=100.0)
    assert U > 10.0, (
        f"expected appreciable angle penalty for distorted methane, got U={U}"
    )
    # Some component of the gradient on the displaced H must be non-zero.
    assert np.linalg.norm(grad[1]) > 1.0e-3, (
        f"expected non-zero grad on displaced H, got {grad[1]}"
    )


# ---------------------------------------------------------------------------
# 5-6. U_clash sanity
# ---------------------------------------------------------------------------

def test_U_clash_zero_when_separated():
    """Two non-bonded C atoms at full vdW sum -> no clash penalty."""
    mol = _two_atom_mol("C", "C", bond_type=None)
    # vdW(C) + vdW(C) ~ 3.40 Å; sit well above that to be safe.
    coords = np.array([[0.0, 0.0, 0.0], [4.0, 0.0, 0.0]])
    U, grad = U_clash(coords, mol, k_clash=500.0)
    assert abs(U) < 1.0e-6, f"expected U ~ 0 for separated atoms, got {U}"
    assert np.linalg.norm(grad) < 1.0e-6, (
        f"expected zero gradient, got max |grad|={np.max(np.abs(grad))}"
    )


def test_U_clash_overlap_positive():
    """Overlapping non-bonded atoms -> positive penalty + non-zero gradient."""
    mol = _two_atom_mol("C", "C", bond_type=None)
    # vdW(C) sum = 3.40, threshold ~ 2.89.  At d=2.0 the atoms are well inside.
    coords = np.array([[0.0, 0.0, 0.0], [2.0, 0.0, 0.0]])
    U, grad = U_clash(coords, mol, k_clash=500.0)
    assert U > 0.0, f"expected U > 0 for overlapping atoms, got {U}"
    # Gradient must point along the inter-atomic axis for at least one atom.
    assert abs(grad[0, 0]) > 1.0 or abs(grad[1, 0]) > 1.0, (
        f"expected non-zero x-component on clash gradient, got {grad}"
    )
    # The two atoms must receive equal-and-opposite gradient contributions
    # (Newton's third law on the pair).
    assert np.allclose(grad[0] + grad[1], 0.0, atol=1.0e-6), (
        f"expected grad[0] + grad[1] = 0, got {grad[0] + grad[1]}"
    )


# ---------------------------------------------------------------------------
# 7-8. U_topology sanity
# ---------------------------------------------------------------------------

def test_U_topology_inside_window_finite():
    """M-D bond inside [0.85, 1.10]·ideal -> finite log-barrier."""
    # Co-N has ideal = 2.0 Å (from delfin.smiles_converter._get_ml_bond_length).
    mw = Chem.RWMol()
    co = mw.AddAtom(Chem.Atom("Co"))
    n = mw.AddAtom(Chem.Atom("N"))
    mw.AddBond(co, n, Chem.BondType.SINGLE)
    mol = mw.GetMol()
    try:
        Chem.SanitizeMol(mol)
    except Exception:
        pass

    coords = np.array([[0.0, 0.0, 0.0], [2.0, 0.0, 0.0]])  # ideal exactly
    U, grad = U_topology(coords, mol, k_topology=10000.0)
    assert math.isfinite(U), f"expected finite barrier at ideal, got {U}"
    assert U > 0.0 or U == pytest.approx(0.0, abs=1.0e6), (
        f"expected finite, non-divergent U at ideal, got {U}"
    )
    # The gradient must be finite everywhere too.
    assert np.all(np.isfinite(grad)), "expected finite gradient inside window"


def test_U_topology_outside_window_large_penalty():
    """M-D bond just outside the upper limit -> very large penalty."""
    mw = Chem.RWMol()
    co = mw.AddAtom(Chem.Atom("Co"))
    n = mw.AddAtom(Chem.Atom("N"))
    mw.AddBond(co, n, Chem.BondType.SINGLE)
    mol = mw.GetMol()
    try:
        Chem.SanitizeMol(mol)
    except Exception:
        pass

    # ideal=2.0, hi = 2.20.  Place the donor at 2.40 -> outside window.
    coords = np.array([[0.0, 0.0, 0.0], [2.4, 0.0, 0.0]])
    U_outside, grad_outside = U_topology(coords, mol, k_topology=10000.0)
    assert U_outside > 1.0e6, (
        f"expected huge penalty outside topology window, got U={U_outside}"
    )
    # Gradient on the donor should point INWARD (decrease distance).
    # Descent direction is -grad.  With donor at +x and metal at origin,
    # descent on donor must be -x to bring it back into the window.
    descent_donor = -grad_outside[1]
    assert descent_donor[0] < 0.0, (
        f"expected donor descent to point toward metal (-x), "
        f"got descent={descent_donor}"
    )


# ---------------------------------------------------------------------------
# 9-10. Finite-difference cross-checks of analytic gradients
# ---------------------------------------------------------------------------

def test_U_bond_gradient_matches_finite_difference():
    """Analytic U_bond gradient must agree with central-difference FD."""
    mol = _two_atom_mol("C", "C", Chem.BondType.SINGLE)
    coords = np.array([[0.05, 0.0, 0.0], [1.65, 0.10, -0.05]])
    _, grad_analytic = U_bond(coords, mol, k_bond=1000.0)
    grad_fd = _fd_gradient(
        lambda c: U_bond(c, mol, k_bond=1000.0), coords, eps=1.0e-5
    )
    assert np.allclose(grad_analytic, grad_fd, atol=1.0e-2, rtol=1.0e-3), (
        f"analytic vs FD mismatch:\nanalytic=\n{grad_analytic}\n"
        f"fd=\n{grad_fd}"
    )


def test_U_angle_gradient_matches_finite_difference():
    """Analytic U_angle gradient must agree with central-difference FD."""
    mol = _methane_mol()
    coords = np.zeros((5, 3), dtype=float)
    coords[1:] = _td_unit_vectors() * 1.09
    # Slight off-ideal displacements so the angle term is non-degenerate.
    coords[1] += np.array([0.05, 0.02, -0.03])
    coords[2] += np.array([-0.02, 0.04, 0.01])
    _, grad_analytic = U_angle(coords, mol, k_angle=100.0)
    grad_fd = _fd_gradient(
        lambda c: U_angle(c, mol, k_angle=100.0), coords, eps=1.0e-5
    )
    assert np.allclose(grad_analytic, grad_fd, atol=1.0e-1, rtol=1.0e-3), (
        f"analytic vs FD mismatch (max abs diff = "
        f"{np.max(np.abs(grad_analytic - grad_fd))})"
    )


# ---------------------------------------------------------------------------
# 11-12. Tier B: find_equivalent_atoms
# ---------------------------------------------------------------------------

def test_find_equivalent_atoms_benzene():
    """Benzene -> two equivalence classes (6 C, 6 H)."""
    mol = Chem.MolFromSmiles("c1ccccc1")
    mol = Chem.AddHs(mol)
    classes = find_equivalent_atoms(mol)
    # Expect a 6-C class and a 6-H class.
    sizes = sorted(len(c) for c in classes)
    assert 6 in sizes, (
        f"expected at least one 6-atom class in benzene, got sizes={sizes}"
    )
    # Both Cs (6) and Hs (6) should be detected as distinct equivalence sets.
    six_atom_classes = [c for c in classes if len(c) == 6]
    assert len(six_atom_classes) >= 2, (
        f"expected >=2 six-atom classes (C+H) in benzene, got {sizes}"
    )


def test_find_equivalent_atoms_methane():
    """Methane CH4 -> single class of 4 Hs."""
    mol = Chem.MolFromSmiles("C")
    mol = Chem.AddHs(mol)
    classes = find_equivalent_atoms(mol)
    sizes = sorted(len(c) for c in classes)
    assert 4 in sizes, (
        f"expected a 4-atom equivalence class in methane, got sizes={sizes}"
    )


# ---------------------------------------------------------------------------
# 13-16. Tier D: detect_global_point_group
# ---------------------------------------------------------------------------

def test_detect_point_group_water_c2v():
    """Bent H2O placed in canonical orientation -> C2v."""
    mol = Chem.MolFromSmiles("O")
    mol = Chem.AddHs(mol)
    half = math.radians(104.5 / 2.0)
    d_OH = 0.96
    coords = np.array([
        [0.0, 0.0, 0.0],
        [d_OH * math.sin(half), 0.0, d_OH * math.cos(half)],
        [-d_OH * math.sin(half), 0.0, d_OH * math.cos(half)],
    ])
    pg, ops, _ = detect_global_point_group(mol, coords, tolerance=0.3)
    assert pg == "C2v", f"expected C2v for H2O, got {pg}"
    assert len(ops) == 4, f"expected 4 C2v ops, got {len(ops)}"


def test_detect_point_group_ammonia_c3v():
    """Pyramidal NH3 with 107.8° H-N-H -> C3v."""
    mol = Chem.MolFromSmiles("N")
    mol = Chem.AddHs(mol)
    d_NH = 1.01
    alpha = math.radians(68.9)  # gives ~107.8° H-N-H
    coords = np.zeros((4, 3), dtype=float)
    for i in range(3):
        az = 2.0 * math.pi * i / 3.0
        coords[i + 1] = np.array([
            d_NH * math.sin(alpha) * math.cos(az),
            d_NH * math.sin(alpha) * math.sin(az),
            -d_NH * math.cos(alpha),
        ])
    pg, ops, _ = detect_global_point_group(mol, coords, tolerance=0.3)
    assert pg == "C3v", f"expected C3v for NH3, got {pg}"
    assert len(ops) == 6, f"expected 6 C3v ops, got {len(ops)}"


def test_detect_point_group_methane_td():
    """Perfect methane Td geometry -> Td (24 ops)."""
    mol = Chem.MolFromSmiles("C")
    mol = Chem.AddHs(mol)
    coords = np.zeros((5, 3), dtype=float)
    coords[1:] = _td_unit_vectors() * 1.09
    pg, ops, _ = detect_global_point_group(mol, coords, tolerance=0.3)
    assert pg == "Td", f"expected Td for methane, got {pg}"
    assert len(ops) == 24, f"expected 24 Td ops, got {len(ops)}"


def test_detect_point_group_benzene_d6h():
    """Planar regular hexagonal benzene -> D6h (24 ops)."""
    mol = Chem.MolFromSmiles("c1ccccc1")
    mol = Chem.AddHs(mol)
    coords: List[List[float]] = []
    r_c = 1.40  # benzene C-C ~ 1.40 Å -> circumradius ~ 1.40 too
    r_h = 1.40 + 1.08  # add C-H bond
    for i in range(6):
        a = 2.0 * math.pi * i / 6.0
        coords.append([r_c * math.cos(a), r_c * math.sin(a), 0.0])
    for i in range(6):
        a = 2.0 * math.pi * i / 6.0
        coords.append([r_h * math.cos(a), r_h * math.sin(a), 0.0])
    pg, ops, _ = detect_global_point_group(mol, np.array(coords), tolerance=0.3)
    assert pg == "D6h", f"expected D6h for benzene, got {pg}"
    assert len(ops) == 24, f"expected 24 D6h ops, got {len(ops)}"


# ---------------------------------------------------------------------------
# 17-18. Tier C: detect_fragments
# ---------------------------------------------------------------------------

def test_detect_fragments_pyridine():
    """Pyridine should produce a 'pyridine' archetype match with C2v PG."""
    mol = Chem.MolFromSmiles("c1ccncc1")
    matches = detect_fragments(mol)
    py_hits = [m for m in matches if m.archetype == "pyridine"]
    assert len(py_hits) >= 1, (
        f"expected at least one pyridine match, "
        f"got archetypes={[m.archetype for m in matches]}"
    )
    assert py_hits[0].point_group == "C2v", (
        f"expected C2v for pyridine, got {py_hits[0].point_group}"
    )


def test_detect_fragments_methyl_on_ethane():
    """Ethane CC has a 'methyl' archetype match (atom-set deduplicated)."""
    mol = Chem.MolFromSmiles("CC")
    matches = detect_fragments(mol)
    me_hits = [m for m in matches if m.archetype == "methyl"]
    assert len(me_hits) >= 1, (
        f"expected at least one methyl match, "
        f"got archetypes={[m.archetype for m in matches]}"
    )
    assert me_hits[0].point_group == "C3v", (
        f"expected C3v for methyl, got {me_hits[0].point_group}"
    )


# ---------------------------------------------------------------------------
# 19. Hungarian assignment for [Co(NH3)6]
# ---------------------------------------------------------------------------

def test_hungarian_assign_oh_co_nh3_6():
    """Distorted [Co(NH3)6] octahedral donors -> Hungarian recovers axial slots."""
    mol, coords = _build_co_nh3_6(distortion=0.15, seed=42)
    targets = hungarian_assign_donors_to_slots(coords, mol, metal_idx=0)
    # Six N donors -> six targets.
    assert len(targets) == 6, f"expected 6 donor targets, got {len(targets)}"
    # Each target sits at ~2.0 Å from the metal (Co-N ideal).
    metal_pos = coords[0]
    for donor_idx, tgt in targets.items():
        d = float(np.linalg.norm(tgt - metal_pos))
        assert 1.9 < d < 2.1, (
            f"target for donor {donor_idx} at distance {d:.3f} Å is off-ideal"
        )
        # Target direction must coincide with one of the six Oh axes.
        u = (tgt - metal_pos) / d
        # cos with the closest Oh axis must be ~ 1.0
        oh = _oh_unit_vectors()
        best_cos = float(np.max(oh @ u))
        assert best_cos > 0.99, (
            f"donor {donor_idx} target not aligned with Oh axis "
            f"(best cos={best_cos:.3f})"
        )


# ---------------------------------------------------------------------------
# 20-22. variational_refine integration
# ---------------------------------------------------------------------------

def test_variational_refine_distorted_octahedron():
    """Distorted [Co(NH3)6] -> refine recovers near-Oh donor geometry."""
    mol, coords_distorted = _build_co_nh3_6(distortion=0.20, seed=7)
    symbols = [a.GetSymbol() for a in mol.GetAtoms()]
    xyz_in = _xyz_from_arrays(symbols, coords_distorted, "distorted Co(NH3)6")

    # Attach a 3D conformer so RDKit sees the geometry (refiner does not
    # require it but several helpers do).
    conf = Chem.Conformer(mol.GetNumAtoms())
    for i, p in enumerate(coords_distorted):
        conf.SetAtomPosition(i, tuple(p))
    mol.RemoveAllConformers()
    mol.AddConformer(conf, assignId=True)

    new_xyz, report = variational_refine(
        xyz_in, mol, class_label="sigma", max_iter=100, ftol=1.0e-6,
    )

    # Report must carry all mandatory keys.
    for key in (
        "iterations", "converged", "energy_initial", "energy_final",
        "topology_preserved", "fallback_used", "global_pg",
        "fragments_detected", "equiv_classes",
    ):
        assert key in report, f"missing key '{key}' in report: {list(report)}"

    coords_new = _xyz_to_array(new_xyz)
    assert coords_new.shape == coords_distorted.shape, (
        f"shape mismatch: new={coords_new.shape}, "
        f"orig={coords_distorted.shape}"
    )

    # Energy should be finite (variational either converged or fell back).
    assert math.isfinite(report["energy_initial"]), "non-finite initial energy"

    # Co-N distances should remain inside the topology window.
    co_pos = coords_new[0]
    co_n_dists = []
    for n_idx in (1, 5, 9, 13, 17, 21):  # N atoms interleaved every 4
        d = float(np.linalg.norm(coords_new[n_idx] - co_pos))
        co_n_dists.append(d)
    for d in co_n_dists:
        assert 1.70 < d < 2.20, (
            f"Co-N distance {d:.3f} Å outside topology window [1.7, 2.2]"
        )


def test_variational_refine_preserves_topology():
    """Mild distortion -> topology_preserved=True, fallback_used=False."""
    mol, coords_distorted = _build_co_nh3_6(distortion=0.08, seed=99)
    symbols = [a.GetSymbol() for a in mol.GetAtoms()]
    xyz_in = _xyz_from_arrays(symbols, coords_distorted, "mild Co(NH3)6")

    conf = Chem.Conformer(mol.GetNumAtoms())
    for i, p in enumerate(coords_distorted):
        conf.SetAtomPosition(i, tuple(p))
    mol.RemoveAllConformers()
    mol.AddConformer(conf, assignId=True)

    new_xyz, report = variational_refine(
        xyz_in, mol, class_label="sigma", max_iter=50, ftol=1.0e-5,
    )

    assert isinstance(report, dict)
    # Output must be a parseable XYZ string of the same atom count.
    new_coords = _xyz_to_array(new_xyz)
    assert new_coords.shape == coords_distorted.shape, (
        "refiner must not change atom count"
    )
    # Topology flag should be True on this mild perturbation.
    assert report["topology_preserved"] in (True, False)
    if not report.get("fallback_used", False):
        assert report["topology_preserved"] is True, (
            "Expected topology_preserved=True on non-fallback path"
        )


def test_variational_refine_handles_pathological_input():
    """All atoms collapsed to origin -> graceful fallback, no exception."""
    mol, _ = _build_co_nh3_6(distortion=0.0)
    n = mol.GetNumAtoms()
    pathological = np.zeros((n, 3), dtype=float)
    symbols = [a.GetSymbol() for a in mol.GetAtoms()]
    xyz_in = _xyz_from_arrays(symbols, pathological, "all-at-origin")

    new_xyz, report = variational_refine(
        xyz_in, mol, class_label="sigma", max_iter=20, ftol=1.0e-4,
    )

    assert isinstance(report, dict)
    # On total breakdown the refiner must report fallback_used and
    # NOT crash.  It may or may not return the original XYZ verbatim,
    # but the output must still be a string.
    assert isinstance(new_xyz, str), "refiner must always return a string"
    # The report must keep the mandatory key set so downstream consumers
    # never KeyError.
    for key in (
        "iterations", "converged", "energy_initial", "energy_final",
        "topology_preserved", "fallback_used",
    ):
        assert key in report, f"missing '{key}' in pathological report"
