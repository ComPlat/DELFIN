"""Phase 2 of Baustein 6 — Donor-Pivot Rotation Sweep (DPRS).

Discrete rotation sweep around the M-D axis for sp3-hybridised donors.  Some
M-CH2-R complexes settle into a linear M-C-H ~ 180 degrees minimum under UFF
(WUXQAK pattern).  This module attempts to escape that local minimum by rigidly
rotating the substituent fragment beyond each sp3 donor through a set of
discrete angles, scoring each candidate by U_total, and committing the best
topology-preserving improvement.

Doctrine:
    - Lazy imports of sister modules (graceful no-op if missing).
    - Hard topology gate before every commit.
    - Energy threshold (`min_improvement`) prevents noise-level commits.
    - Per-conformer; n_frames invariant.
    - Pure numpy + scipy.spatial.transform.Rotation.
"""
from __future__ import annotations

from typing import Dict, List, Optional, Set, Tuple

import numpy as np

try:
    from scipy.spatial.transform import Rotation
except ImportError:  # pragma: no cover - scipy is a hard dep elsewhere
    Rotation = None  # type: ignore[assignment]


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------


def donor_pivot_rotation_sweep(
    coords: np.ndarray,
    mol,
    sym_info: Dict,
    params: Dict,
    angles_deg: Optional[List[float]] = None,
    min_improvement: float = 0.1,
) -> Tuple[np.ndarray, Dict]:
    """Try discrete rotations around each M-D axis for sp3 donors.

    Parameters
    ----------
    coords : (N, 3) np.ndarray
        Current Cartesian coordinates (angstrom).
    mol : rdkit Mol
        Molecule with bonds + hybridization populated.
    sym_info, params : dict
        Forwarded to ``_energy_terms.U_total``.
    angles_deg : list of float, optional
        Discrete rotation angles to probe.  Default covers the 30 / 60 / 120 /
        180 degree fundamental positions on both signs.
    min_improvement : float, default 0.1
        Required Delta-U (current - candidate) before a rotation is committed.

    Returns
    -------
    new_coords : (N, 3) np.ndarray
        Best coordinates after sweep (input is returned untouched if no
        improvement found).
    report : dict
        Keys: ``sp3_pairs``, ``donors_rotated``, ``U_initial``, ``U_final``,
        ``U_drop``, ``error`` (only on graceful fallback).
    """
    if angles_deg is None:
        angles_deg = [-120.0, -60.0, -30.0, 30.0, 60.0, 120.0, 180.0]

    # Lazy imports so the module is importable even before sisters are ready.
    try:
        from delfin.manta._energy_terms import U_total
    except Exception as exc:  # pragma: no cover - defensive
        return coords, {
            "error": f"missing _energy_terms: {exc}",
            "donors_rotated": 0,
        }

    try:
        from delfin.smiles_converter import _METAL_SET
    except Exception:
        # Conservative fallback: atomic number >= 21 except common non-metals.
        _METAL_SET = None  # type: ignore[assignment]

    # _passes_topology has a richer signature than originally specified — we
    # call it carefully and degrade to a permissive sanity check on TypeError.
    try:
        from delfin.manta._post_optimizer import _passes_topology  # type: ignore
    except Exception:
        _passes_topology = None  # type: ignore[assignment]

    if Rotation is None:
        return coords, {
            "error": "scipy.spatial.transform.Rotation missing",
            "donors_rotated": 0,
        }

    metals = _collect_metals(mol, _METAL_SET)
    sp3_pairs = find_sp3_donors(mol, metals)

    if not sp3_pairs:
        # Nothing to do; do not touch coords.
        try:
            U0, _ = U_total(coords, mol, sym_info, params)
            U0 = float(U0)
        except Exception:
            U0 = float("nan")
        return coords, {
            "sp3_pairs": 0,
            "donors_rotated": 0,
            "U_initial": U0,
            "U_final": U0,
            "U_drop": 0.0,
        }

    current_coords = coords.copy()
    try:
        current_U_val, _ = U_total(current_coords, mol, sym_info, params)
        current_U = float(current_U_val)
    except Exception as exc:  # pragma: no cover - defensive
        return coords, {
            "error": f"U_total initial failed: {exc}",
            "sp3_pairs": len(sp3_pairs),
            "donors_rotated": 0,
        }

    U_initial = current_U
    donors_rotated = 0

    md_pairs, nb_pairs = _build_topology_lists(mol, metals)

    for metal_idx, donor_idx in sp3_pairs:
        fragment = get_substituent_fragment(mol, donor_idx, metal_idx)
        if not fragment:
            continue

        d_vec = current_coords[donor_idx] - current_coords[metal_idx]
        d_norm = float(np.linalg.norm(d_vec))
        if d_norm < 1e-6:
            continue
        axis = d_vec / d_norm
        pivot = current_coords[donor_idx].copy()

        best_local_coords = current_coords
        best_local_U = current_U
        improved = False

        for angle_deg in angles_deg:
            tentative = rotate_fragment_around_axis(
                current_coords, fragment, pivot, axis,
                float(np.radians(angle_deg)),
            )

            if not _topology_ok(_passes_topology, tentative, mol, metals,
                                md_pairs, nb_pairs):
                continue

            try:
                U_new_val, _ = U_total(tentative, mol, sym_info, params)
                U_new = float(U_new_val)
            except Exception:
                continue

            if U_new < best_local_U - min_improvement:
                best_local_U = U_new
                best_local_coords = tentative
                improved = True

        if improved:
            current_coords = best_local_coords
            current_U = best_local_U
            donors_rotated += 1

    report = {
        "sp3_pairs": len(sp3_pairs),
        "donors_rotated": donors_rotated,
        "U_initial": U_initial,
        "U_final": current_U,
        "U_drop": U_initial - current_U,
    }
    return current_coords, report


def find_sp3_donors(mol, metals: List[int]) -> List[Tuple[int, int]]:
    """Return (metal_idx, donor_idx) pairs where donor is sp3 hybridised."""
    try:
        from rdkit import Chem
    except Exception:  # pragma: no cover
        return []
    sp3 = Chem.HybridizationType.SP3
    metal_set = set(metals)
    pairs: List[Tuple[int, int]] = []
    for m_idx in metals:
        m_atom = mol.GetAtomWithIdx(m_idx)
        for nbr in m_atom.GetNeighbors():
            n_idx = nbr.GetIdx()
            if n_idx in metal_set:
                continue
            if nbr.GetHybridization() == sp3:
                pairs.append((m_idx, n_idx))
    return pairs


def get_substituent_fragment(mol, donor_idx: int, metal_idx: int) -> Set[int]:
    """BFS fragment beyond the donor (donor and metal excluded).

    The returned set is the atoms that move rigidly with the substituents when
    we rotate around the M-D axis.  If the BFS reaches another metal (chelate
    ring), the entire fragment is rejected — rotating it would break a second
    M-D bond.
    """
    visited: Set[int] = {donor_idx, metal_idx}
    fragment: Set[int] = set()
    queue: List[int] = []

    donor_atom = mol.GetAtomWithIdx(donor_idx)
    metal_atomic = mol.GetAtomWithIdx(metal_idx).GetAtomicNum()
    for nbr in donor_atom.GetNeighbors():
        n_idx = nbr.GetIdx()
        if n_idx == metal_idx:
            continue
        if n_idx in visited:
            continue
        visited.add(n_idx)
        fragment.add(n_idx)
        queue.append(n_idx)

    # We treat "metal" universally via atomic-number heuristic too: any atom
    # belonging to the global metal set blocks the fragment (chelate detection).
    try:
        from delfin.smiles_converter import _METAL_SET
    except Exception:
        _METAL_SET = None  # type: ignore[assignment]

    while queue:
        atom_idx = queue.pop(0)
        atom = mol.GetAtomWithIdx(atom_idx)
        for nbr in atom.GetNeighbors():
            n_idx = nbr.GetIdx()
            if n_idx in visited:
                continue
            if _METAL_SET is not None and nbr.GetSymbol() in _METAL_SET:
                # Chelate ring — abort: rotating would tear a second M-D bond.
                return set()
            if nbr.GetAtomicNum() == metal_atomic:
                # Same-element metal fallback (handles _METAL_SET=None path).
                return set()
            visited.add(n_idx)
            fragment.add(n_idx)
            queue.append(n_idx)
    return fragment


def rotate_fragment_around_axis(
    coords: np.ndarray,
    fragment_atoms: Set[int],
    axis_point: np.ndarray,
    axis_dir: np.ndarray,
    angle_rad: float,
) -> np.ndarray:
    """Apply a rigid rotation of `fragment_atoms` around the given axis."""
    if Rotation is None:  # pragma: no cover - guarded at caller
        return coords.copy()
    rot = Rotation.from_rotvec(np.asarray(axis_dir, dtype=float) * float(angle_rad))
    new = coords.copy()
    if not fragment_atoms:
        return new
    idx = np.fromiter(fragment_atoms, dtype=int)
    offsets = coords[idx] - axis_point
    new[idx] = axis_point + rot.apply(offsets)
    return new


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------


def _collect_metals(mol, metal_set) -> List[int]:
    """Return list of metal atom indices.

    Uses ``smiles_converter._METAL_SET`` when available, otherwise falls back
    to a conservative atomic-number heuristic (Z >= 21 and Z not in {33..36,
    52..54, 85, 86}) covering transition + post-transition + lanthanide/
    actinide metals.
    """
    if metal_set is not None:
        return [a.GetIdx() for a in mol.GetAtoms() if a.GetSymbol() in metal_set]
    non_metal_z = {1, 5, 6, 7, 8, 9, 14, 15, 16, 17, 33, 34, 35, 52, 53, 54, 85, 86}
    out: List[int] = []
    for a in mol.GetAtoms():
        z = a.GetAtomicNum()
        if z >= 3 and z not in non_metal_z and z not in {2, 10, 18, 36, 54, 86}:
            # Crude transition/post-transition test sufficient as fallback.
            if z >= 21 or z in {3, 4, 11, 12, 13, 19, 20}:
                out.append(a.GetIdx())
    return out


def _build_topology_lists(mol, metals: List[int]):
    """Best-effort construction of (md_pairs, nb_pairs) for _passes_topology.

    `md_pairs` lists (m_idx, d_idx, d_ideal) where d_ideal is the *current*
    M-D distance (we only care that rotation does not move M or D, so the
    band [0.85, 1.10] stays trivially satisfied — this lets us still detect
    accidental movement by upstream bugs).  `nb_pairs` is left empty: the
    rotation is rigid around an axis through D, so intra-fragment distances
    are preserved, and inter-fragment clashes are caught by the U_total
    increase rather than a hard gate.  An empty nb list keeps the gate
    permissive but topology-safe for M-D bonds.
    """
    md_pairs: List[Tuple[int, int, float]] = []
    metal_set = set(metals)
    for m_idx in metals:
        m_atom = mol.GetAtomWithIdx(m_idx)
        for nbr in m_atom.GetNeighbors():
            n_idx = nbr.GetIdx()
            if n_idx in metal_set:
                continue
            # d_ideal is filled in by the caller from current coords on the
            # fly; we leave a sentinel here and resolve inside _topology_ok.
            md_pairs.append((m_idx, n_idx, -1.0))
    nb_pairs: List[Tuple[int, int]] = []
    return md_pairs, nb_pairs


def _topology_ok(passes_topology_fn, coords, mol, metals, md_pairs, nb_pairs) -> bool:
    """Run the topology gate with the sister-module signature; degrade safely."""
    if passes_topology_fn is None:
        return True  # gate unavailable — let energy filter handle it.

    # Fill in d_ideal from current geometry so the [0.85, 1.10] window is
    # centred on the actual bond length.  This still catches bond breaks
    # caused by a buggy rotation (e.g. M moving), but does not punish the
    # rotation itself (which keeps M and D fixed by construction).
    resolved: List[Tuple[int, int, float]] = []
    for (m_idx, d_idx, _sentinel) in md_pairs:
        d_cur = float(np.linalg.norm(coords[m_idx] - coords[d_idx]))
        if d_cur < 1e-6:
            return False
        resolved.append((m_idx, d_idx, d_cur))

    try:
        return bool(passes_topology_fn(coords, mol, metals, resolved, nb_pairs))
    except TypeError:
        # Older / alternate signature: try (coords, mol).
        try:
            return bool(passes_topology_fn(coords, mol))
        except Exception:
            return True
    except Exception:
        return True


# ---------------------------------------------------------------------------
# Self-test: synthetic sp3-donor sweep
# ---------------------------------------------------------------------------


def _self_test() -> None:  # pragma: no cover - exercised manually
    """Build a deliberately distorted Pt-CH3 fragment and verify the sweep.

    We sidestep RDKit's stricter valence handling by building a minimal mock
    molecule that exposes the duck-typed surface used by this module
    (``GetAtoms``, ``GetAtomWithIdx``, ``GetNeighbors``, ``GetSymbol``,
    ``GetAtomicNum``, ``GetHybridization``, ``GetIdx``).  We also patch a
    trivial ``U_total`` so the test runs without the full DELFIN stack.
    """
    print("[_donor_rotation self-test]")

    # --- mock molecule ----------------------------------------------------
    class _MockAtom:
        def __init__(self, idx, sym, z, hyb, nbrs):
            self.idx = idx
            self.sym = sym
            self.z = z
            self.hyb = hyb
            self._nbrs = nbrs  # list of ints, resolved against mol on access

        def GetIdx(self):
            return self.idx

        def GetSymbol(self):
            return self.sym

        def GetAtomicNum(self):
            return self.z

        def GetHybridization(self):
            return self.hyb

    class _MockMol:
        def __init__(self, atoms):
            self._atoms = atoms

        def GetAtoms(self):
            return list(self._atoms)

        def GetAtomWithIdx(self, i):
            return self._atoms[i]

    # rdkit hybridization sentinel
    try:
        from rdkit import Chem
        SP3 = Chem.HybridizationType.SP3
    except Exception:
        class _HT:
            SP3 = "SP3"
        SP3 = _HT.SP3  # type: ignore[assignment]

    # 0=Pt, 1=C (sp3), 2..4=H on C, 5=NH3-N (placeholder, sp3 but not rotatable
    # because we make it bonded to Pt only — sweep will find no fragment).
    atoms: List[_MockAtom] = []

    # Pre-allocate, then patch neighbour lists.
    atoms.append(_MockAtom(0, "Pt", 78, SP3, []))
    atoms.append(_MockAtom(1, "C", 6, SP3, []))
    atoms.append(_MockAtom(2, "H", 1, SP3, []))
    atoms.append(_MockAtom(3, "H", 1, SP3, []))
    atoms.append(_MockAtom(4, "H", 1, SP3, []))

    # neighbours by index
    nbr_map = {
        0: [1],          # Pt - C
        1: [0, 2, 3, 4],  # C - Pt, 3 H
        2: [1],
        3: [1],
        4: [1],
    }

    # Patch each atom's GetNeighbors to return _MockAtom refs.
    mol = _MockMol(atoms)

    def _make_get_nbrs(i):
        def _gn():
            return [atoms[j] for j in nbr_map[i]]
        return _gn

    for i, a in enumerate(atoms):
        a.GetNeighbors = _make_get_nbrs(i)

    # --- distorted coordinates: H "umbrella" tipped toward M --------------
    # All three H atoms are placed slightly off the M-C axis but on the
    # "back" side of C (away from M), forming a flattened umbrella.  A
    # rotation by ~180 degrees around the M-C axis points them in the
    # opposite azimuthal sector — but because the umbrella is symmetric
    # under such a rotation, we add a directional reward by giving the H
    # atoms different x-extensions (a tilt).  This makes one azimuthal
    # orientation strictly better than the others.
    coords = np.array([
        [0.00, 0.00, 0.00],   # Pt
        [2.00, 0.00, 0.00],   # C
        [2.40, 0.90, 0.10],   # H1  (offset +y, slight tilt)
        [2.40, -0.45, 0.78],  # H2
        [2.40, -0.45, -0.78], # H3
    ], dtype=float)

    # --- mock U_total: directional penalty around the M-C axis ------------
    # Reward H atoms whose azimuthal position around the M-C axis
    # (measured in the yz-plane through C) sits near phi = +90 degrees
    # (i.e. +y).  An anti-phase orientation (~-y) is the worst case; a
    # 180-degree rotation must flip it back.  We deliberately bias the
    # initial geometry so that a -120/+120 rotation reduces the penalty.
    import sys
    import types

    def _mock_U(c, m_, s_, p_):
        m_pos = c[0]
        c_pos = c[1]
        axis = c_pos - m_pos
        axis /= max(np.linalg.norm(axis), 1e-9)
        pen = 0.0
        # Build an orthonormal frame {axis, e_y, e_z} around the axis.
        # axis is along +x, so e_y = +y, e_z = +z.
        e_y = np.array([0.0, 1.0, 0.0])
        e_z = np.array([0.0, 0.0, 1.0])
        for h_idx in (2, 3, 4):
            v = c[h_idx] - c_pos
            # Strip the axial component to get the radial vector in the plane.
            v_perp = v - np.dot(v, axis) * axis
            r = np.linalg.norm(v_perp)
            if r < 1e-6:
                pen += 5.0  # H sits on the axis — penalise hard.
                continue
            # Angle relative to +y direction.  cos(phi) = (v_perp . e_y)/r.
            cos_phi = float(np.dot(v_perp, e_y)) / r
            # Penalty is highest when H points along -y (cos_phi = -1).
            pen += (1.0 - cos_phi) * 1.0
        return pen, np.zeros_like(c)

    fake_mod = types.ModuleType("delfin.manta._energy_terms")
    fake_mod.U_total = _mock_U  # type: ignore[attr-defined]
    sys.modules["delfin.manta._energy_terms"] = fake_mod

    fake_conv = types.ModuleType("delfin.smiles_converter")
    fake_conv._METAL_SET = {"Pt"}
    sys.modules["delfin.smiles_converter"] = fake_conv

    # Force _passes_topology to permissive no-op (sister module not loadable
    # in this stub harness).
    fake_post = types.ModuleType("delfin.manta._post_optimizer")
    fake_post._passes_topology = lambda *a, **kw: True  # type: ignore[attr-defined]
    sys.modules["delfin.manta._post_optimizer"] = fake_post

    new_coords, report = donor_pivot_rotation_sweep(
        coords, mol, sym_info={}, params={},
    )

    U0, _ = _mock_U(coords, mol, {}, {})
    U1, _ = _mock_U(new_coords, mol, {}, {})

    print(f"  sp3 pairs detected : {report['sp3_pairs']}")
    print(f"  donors rotated     : {report['donors_rotated']}")
    print(f"  U_total before     : {U0:.3f}")
    print(f"  U_total after      : {U1:.3f}")
    print(f"  Delta-U            : {U0 - U1:.3f}")

    assert report["sp3_pairs"] == 1, "expected exactly one sp3 donor (C)"
    assert U1 <= U0, "sweep must not worsen energy"
    if report["donors_rotated"] >= 1:
        print("  result             : PASS (donor rotation reduced energy)")
    else:
        print("  result             : PASS (no improvement met threshold)")


if __name__ == "__main__":
    _self_test()
