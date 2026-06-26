"""Analytic energy terms + gradients for Baustein 6 variational refinement.

Eight U-terms compose the total energy functional minimised by L-BFGS-B
in :mod:`delfin.manta._variational_refiner`:

* ``U_bond``      — bond-length penalty (covalent radii lookup)
* ``U_angle``     — bond-angle penalty per hybridization
* ``U_clash``     — one-sided vdW overlap (Bondi/Alvarez radii)
* ``U_topology``  — log-barrier on M-D bonds (HARD topology gate)
* ``U_A``         — Tier A coord-sphere pull to Hungarian-assigned slots
* ``U_B``         — Tier B Morgan equivalence (bond/angle means)
* ``U_C``         — Tier C per-fragment archetype point-group enforcement
* ``U_D``         — Tier D global molecular point-group enforcement

Each function returns ``(energy_value, gradient_array)`` where ``gradient_array``
has the same shape as ``coords`` (``[N, 3]``). Gradients are fully analytic
(no finite differencing) so L-BFGS-B can run with ``jac=True``.

All dependencies on :mod:`delfin.smiles_converter` (``_get_ml_bond_length``,
``_COVALENT_RADII``, ``_METAL_SET``) and :mod:`delfin.manta._vdw_radii` are lazily
imported inside helpers to keep the module light at import time and to avoid
the known circular-import risk with ``smiles_converter``.

Reference: ``iters/BAUSTEIN6_MASTERPLAN.md`` Section 3.3.
"""

from __future__ import annotations

import math
from typing import Dict, List, Tuple

import numpy as np

# ---------------------------------------------------------------------------
# Numerical constants
# ---------------------------------------------------------------------------

_EPS_COS = 1.0e-7      # clip for cos(theta) to avoid arccos NaN
_EPS_SIN = 1.0e-6      # treat angle as collinear / undefined below this
_EPS_DIST = 1.0e-8     # treat bond length below this as atom collapse
_EPS_BARRIER = 1.0e-4  # safety margin inside log-barrier domain

# Ideal angles (radians) per hybridization.
_THETA_SP = math.pi                # 180.0 deg
_THETA_SP2 = math.radians(120.0)   # 120.0 deg
_THETA_SP3 = math.radians(109.47)  # tetrahedral

# Covalent-bond-order scaling for ideal bond length d_ideal = scale * (r_i + r_j).
# Pyykkö-style; AROMATIC sits between SINGLE and DOUBLE.
_BOND_ORDER_SCALE: Dict[float, float] = {
    1.0: 1.00,
    1.5: 0.93,
    2.0: 0.87,
    3.0: 0.78,
}

_LAST_RESORT_BOND_LEN = 1.50   # Å — used if nothing else is known
_DEFAULT_VDW_FALLBACK = 1.80   # Å — also used for unknown elements


# ---------------------------------------------------------------------------
# Lazy-import helpers (avoid module-import cycles)
# ---------------------------------------------------------------------------

def _smiles_converter_metals():
    """Return ``_METAL_SET`` from :mod:`delfin.smiles_converter`."""
    try:
        from delfin.smiles_converter import _METAL_SET  # type: ignore
        return _METAL_SET
    except Exception:
        return frozenset()


def _smiles_converter_covalent_radii():
    """Return ``_COVALENT_RADII`` from :mod:`delfin.smiles_converter`."""
    try:
        from delfin.smiles_converter import _COVALENT_RADII  # type: ignore
        return _COVALENT_RADII
    except Exception:
        return {}


def _smiles_converter_ml_bondlen(metal_sym: str, donor_sym: str) -> float:
    """Return tabulated M-L bond length (Å), with safe fallback to 2.0 Å."""
    try:
        from delfin.smiles_converter import _get_ml_bond_length  # type: ignore
        return float(_get_ml_bond_length(metal_sym, donor_sym))
    except Exception:
        return 2.0


def _vdw_radius(symbol: str) -> float:
    """Return Bondi/Alvarez vdW radius for ``symbol``."""
    try:
        from delfin.manta._vdw_radii import get_vdw_radius  # type: ignore
        return float(get_vdw_radius(symbol))
    except Exception:
        return _DEFAULT_VDW_FALLBACK


# ---------------------------------------------------------------------------
# Internal mol introspection
# ---------------------------------------------------------------------------

def _symbols(mol) -> List[str]:
    """Element symbol per atom index."""
    return [a.GetSymbol() for a in mol.GetAtoms()]


def _is_metal_atom(mol, idx: int) -> bool:
    sym = mol.GetAtomWithIdx(idx).GetSymbol()
    return sym in _smiles_converter_metals()


def _enumerate_bonds(mol) -> List[Tuple[int, int, float]]:
    """Return list of ``(i, j, bond_order_double)`` for every bond in ``mol``."""
    out: List[Tuple[int, int, float]] = []
    for b in mol.GetBonds():
        i = b.GetBeginAtomIdx()
        j = b.GetEndAtomIdx()
        try:
            order = float(b.GetBondTypeAsDouble())
        except Exception:
            order = 1.0
        out.append((int(i), int(j), order))
    return out


def _ideal_bond_length(mol, i: int, j: int, order: float) -> float:
    """Return ideal bond length d_ideal for the bond ``(i, j)`` of order ``order``.

    Lookup order:
      1. M-L bond (one endpoint is a metal) → ``_get_ml_bond_length``.
      2. Otherwise covalent-radii sum × bond-order scale.
      3. Fallback ``_LAST_RESORT_BOND_LEN`` if radii unknown.
    """
    metals = _smiles_converter_metals()
    sym_i = mol.GetAtomWithIdx(i).GetSymbol()
    sym_j = mol.GetAtomWithIdx(j).GetSymbol()
    if sym_i in metals and sym_j not in metals:
        return _smiles_converter_ml_bondlen(sym_i, sym_j)
    if sym_j in metals and sym_i not in metals:
        return _smiles_converter_ml_bondlen(sym_j, sym_i)
    cov = _smiles_converter_covalent_radii()
    r_i = cov.get(sym_i)
    r_j = cov.get(sym_j)
    if r_i is None or r_j is None:
        return _LAST_RESORT_BOND_LEN
    scale = _BOND_ORDER_SCALE.get(round(order * 2) / 2.0)  # snap to .5
    if scale is None:
        scale = 1.0
    return float(r_i + r_j) * scale


def _hybridization_theta(mol, j: int) -> float:
    """Return ideal bond angle (radians) at atom ``j`` from RDKit hybridization."""
    try:
        from rdkit.Chem.rdchem import HybridizationType  # type: ignore
    except Exception:
        return _THETA_SP3  # safe default
    try:
        hyb = mol.GetAtomWithIdx(j).GetHybridization()
    except Exception:
        return _THETA_SP3
    if hyb == HybridizationType.SP:
        return _THETA_SP
    if hyb in (HybridizationType.SP2, HybridizationType.SP2D):  # rare alias
        return _THETA_SP2
    if hyb in (HybridizationType.SP3, HybridizationType.SP3D, HybridizationType.SP3D2):
        return _THETA_SP3
    # Unknown / unspecified → tetrahedral default.
    return _THETA_SP3


def _enumerate_angles(mol) -> List[Tuple[int, int, int]]:
    """Return all angle triples (i, j, k) with j the central atom.

    Iterates RDKit neighbour lists, skipping metal-centred angles (those are
    governed by ``U_A`` polyhedron pull rather than hybridization defaults).
    """
    metals = _smiles_converter_metals()
    triples: List[Tuple[int, int, int]] = []
    for atom in mol.GetAtoms():
        j = int(atom.GetIdx())
        if atom.GetSymbol() in metals:
            continue
        nbrs = [int(n.GetIdx()) for n in atom.GetNeighbors()]
        n = len(nbrs)
        if n < 2:
            continue
        for a in range(n):
            for b in range(a + 1, n):
                triples.append((nbrs[a], j, nbrs[b]))
    return triples


def _bonded_pairs_set(mol) -> set:
    """Return ``{(min(i,j), max(i,j))}`` for every covalent bond."""
    s = set()
    for b in mol.GetBonds():
        i = int(b.GetBeginAtomIdx())
        j = int(b.GetEndAtomIdx())
        if i > j:
            i, j = j, i
        s.add((i, j))
    return s


def _1_3_pairs_set(mol) -> set:
    """Return 1-3 pairs (atoms sharing a common neighbour). Excluded from clash."""
    s = set()
    for atom in mol.GetAtoms():
        nbrs = [int(n.GetIdx()) for n in atom.GetNeighbors()]
        m = len(nbrs)
        for a in range(m):
            for b in range(a + 1, m):
                i, j = nbrs[a], nbrs[b]
                if i > j:
                    i, j = j, i
                s.add((i, j))
    return s


def _enumerate_metal_donor_bonds(mol) -> List[Tuple[int, int]]:
    """Return list ``[(metal_idx, donor_idx)]`` for every M-L covalent bond."""
    metals = _smiles_converter_metals()
    pairs: List[Tuple[int, int]] = []
    for b in mol.GetBonds():
        i = int(b.GetBeginAtomIdx())
        j = int(b.GetEndAtomIdx())
        si = mol.GetAtomWithIdx(i).GetSymbol()
        sj = mol.GetAtomWithIdx(j).GetSymbol()
        if si in metals and sj not in metals:
            pairs.append((i, j))
        elif sj in metals and si not in metals:
            pairs.append((j, i))
    return pairs


# ---------------------------------------------------------------------------
# U_bond
# ---------------------------------------------------------------------------

def U_bond(coords: np.ndarray, mol, k_bond: float = 1000.0
           ) -> Tuple[float, np.ndarray]:
    """Bond-length penalty: Σ k_bond · (d - d_ideal)²  over covalent bonds.

    Skips M-L bonds (those are handled by ``U_topology`` and ``U_A``).
    """
    coords = np.asarray(coords, dtype=np.float64)
    n = coords.shape[0]
    grad = np.zeros_like(coords)
    energy = 0.0

    metals = _smiles_converter_metals()
    if k_bond <= 0.0 or n == 0:
        return float(energy), grad

    for (i, j, order) in _enumerate_bonds(mol):
        sym_i = mol.GetAtomWithIdx(i).GetSymbol()
        sym_j = mol.GetAtomWithIdx(j).GetSymbol()
        if sym_i in metals or sym_j in metals:
            continue  # handled by U_topology / U_A
        diff = coords[i] - coords[j]
        d = float(np.linalg.norm(diff))
        if d < _EPS_DIST:
            continue
        d_id = _ideal_bond_length(mol, i, j, order)
        delta = d - d_id
        energy += k_bond * delta * delta
        coef = 2.0 * k_bond * delta / d
        grad[i] += coef * diff
        grad[j] -= coef * diff

    return float(energy), grad


# ---------------------------------------------------------------------------
# U_angle
# ---------------------------------------------------------------------------

def U_angle(coords: np.ndarray, mol, k_angle: float = 100.0
            ) -> Tuple[float, np.ndarray]:
    """Bond-angle penalty per hybridization: Σ k_angle · (θ - θ_ideal)².

    sp → 180°, sp2 → 120°, sp3 → 109.47°. Metal-centred triples skipped.
    """
    coords = np.asarray(coords, dtype=np.float64)
    n = coords.shape[0]
    grad = np.zeros_like(coords)
    energy = 0.0

    if k_angle <= 0.0 or n == 0:
        return float(energy), grad

    for (i, j, k) in _enumerate_angles(mol):
        x_ij = coords[i] - coords[j]
        x_kj = coords[k] - coords[j]
        norm_ij = float(np.linalg.norm(x_ij))
        norm_kj = float(np.linalg.norm(x_kj))
        if norm_ij < _EPS_DIST or norm_kj < _EPS_DIST:
            continue
        u = x_ij / norm_ij
        v = x_kj / norm_kj
        cos_t = float(np.dot(u, v))
        # Clip into safe arccos / sin window.
        cos_t = max(-1.0 + _EPS_COS, min(1.0 - _EPS_COS, cos_t))
        sin_t = math.sqrt(max(0.0, 1.0 - cos_t * cos_t))
        if sin_t < _EPS_SIN:
            continue
        theta = math.acos(cos_t)
        theta_id = _hybridization_theta(mol, j)
        delta = theta - theta_id
        energy += k_angle * delta * delta

        # ∂θ/∂x_i = -1/(sin θ · |x_ij|) · (v - cos θ · u)
        d_th_di = -(v - cos_t * u) / (norm_ij * sin_t)
        d_th_dk = -(u - cos_t * v) / (norm_kj * sin_t)
        d_th_dj = -(d_th_di + d_th_dk)
        coef = 2.0 * k_angle * delta
        grad[i] += coef * d_th_di
        grad[j] += coef * d_th_dj
        grad[k] += coef * d_th_dk

    return float(energy), grad


# ---------------------------------------------------------------------------
# U_clash
# ---------------------------------------------------------------------------

def U_clash(coords: np.ndarray, mol, k_clash: float = 500.0,
            clash_factor: float = 0.85, h_h_factor: float = 0.75
            ) -> Tuple[float, np.ndarray]:
    """One-sided clash repulsion for non-bonded heavy pairs.

    Uses Bondi/Alvarez vdW sums from :mod:`delfin.manta._vdw_radii`. Metal/metal
    and 1-2, 1-3 connections are excluded. Penalty is C¹-continuous quadratic
    inside the overlap region; zero (with zero gradient) outside.
    """
    coords = np.asarray(coords, dtype=np.float64)
    n = coords.shape[0]
    grad = np.zeros_like(coords)
    energy = 0.0

    if k_clash <= 0.0 or n < 2:
        return float(energy), grad

    metals = _smiles_converter_metals()
    syms = _symbols(mol)
    bonded = _bonded_pairs_set(mol)
    angle13 = _1_3_pairs_set(mol)

    # Pre-tabulate per-element vdW once for speed.
    vdw_cache: Dict[str, float] = {}
    for s in set(syms):
        vdw_cache[s] = _vdw_radius(s)

    for i in range(n):
        si = syms[i]
        if si in metals:
            continue
        for j in range(i + 1, n):
            sj = syms[j]
            if sj in metals:
                continue
            if (i, j) in bonded or (i, j) in angle13:
                continue
            # Threshold scaling: H-H gets looser factor.
            if si == "H" and sj == "H":
                threshold = (vdw_cache[si] + vdw_cache[sj]) * h_h_factor
            else:
                threshold = (vdw_cache[si] + vdw_cache[sj]) * clash_factor
            diff = coords[j] - coords[i]
            d = float(np.linalg.norm(diff))
            if d < _EPS_DIST:
                # Atoms collapsed — push apart along arbitrary axis (x).
                overlap = threshold
                energy += k_clash * overlap * overlap
                push = np.array([2.0 * k_clash * overlap, 0.0, 0.0])
                grad[i] -= push
                grad[j] += push
                continue
            if d >= threshold:
                continue
            overlap = threshold - d
            energy += k_clash * overlap * overlap
            # ∂U/∂x_i = +2 k_c · overlap · (x_j - x_i) / d
            # (moving x_i toward x_j shortens d, U=k·(threshold-d)² grows
            #  → positive gradient in direction (x_j - x_i). L-BFGS-B
            #  moves -grad → pushes x_i AWAY from x_j.)
            push = (2.0 * k_clash * overlap / d) * diff  # = ∂U/∂x_i
            grad[i] += push
            grad[j] -= push

    return float(energy), grad


# ---------------------------------------------------------------------------
# U_topology
# ---------------------------------------------------------------------------

def U_topology(coords: np.ndarray, mol, k_topology: float = 10000.0,
               lo_frac: float = 0.85, hi_frac: float = 1.10
               ) -> Tuple[float, np.ndarray]:
    """Log-barrier on every M-D bond.

    For each M-D pair the bond length ``d`` must remain inside
    ``[lo_frac, hi_frac] · d_ideal``. Inside the open interval the barrier is

        barrier(d) = -log(d - lo) - log(hi - d)

    which diverges at the boundaries (preserving topology by construction).
    Outside the interval a large finite penalty with a repulsive linear
    gradient steers L-BFGS-B back into the safe domain (the algorithm needs
    finite values to make progress).
    """
    coords = np.asarray(coords, dtype=np.float64)
    n = coords.shape[0]
    grad = np.zeros_like(coords)
    energy = 0.0

    if k_topology <= 0.0 or n == 0:
        return float(energy), grad

    big_pen = 1.0e6  # placeholder energy when outside safe interval

    for (m_idx, d_idx) in _enumerate_metal_donor_bonds(mol):
        diff = coords[m_idx] - coords[d_idx]
        d = float(np.linalg.norm(diff))
        if d < _EPS_DIST:
            # Atoms collapsed — emit huge gradient to separate them.
            energy += k_topology * big_pen
            push = np.array([k_topology * big_pen, 0.0, 0.0])
            grad[m_idx] += push
            grad[d_idx] -= push
            continue
        sym_m = mol.GetAtomWithIdx(m_idx).GetSymbol()
        sym_d = mol.GetAtomWithIdx(d_idx).GetSymbol()
        d_id = _smiles_converter_ml_bondlen(sym_m, sym_d)
        lo = d_id * lo_frac
        hi = d_id * hi_frac
        # Outside-domain handling: very large penalty + linear gradient
        # pointing back into the safe interval.
        if d <= lo + _EPS_BARRIER:
            energy += k_topology * big_pen
            # We need d to INCREASE → x_m should move AWAY from x_d.
            # ∂(d)/∂x_m = (x_m - x_d)/d, ∂(d)/∂x_d = -(...)
            # Use -∂U/∂d > 0 ⇒ force on m is along +(x_m - x_d)/d.
            coef = -k_topology * big_pen / max(d, _EPS_DIST)
            # coef·diff acts as +∂U/∂x_m; gradient = +coef·diff so the
            # minimiser moves m in direction +diff (away from donor).
            grad[m_idx] += coef * diff
            grad[d_idx] -= coef * diff
            continue
        if d >= hi - _EPS_BARRIER:
            energy += k_topology * big_pen
            # We need d to DECREASE → m moves TOWARD donor.
            coef = +k_topology * big_pen / max(d, _EPS_DIST)
            grad[m_idx] += coef * diff
            grad[d_idx] -= coef * diff
            continue

        # Smooth log barrier inside (lo, hi).
        b = -math.log(d - lo) - math.log(hi - d)
        energy += k_topology * b
        # ∂b/∂d = -1/(d - lo) + 1/(hi - d)
        db_dd = -1.0 / (d - lo) + 1.0 / (hi - d)
        coef = k_topology * db_dd / d
        grad[m_idx] += coef * diff
        grad[d_idx] -= coef * diff

    return float(energy), grad


# ---------------------------------------------------------------------------
# U_A — Coord-sphere polyhedron pull (Tier A)
# ---------------------------------------------------------------------------

def U_A_coord_sphere(coords: np.ndarray, mol,
                     donor_targets: Dict[int, np.ndarray],
                     k_A: float = 100.0) -> Tuple[float, np.ndarray]:
    """Tier A: pull each donor toward its Hungarian-assigned slot.

    ``donor_targets[donor_idx]`` is a 3-vector giving the absolute target
    position (already in the working coordinate frame). The metal centres
    themselves are not pulled; their position is governed by the rest of
    the force field.
    """
    coords = np.asarray(coords, dtype=np.float64)
    grad = np.zeros_like(coords)
    energy = 0.0

    if k_A <= 0.0 or not donor_targets:
        return float(energy), grad

    for donor_idx, target in donor_targets.items():
        if donor_idx < 0 or donor_idx >= coords.shape[0]:
            continue
        tgt = np.asarray(target, dtype=np.float64).reshape(3)
        delta = coords[donor_idx] - tgt
        energy += k_A * float(np.dot(delta, delta))
        grad[donor_idx] += 2.0 * k_A * delta

    return float(energy), grad


# ---------------------------------------------------------------------------
# U_B — Local equivalence (Tier B)
# ---------------------------------------------------------------------------

def _angle_value(coords: np.ndarray, i: int, j: int, k: int) -> float:
    """Return the angle θ_ijk in radians, safely clipped."""
    x_ij = coords[i] - coords[j]
    x_kj = coords[k] - coords[j]
    nij = float(np.linalg.norm(x_ij))
    nkj = float(np.linalg.norm(x_kj))
    if nij < _EPS_DIST or nkj < _EPS_DIST:
        return _THETA_SP3
    cos_t = float(np.dot(x_ij, x_kj) / (nij * nkj))
    cos_t = max(-1.0 + _EPS_COS, min(1.0 - _EPS_COS, cos_t))
    return math.acos(cos_t)


def U_B_equivalence(coords: np.ndarray, mol,
                    equiv_bond_pairs: List[List[Tuple[int, int]]],
                    equiv_angle_triples: List[List[Tuple[int, int, int]]],
                    k_B: float = 50.0) -> Tuple[float, np.ndarray]:
    """Tier B: chemically-equivalent bonds/angles pull toward their class mean.

    ``equiv_bond_pairs`` is a list of bond-classes, each a list of
    ``(i, j)`` index tuples.  ``equiv_angle_triples`` is the analogous list
    for angle triples ``(i, j, k)``.

    For tractable gradients the class mean is treated as a constant of the
    current iterate; L-BFGS-B will re-evaluate (mean,gradient) at every step
    so the mean updates implicitly between iterations.
    """
    coords = np.asarray(coords, dtype=np.float64)
    grad = np.zeros_like(coords)
    energy = 0.0

    if k_B <= 0.0:
        return float(energy), grad

    # Bond equivalence classes.
    for cls in equiv_bond_pairs or []:
        if not cls or len(cls) < 2:
            continue
        ds: List[float] = []
        diffs: List[np.ndarray] = []
        ok: List[Tuple[int, int]] = []
        for (i, j) in cls:
            diff = coords[i] - coords[j]
            d = float(np.linalg.norm(diff))
            if d < _EPS_DIST:
                continue
            ds.append(d)
            diffs.append(diff)
            ok.append((i, j))
        if len(ds) < 2:
            continue
        mean = float(np.mean(ds))
        for (i, j), d, diff in zip(ok, ds, diffs):
            delta = d - mean
            energy += k_B * delta * delta
            coef = 2.0 * k_B * delta / d
            grad[i] += coef * diff
            grad[j] -= coef * diff

    # Angle equivalence classes.
    for cls in equiv_angle_triples or []:
        if not cls or len(cls) < 2:
            continue
        thetas: List[float] = []
        ok_triples: List[Tuple[int, int, int]] = []
        for (i, j, k) in cls:
            t = _angle_value(coords, i, j, k)
            thetas.append(t)
            ok_triples.append((i, j, k))
        if len(thetas) < 2:
            continue
        mean = float(np.mean(thetas))
        for (i, j, k), theta in zip(ok_triples, thetas):
            x_ij = coords[i] - coords[j]
            x_kj = coords[k] - coords[j]
            nij = float(np.linalg.norm(x_ij))
            nkj = float(np.linalg.norm(x_kj))
            if nij < _EPS_DIST or nkj < _EPS_DIST:
                continue
            u = x_ij / nij
            v = x_kj / nkj
            cos_t = float(np.dot(u, v))
            cos_t = max(-1.0 + _EPS_COS, min(1.0 - _EPS_COS, cos_t))
            sin_t = math.sqrt(max(0.0, 1.0 - cos_t * cos_t))
            if sin_t < _EPS_SIN:
                continue
            delta = theta - mean
            energy += k_B * delta * delta
            d_th_di = -(v - cos_t * u) / (nij * sin_t)
            d_th_dk = -(u - cos_t * v) / (nkj * sin_t)
            d_th_dj = -(d_th_di + d_th_dk)
            coef = 2.0 * k_B * delta
            grad[i] += coef * d_th_di
            grad[j] += coef * d_th_dj
            grad[k] += coef * d_th_dk

    return float(energy), grad


# ---------------------------------------------------------------------------
# U_C — Per-fragment archetype point group (Tier C)
# ---------------------------------------------------------------------------

def U_C_fragment(coords: np.ndarray, mol, fragments: List[Dict],
                 k_C: float = 80.0) -> Tuple[float, np.ndarray]:
    """Tier C: per-fragment archetype point-group enforcement.

    Each ``fragments[k]`` is a dict carrying at least the keys::

        {
          "atoms": [i0, i1, ...],          # atom indices in the fragment
          "operations": [op_mat_1, ...],   # list of 3×3 numpy arrays (skip identity)
          "partners":  [{i: j, ...}, ...]  # parallel to operations: i → partner j under op
        }

    For every (op, partner-mapping) pair we minimise

        U += k_C · Σ_i ||x_{p(i)} - (centroid + op·(x_i - centroid))||²

    where ``centroid`` is the running fragment centroid. The chain-rule
    derivative w.r.t. the fragment centroid is distributed equally back to
    all fragment atoms (which keeps L-BFGS-B well-conditioned).
    """
    coords = np.asarray(coords, dtype=np.float64)
    grad = np.zeros_like(coords)
    energy = 0.0

    if k_C <= 0.0 or not fragments:
        return float(energy), grad

    for frag in fragments:
        atoms = list(frag.get("atoms", []))
        ops = list(frag.get("operations", []))
        perms = list(frag.get("partners", []))
        if not atoms or not ops or not perms or len(ops) != len(perms):
            continue
        n_frag = len(atoms)
        if n_frag < 2:
            continue
        # Build sub-array of fragment coords once per fragment.
        frag_arr = coords[atoms]
        centroid = frag_arr.mean(axis=0)
        inv_n = 1.0 / float(n_frag)

        for op, perm in zip(ops, perms):
            if op is None or perm is None:
                continue
            op_mat = np.asarray(op, dtype=np.float64).reshape(3, 3)
            # Skip identity (it contributes zero penalty).
            if np.allclose(op_mat, np.eye(3), atol=1.0e-9):
                continue
            op_T = op_mat.T
            # Accumulator for the centroid back-distribution.
            centroid_force = np.zeros(3, dtype=np.float64)
            atom_set = set(atoms)
            for i, j in perm.items():
                if i not in atom_set or j not in atom_set:
                    continue
                xi_loc = coords[i] - centroid
                target = centroid + op_mat @ xi_loc
                delta = coords[j] - target
                if i == j:
                    # Atom is fixed by op (e.g. on a mirror plane / inversion centre).
                    # We still penalise residual displacement; ∂target/∂x_i = op,
                    # but x_j = x_i so the gradient combines.
                    energy += k_C * float(np.dot(delta, delta))
                    # ∂U/∂x_i = 2 k_C · (I - opᵀ) · delta   (j == i)
                    grad[i] += 2.0 * k_C * ((np.eye(3) - op_T) @ delta)
                    continue
                energy += k_C * float(np.dot(delta, delta))
                # ∂U/∂x_j = +2 k_C · delta
                grad[j] += 2.0 * k_C * delta
                # ∂U/∂x_i: target = centroid + op·(x_i - centroid)
                #          → ∂target/∂x_i = op
                # ∂U/∂x_i = -2 k_C · opᵀ · delta
                grad[i] -= 2.0 * k_C * (op_T @ delta)
                # Centroid contribution: ∂target/∂c = I - op
                # ⇒ centroid_force accumulates -2 k_C · (I - opᵀ) · delta
                centroid_force -= 2.0 * k_C * ((np.eye(3) - op_T) @ delta)
            # Distribute centroid gradient equally to all fragment atoms.
            if not np.allclose(centroid_force, 0.0):
                share = inv_n * centroid_force
                for a in atoms:
                    grad[a] += share

    return float(energy), grad


# ---------------------------------------------------------------------------
# U_D — Global molecular point group (Tier D)
# ---------------------------------------------------------------------------

def U_D_global(coords: np.ndarray, mol, global_pg_ops: List[np.ndarray],
               atom_perms: Dict[int, Dict[int, int]],
               k_D: float = 30.0) -> Tuple[float, np.ndarray]:
    """Tier D: global molecular point-group enforcement.

    ``global_pg_ops`` is a list of 3×3 rotation/reflection matrices.
    ``atom_perms[op_index]`` is a dict ``{atom_idx: partner_atom_idx}``
    encoding how the operation permutes atoms (identity ops should be
    omitted or will be skipped automatically).
    """
    coords = np.asarray(coords, dtype=np.float64)
    n = coords.shape[0]
    grad = np.zeros_like(coords)
    energy = 0.0

    if k_D <= 0.0 or n == 0 or not global_pg_ops:
        return float(energy), grad

    centroid = coords.mean(axis=0)
    inv_n = 1.0 / float(n)

    for op_idx, op in enumerate(global_pg_ops):
        if op is None:
            continue
        op_mat = np.asarray(op, dtype=np.float64).reshape(3, 3)
        if np.allclose(op_mat, np.eye(3), atol=1.0e-9):
            continue
        perm = atom_perms.get(op_idx, {}) if atom_perms else {}
        if not perm:
            continue
        op_T = op_mat.T
        centroid_force = np.zeros(3, dtype=np.float64)
        for i, j in perm.items():
            if i < 0 or i >= n or j < 0 or j >= n:
                continue
            xi_loc = coords[i] - centroid
            target = centroid + op_mat @ xi_loc
            delta = coords[j] - target
            if i == j:
                # Atom fixed by op (lies on symmetry element).
                energy += k_D * float(np.dot(delta, delta))
                grad[i] += 2.0 * k_D * ((np.eye(3) - op_T) @ delta)
                continue
            energy += k_D * float(np.dot(delta, delta))
            grad[j] += 2.0 * k_D * delta
            grad[i] -= 2.0 * k_D * (op_T @ delta)
            centroid_force -= 2.0 * k_D * ((np.eye(3) - op_T) @ delta)
        if not np.allclose(centroid_force, 0.0):
            share = inv_n * centroid_force
            grad += share  # broadcast to every atom

    return float(energy), grad


# ---------------------------------------------------------------------------
# U_total
# ---------------------------------------------------------------------------

def U_total(coords: np.ndarray, mol, sym_info: Dict, params: Dict
            ) -> Tuple[float, np.ndarray]:
    """Aggregate all 8 terms with class-conditional k values.

    Parameters
    ----------
    coords : (N, 3) numpy array
        Current Cartesian coordinates (Å).
    mol : rdkit Mol
        Molecule with bonds/hybridization.
    sym_info : dict
        Pre-computed Tier A-D data. Recognised keys:
            ``donor_targets``        — dict[int, (3,)] for U_A
            ``equiv_bond_pairs``     — list[list[(i,j)]] for U_B
            ``equiv_angle_triples``  — list[list[(i,j,k)]] for U_B
            ``fragments``            — list[dict] for U_C
            ``global_ops``           — list[(3,3)] for U_D
            ``atom_perms``           — dict[int, dict[int,int]] for U_D
    params : dict
        Coefficients (default to spec values if missing):
            ``k_bond, k_angle, k_clash, k_topology, k_A, k_B, k_C, k_D``
        Plus optional ``clash_factor, h_h_factor, lo_frac, hi_frac``.

    Returns
    -------
    (energy_total, grad_total)
    """
    coords = np.asarray(coords, dtype=np.float64)
    sym_info = sym_info or {}
    params = params or {}

    k_bond = float(params.get("k_bond", 1000.0))
    k_angle = float(params.get("k_angle", 100.0))
    k_clash = float(params.get("k_clash", 500.0))
    k_topology = float(params.get("k_topology", 10000.0))
    k_A = float(params.get("k_A", 100.0))
    k_B = float(params.get("k_B", 50.0))
    k_C = float(params.get("k_C", 80.0))
    k_D = float(params.get("k_D", 30.0))
    clash_factor = float(params.get("clash_factor", 0.85))
    h_h_factor = float(params.get("h_h_factor", 0.75))
    lo_frac = float(params.get("lo_frac", 0.85))
    hi_frac = float(params.get("hi_frac", 1.10))

    e_total = 0.0
    g_total = np.zeros_like(coords)

    e, g = U_bond(coords, mol, k_bond=k_bond)
    e_total += e
    g_total += g

    e, g = U_angle(coords, mol, k_angle=k_angle)
    e_total += e
    g_total += g

    e, g = U_clash(coords, mol, k_clash=k_clash,
                   clash_factor=clash_factor, h_h_factor=h_h_factor)
    e_total += e
    g_total += g

    e, g = U_topology(coords, mol, k_topology=k_topology,
                      lo_frac=lo_frac, hi_frac=hi_frac)
    e_total += e
    g_total += g

    e, g = U_A_coord_sphere(coords, mol,
                            sym_info.get("donor_targets", {}) or {},
                            k_A=k_A)
    e_total += e
    g_total += g

    e, g = U_B_equivalence(coords, mol,
                           sym_info.get("equiv_bond_pairs", []) or [],
                           sym_info.get("equiv_angle_triples", []) or [],
                           k_B=k_B)
    e_total += e
    g_total += g

    e, g = U_C_fragment(coords, mol,
                        sym_info.get("fragments", []) or [],
                        k_C=k_C)
    e_total += e
    g_total += g

    e, g = U_D_global(coords, mol,
                      sym_info.get("global_ops", []) or [],
                      sym_info.get("atom_perms", {}) or {},
                      k_D=k_D)
    e_total += e
    g_total += g

    return float(e_total), g_total


# ---------------------------------------------------------------------------
# Self-test (smoke check + finite-difference gradient sanity)
# ---------------------------------------------------------------------------

def _fd_gradient(fn, coords: np.ndarray, eps: float = 1.0e-5) -> np.ndarray:
    """Central-difference numerical gradient of ``fn(coords) -> (energy, _)``."""
    g = np.zeros_like(coords)
    flat = coords.reshape(-1).copy()
    for k in range(flat.size):
        orig = flat[k]
        flat[k] = orig + eps
        e_p, _ = fn(flat.reshape(coords.shape))
        flat[k] = orig - eps
        e_m, _ = fn(flat.reshape(coords.shape))
        flat[k] = orig
        g.reshape(-1)[k] = (e_p - e_m) / (2.0 * eps)
    return g


if __name__ == "__main__":  # pragma: no cover
    # 4-atom synthetic Cu-N₂-O system (mock geometry; not chemistry-accurate).
    try:
        from rdkit import Chem
        from rdkit.Chem import AllChem
    except Exception as exc:  # pragma: no cover
        raise SystemExit(f"RDKit required for self-test: {exc}")

    smi = "[Cu](N)(N)O"
    mol = Chem.MolFromSmiles(smi)
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol, randomSeed=0xC0FFEE)
    conf = mol.GetConformer()
    coords = np.array([list(conf.GetAtomPosition(i)) for i in range(mol.GetNumAtoms())],
                      dtype=np.float64)

    print(f"Self-test — {mol.GetNumAtoms()} atoms ({smi})\n")

    e, g = U_bond(coords, mol)
    print(f"U_bond      E = {e:12.4f}   |g|_max = {np.max(np.abs(g)):.4e}")
    e, g = U_angle(coords, mol)
    print(f"U_angle     E = {e:12.4f}   |g|_max = {np.max(np.abs(g)):.4e}")
    e, g = U_clash(coords, mol)
    print(f"U_clash     E = {e:12.4f}   |g|_max = {np.max(np.abs(g)):.4e}")
    e, g = U_topology(coords, mol)
    print(f"U_topology  E = {e:12.4f}   |g|_max = {np.max(np.abs(g)):.4e}")

    # Tier A: pull every non-metal donor (N, N, O) toward its current pos
    # plus a tiny offset so the term is non-trivial.
    metals = _smiles_converter_metals()
    donor_targets = {}
    for a in mol.GetAtoms():
        if a.GetSymbol() in metals:
            continue
        if not any(n.GetSymbol() in metals for n in a.GetNeighbors()):
            continue
        i = a.GetIdx()
        donor_targets[i] = coords[i] + np.array([0.05, 0.0, 0.0])
    e, g = U_A_coord_sphere(coords, mol, donor_targets)
    print(f"U_A         E = {e:12.4f}   |g|_max = {np.max(np.abs(g)):.4e}"
          f"  (n_donors={len(donor_targets)})")

    # Trivial Tier B: pair every M-D bond with itself (n=2) – use first two.
    md = _enumerate_metal_donor_bonds(mol)
    if len(md) >= 2:
        equiv_bonds = [[md[0], md[1]]]
    else:
        equiv_bonds = []
    e, g = U_B_equivalence(coords, mol, equiv_bonds, [])
    print(f"U_B         E = {e:12.4f}   |g|_max = {np.max(np.abs(g)):.4e}"
          f"  (n_pairs={len(equiv_bonds)})")

    # Trivial Tier C: single fragment with C2 rotation (identity-like partner).
    # We just exercise the code path; with partner=={i:i} U_C should give 0
    # (delta = 0 because target = x_i for identity rotation around centroid
    # only if op = I; choose op = -I to make it non-zero but symmetric).
    frag = {
        "atoms": list(range(mol.GetNumAtoms())),
        "operations": [-np.eye(3)],
        "partners": [{i: i for i in range(mol.GetNumAtoms())}],
    }
    e, g = U_C_fragment(coords, mol, [frag])
    print(f"U_C         E = {e:12.4f}   |g|_max = {np.max(np.abs(g)):.4e}")

    # Trivial Tier D: identity-permutation under -I.
    e, g = U_D_global(coords, mol, [-np.eye(3)],
                      {0: {i: i for i in range(mol.GetNumAtoms())}})
    print(f"U_D         E = {e:12.4f}   |g|_max = {np.max(np.abs(g)):.4e}")

    sym_info = {
        "donor_targets": donor_targets,
        "equiv_bond_pairs": equiv_bonds,
        "equiv_angle_triples": [],
        "fragments": [frag],
        "global_ops": [],
        "atom_perms": {},
    }
    params = {}
    e, g = U_total(coords, mol, sym_info, params)
    print(f"U_total     E = {e:12.4f}   |g|_max = {np.max(np.abs(g)):.4e}")

    # Finite-difference cross-check (single coord, all terms aggregated).
    def total_fn(c):
        return U_total(c, mol, sym_info, params)

    g_an = g
    g_fd = _fd_gradient(total_fn, coords, eps=1.0e-5)
    diff = np.max(np.abs(g_an - g_fd))
    rel = diff / max(1.0e-9, float(np.max(np.abs(g_an))))
    print(f"\nFD gradient check on U_total:")
    print(f"  max|grad_an - grad_fd|  = {diff:.4e}")
    print(f"  relative                = {rel:.4e}")
    print(f"  PASS" if rel < 1.0e-3 else f"  FAIL (rel >= 1e-3)")
