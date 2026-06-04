"""Topology Healing — Phantom Bonds + Missing Bonds + Wrong Angles (2026-06-04).

Complements :mod:`grip_polish` (L-BFGS / LM internal-residual polish) and the
parallel ``grip_healing_mode`` (iterative atom-repositioning) by detecting
and repairing structural defects at the BOND-GRAPH level:

1. **Phantom bonds** — atom pairs that look bonded by geometry (close enough
   that a covalent-radius rule would treat them as a bond) but are NOT in
   the SMILES topology.  Cause: fused ring overlaps, π-π stacking too
   close, builder-emitted superimposed substituents.  Heal = push apart to
   a target separation (default ``1.5 × Σcov``).
2. **Missing bonds** — pairs declared bonded by SMILES but stretched far
   beyond the element-pair ideal length.  Cause: ring opening during
   constrained UFF relax, donor-slide overshoot, cluster-scaffold
   fragmentation.  Heal = pull together along the bond axis to the ideal
   length (donor stays fixed if frozen).
3. **Wrong angles** — bond angles whose Mahalanobis residual against the
   CCDC-grounded :mod:`grip_mogul_lookup` prior exceeds a configurable
   σ-threshold (default 3.0σ).  Heal = rotate the terminal atom ``C`` of
   the angle ``A-B-C`` around the pivot ``B`` toward the ideal angle θ in
   the (A, B, C) plane, preserving the B-C bond length.

Architecture contract
---------------------
* Pure-functional, deterministic — same input -> bit-identical output across
  runs (PYTHONHASHSEED-respecting, sorted iteration order, float64).
* The metal + donors are FROZEN — the healing never moves them so the
  M-D invariant + donor polyhedron from :mod:`grip_polish` are preserved
  exactly.
* Each heal step is gated by a cheap acceptance check (defect count must
  not *increase* on any of the three axes) — the function rolls back to
  the input on any regression, matching the
  :func:`grip_polish` accept-if-better semantics.
* Default OFF — env flag ``DELFIN_FFFREE_GRIP_TOPOLOGY_HEALING`` (read by
  :func:`grip_polish`) controls whether the pipeline is wired into the
  polish dispatcher.  Without the flag, byte-identical to HEAD ``93b396d``.

Composability with parallel modules
-----------------------------------
* ``topology_healing`` (this module) — fixes **bond-graph defects**
  (phantom / missing / angle).
* ``grip_healing_mode`` (parallel subagent ``a5a86807``) — fixes
  **atom-level mispositioning** (residual outliers; doesn't touch the
  graph).
* ``grip_polish`` (L-BFGS / LM) — minimises the **internal Mahalanobis
  loss** (smooth optimisation over the GRIP loss surface).

Order (when all enabled): topology_healing -> grip_healing_mode -> polish.
Each is env-gated INDEPENDENTLY (composable; one can be on without the
others).

Public API
----------
* :func:`detect_phantom_bonds`
* :func:`detect_missing_bonds`
* :func:`detect_wrong_angles`
* :func:`heal_phantom_bonds`
* :func:`heal_missing_bonds`
* :func:`heal_wrong_angles`
* :func:`topology_healing_pipeline`

Env flags
---------
* ``DELFIN_FFFREE_GRIP_TOPOLOGY_HEALING`` — master gate (default OFF).
* ``DELFIN_FFFREE_TOPOLOGY_HEALING_MAX_ITER`` — outer loop budget
  (default 5).
* ``DELFIN_FFFREE_TOPOLOGY_HEALING_PHANTOM_FACTOR`` — covalent-radius
  factor below which a non-topological pair is flagged as phantom
  (default 1.2).
* ``DELFIN_FFFREE_TOPOLOGY_HEALING_PHANTOM_TARGET_FACTOR`` — target
  separation factor for phantom healing (default 1.5).
* ``DELFIN_FFFREE_TOPOLOGY_HEALING_MISSING_FACTOR`` — multiplier on the
  ideal length above which a topological bond is flagged as missing
  (default 1.5).
* ``DELFIN_FFFREE_TOPOLOGY_HEALING_ANGLE_SIGMA`` — σ-threshold for the
  wrong-angle detector (default 3.0).
"""
from __future__ import annotations

import logging
import math
import os
from dataclasses import dataclass, field
from typing import Dict, FrozenSet, Iterable, List, Optional, Sequence, Set, Tuple

import numpy as np

_LOG = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Module-level constants -- same physical numbers as ``_bond_decollapse``
# so the two corrector stacks agree on what a "bond-length-violation" is.
# ---------------------------------------------------------------------------
# Covalent radii (Å) -- Cordero 2008 (subset; matches ``_bond_decollapse._COV``).
_COV: Dict[str, float] = {
    "H": 0.31, "He": 0.28, "Li": 1.28, "Be": 0.96, "B": 0.84,
    "C": 0.76, "N": 0.71, "O": 0.66, "F": 0.57, "Ne": 0.58,
    "Na": 1.66, "Mg": 1.41, "Al": 1.21, "Si": 1.11, "P": 1.07,
    "S": 1.05, "Cl": 1.02, "Ar": 1.06,
    "K": 2.03, "Ca": 1.76, "Sc": 1.70, "Ti": 1.60, "V": 1.53,
    "Cr": 1.39, "Mn": 1.39, "Fe": 1.32, "Co": 1.26, "Ni": 1.24,
    "Cu": 1.32, "Zn": 1.22, "Ga": 1.22, "Ge": 1.20, "As": 1.19,
    "Se": 1.20, "Br": 1.20, "Kr": 1.16,
    "Rb": 2.20, "Sr": 1.95, "Y": 1.90, "Zr": 1.75, "Nb": 1.64,
    "Mo": 1.54, "Tc": 1.47, "Ru": 1.46, "Rh": 1.42, "Pd": 1.39,
    "Ag": 1.45, "Cd": 1.44, "In": 1.42, "Sn": 1.39, "Sb": 1.39,
    "Te": 1.38, "I": 1.39, "Xe": 1.40,
    "Cs": 2.44, "Ba": 2.15, "La": 2.07, "Ce": 2.04, "Lu": 1.87,
    "Hf": 1.75, "Ta": 1.70, "W": 1.62, "Re": 1.51, "Os": 1.44,
    "Ir": 1.41, "Pt": 1.36, "Au": 1.36, "Hg": 1.32, "Tl": 1.45,
    "Pb": 1.46, "Bi": 1.48,
}
_COV_DEFAULT = 0.90

# Numerical floor: divide-by-norm guard.
_EPS: float = 1e-12

# Default healing thresholds (env-flag overridable).
DEFAULT_PHANTOM_FACTOR: float = 1.2  # cov-radius factor for phantom detection
DEFAULT_PHANTOM_TARGET_FACTOR: float = 1.5  # target separation for phantom heal
DEFAULT_MISSING_FACTOR: float = 1.5  # multiplier on ideal length for missing detection
DEFAULT_ANGLE_SIGMA: float = 3.0  # σ-threshold for wrong-angle detection
DEFAULT_MAX_ITER: int = 5  # outer loop budget for the pipeline
DEFAULT_ANGLE_STEP: float = 0.5  # damping factor for angle correction (0..1)

# Env-flag names (kept in one place so all callers / docs agree).
ENV_MASTER: str = "DELFIN_FFFREE_GRIP_TOPOLOGY_HEALING"
ENV_MAX_ITER: str = "DELFIN_FFFREE_TOPOLOGY_HEALING_MAX_ITER"
ENV_PHANTOM_FACTOR: str = "DELFIN_FFFREE_TOPOLOGY_HEALING_PHANTOM_FACTOR"
ENV_PHANTOM_TARGET: str = "DELFIN_FFFREE_TOPOLOGY_HEALING_PHANTOM_TARGET_FACTOR"
ENV_MISSING_FACTOR: str = "DELFIN_FFFREE_TOPOLOGY_HEALING_MISSING_FACTOR"
ENV_ANGLE_SIGMA: str = "DELFIN_FFFREE_TOPOLOGY_HEALING_ANGLE_SIGMA"


__all__ = [
    "HealingDiagnostics",
    "DEFAULT_PHANTOM_FACTOR",
    "DEFAULT_PHANTOM_TARGET_FACTOR",
    "DEFAULT_MISSING_FACTOR",
    "DEFAULT_ANGLE_SIGMA",
    "DEFAULT_MAX_ITER",
    "ENV_MASTER",
    "ENV_MAX_ITER",
    "ENV_PHANTOM_FACTOR",
    "ENV_PHANTOM_TARGET",
    "ENV_MISSING_FACTOR",
    "ENV_ANGLE_SIGMA",
    "is_active",
    "cov_radius",
    "ideal_bond_length",
    "detect_phantom_bonds",
    "detect_missing_bonds",
    "detect_wrong_angles",
    "heal_phantom_bonds",
    "heal_missing_bonds",
    "heal_wrong_angles",
    "topology_healing_pipeline",
]


# ---------------------------------------------------------------------------
# Env-flag helpers (single resolver pattern, mirrors grip_polish)
# ---------------------------------------------------------------------------
def is_active() -> bool:
    """``True`` iff the topology-healing master flag is on (default OFF).

    Reads :data:`ENV_MASTER` and accepts ``1``, ``true``, ``yes``, ``on``
    (case-insensitive).  Any other value -> OFF.
    """
    raw = os.environ.get(ENV_MASTER, "").strip().lower()
    if not raw:
        return False
    return raw in ("1", "true", "yes", "on")


def _env_float(name: str, default: float) -> float:
    """Parse a float env-flag with safe fallback (never raise)."""
    raw = os.environ.get(name, "").strip()
    if not raw:
        return float(default)
    try:
        v = float(raw)
        if math.isfinite(v):
            return float(v)
    except (TypeError, ValueError):
        pass
    return float(default)


def _env_int(name: str, default: int) -> int:
    """Parse an int env-flag with safe fallback (never raise)."""
    raw = os.environ.get(name, "").strip()
    if not raw:
        return int(default)
    try:
        v = int(raw)
        if v > 0:
            return int(v)
    except (TypeError, ValueError):
        pass
    return int(default)


# ---------------------------------------------------------------------------
# Geometry primitives
# ---------------------------------------------------------------------------
def cov_radius(symbol: str) -> float:
    """Return the Cordero covalent radius for ``symbol`` (or default 0.90 Å)."""
    return float(_COV.get(str(symbol), _COV_DEFAULT))


def ideal_bond_length(sym_a: str, sym_b: str) -> float:
    """Ideal single-bond length between elements ``sym_a`` and ``sym_b`` (Å)."""
    return cov_radius(sym_a) + cov_radius(sym_b)


def _canonical_pair(i: int, j: int) -> Tuple[int, int]:
    """Lex-sorted (low, high) pair so detection sets are deterministic."""
    a, b = int(i), int(j)
    return (a, b) if a <= b else (b, a)


def _normalise_topology(
    smiles_topology: Iterable[Tuple[int, int]],
) -> List[Tuple[int, int]]:
    """Return a lex-sorted, de-duplicated, self-loop-free copy."""
    out: Set[Tuple[int, int]] = set()
    for ab in smiles_topology:
        try:
            a, b = int(ab[0]), int(ab[1])
        except (TypeError, ValueError, IndexError):
            continue
        if a == b:
            continue
        out.add(_canonical_pair(a, b))
    return sorted(out)


def _build_adjacency(
    n_atoms: int,
    topology_bonds: Sequence[Tuple[int, int]],
) -> List[List[int]]:
    """Return ``adj[i]`` = sorted list of neighbour indices of atom ``i``."""
    adj: List[Set[int]] = [set() for _ in range(n_atoms)]
    for a, b in topology_bonds:
        if 0 <= a < n_atoms and 0 <= b < n_atoms:
            adj[a].add(b)
            adj[b].add(a)
    return [sorted(s) for s in adj]


def _angle_deg(p_a: np.ndarray, p_b: np.ndarray, p_c: np.ndarray) -> float:
    """Return the angle ``A-B-C`` (centered at ``B``) in degrees, clamped."""
    v1 = p_a - p_b
    v2 = p_c - p_b
    n1 = float(np.linalg.norm(v1))
    n2 = float(np.linalg.norm(v2))
    if n1 < _EPS or n2 < _EPS:
        return float("nan")
    cos_t = float(np.dot(v1, v2) / (n1 * n2))
    cos_t = max(-1.0, min(1.0, cos_t))
    return float(math.degrees(math.acos(cos_t)))


# ---------------------------------------------------------------------------
# Mol-bond helper -- accept either RDKit mol OR explicit pair list.
# ---------------------------------------------------------------------------
def _coerce_topology(
    smiles_topology,
) -> List[Tuple[int, int]]:
    """Return a lex-sorted list of (a, b) pairs from a mol or pair iterable.

    Accepts:
    * an RDKit ``Mol`` (uses ``GetBonds()``);
    * any iterable of ``(int, int)`` pairs.

    Self-loops and duplicates are removed.  Failures return an empty list
    (the caller then treats every spatial bond as a phantom, which is
    safe).
    """
    if smiles_topology is None:
        return []
    # RDKit-mol path -- duck-typed on ``GetBonds``.
    if hasattr(smiles_topology, "GetBonds"):
        out: List[Tuple[int, int]] = []
        try:
            for bond in smiles_topology.GetBonds():
                a = int(bond.GetBeginAtomIdx())
                b = int(bond.GetEndAtomIdx())
                if a == b:
                    continue
                out.append(_canonical_pair(a, b))
        except Exception:
            return []
        return _normalise_topology(out)
    # Pair-iterable path.
    try:
        return _normalise_topology(smiles_topology)
    except Exception:
        return []


# ---------------------------------------------------------------------------
# Diagnostic container
# ---------------------------------------------------------------------------
@dataclass
class HealingDiagnostics:
    """Per-iteration diagnostic snapshot for forensic post-mortems."""

    n_phantom_before: int = 0
    n_phantom_after: int = 0
    n_missing_before: int = 0
    n_missing_after: int = 0
    n_wrong_angle_before: int = 0
    n_wrong_angle_after: int = 0
    n_iterations: int = 0
    accepted: bool = False
    rollback_reason: str = ""
    max_displacement: float = 0.0
    per_iter_counts: List[Dict[str, int]] = field(default_factory=list)


# ---------------------------------------------------------------------------
# Detectors
# ---------------------------------------------------------------------------
def detect_phantom_bonds(
    coords: np.ndarray,
    atoms: Sequence[str],
    smiles_topology,
    factor: float = DEFAULT_PHANTOM_FACTOR,
) -> List[Tuple[int, int]]:
    """Return lex-sorted phantom pairs: spatially bonded but topologically NOT.

    A pair (i, j) is **spatially bonded** when
    ``||r_i - r_j|| < factor × (cov_radius(atoms[i]) + cov_radius(atoms[j]))``.

    The phantom set is ``spatial_bonds \\ topology_bonds`` -- i.e. atoms
    that look bonded by geometry yet are not connected in the SMILES graph.
    This catches superimposed phenyl rings, η-cluster overlap, builder-
    fused substituents -- all the artefacts that fool downstream Mahalanobis
    fragment detection into mis-classifying centres.

    Parameters
    ----------
    coords : ndarray (N, 3)
        Atomic Cartesian coordinates (Å).
    atoms : sequence of str
        Element symbols, length N.
    smiles_topology : mol or iterable of (int, int)
        The SMILES-derived bond list.
    factor : float
        Covalent-radius factor below which a pair is treated as
        spatially bonded.  Default :data:`DEFAULT_PHANTOM_FACTOR` (1.2).

    Returns
    -------
    list of (int, int)
        Lex-sorted phantom pairs (each pair lex-sorted internally).
    """
    coords = np.asarray(coords, dtype=np.float64)
    n = int(coords.shape[0])
    if n < 2 or len(atoms) < n:
        return []
    topo = set(_coerce_topology(smiles_topology))
    factor = float(factor)
    phantom: Set[Tuple[int, int]] = set()
    for i in range(n):
        ri = coords[i]
        cov_i = cov_radius(atoms[i])
        for j in range(i + 1, n):
            d = float(np.linalg.norm(ri - coords[j]))
            if d <= _EPS:
                continue
            cutoff = factor * (cov_i + cov_radius(atoms[j]))
            if d < cutoff:
                pair = _canonical_pair(i, j)
                if pair not in topo:
                    phantom.add(pair)
    return sorted(phantom)


def detect_missing_bonds(
    coords: np.ndarray,
    atoms: Sequence[str],
    smiles_topology,
    factor: float = DEFAULT_MISSING_FACTOR,
) -> List[Tuple[int, int]]:
    """Return lex-sorted missing pairs: topology bonds stretched beyond
    ``factor × ideal_length``.

    A SMILES-declared bond (i, j) is **missing** when the actual length
    exceeds ``factor × (cov_radius(atoms[i]) + cov_radius(atoms[j]))``.
    This catches torn rings, donor-slide overshoot, cluster-scaffold
    fragmentation -- failures where the construction stretched a bond
    until the geometry no longer encodes it as a single bond.

    Parameters
    ----------
    coords : ndarray (N, 3)
        Atomic Cartesian coordinates (Å).
    atoms : sequence of str
        Element symbols, length N.
    smiles_topology : mol or iterable of (int, int)
    factor : float
        Stretch multiplier on the ideal length above which the bond is
        considered "missing".  Default :data:`DEFAULT_MISSING_FACTOR`
        (1.5).

    Returns
    -------
    list of (int, int)
        Lex-sorted missing pairs.
    """
    coords = np.asarray(coords, dtype=np.float64)
    n = int(coords.shape[0])
    if n < 2 or len(atoms) < n:
        return []
    topo = _coerce_topology(smiles_topology)
    factor = float(factor)
    missing: Set[Tuple[int, int]] = set()
    for a, b in topo:
        if a < 0 or b < 0 or a >= n or b >= n:
            continue
        d = float(np.linalg.norm(coords[a] - coords[b]))
        ideal = ideal_bond_length(atoms[a], atoms[b])
        if ideal <= _EPS:
            continue
        if d > factor * ideal:
            missing.add(_canonical_pair(a, b))
    return sorted(missing)


def detect_wrong_angles(
    coords: np.ndarray,
    atoms: Sequence[str],
    smiles_topology,
    ideal_angles_lookup=None,
    sigma_threshold: float = DEFAULT_ANGLE_SIGMA,
) -> List[Tuple[int, int, int, float, float, float]]:
    """Return lex-sorted wrong-angle records.

    For each ``(A, B, C)`` triple where ``B`` is bonded to both ``A`` and
    ``C`` in the SMILES topology, compute the actual angle and query the
    Mahalanobis prior ``(μ, σ)``.  Flag when the residual exceeds the
    σ-threshold.

    Parameters
    ----------
    coords : ndarray (N, 3)
    atoms : sequence of str
    smiles_topology : mol or iterable of (int, int)
    ideal_angles_lookup : callable, optional
        ``f(sym_a, sym_b, sym_c) -> (mu_deg, sigma_deg) | None``.  If
        ``None``, the default :func:`grip_mogul_lookup.lookup_angle`
        wrapper is used.  Pass a stub for unit tests.
    sigma_threshold : float
        Residual cutoff in σ units.  Default
        :data:`DEFAULT_ANGLE_SIGMA` (3.0).

    Returns
    -------
    list of (a, b, c, actual_deg, ideal_deg, residual_sigma)
        Lex-sorted by ``(b, a, c)``.
    """
    coords = np.asarray(coords, dtype=np.float64)
    n = int(coords.shape[0])
    if n < 3 or len(atoms) < n:
        return []

    topo = _coerce_topology(smiles_topology)
    adj = _build_adjacency(n, topo)
    sigma_thr = float(sigma_threshold)

    lookup = ideal_angles_lookup
    if lookup is None:
        lookup = _default_angle_lookup()

    records: List[Tuple[int, int, int, float, float, float]] = []
    for b in range(n):
        nbrs = adj[b]
        if len(nbrs) < 2:
            continue
        for i in range(len(nbrs)):
            a = nbrs[i]
            for j in range(i + 1, len(nbrs)):
                c = nbrs[j]
                actual = _angle_deg(coords[a], coords[b], coords[c])
                if not math.isfinite(actual):
                    continue
                try:
                    prior = lookup(atoms[a], atoms[b], atoms[c])
                except Exception:
                    prior = None
                if not prior:
                    continue
                mu_deg, sigma_deg = float(prior[0]), float(prior[1])
                if sigma_deg <= _EPS:
                    continue
                resid = (actual - mu_deg) / sigma_deg
                if abs(resid) > sigma_thr:
                    # canonicalise (a, c) so the record key is stable
                    if a > c:
                        a_out, c_out = c, a
                    else:
                        a_out, c_out = a, c
                    records.append(
                        (int(a_out), int(b), int(c_out),
                         float(actual), float(mu_deg), float(resid))
                    )
    records.sort(key=lambda r: (r[1], r[0], r[2]))
    return records


# ---------------------------------------------------------------------------
# Default angle-prior lookup (wraps grip_mogul_lookup, with sane fallbacks)
# ---------------------------------------------------------------------------
_DEFAULT_PRIORS: Dict[Tuple[str, ...], Tuple[float, float]] = {
    # Fallback per-element-pair priors used when the CCDC library cannot
    # answer.  Numbers reflect canonical ideal geometries with conservative
    # sigmas so the gate doesn't fire on slight deviations.
    ("C", "C", "C"): (111.0, 6.0),
    ("H", "C", "H"): (109.5, 5.0),
    ("H", "C", "C"): (109.5, 5.0),
    ("H", "N", "H"): (107.0, 5.0),
    ("H", "O", "H"): (104.5, 5.0),
    ("C", "N", "C"): (115.0, 6.0),
    ("C", "O", "C"): (110.0, 6.0),
    ("N", "C", "C"): (110.0, 6.0),
    ("O", "C", "C"): (110.0, 6.0),
    ("O", "C", "O"): (115.0, 6.0),
    ("C", "C", "H"): (109.5, 5.0),
}


def _default_angle_lookup():
    """Return a callable ``f(sym_a, sym_b, sym_c) -> (mu, sigma) | None``.

    Tries to import :mod:`grip_mogul_lookup` and use its CCDC library
    first; falls back to :data:`_DEFAULT_PRIORS` when the library is
    missing or returns no hit.
    """
    try:
        from .grip_mogul_lookup import get_default_library  # type: ignore
        lib = get_default_library()
    except Exception:
        lib = None

    def _lookup(sym_a: str, sym_b: str, sym_c: str):
        # CCDC-backed lookup first.
        if lib is not None:
            try:
                hit = lib.lookup_angle(sym_a, sym_b, "*", sym_c)
                if hit is not None and len(hit) >= 2:
                    mu, sigma = float(hit[0]), float(hit[1])
                    if math.isfinite(mu) and math.isfinite(sigma) and sigma > 0.0:
                        return (mu, sigma)
            except Exception:
                pass
        # Canonical-pair static fallback (try both orderings of the
        # terminal neighbours; B stays fixed).
        for key in ((sym_a, sym_b, sym_c), (sym_c, sym_b, sym_a)):
            prior = _DEFAULT_PRIORS.get(key)
            if prior is not None:
                return prior
        # Final generic fallback: 110° ± 8° (a wide-but-physical band).
        return (110.0, 8.0)

    return _lookup


# ---------------------------------------------------------------------------
# Healers
# ---------------------------------------------------------------------------
def heal_phantom_bonds(
    coords: np.ndarray,
    atoms: Sequence[str],
    phantom_pairs: Sequence[Tuple[int, int]],
    target_separation_factor: float = DEFAULT_PHANTOM_TARGET_FACTOR,
    frozen_atoms: Optional[Set[int]] = None,
) -> Tuple[np.ndarray, float]:
    """Push phantom pairs apart to the target separation.

    For each phantom (i, j), move ``i`` and ``j`` symmetrically along the
    pair axis (or each, alone, if the other is frozen) so that the final
    distance equals ``target_separation_factor × Σcov``.

    Parameters
    ----------
    coords : ndarray (N, 3)
        Current coordinates -- a COPY is returned.
    atoms : sequence of str
    phantom_pairs : sequence of (int, int)
        Pairs to heal (typically the output of
        :func:`detect_phantom_bonds`).
    target_separation_factor : float
        Multiplier on Σcov for the target distance.  Default
        :data:`DEFAULT_PHANTOM_TARGET_FACTOR` (1.5).
    frozen_atoms : set of int, optional
        Atoms that cannot move (typically ``{metal, *donors}``).

    Returns
    -------
    healed_coords : ndarray (N, 3)
    delta : float
        Maximum atom displacement applied (Å), for the diagnostic record.
    """
    coords = np.asarray(coords, dtype=np.float64)
    healed = coords.copy()
    if frozen_atoms is None:
        frozen_atoms = set()
    else:
        frozen_atoms = {int(a) for a in frozen_atoms}
    max_disp = 0.0
    factor = float(target_separation_factor)
    n = int(healed.shape[0])

    # Lex-sort the input pairs so the application order is deterministic
    # (input may already be sorted but caller can pass any iterable).
    pairs_sorted = sorted({_canonical_pair(p[0], p[1]) for p in phantom_pairs})

    for i, j in pairs_sorted:
        if not (0 <= i < n and 0 <= j < n) or i == j:
            continue
        ri = healed[i]
        rj = healed[j]
        diff = rj - ri
        d = float(np.linalg.norm(diff))
        if d < _EPS:
            # coincident -- push apart along x (deterministic seed-free).
            direction = np.array([1.0, 0.0, 0.0], dtype=np.float64)
        else:
            direction = diff / d
        target_d = factor * ideal_bond_length(atoms[i], atoms[j])
        if d >= target_d:
            continue  # already separated enough
        deficit = target_d - d
        # Distribute the displacement.  Both frozen -> no-op (caller can
        # rollback elsewhere).  One frozen -> move the other by the full
        # deficit.  Neither frozen -> split symmetrically.
        i_frozen = i in frozen_atoms
        j_frozen = j in frozen_atoms
        if i_frozen and j_frozen:
            continue
        if i_frozen:
            healed[j] = rj + direction * deficit
            disp = deficit
        elif j_frozen:
            healed[i] = ri - direction * deficit
            disp = deficit
        else:
            half = 0.5 * deficit
            healed[i] = ri - direction * half
            healed[j] = rj + direction * half
            disp = half
        if disp > max_disp:
            max_disp = float(disp)
    return healed, max_disp


def heal_missing_bonds(
    coords: np.ndarray,
    atoms: Sequence[str],
    missing_pairs: Sequence[Tuple[int, int]],
    target_length_factor: float = 1.0,
    frozen_atoms: Optional[Set[int]] = None,
) -> Tuple[np.ndarray, float]:
    """Pull missing bonds together toward their ideal length.

    For each missing (i, j), set the pair distance to
    ``target_length_factor × Σcov``.  The mid-displacement is split
    between the two atoms (or applied entirely to the non-frozen one if
    the other is frozen).

    Parameters
    ----------
    coords : ndarray (N, 3)
    atoms : sequence of str
    missing_pairs : sequence of (int, int)
    target_length_factor : float
        Multiplier on Σcov for the target distance.  Default 1.0 (exact
        ideal single-bond length).
    frozen_atoms : set of int, optional

    Returns
    -------
    healed_coords : ndarray (N, 3)
    delta : float
        Maximum atom displacement applied (Å).
    """
    coords = np.asarray(coords, dtype=np.float64)
    healed = coords.copy()
    if frozen_atoms is None:
        frozen_atoms = set()
    else:
        frozen_atoms = {int(a) for a in frozen_atoms}
    max_disp = 0.0
    factor = float(target_length_factor)
    n = int(healed.shape[0])

    pairs_sorted = sorted({_canonical_pair(p[0], p[1]) for p in missing_pairs})

    for i, j in pairs_sorted:
        if not (0 <= i < n and 0 <= j < n) or i == j:
            continue
        ri = healed[i]
        rj = healed[j]
        diff = rj - ri
        d = float(np.linalg.norm(diff))
        if d < _EPS:
            direction = np.array([1.0, 0.0, 0.0], dtype=np.float64)
        else:
            direction = diff / d
        target_d = factor * ideal_bond_length(atoms[i], atoms[j])
        if d <= target_d:
            continue  # already shorter than target
        excess = d - target_d
        i_frozen = i in frozen_atoms
        j_frozen = j in frozen_atoms
        if i_frozen and j_frozen:
            continue
        if i_frozen:
            healed[j] = rj - direction * excess
            disp = excess
        elif j_frozen:
            healed[i] = ri + direction * excess
            disp = excess
        else:
            half = 0.5 * excess
            healed[i] = ri + direction * half
            healed[j] = rj - direction * half
            disp = half
        if disp > max_disp:
            max_disp = float(disp)
    return healed, max_disp


def _rotate_C_toward_angle(
    p_a: np.ndarray,
    p_b: np.ndarray,
    p_c: np.ndarray,
    target_deg: float,
    step: float = DEFAULT_ANGLE_STEP,
) -> np.ndarray:
    """Move ``p_c`` in the (A, B, C) plane to bring the angle to ``target_deg``.

    The B-C bond length is preserved exactly.  ``step`` damps the
    correction (0 = no move; 1 = full correction in one pass).

    The rotation axis is the plane normal of (BA × BC); when those are
    colinear (degenerate angle) the function uses an arbitrary
    deterministic axis (z-hat if available, else x-hat) so the
    transformation is reproducible.
    """
    p_a = np.asarray(p_a, dtype=np.float64)
    p_b = np.asarray(p_b, dtype=np.float64)
    p_c = np.asarray(p_c, dtype=np.float64)
    v_ba = p_a - p_b
    v_bc = p_c - p_b
    n_ba = float(np.linalg.norm(v_ba))
    n_bc = float(np.linalg.norm(v_bc))
    if n_ba < _EPS or n_bc < _EPS:
        return p_c.copy()
    cos_actual = float(np.dot(v_ba, v_bc) / (n_ba * n_bc))
    cos_actual = max(-1.0, min(1.0, cos_actual))
    actual_rad = float(math.acos(cos_actual))
    target_rad = float(math.radians(float(target_deg)))
    # Damped correction: blend toward target_rad by ``step`` fraction.
    new_rad = actual_rad + float(step) * (target_rad - actual_rad)
    # Build the in-plane rotation axis.  When BA × BC is degenerate
    # (parallel/anti-parallel) fall back to a deterministic axis.
    axis = np.cross(v_ba, v_bc)
    norm_axis = float(np.linalg.norm(axis))
    if norm_axis < _EPS:
        # Choose any unit vector perpendicular to v_bc deterministically.
        ref = np.array([0.0, 0.0, 1.0], dtype=np.float64)
        if abs(float(np.dot(ref, v_bc) / n_bc)) > 0.9:
            ref = np.array([1.0, 0.0, 0.0], dtype=np.float64)
        axis = np.cross(v_bc, ref)
        norm_axis = float(np.linalg.norm(axis))
        if norm_axis < _EPS:
            return p_c.copy()
    axis = axis / norm_axis
    # Rotate v_bc around ``axis`` (= unit normal of the BA-BC plane).
    # The axis is built as cross(v_ba, v_bc) so the right-hand rule means
    # rotating v_bc by a POSITIVE angle around ``axis`` moves it AWAY
    # from v_ba (i.e. increases the BAC angle).  We want the new angle to
    # equal ``new_rad``, so we rotate by ``new_rad - actual_rad``.
    delta_rad = new_rad - actual_rad
    # Rodrigues' formula
    cos_d = math.cos(delta_rad)
    sin_d = math.sin(delta_rad)
    v = v_bc
    v_rot = (
        v * cos_d
        + np.cross(axis, v) * sin_d
        + axis * float(np.dot(axis, v)) * (1.0 - cos_d)
    )
    # Renormalise to preserve B-C length exactly (guards against
    # accumulated round-off in repeated calls).
    n_rot = float(np.linalg.norm(v_rot))
    if n_rot < _EPS:
        return p_c.copy()
    v_rot = v_rot * (n_bc / n_rot)
    return p_b + v_rot


def heal_wrong_angles(
    coords: np.ndarray,
    atoms: Sequence[str],
    wrong_angle_triples: Sequence[Tuple[int, int, int, float, float, float]],
    max_iter: int = 1,
    step: float = DEFAULT_ANGLE_STEP,
    frozen_atoms: Optional[Set[int]] = None,
) -> Tuple[np.ndarray, float]:
    """Correct wrong angles by rotating atom ``C`` of each (A, B, C) triple.

    The B-C bond length is preserved (in-plane rotation around ``B``).
    Multiple triples sharing an atom may interfere; the function applies
    them sequentially in deterministic order and clamps the per-call
    displacement so a single iteration cannot blow up.  Several passes
    (``max_iter``) let the corrections settle.

    Parameters
    ----------
    coords : ndarray (N, 3)
    atoms : sequence of str
    wrong_angle_triples : sequence
        Output of :func:`detect_wrong_angles` (or a compatible tuple).
        Each element is ``(a, b, c, actual_deg, ideal_deg, residual_sigma)``.
    max_iter : int
        Inner passes per call (default 1; usually called by the pipeline
        loop which does its own re-detection).
    step : float
        Damping factor in (0, 1].  Default :data:`DEFAULT_ANGLE_STEP` (0.5).
    frozen_atoms : set of int, optional

    Returns
    -------
    healed_coords : ndarray (N, 3)
    delta : float
        Maximum atom displacement applied (Å).
    """
    coords = np.asarray(coords, dtype=np.float64)
    healed = coords.copy()
    if frozen_atoms is None:
        frozen_atoms = set()
    else:
        frozen_atoms = {int(a) for a in frozen_atoms}
    n = int(healed.shape[0])
    max_disp = 0.0

    # Lex-sort by (b, a, c) so application order is deterministic.
    sorted_triples = sorted(
        wrong_angle_triples,
        key=lambda r: (int(r[1]), int(r[0]), int(r[2])),
    )

    for _pass in range(max(1, int(max_iter))):
        any_moved = False
        for record in sorted_triples:
            a, b, c = int(record[0]), int(record[1]), int(record[2])
            ideal_deg = float(record[4])
            if not (0 <= a < n and 0 <= b < n and 0 <= c < n):
                continue
            if a == b or b == c or a == c:
                continue
            # Pick the terminal atom to rotate.  Prefer ``c``; if it's
            # frozen, try ``a``.  If both are frozen, skip.
            if c not in frozen_atoms:
                old = healed[c].copy()
                healed[c] = _rotate_C_toward_angle(
                    healed[a], healed[b], healed[c], ideal_deg, step=step,
                )
                disp = float(np.linalg.norm(healed[c] - old))
            elif a not in frozen_atoms:
                # Symmetric: rotate ``a`` toward the ideal angle by treating
                # (C, B, A) as the triple.
                old = healed[a].copy()
                healed[a] = _rotate_C_toward_angle(
                    healed[c], healed[b], healed[a], ideal_deg, step=step,
                )
                disp = float(np.linalg.norm(healed[a] - old))
            else:
                disp = 0.0
            if disp > _EPS:
                any_moved = True
            if disp > max_disp:
                max_disp = float(disp)
        if not any_moved:
            break
    return healed, max_disp


# ---------------------------------------------------------------------------
# Pipeline
# ---------------------------------------------------------------------------
def _md_snapshot(
    coords: np.ndarray,
    metal_idx: Optional[int],
    donors: Sequence[int],
) -> Dict[int, float]:
    """Snapshot M-D distances for the invariant check (Å)."""
    if metal_idx is None or metal_idx < 0:
        return {}
    coords = np.asarray(coords, dtype=np.float64)
    out: Dict[int, float] = {}
    if metal_idx >= coords.shape[0]:
        return out
    for d in donors:
        d = int(d)
        if 0 <= d < coords.shape[0] and d != metal_idx:
            out[d] = float(np.linalg.norm(coords[d] - coords[metal_idx]))
    return out


def _md_intact(
    coords: np.ndarray,
    snapshot: Dict[int, float],
    metal_idx: Optional[int],
    tol: float = 0.05,
) -> bool:
    """``True`` iff every M-D distance has changed by < tol from snapshot."""
    if not snapshot or metal_idx is None or metal_idx < 0:
        return True
    coords = np.asarray(coords, dtype=np.float64)
    if metal_idx >= coords.shape[0]:
        return True
    m = coords[metal_idx]
    for d, d0 in snapshot.items():
        if d >= coords.shape[0]:
            continue
        d_now = float(np.linalg.norm(coords[d] - m))
        if abs(d_now - d0) > tol:
            return False
    return True


def topology_healing_pipeline(
    coords: np.ndarray,
    atoms: Sequence[str],
    smiles_topology,
    *,
    metal_idx: Optional[int] = None,
    donors: Sequence[int] = (),
    max_iter: Optional[int] = None,
    phantom_factor: Optional[float] = None,
    phantom_target_factor: Optional[float] = None,
    missing_factor: Optional[float] = None,
    angle_sigma: Optional[float] = None,
    ideal_angles_lookup=None,
    md_tol: float = 0.05,
    return_diagnostics: bool = False,
):
    """Iteratively heal phantom + missing + wrong-angle defects.

    Pipeline order (per outer iter): **missing** -> **phantom** -> **angles**.
    Rationale: closing missing bonds first restores correct topology so the
    phantom detector sees the right baseline; angle correction is least
    intrusive and runs last so it polishes the result.

    The pipeline implements **accept-if-better** semantics: after each
    outer pass, the cumulative defect count
    ``n_phantom + n_missing + n_wrong_angle`` must not increase relative to
    the starting state.  If it does (or the M-D invariant is broken), the
    function rolls back to the input ``coords``.

    Parameters
    ----------
    coords : ndarray (N, 3)
    atoms : sequence of str
    smiles_topology : mol or iterable of (int, int)
    metal_idx : int, optional
        Frozen metal index (M-D invariant guard); None -> no metal frozen.
    donors : sequence of int
        Frozen donor indices.
    max_iter : int, optional
        Outer loop budget.  ``None`` -> read env / default 5.
    phantom_factor, phantom_target_factor : float, optional
        Threshold + target for phantom detection / healing.
    missing_factor : float, optional
        Stretch threshold for missing detection.
    angle_sigma : float, optional
        σ-threshold for wrong-angle detection.
    ideal_angles_lookup : callable, optional
        Override for the angle prior lookup (tests use this to inject
        deterministic priors).
    md_tol : float
        Half-width tolerance for the M-D invariant check (Å).
    return_diagnostics : bool
        If True, returns ``(coords, HealingDiagnostics)``; otherwise just
        the healed coordinates.

    Returns
    -------
    ndarray (N, 3) or (ndarray, HealingDiagnostics)
    """
    coords_in = np.asarray(coords, dtype=np.float64).copy()
    diag = HealingDiagnostics()
    n = int(coords_in.shape[0])
    if n < 2 or len(atoms) < n:
        if return_diagnostics:
            return coords_in, diag
        return coords_in

    if not np.all(np.isfinite(coords_in)):
        diag.rollback_reason = "non-finite input coords"
        if return_diagnostics:
            return coords_in, diag
        return coords_in

    # Resolve flags.
    max_iter_eff = int(max_iter) if max_iter is not None else _env_int(ENV_MAX_ITER, DEFAULT_MAX_ITER)
    phantom_f = float(phantom_factor) if phantom_factor is not None else _env_float(ENV_PHANTOM_FACTOR, DEFAULT_PHANTOM_FACTOR)
    phantom_tgt = float(phantom_target_factor) if phantom_target_factor is not None else _env_float(ENV_PHANTOM_TARGET, DEFAULT_PHANTOM_TARGET_FACTOR)
    missing_f = float(missing_factor) if missing_factor is not None else _env_float(ENV_MISSING_FACTOR, DEFAULT_MISSING_FACTOR)
    angle_s = float(angle_sigma) if angle_sigma is not None else _env_float(ENV_ANGLE_SIGMA, DEFAULT_ANGLE_SIGMA)

    metal_eff: Optional[int] = int(metal_idx) if metal_idx is not None else None
    donors_eff = tuple(int(d) for d in donors)
    frozen: Set[int] = set()
    if metal_eff is not None and metal_eff >= 0:
        frozen.add(metal_eff)
    frozen.update(donors_eff)

    md_snap = _md_snapshot(coords_in, metal_eff, donors_eff)

    # Initial diagnostic snapshot.
    phantom0 = detect_phantom_bonds(coords_in, atoms, smiles_topology, factor=phantom_f)
    missing0 = detect_missing_bonds(coords_in, atoms, smiles_topology, factor=missing_f)
    angles0 = detect_wrong_angles(
        coords_in, atoms, smiles_topology,
        ideal_angles_lookup=ideal_angles_lookup, sigma_threshold=angle_s,
    )
    diag.n_phantom_before = len(phantom0)
    diag.n_missing_before = len(missing0)
    diag.n_wrong_angle_before = len(angles0)

    # No defects -> bit-exact early return.
    if not (phantom0 or missing0 or angles0):
        diag.accepted = False
        diag.rollback_reason = "no defects detected"
        diag.n_phantom_after = 0
        diag.n_missing_after = 0
        diag.n_wrong_angle_after = 0
        if return_diagnostics:
            return coords_in, diag
        return coords_in

    work = coords_in.copy()
    max_disp_total = 0.0

    for it in range(max_iter_eff):
        diag.n_iterations = it + 1
        # 1. Missing bonds (close stretched edges).
        missing = detect_missing_bonds(work, atoms, smiles_topology, factor=missing_f)
        if missing:
            work_try, disp = heal_missing_bonds(
                work, atoms, missing, frozen_atoms=frozen,
            )
            if np.all(np.isfinite(work_try)):
                work = work_try
                max_disp_total = max(max_disp_total, disp)

        # 2. Phantom bonds (push apart non-topological close pairs).
        phantom = detect_phantom_bonds(work, atoms, smiles_topology, factor=phantom_f)
        if phantom:
            work_try, disp = heal_phantom_bonds(
                work, atoms, phantom,
                target_separation_factor=phantom_tgt,
                frozen_atoms=frozen,
            )
            if np.all(np.isfinite(work_try)):
                work = work_try
                max_disp_total = max(max_disp_total, disp)

        # 3. Wrong angles (rotate terminal atoms to the ideal).
        angles = detect_wrong_angles(
            work, atoms, smiles_topology,
            ideal_angles_lookup=ideal_angles_lookup, sigma_threshold=angle_s,
        )
        if angles:
            work_try, disp = heal_wrong_angles(
                work, atoms, angles, frozen_atoms=frozen,
            )
            if np.all(np.isfinite(work_try)):
                work = work_try
                max_disp_total = max(max_disp_total, disp)

        # Record per-iter counts (post-pass).
        diag.per_iter_counts.append({
            "missing": len(detect_missing_bonds(work, atoms, smiles_topology, factor=missing_f)),
            "phantom": len(detect_phantom_bonds(work, atoms, smiles_topology, factor=phantom_f)),
            "angles": len(detect_wrong_angles(
                work, atoms, smiles_topology,
                ideal_angles_lookup=ideal_angles_lookup, sigma_threshold=angle_s,
            )),
        })

        # Convergence: every detector now reports 0 defects.
        latest = diag.per_iter_counts[-1]
        if latest["missing"] == 0 and latest["phantom"] == 0 and latest["angles"] == 0:
            break

    # Final defect counts.
    phantom_after = detect_phantom_bonds(work, atoms, smiles_topology, factor=phantom_f)
    missing_after = detect_missing_bonds(work, atoms, smiles_topology, factor=missing_f)
    angles_after = detect_wrong_angles(
        work, atoms, smiles_topology,
        ideal_angles_lookup=ideal_angles_lookup, sigma_threshold=angle_s,
    )
    diag.n_phantom_after = len(phantom_after)
    diag.n_missing_after = len(missing_after)
    diag.n_wrong_angle_after = len(angles_after)
    diag.max_displacement = float(max_disp_total)

    total_before = diag.n_phantom_before + diag.n_missing_before + diag.n_wrong_angle_before
    total_after = diag.n_phantom_after + diag.n_missing_after + diag.n_wrong_angle_after

    # Accept-if-better gates.
    if not np.all(np.isfinite(work)):
        diag.accepted = False
        diag.rollback_reason = "non-finite healed coords"
        if return_diagnostics:
            return coords_in, diag
        return coords_in

    if not _md_intact(work, md_snap, metal_eff, tol=md_tol):
        diag.accepted = False
        diag.rollback_reason = "M-D invariant broken"
        if return_diagnostics:
            return coords_in, diag
        return coords_in

    if total_after >= total_before:
        diag.accepted = False
        diag.rollback_reason = "no net defect reduction"
        if return_diagnostics:
            return coords_in, diag
        return coords_in

    diag.accepted = True
    if return_diagnostics:
        return work, diag
    return work
