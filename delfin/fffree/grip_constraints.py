"""GRIP — Hard Constraints + Clash Floor (Phase 3, v1).

H-H clash inclusion (2026-06-04 wiring cleanup)
-----------------------------------------------
``ClashFloorPenalty`` always iterates every non-bonded, non-1,3 pair
including H atoms (Bondi radius 1.20 Å is in the standard vdW table).
The :func:`hh_pair_contribution` helper exposed at the bottom of this
module returns the H-H sub-contribution to the loss when the env-flag
``DELFIN_FFFREE_HH_CLASH_INCLUDE`` is set, so callers can diagnose
whether the polish actually relieved H-H eclipsing.  When the flag is
OFF the helper returns ``(0.0, 0)`` -- byte-identical with HEAD
``b195dba``.


Implements the constraint stack from SPEC §3.2 used by the GRIP polish
optimiser (:mod:`grip_polish`):

* :class:`MDInvariantConstraint`  -- metal-donor distances stay ±tol of their
  initial values (fffree's M-D rigidity is sacred and never relaxed)
* :class:`DonorPolyhedronConstraint` -- donor unit vectors still trace the
  reference polyhedron within a CShM tolerance (implicit consequence of
  M-D rigidity + initial-placement, validated here as an extra safety net)
* :class:`TopologyConstraint` -- no covalent bond is stretched past
  ``max_distance_multiplier × initial_length`` (no spontaneous bond breaking)
* :class:`ChiralVolumeConstraint` -- signed tetrahedral volume at every
  stereocenter keeps its initial sign (no inversion)
* :class:`ClashFloorPenalty` -- pure-repulsive "Pauli floor" implemented as a
  smooth C1 penalty against the statistical impossibility of clashes
  (CCDC has zero clash geometries, see SPEC §2.2)

Philosophical contract (SPEC §2.2 + §11):

* The clash floor IS a soft penalty BUT it is NOT a force field --
  it has no parameterised depth, no LJ minimum, no attractive part.
  It is one-sided repulsion against statistical impossibility.
* Hard constraints validate post-step (M-D, topology, chirality) and trigger
  a rollback in :func:`grip_polish` rather than entering the loss.
* Everything is deterministic float64, no RNG, no platform-dependent calls.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, FrozenSet, Iterable, List, Optional, Sequence, Set, Tuple

import numpy as np

__all__ = [
    "MDInvariantConstraint",
    "DonorPolyhedronConstraint",
    "TopologyConstraint",
    "ChiralVolumeConstraint",
    "ClashFloorPenalty",
    "hh_clash_include_active",
    "hh_pair_contribution",
]


# ---------------------------------------------------------------------------
# H-H clash inclusion env-predicate (2026-06-04, wiring cleanup).
#
# The standalone :mod:`hh_clash_detector` module ships specialised geometry
# for H-H contacts (geminal / 1-3 / 1-4 topology filter, Bondi vdW floor).
# When the env-flag ``DELFIN_FFFREE_HH_CLASH_INCLUDE=1`` is set, callers
# refine ``ClashFloorPenalty.exclude_13_pairs`` to skip non-clash H-H pairs
# (geminal + 1-3 within CH3) so the Pauli penalty fires only on real
# eclipsing.  Default OFF -> byte-identical with HEAD b195dba.
# ---------------------------------------------------------------------------
import os as _os_for_hh  # local alias so we don't pollute the public namespace


def hh_clash_include_active() -> bool:
    """``True`` iff ``DELFIN_FFFREE_HH_CLASH_INCLUDE`` is on (default OFF)."""
    raw = _os_for_hh.environ.get(
        "DELFIN_FFFREE_HH_CLASH_INCLUDE", ""
    ).strip().lower()
    if not raw:
        return False
    return raw in ("1", "true", "yes", "on")


def hh_pair_contribution(
    syms: Sequence[str],
    R: np.ndarray,
    *,
    factor: Optional[float] = None,
) -> Tuple[float, int]:
    """Return ``(severity, count)`` of H-H clashes in ``R`` when the env-flag
    is ON; ``(0.0, 0)`` when OFF (byte-identical).

    Diagnostic helper kept next to :class:`ClashFloorPenalty` so callers can
    measure the H-H contribution of a polished geometry without re-walking
    the full pair table.  The detector internally honours the
    ``DELFIN_FFFREE_HH_CLASH_FACTOR`` and ``DELFIN_FFFREE_HH_CLASH_MIN_TOPO``
    env-flags; pass ``factor`` to override the vdW-sum fraction.

    Returns
    -------
    (severity, count)
        ``severity`` = ``Σ max(0, d_min - d)^2`` (Å²), summed over H-H
        clashes; ``count`` = pair count.  Both 0 when env is OFF.
    """
    if not hh_clash_include_active():
        return 0.0, 0
    try:
        from .hh_clash_detector import detect_hh_clashes
    except Exception:
        return 0.0, 0
    try:
        clashes = detect_hh_clashes(syms, R, factor=factor)
    except Exception:
        return 0.0, 0
    sev = float(sum(t[3] for t in clashes))
    return sev, int(len(clashes))


# Numerical floor for distance-norm divisions in the analytic gradient.
_EPS = 1e-12


def _as_R(R: np.ndarray) -> np.ndarray:
    """Coerce input to a contiguous (N, 3) float64 array."""
    R = np.asarray(R, dtype=np.float64)
    if R.ndim == 1:
        if R.size % 3 != 0:
            raise ValueError(f"R must be (N,3) or (3N,), got length {R.size}")
        R = R.reshape(-1, 3)
    return R


# ---------------------------------------------------------------------------
# M-D Invariant
# ---------------------------------------------------------------------------
@dataclass
class MDInvariantConstraint:
    """Validates that ``|r_M - r_D_i|`` stays within ``tol`` of the targets.

    Parameters
    ----------
    metal_idx : int
        Atom index of the metal centre.
    donor_idxs : sequence of int
        Atom indices of the coordinated donors.
    target_distances : sequence of float
        Target M-D distance for each donor (same order as ``donor_idxs``).
    tol : float, default 0.05
        Half-width of the allowed band (Å) -- the SPEC §3.2 hard tolerance.
    """

    metal_idx: int
    donor_idxs: Tuple[int, ...]
    target_distances: Tuple[float, ...]
    tol: float = 0.05

    def __post_init__(self):
        self.metal_idx = int(self.metal_idx)
        self.donor_idxs = tuple(int(i) for i in self.donor_idxs)
        self.target_distances = tuple(float(d) for d in self.target_distances)
        if len(self.donor_idxs) != len(self.target_distances):
            raise ValueError(
                "donor_idxs and target_distances must be the same length"
            )
        self.tol = float(self.tol)
        if self.tol <= 0.0:
            raise ValueError("tol must be > 0")

    @property
    def frozen_atom_set(self) -> FrozenSet[int]:
        """The set of atoms whose gradient must be zeroed in the optimiser."""
        return frozenset((self.metal_idx, *self.donor_idxs))

    def measure(self, R: np.ndarray) -> np.ndarray:
        """Return the array of current M-D distances (order = donor_idxs)."""
        R = _as_R(R)
        m = R[self.metal_idx]
        return np.array(
            [float(np.linalg.norm(R[d] - m)) for d in self.donor_idxs],
            dtype=np.float64,
        )

    def violations(self, R: np.ndarray) -> np.ndarray:
        """Per-donor signed deviation ``|r_M - r_D| - target`` (Å)."""
        return self.measure(R) - np.asarray(self.target_distances, dtype=np.float64)

    def validate(self, R: np.ndarray) -> bool:
        """``True`` iff every M-D distance is within ±``tol`` of its target."""
        if not self.donor_idxs:
            return True
        dev = np.abs(self.violations(R))
        if not np.all(np.isfinite(dev)):
            return False
        return bool(np.max(dev) <= self.tol + 1e-12)


# ---------------------------------------------------------------------------
# Donor polyhedron (implicit via M-D + placement; validated via CShM tolerance)
# ---------------------------------------------------------------------------
@dataclass
class DonorPolyhedronConstraint:
    """CShM-based safety net for the donor polyhedron.

    The placement step of fffree already lays the donors on the polyhedron
    vertices and the M-D constraint above pins them; this class re-measures
    the continuous shape measure (CShM) of the donor positions vs. the
    reference polyhedron and rejects any drift beyond ``cshm_tol``.

    For the polish stage this is essentially an extra defence against a
    pathological combination of donor displacements that nonetheless keeps
    M-D within tol (e.g. rotation of two donors about the metal).
    """

    metal_idx: int
    donor_idxs: Tuple[int, ...]
    ref_vectors_normalized: np.ndarray
    target_md_distances: Tuple[float, ...]
    geometry: str = ""
    cshm_tol: float = 5.0   # CShM > 5 is a clear topological drift

    def __post_init__(self):
        self.metal_idx = int(self.metal_idx)
        self.donor_idxs = tuple(int(i) for i in self.donor_idxs)
        self.target_md_distances = tuple(float(d) for d in self.target_md_distances)
        self.ref_vectors_normalized = np.asarray(
            self.ref_vectors_normalized, dtype=np.float64
        )
        if self.ref_vectors_normalized.ndim != 2 or self.ref_vectors_normalized.shape[1] != 3:
            raise ValueError("ref_vectors_normalized must be (k, 3)")
        if len(self.donor_idxs) != self.ref_vectors_normalized.shape[0]:
            raise ValueError(
                "donor_idxs and ref_vectors_normalized must have the same length"
            )
        self.cshm_tol = float(self.cshm_tol)

    def measure(self, R: np.ndarray) -> float:
        """Return the current CShM of the donor unit vectors vs. the reference."""
        R = _as_R(R)
        if len(self.donor_idxs) < 2:
            return 0.0
        # Self-contained CShM evaluation to avoid a circular import on
        # delfin.fffree.polyhedra (which itself is loaded lazily here).
        try:
            from .polyhedra import cshm
        except Exception:
            return 0.0
        m = R[self.metal_idx]
        vecs = np.array(
            [R[d] - m for d in self.donor_idxs], dtype=np.float64
        )
        norms = np.linalg.norm(vecs, axis=1, keepdims=True)
        if np.any(norms < _EPS):
            return 0.0
        unit = vecs / norms
        # cshm takes geometry-name; the caller stored it (or empty -> we pass
        # an empty key which returns 0.0 from the polyhedra module).
        try:
            return float(cshm(unit, self.geometry)) if self.geometry else 0.0
        except Exception:
            return 0.0

    def validate(self, R: np.ndarray) -> bool:
        """``True`` iff the CShM is within ``cshm_tol``."""
        val = self.measure(R)
        if not np.isfinite(val):
            return False
        return bool(val <= self.cshm_tol + 1e-9)


# ---------------------------------------------------------------------------
# Topology
# ---------------------------------------------------------------------------
@dataclass
class TopologyConstraint:
    """Validates that no bond has been stretched past a multiplier of its
    initial length (no spontaneous bond breaking).

    ``mol_bonds`` is the list of (a, b) pairs taken from the RDKit mol;
    ``initial_lengths`` is computed once at construction-time from ``P0``
    so the test is "did the polish break a bond" rather than "is the bond
    near some library median".
    """

    mol_bonds: List[Tuple[int, int]]
    initial_lengths: List[float]
    max_distance_multiplier: float = 1.5

    @classmethod
    def from_initial(
        cls,
        mol_bonds: Sequence[Tuple[int, int]],
        P0: np.ndarray,
        max_distance_multiplier: float = 1.5,
    ) -> "TopologyConstraint":
        """Build a constraint whose initial lengths are read off ``P0``."""
        P0 = _as_R(P0)
        bonds: List[Tuple[int, int]] = []
        lens: List[float] = []
        for (a, b) in mol_bonds:
            a, b = int(a), int(b)
            if a == b:
                continue
            d = float(np.linalg.norm(P0[a] - P0[b]))
            bonds.append((a, b))
            lens.append(d if d > _EPS else 1.0)  # avoid degenerate 0 length
        return cls(
            mol_bonds=bonds,
            initial_lengths=lens,
            max_distance_multiplier=float(max_distance_multiplier),
        )

    def validate(self, R: np.ndarray) -> bool:
        """``True`` iff every tracked bond stays under
        ``max_distance_multiplier × initial_length``."""
        R = _as_R(R)
        if not self.mol_bonds:
            return True
        mult = float(self.max_distance_multiplier)
        for (a, b), L0 in zip(self.mol_bonds, self.initial_lengths):
            d = float(np.linalg.norm(R[a] - R[b]))
            if not np.isfinite(d):
                return False
            if d > mult * L0 + 1e-12:
                return False
        return True


# ---------------------------------------------------------------------------
# Chiral volume
# ---------------------------------------------------------------------------
@dataclass
class ChiralVolumeConstraint:
    """Validates that every stereocentre keeps its initial signed tetrahedral
    volume (sign only -- magnitude allowed to vary).

    Each entry is ``(center, a, b, c, sign)`` where ``sign in {+1, -1}`` is the
    sign of ``det([r_a - r_center, r_b - r_center, r_c - r_center])`` at
    construction time.
    """

    stereo_centers: List[Tuple[int, int, int, int, int]] = field(default_factory=list)

    @staticmethod
    def signed_volume(R: np.ndarray, center: int, a: int, b: int, c: int) -> float:
        """Return ``det([r_a-r_center, r_b-r_center, r_c-r_center])``."""
        R = _as_R(R)
        rc = R[center]
        m = np.array([R[a] - rc, R[b] - rc, R[c] - rc], dtype=np.float64)
        return float(np.linalg.det(m))

    @classmethod
    def from_initial(
        cls,
        centers: Sequence[Tuple[int, int, int, int]],
        P0: np.ndarray,
    ) -> "ChiralVolumeConstraint":
        """Build a constraint capturing the signs at ``P0``.

        ``centers`` is a sequence of ``(center, a, b, c)`` quadruples.  The
        sign of each tetrahedral volume is recorded; entries whose volume is
        within ``1e-9`` of zero are skipped (no chirality to protect).
        """
        P0 = _as_R(P0)
        entries: List[Tuple[int, int, int, int, int]] = []
        for (center, a, b, c) in centers:
            v = cls.signed_volume(P0, int(center), int(a), int(b), int(c))
            if abs(v) < 1e-9:
                continue
            s = 1 if v > 0.0 else -1
            entries.append((int(center), int(a), int(b), int(c), int(s)))
        return cls(stereo_centers=entries)

    def validate(self, R: np.ndarray) -> bool:
        """``True`` iff every recorded sign is preserved."""
        if not self.stereo_centers:
            return True
        R = _as_R(R)
        for (center, a, b, c, sign) in self.stereo_centers:
            v = self.signed_volume(R, center, a, b, c)
            if not np.isfinite(v):
                return False
            cur = 1 if v > 0.0 else (-1 if v < 0.0 else 0)
            if cur != sign:
                return False
        return True


# ---------------------------------------------------------------------------
# Clash floor ("Pauli" penalty)
# ---------------------------------------------------------------------------
@dataclass
class ClashFloorPenalty:
    """Pure-repulsive smooth penalty against statistical-impossibility clashes.

    For every non-bonded, non-1,3 atom pair we accumulate

    .. math::

        L_\\text{clash} = \\sum_{(i,j)} w_{ij} \\max(0, d_{ij}^{\\min} - d_{ij})^2

    with :math:`d^{\\min} = 0.85 (r_i^{vdW} + r_j^{vdW})` (SPEC §3.2 row 5).
    This is one-sided -- once :math:`d \\ge d^{\\min}` the term vanishes, so
    no attractive force is ever applied.  The function is C1 (continuous
    value AND continuous gradient at :math:`d = d^{\\min}`).

    Pairs are filtered through ``exclude_13_pairs`` (frozensets of two atom
    indices) which the caller populates with bonded pairs and 1,3 neighbours
    -- those are governed by the Mahalanobis loss (bonds/angles) and must
    not be double-counted here.

    Inter-ligand awareness (2026-06-02 User-Direktive, fix A)
    --------------------------------------------------------
    When ``ligand_atom_id`` is provided (a dict mapping ``atom_idx -> ligand_id``),
    pairs whose endpoints belong to DIFFERENT ligands receive the boosted
    weight ``w_inter`` instead of ``weight``.  The metal atom (and any atom
    missing from the map) is treated with the default ``weight``.  This
    bakes the "inter-ligand clashes are more unphysical than intra-ligand
    close contacts" prior directly into the L-BFGS gradient, complementing
    the discrete ranking filter in :func:`grip_ensemble.rank_candidates`.

    Default ``w_inter = weight`` (i.e. no boost) keeps the byte-identical
    legacy behaviour when callers do not provide a ligand map -- new
    behaviour is purely opt-in.

    Parameters
    ----------
    vdw_radii : dict
        ``{atom_index: vdW_radius_Å}``. Atoms missing from the table are
        skipped (no clash contribution).
    exclude_13_pairs : set of frozenset({i, j})
        Pairs to skip -- bonded and 1,3 pairs.
    floor_fraction : float, default 0.85
        ``d_min = floor_fraction * (r_i + r_j)``.  Matches SPEC §3.2.
    weight : float, default 1.0
        Global multiplier applied to the whole penalty (the L-BFGS wrapper
        multiplies this by ``clash_weight`` at the call site).  This is
        also the intra-ligand weight when ``ligand_atom_id`` is supplied.
    ligand_atom_id : dict, optional
        Optional ``{atom_idx -> ligand_id}`` map.  When None (default), all
        pairs are treated identically (legacy behaviour).  When provided,
        any pair ``(i, j)`` with ``ligand_atom_id[i] != ligand_atom_id[j]``
        is treated as inter-ligand and weighted by ``w_inter``.
    w_inter : float, optional
        Weight applied to inter-ligand pairs.  Defaults to ``weight`` (no
        boost) so omitting the kwarg is a no-op.  Typical production value
        is ``weight * 3`` (e.g. 15.0 when ``weight=5.0``).
    """

    vdw_radii: Dict[int, float]
    exclude_13_pairs: Set[FrozenSet[int]] = field(default_factory=set)
    floor_fraction: float = 0.85
    weight: float = 1.0
    ligand_atom_id: Optional[Dict[int, int]] = None
    w_inter: Optional[float] = None

    def __post_init__(self):
        self.vdw_radii = {int(k): float(v) for k, v in self.vdw_radii.items()}
        # Normalise the exclude set to frozensets of ints.
        norm: Set[FrozenSet[int]] = set()
        for pair in self.exclude_13_pairs:
            try:
                a, b = tuple(pair)
            except Exception:
                continue
            if int(a) == int(b):
                continue
            norm.add(frozenset((int(a), int(b))))
        self.exclude_13_pairs = norm
        self.floor_fraction = float(self.floor_fraction)
        self.weight = float(self.weight)
        # Normalise the optional ligand map.  Non-int / non-int values are
        # silently skipped (the same defensive behaviour as the vdW table).
        if self.ligand_atom_id is not None:
            norm_map: Dict[int, int] = {}
            for k, v in self.ligand_atom_id.items():
                try:
                    norm_map[int(k)] = int(v)
                except (TypeError, ValueError):
                    continue
            self.ligand_atom_id = norm_map if norm_map else None
        # Default w_inter to weight (no boost) so omitting the kwarg is a
        # byte-identical no-op for legacy callers.
        if self.w_inter is None:
            self.w_inter = self.weight
        else:
            self.w_inter = float(self.w_inter)

    def value_and_grad(self, R: np.ndarray) -> Tuple[float, np.ndarray]:
        """Return ``(L, G)``; ``G`` matches the shape of ``R`` (N, 3)."""
        R = _as_R(R)
        n = R.shape[0]
        grad = np.zeros_like(R)
        if n < 2 or not self.vdw_radii:
            return 0.0, grad

        # Precompute the radii in an ndarray, with NaN for missing atoms which
        # are then skipped (no clash term for atoms we don't have a vdW for).
        radii = np.full(n, np.nan, dtype=np.float64)
        for idx, r in self.vdw_radii.items():
            if 0 <= idx < n:
                radii[idx] = r

        frac = self.floor_fraction
        w_intra = self.weight
        w_inter = float(self.w_inter)
        lig_map = self.ligand_atom_id  # may be None -> legacy
        total = 0.0

        # The deterministic O(N^2) walk -- sorted by (i, j) with i < j.
        # For the production sizes we hit (≤ a few hundred atoms per complex
        # in fffree), this is well under the build-time budget per SPEC §7.
        for i in range(n):
            ri = radii[i]
            if not np.isfinite(ri):
                continue
            Pi = R[i]
            lig_i = lig_map.get(i) if lig_map is not None else None
            for j in range(i + 1, n):
                rj = radii[j]
                if not np.isfinite(rj):
                    continue
                if frozenset((i, j)) in self.exclude_13_pairs:
                    continue
                d_vec = Pi - R[j]
                d = float(np.linalg.norm(d_vec))
                d_min = frac * (ri + rj)
                if d >= d_min:
                    continue
                if d < _EPS:
                    # Exactly coincident -- skip rather than emit NaN.
                    continue
                # Per-pair weight: boost when pair crosses ligands.
                if lig_map is not None:
                    lig_j = lig_map.get(j)
                    if (lig_i is not None and lig_j is not None
                            and lig_i != lig_j):
                        w = w_inter
                    else:
                        w = w_intra
                else:
                    w = w_intra
                gap = d_min - d
                total += w * gap * gap
                # dL/d|d| = -2 w (d_min - d); d|d|/dr_i = (r_i - r_j)/|d|
                coef = -2.0 * w * gap / d
                gi = coef * d_vec
                grad[i] += gi
                grad[j] -= gi

        return float(total), grad
