"""delfin.fffree.build_time_clash_gate — Per-candidate collapse / clash gate.

Construction-side fix #3 (User mandate 2026-06-03): structqual_xh_collapse
+154% and heavy_collapse +/- 0 still fire in voll-pool because the
construction phase emits candidates that already contain X-H or heavy-heavy
collapses (ETKDG numerics + degenerate metallacycle embeds), and the existing
`_has_collapsed_heavy_bonds` guard is env-gated default-OFF.

This module exposes a UNIVERSAL build-time check that mirrors the structqual
criteria EXACTLY (so a candidate that passes the gate cannot fail structqual):

  1. SUPERPOSITION   any-pair distance < SUPERPOS_TOL (0.5 A).
  2. HEAVY_COLLAPSE  any heavy-heavy bonded pair < HH_COLLAPSE_RATIO (0.70) x
                     element-pair ideal length.
  3. XH_COLLAPSE     any X-H bonded pair < XH_COLLAPSE_TOL (0.70 A).

Geometry-only, deterministic, no SMILES.  Env-gated default-OFF byte-identical.

Surgical fix F24-interlig (User 2026-06-03): the original 3 criteria are
intra-candidate.  Inter-fragment (inter-ligand) van-der-Waals clashes (CO-CO
adjacent axes, PMe3 methyl-methyl, π-π too-close) are NOT covered and
re-emerge in the voll-pool interlig metric.  ``compute_interlig_vdw_penalty``
adds a non-bonded inter-fragment vdW penalty (sum of squared overlaps) so
the candidate-selection loop deprioritises collided configurations.  Env-
gated default-OFF byte-identical (``DELFIN_FFFREE_BUILD_INTERLIG_PENALTY``).

API:
    has_collapse(syms, P) -> bool
    collapse_count(syms, P) -> int       # total violation count for ranking
    collapse_details(syms, P) -> list[dict]   # per-violation records
    apply_build_gate(syms, P, candidates_iter) -> chosen_P | None
    compute_interlig_vdw_penalty(P, frag_to_atoms, vdw_radii, factor=0.85,
        return_count=False) -> float | (float, int)
    interlig_penalty_active() -> bool

Env flag:
    DELFIN_FFFREE_BUILD_CLASH_GATE=1    (per-fix activation)
    DELFIN_FFFREE_BUILD_INTERLIG_PENALTY=1  (surgical fix 1 activation)
    DELFIN_FFFREE_BUILD_INTERLIG_PENALTY_W=15.0  (penalty weight, optional)
    DELFIN_FFFREE_BUILD_INTERLIG_PENALTY_FACTOR=0.85  (vdW fraction, optional)
    DELFIN_FFFREE_F24_INTERLIG_FIX_ALL=1   (F24-interlig surgical master)
    DELFIN_FFFREE_CONSTRUCTION_FIX_ALL=1   (master activation)
"""
from __future__ import annotations

import os
from typing import Dict, Iterable, List, Mapping, Optional, Sequence, Tuple

import numpy as np

import delfin._bond_decollapse as _bd


# Calibrated thresholds: identical to structure_quality.py so the build-gate
# fires on EXACTLY the cases structqual would fail.
SUPERPOS_TOL = 0.5
HH_COLLAPSE_RATIO = 0.70
XH_COLLAPSE_TOL = 0.70
_BOND_FACTOR = 1.30

# vdW radii (Bondi 1964 + Alvarez 2013 for TMs).  Mirrors
# ``delfin.fffree.grip_polish.DEFAULT_VDW_RADII`` so the build-time penalty
# and the GRIP / detector both see the same vdW geometry.
INTERLIG_VDW_RADII: Dict[str, float] = {
    "H": 1.20, "He": 1.40, "Li": 1.82, "Be": 1.53, "B": 1.92,
    "C": 1.70, "N": 1.55, "O": 1.52, "F": 1.47, "Ne": 1.54,
    "Na": 2.27, "Mg": 1.73, "Al": 1.84, "Si": 2.10, "P": 1.80,
    "S": 1.80, "Cl": 1.75, "Ar": 1.88,
    "K": 2.75, "Ca": 2.31, "Sc": 2.11, "Ti": 1.95, "V": 2.06,
    "Cr": 1.89, "Mn": 1.97, "Fe": 1.94, "Co": 1.92, "Ni": 1.84,
    "Cu": 1.86, "Zn": 2.10, "Ga": 1.87, "Ge": 2.11, "As": 1.85,
    "Se": 1.90, "Br": 1.85, "Kr": 2.02,
    "Rb": 3.03, "Sr": 2.49, "Y": 2.11, "Zr": 1.86, "Nb": 2.07,
    "Mo": 2.09, "Ru": 1.97, "Rh": 1.95, "Pd": 2.02, "Ag": 2.03,
    "Cd": 2.30, "In": 1.93, "Sn": 2.17, "Sb": 2.06, "Te": 2.06,
    "I": 1.98, "Xe": 2.16,
    "Cs": 3.43, "Ba": 2.68, "La": 2.43,
    "Hf": 2.12, "Ta": 2.17, "W": 2.10, "Re": 2.17, "Os": 2.16,
    "Ir": 2.13, "Pt": 2.13, "Au": 2.14, "Hg": 2.23, "Tl": 1.96,
    "Pb": 2.02, "Bi": 2.07,
}

# Surgical-fix-1 defaults.  ``_INTERLIG_FACTOR`` (default 0.85) is the vdW-sum
# fraction below which a non-bonded inter-fragment pair counts as a clash —
# same threshold as
# :func:`delfin.fffree.grip_ensemble.count_inter_ligand_clashes`.
# ``_INTERLIG_PENALTY_W`` (default 15.0) is the quadratic-overlap weight
# wired into the candidate selection so a single 1.85 A C-O contact
# dominates over a clean alternative.
_INTERLIG_FACTOR_DEFAULT = 0.85
_INTERLIG_PENALTY_W_DEFAULT = 15.0


def _flag_active() -> bool:
    if os.environ.get("DELFIN_FFFREE_CONSTRUCTION_FIX_ALL", "0") == "1":
        return True
    if os.environ.get("DELFIN_FFFREE_F24_INTERLIG_FIX_ALL", "0") == "1":
        return True
    return os.environ.get("DELFIN_FFFREE_BUILD_CLASH_GATE", "0") == "1"


def interlig_penalty_active() -> bool:
    """True iff the inter-ligand vdW penalty (surgical fix 1) is enabled.

    Activated by either the per-fix flag
    ``DELFIN_FFFREE_BUILD_INTERLIG_PENALTY=1``, the F24-interlig surgical
    master ``DELFIN_FFFREE_F24_INTERLIG_FIX_ALL=1``, or the global
    construction master ``DELFIN_FFFREE_CONSTRUCTION_FIX_ALL=1``.  Default
    OFF -> byte-identical to pre-fix HEAD.
    """
    if os.environ.get("DELFIN_FFFREE_CONSTRUCTION_FIX_ALL", "0") == "1":
        return True
    if os.environ.get("DELFIN_FFFREE_F24_INTERLIG_FIX_ALL", "0") == "1":
        return True
    return os.environ.get("DELFIN_FFFREE_BUILD_INTERLIG_PENALTY", "0") == "1"


def _interlig_weight() -> float:
    """Read ``DELFIN_FFFREE_BUILD_INTERLIG_PENALTY_W`` (default 15.0)."""
    try:
        return float(os.environ.get(
            "DELFIN_FFFREE_BUILD_INTERLIG_PENALTY_W",
            str(_INTERLIG_PENALTY_W_DEFAULT),
        ))
    except (TypeError, ValueError):
        return _INTERLIG_PENALTY_W_DEFAULT


def _interlig_factor() -> float:
    """Read ``DELFIN_FFFREE_BUILD_INTERLIG_PENALTY_FACTOR`` (default 0.85)."""
    try:
        return float(os.environ.get(
            "DELFIN_FFFREE_BUILD_INTERLIG_PENALTY_FACTOR",
            str(_INTERLIG_FACTOR_DEFAULT),
        ))
    except (TypeError, ValueError):
        return _INTERLIG_FACTOR_DEFAULT


def _min_pair_distance(P: np.ndarray) -> float:
    """Minimum atom-atom distance in the structure (any pair)."""
    n = len(P)
    if n < 2:
        return float("inf")
    mn = float("inf")
    for i in range(n):
        if i + 1 >= n:
            break
        d = np.linalg.norm(P[i + 1:] - P[i], axis=1)
        if d.size:
            md = float(d.min())
            if md < mn:
                mn = md
    return mn


def _geometric_bonds(syms: Sequence[str], P: np.ndarray
                     ) -> List[Tuple[int, int]]:
    """All geometric bonds (including X-H) via covalent-radii sum cutoff.

    Excludes metal bonds (matches structure_quality behaviour).
    """
    n = len(syms)
    out: List[Tuple[int, int]] = []
    for i in range(n):
        if _bd._is_metal(syms[i]):
            continue
        for j in range(i + 1, n):
            if _bd._is_metal(syms[j]):
                continue
            d = float(np.linalg.norm(P[i] - P[j]))
            if d < _BOND_FACTOR * _bd._ideal_bond(syms[i], syms[j]):
                out.append((i, j))
    return out


def has_collapse(syms: Sequence[str], P: np.ndarray) -> bool:
    """True iff ANY of the 3 structqual collapse criteria fires.

    Env-gated default-OFF: returns False when no flag set (no rejection ->
    byte-identical behaviour).

    When ``DELFIN_FFFREE_HH_CLASH_INCLUDE=1`` AND the build-clash-gate
    master flag is also set, additionally returns True when ANY
    non-bonded H-H pair is below the Bondi vdW floor (geminal/1-3
    excluded).  Independent of the H-H detector activation, the
    existing 3 criteria are unchanged.
    """
    if not _flag_active():
        return False
    if len(syms) == 0:
        return False
    P = np.asarray(P, dtype=float)
    # 1. superposition
    if _min_pair_distance(P) < SUPERPOS_TOL:
        return True
    # 2 + 3 via geometric bonds
    for i, j in _geometric_bonds(syms, P):
        d = float(np.linalg.norm(P[i] - P[j]))
        if syms[i] == "H" or syms[j] == "H":
            # X-H bond
            if d < XH_COLLAPSE_TOL:
                return True
        else:
            ideal = _bd._ideal_bond(syms[i], syms[j])
            if d < HH_COLLAPSE_RATIO * ideal:
                return True
    # 4. (opt-in) non-bonded H-H below Bondi vdW floor.
    try:
        from delfin.fffree.hh_clash_detector import (
            hh_clash_active as _hh_active,
            count_hh_clashes as _hh_count,
        )
        if _hh_active() and _hh_count(syms, P) > 0:
            return True
    except Exception:
        pass
    return False


def collapse_count(syms: Sequence[str], P: np.ndarray) -> int:
    """Total number of collapse-violation atom-pairs (used for ranking).

    Always counts (no env flag) -- this is read-only.  Callers gate on
    `has_collapse` and use `collapse_count` to pick the candidate with the
    FEWEST violations when none is fully clean.

    When ``DELFIN_FFFREE_HH_CLASH_INCLUDE=1`` is set, additionally
    accumulates non-bonded H-H pairs below the vdW floor (Bondi 0.85 ×
    2 r_H = 2.04 Å, geminal/1-3 H-H within CH3 are skipped via topology).
    Default OFF -> byte-identical to HEAD ``93b396d``.
    """
    if len(syms) == 0:
        return 0
    P = np.asarray(P, dtype=float)
    n_pairs = 0
    # 1. superposition
    if _min_pair_distance(P) < SUPERPOS_TOL:
        # count all overlapping pairs
        n = len(P)
        for i in range(n):
            for j in range(i + 1, n):
                if float(np.linalg.norm(P[i] - P[j])) < SUPERPOS_TOL:
                    n_pairs += 1
    for i, j in _geometric_bonds(syms, P):
        d = float(np.linalg.norm(P[i] - P[j]))
        if syms[i] == "H" or syms[j] == "H":
            if d < XH_COLLAPSE_TOL:
                n_pairs += 1
        else:
            ideal = _bd._ideal_bond(syms[i], syms[j])
            if d < HH_COLLAPSE_RATIO * ideal:
                n_pairs += 1
    # Opt-in: H-H non-bonded clashes (default OFF byte-identical).
    try:
        from delfin.fffree.hh_clash_detector import (
            hh_clash_active as _hh_active,
            count_hh_clashes as _hh_count,
        )
        if _hh_active():
            n_pairs += int(_hh_count(syms, P))
    except Exception:
        pass
    return n_pairs


def collapse_details(syms: Sequence[str], P: np.ndarray) -> List[dict]:
    """Per-violation records (used for tests/debug)."""
    if len(syms) == 0:
        return []
    P = np.asarray(P, dtype=float)
    out: List[dict] = []
    # 1. superposition: only any-pair distances under tol
    n = len(P)
    for i in range(n):
        for j in range(i + 1, n):
            d = float(np.linalg.norm(P[i] - P[j]))
            if d < SUPERPOS_TOL:
                out.append({"mode": "superposition", "i": i, "j": j,
                            "a": syms[i], "b": syms[j], "d": round(d, 4)})
    for i, j in _geometric_bonds(syms, P):
        d = float(np.linalg.norm(P[i] - P[j]))
        if syms[i] == "H" or syms[j] == "H":
            if d < XH_COLLAPSE_TOL:
                out.append({"mode": "xh_collapse", "i": i, "j": j,
                            "a": syms[i], "b": syms[j], "d": round(d, 4)})
        else:
            ideal = _bd._ideal_bond(syms[i], syms[j])
            if d < HH_COLLAPSE_RATIO * ideal:
                out.append({"mode": "heavy_collapse", "i": i, "j": j,
                            "a": syms[i], "b": syms[j], "d": round(d, 4),
                            "ratio": round(d / max(ideal, 1e-9), 4)})
    return out


def select_clean_candidate(
    syms_per_candidate: Sequence[Sequence[str]],
    P_per_candidate: Sequence[np.ndarray],
) -> Optional[int]:
    """Pick the index of the first FULLY-CLEAN candidate (collapse_count==0),
    falling back to the candidate with the FEWEST collapse violations.

    Returns None only when there are no candidates.

    Env-gated default-OFF: when no flag, returns 0 (first candidate)
    unconditionally -> byte-identical to picking the first candidate.
    """
    if not syms_per_candidate or not P_per_candidate:
        return None
    if not _flag_active():
        return 0
    best_idx = 0
    best_count = collapse_count(syms_per_candidate[0], P_per_candidate[0])
    if best_count == 0:
        return 0
    for k in range(1, len(syms_per_candidate)):
        c = collapse_count(syms_per_candidate[k], P_per_candidate[k])
        if c < best_count:
            best_count = c
            best_idx = k
            if best_count == 0:
                return best_idx
    return best_idx


def compute_interlig_vdw_penalty(
    coords: np.ndarray,
    fragment_to_atoms: Mapping[int, Sequence[int]],
    vdw_radii: Sequence[float],
    factor: float = _INTERLIG_FACTOR_DEFAULT,
    return_count: bool = False,
):
    """Sum of squared vdW-overlaps between atom pairs in DIFFERENT fragments.

    Implements surgical fix 1 (F24/interlig deep forensik, User 2026-06-03):
    the existing intra-Q collapse criteria do not see inter-ligand
    non-bonded contacts (CO-CO 1.85 A, PMe3 methyl-methyl 2.12 A,
    aromatic-aromatic 2.01 A).  This routine adds a quadratic vdW-overlap
    penalty between fragments so the selection loop deprioritises collided
    configurations.

    Parameters
    ----------
    coords : (N, 3) ndarray
        Cartesian atom coordinates (Å).
    fragment_to_atoms : mapping of int -> sequence[int]
        Fragment-id -> atom-index list.  Atom indices that appear in NO
        fragment are excluded (e.g. the metal atom when M is at index 0
        and only ligand atoms are passed in).  Each atom index must
        appear in at most ONE fragment (deterministic; duplicates ignored).
    vdw_radii : sequence[float]
        Per-atom vdW radius (Å).  Use ``np.nan`` to skip an atom (e.g.
        unknown element).
    factor : float, default 0.85
        vdW-sum fraction below which a pair counts.  Penalty per pair is
        ``max(0, factor * (r_i + r_j) - d_ij)^2`` so a clean (above-target)
        pair contributes 0.0.  Matches the
        :mod:`grip_ensemble.count_inter_ligand_clashes` threshold.
    return_count : bool, default False
        When True, also return the integer clash-pair count.

    Returns
    -------
    float  (or  (float, int)  when ``return_count=True``)
        Sum of squared overlaps (Å^2).  Always >= 0.0.  Atom-symbol
        lookup is the caller's responsibility (use
        ``[INTERLIG_VDW_RADII.get(s, np.nan) for s in syms]``).

    Notes
    -----
    Deterministic: iterates fragments in sorted key order and atom indices
    in ascending order.  Geometry-only: no SMILES, no chemistry tags.
    Cost is O(N_inter_pairs) where N_inter_pairs is the number of
    cross-fragment atom pairs.
    """
    P = np.asarray(coords, dtype=float)
    if P.ndim == 1:
        if P.size % 3 != 0:
            raise ValueError("coords must reshape to (N, 3)")
        P = P.reshape(-1, 3)
    n = P.shape[0]
    radii = np.asarray(vdw_radii, dtype=float)
    if radii.shape[0] != n:
        raise ValueError(
            f"vdw_radii length {radii.shape[0]} != coords rows {n}"
        )
    # Atom -> fragment-id map (deterministic; first-fragment wins on dup).
    atom_to_frag: Dict[int, int] = {}
    frag_ids = sorted(fragment_to_atoms.keys())
    for fid in frag_ids:
        for a in fragment_to_atoms[fid]:
            ai = int(a)
            if 0 <= ai < n and ai not in atom_to_frag:
                atom_to_frag[ai] = fid

    penalty = 0.0
    clash_count = 0
    for i in range(n):
        if i not in atom_to_frag:
            continue
        ri = radii[i]
        if not np.isfinite(ri):
            continue
        fi = atom_to_frag[i]
        Pi = P[i]
        for j in range(i + 1, n):
            if j not in atom_to_frag:
                continue
            fj = atom_to_frag[j]
            if fi == fj:
                continue  # intra-fragment -- handled by intra-Q gates
            rj = radii[j]
            if not np.isfinite(rj):
                continue
            d = float(np.linalg.norm(Pi - P[j]))
            target = float(factor) * (float(ri) + float(rj))
            overlap = target - d
            if overlap > 0.0:
                penalty += overlap * overlap
                clash_count += 1
    if return_count:
        return penalty, clash_count
    return penalty


def _interlig_penalty_for_pair(
    syms_a: Sequence[str], P_a: np.ndarray,
    syms_b: Sequence[str], P_b: np.ndarray,
    factor: Optional[float] = None,
) -> Tuple[float, int]:
    """Convenience: vdW-overlap penalty between two atom blocks ``a`` and ``b``.

    Both blocks are treated as one fragment each (cross-block pairs only).
    Returns (penalty_squared_A^2, clash_pair_count).  Used by the
    assemble-complex hook to score a candidate ``Q`` against the
    already-placed atoms ``placed`` without changing the call shape.
    """
    if factor is None:
        factor = _interlig_factor()
    Pa = np.asarray(P_a, dtype=float).reshape(-1, 3)
    Pb = np.asarray(P_b, dtype=float).reshape(-1, 3)
    na, nb = Pa.shape[0], Pb.shape[0]
    if na == 0 or nb == 0:
        return 0.0, 0
    ra = np.array(
        [INTERLIG_VDW_RADII.get(str(s), np.nan) for s in syms_a],
        dtype=float,
    )
    rb = np.array(
        [INTERLIG_VDW_RADII.get(str(s), np.nan) for s in syms_b],
        dtype=float,
    )
    penalty = 0.0
    clash_count = 0
    for i in range(na):
        ri = ra[i]
        if not np.isfinite(ri):
            continue
        Pi = Pa[i]
        for j in range(nb):
            rj = rb[j]
            if not np.isfinite(rj):
                continue
            d = float(np.linalg.norm(Pi - Pb[j]))
            target = float(factor) * (float(ri) + float(rj))
            overlap = target - d
            if overlap > 0.0:
                penalty += overlap * overlap
                clash_count += 1
    return penalty, clash_count


__all__ = [
    "has_collapse", "collapse_count", "collapse_details",
    "select_clean_candidate", "_flag_active",
    "compute_interlig_vdw_penalty", "interlig_penalty_active",
    "_interlig_penalty_for_pair",
    "INTERLIG_VDW_RADII",
    "SUPERPOS_TOL", "HH_COLLAPSE_RATIO", "XH_COLLAPSE_TOL",
]
