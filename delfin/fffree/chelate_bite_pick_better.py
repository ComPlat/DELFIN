"""chelate_bite_pick_better.py — per-structure pick-better gate for the
σ-tail CHELATE_BITE constrained per-chelate embed.

Background (memory: ``feedback_sigma_tail_cn6_chelate_distortion`` 2026-05-27)
============================================================================
σ-tail forensik on Iter-27 full-pool: CN6 chelates distort the octahedron
(27% tail vs 8% monodentate).  Pinning the donor-donor distance to the
polyhedron's ideal vertex spacing (DG bounds on the per-chelate metallacycle
embed) pulls the bite angle to 90° and DOES drive the median CShM down
(ABUMET 11.97 → 0.055, etc.), but on rigid/infeasible chelates it forces a
geometry that the ring cannot accommodate → catastrophic regressions (D-YOMHIR
0.6 → 36.1, D-DIMXUS 0.9 → 24.8, D-JAHZUQ 0.7 → 23.3).  Median improved,
mean got worse.  The fix was VALIDATED but smoke MIXED, so it was NOT
shipped: it stayed env-gated default-OFF.

This module implements the per-structure pick-better gate the memory note
prescribes: when ``DELFIN_FFFREE_CHELATE_BITE_PICK_BETTER=1`` AND
``DELFIN_FFFREE_CHELATE_BITE=1`` are BOTH set, the caller embeds the
chelate BOTH ways (constrained + unconstrained), measures the local
polyhedron CShM of the placed donors vs the assigned vertex set, and
keeps the lower-CShM result.  The catastrophic regressions are thus
filtered out STRUCTURALLY (each chelate that gets worse falls back to
the unconstrained behaviour) while the wins are kept.

Default-OFF byte-identical: when either env flag is unset, this module
is not consulted, so the build path is identical to HEAD ec7fb0d.
"""
from __future__ import annotations
import os
from typing import Optional, Sequence, Tuple

import numpy as np


_ENV_BITE = "DELFIN_FFFREE_CHELATE_BITE"
_ENV_PICK = "DELFIN_FFFREE_CHELATE_BITE_PICK_BETTER"


def pick_better_active() -> bool:
    """Return True iff BOTH env-gates are set.  When False, callers must
    take the legacy path (byte-identical to HEAD)."""
    return (os.environ.get(_ENV_BITE, "0") == "1"
            and os.environ.get(_ENV_PICK, "0") == "1")


def _donor_unit_vectors(P: np.ndarray, donor_idxs: Sequence[int],
                        metal_pos: Optional[np.ndarray] = None) -> Optional[np.ndarray]:
    """Return unit vectors metal -> donor (shape (k, 3)), or None on
    degenerate input.  ``P`` are the coords AFTER orient (metal-frame:
    metal at origin) unless ``metal_pos`` is given."""
    if P is None:
        return None
    P = np.asarray(P, dtype=float)
    if P.ndim != 2 or P.shape[1] != 3:
        return None
    if metal_pos is None:
        m = np.zeros(3, dtype=float)
    else:
        m = np.asarray(metal_pos, dtype=float)
    try:
        vecs = P[list(donor_idxs)] - m
    except Exception:
        return None
    norms = np.linalg.norm(vecs, axis=1)
    if np.any(norms < 1e-9):
        return None
    return vecs / norms[:, None]


def cshm_local(P: np.ndarray, donor_idxs: Sequence[int],
               geometry: str, metal_pos: Optional[np.ndarray] = None) -> float:
    """Continuous shape measure of the placed donors vs the ``geometry``
    polyhedron.  Returns ``+inf`` on failure so a failing build always
    LOSES the pick.  Delegates to :func:`delfin.fffree.polyhedra.cshm`.
    """
    if not geometry or donor_idxs is None or len(donor_idxs) < 2:
        return float("inf")
    obs = _donor_unit_vectors(P, donor_idxs, metal_pos=metal_pos)
    if obs is None:
        return float("inf")
    try:
        from delfin.fffree import polyhedra as PLY
        return float(PLY.cshm(obs, geometry))
    except Exception:
        return float("inf")


def pick_better(coords_constrained: Optional[np.ndarray],
                coords_unconstrained: Optional[np.ndarray],
                donor_idxs: Sequence[int],
                geometry: str,
                metal_pos: Optional[np.ndarray] = None,
                ) -> Tuple[Optional[np.ndarray], str, float, float]:
    """Compute CShM for both placements and pick the lower-CShM one.

    Parameters
    ----------
    coords_constrained, coords_unconstrained
        Two candidate placements of the chelate's atoms (same atom order,
        usually in metal frame: metal at origin).  Either may be ``None``
        (embed failure) — the other wins by default.
    donor_idxs
        Local atom indices of the donors in the chelate block.
    geometry
        Polyhedron name (e.g. ``"OC-6 octahedron"``) for the CShM gauge.
    metal_pos
        Optional metal position; default ``[0,0,0]`` (metal-frame).

    Returns
    -------
    (coords, source, cshm_constrained, cshm_unconstrained)
        ``coords`` = winner; ``source`` in ``{"constrained", "unconstrained",
        "none"}``; both CShMs returned for diagnostics / smoke gating.

    Tie-breaking
    ------------
    On exact tie (cshm_c == cshm_u) the CONSTRAINED placement wins (it was
    requested by the env flag).  This is a deterministic, well-ordered rule.
    """
    cshm_c = cshm_local(coords_constrained, donor_idxs, geometry, metal_pos=metal_pos) \
        if coords_constrained is not None else float("inf")
    cshm_u = cshm_local(coords_unconstrained, donor_idxs, geometry, metal_pos=metal_pos) \
        if coords_unconstrained is not None else float("inf")
    # Both failed
    if not np.isfinite(cshm_c) and not np.isfinite(cshm_u):
        return None, "none", cshm_c, cshm_u
    # One failed
    if not np.isfinite(cshm_u):
        return coords_constrained, "constrained", cshm_c, cshm_u
    if not np.isfinite(cshm_c):
        return coords_unconstrained, "unconstrained", cshm_c, cshm_u
    # Both finite -> pick lower; tie -> constrained
    if cshm_c <= cshm_u:
        return coords_constrained, "constrained", cshm_c, cshm_u
    return coords_unconstrained, "unconstrained", cshm_c, cshm_u


def md_invariant_preserved(coords_before: np.ndarray, coords_after: np.ndarray,
                           donor_idxs: Sequence[int],
                           metal_pos: Optional[np.ndarray] = None,
                           tol: float = 0.05) -> bool:
    """Defence-in-depth: assert ‖metal-donor‖ is preserved (within ``tol``)
    by the pick.  Used by the assemble_from_config wiring as a sanity guard
    before swapping in the picked coords."""
    if coords_before is None or coords_after is None:
        return False
    a = np.asarray(coords_before, dtype=float)
    b = np.asarray(coords_after, dtype=float)
    if a.shape != b.shape:
        return False
    if metal_pos is None:
        m = np.zeros(3, dtype=float)
    else:
        m = np.asarray(metal_pos, dtype=float)
    for d in donor_idxs:
        d_old = float(np.linalg.norm(a[d] - m))
        d_new = float(np.linalg.norm(b[d] - m))
        if abs(d_old - d_new) > tol:
            return False
    return True
