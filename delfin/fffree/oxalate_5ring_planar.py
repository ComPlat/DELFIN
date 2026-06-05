"""oxalate_5ring_planar.py — post-orient projector for the M-O-C-C-O 5-ring
oxalate chelate.

Background (Iter-32e YILNUF class, 2026-06-05)
==============================================
The DG metallacycle embed (``_embed_metallacycle``) handles generic chelate
ring topology, but oxalate (O=C-C(=O)-O bidentate via both deprotonated Os
binding the same metal) is geometrically peculiar:

* The 5-ring is **planar** because both carboxylate carbons are sp²: the
  whole O-C(=O)-C(=O)-O system is conjugated.
* The bite angle O-M-O is narrow (**~78-82°** in CSD / literature) because
  the C-C bond pulls the two Os toward each other.
* The dihedral O-C-C-O across the ring should be **~0°** (eclipsed planar)
  — the exocyclic C=O groups likewise stay coplanar with the ring.

Forensik on the 2792332-aromatic-symmetry-VOLLPOOL (9868 frames, 38 files
contain an oxalate 5-ring → 96 chelates total) measured the current state:

==========================  ===========  =========
Metric                      Mean ± Std   Anomaly%
==========================  ===========  =========
bite O-M-O (deg)            59 ± 19      73% (outside 70-85°)
ring planarity (Å, max dev) 0.25 ± 0.27  54% (>0.10 Å)
dihedral O-C-C-O (deg)      42 ± 54      49% (>15°)
C-C distance (Å)            1.53 ± 0.24  83% (outside 1.45-1.65 Å)
==========================  ===========  =========

The ring is collapsing into non-planar / twisted geometries because the
generic chelate embed doesn't know oxalate's sp² constraint.

Fix (this module)
=================
**Detection** (universal, geometry-only, SMARTS-free):

    1. Find oxalate fragments in a ligand block by graph topology:
       two sp² C atoms bonded to each other, each carrying ≥2 O bonds, and
       two of those Os bond the same metal at chelating distance.

**Projection** (in-place, deterministic):

    2. Choose a "ring plane" by least-squares fit through M + 4 ring atoms
       (O1, C1, C2, O2). Project all 5 ring atoms onto that plane.
    3. Re-snap intra-ring distances to literature ideals
       (C-C = 1.55 Å, C-O endo = 1.27 Å, O-M from MSB table,
       O-M-O bite = 78° → O...O target derived from M-O distance and bite).
       Solved as a 2D linear-algebra problem in the chosen plane.
    4. Project the exocyclic C=O (the "other" oxygen on each carboxylate C)
       onto the ring plane so the entire OC-CO π system stays planar.
       Substituents on those exocyclic oxygens (e.g. H, alkyl) are
       translated with their parent oxygen — internal bonds preserved.

**Env-gate**: ``DELFIN_FFFREE_OXALATE_5RING_PLANAR=1`` activates the fix.
Default OFF → this module is a no-op and ``apply()`` returns the input
coords unchanged (byte-identical).
"""
from __future__ import annotations
import os
from typing import List, Optional, Sequence, Tuple

import numpy as np

import delfin._bond_decollapse as _bd

_ENV = "DELFIN_FFFREE_OXALATE_5RING_PLANAR"

# Literature targets (Mogul-COD on oxalato-metal fragments; conservative).
_TARGET_BITE_DEG = 78.0          # O-M-O bite angle (median ~77-82° across TMs)
_TARGET_CC = 1.55                # C-C bond, sp²-sp² conjugated
_TARGET_CO_ENDO = 1.275          # C-O ring-side (delocalised carboxylate)
_TARGET_CO_EXO = 1.225           # C=O exocyclic
# Allowable ring plane deviation BEFORE we touch anything (already-flat rings
# are left untouched to stay byte-identical for the wins).
_PLANE_DEV_TRIGGER_A = 0.05


def flag_active() -> bool:
    """True iff the env-gate is set.  When False, callers MUST short-circuit
    to byte-identical behaviour."""
    return os.environ.get(_ENV, "0") == "1"


# ---------------------------------------------------------------------------
# Detection
# ---------------------------------------------------------------------------

def _bonded(syms: Sequence[str], P: np.ndarray, i: int, j: int,
            factor: float = 1.30) -> bool:
    """Heavy-bond proximity test (covalent-sum × factor)."""
    if i == j:
        return False
    d = float(np.linalg.norm(P[i] - P[j]))
    ideal = _bd._ideal_bond(syms[i], syms[j])
    return d < factor * ideal


def detect_oxalate_5rings(syms: Sequence[str], P: np.ndarray,
                          donor_idxs: Sequence[int],
                          metal_pos: np.ndarray,
                          metal_sym: Optional[str] = None) -> List[dict]:
    """Find all oxalate 5-ring chelates with donors in ``donor_idxs``.

    ``P`` is the LIGAND-ONLY coordinate matrix (metal is NOT a row in P).
    ``metal_pos`` is the metal's position in the same frame as ``P``
    (typically ``np.zeros(3)`` because the chelate is in metal-frame).
    ``metal_sym`` is optional — used only to choose a covalent-radius
    threshold for M-O coordination; defaults to a generous 2.6 Å gate.

    A chelate qualifies iff:
      - two carbons C1, C2 in the ligand block are bonded to each other,
      - each Ci has ≥ 2 oxygen neighbours,
      - the two donor Os in ``donor_idxs`` bond C1 / C2 respectively AND
        both sit at chelating distance from ``metal_pos``.

    Returns list of dicts (empty if no fragment qualifies):
        ``{"O1": int, "C1": int, "C2": int, "O2": int,
           "O1_exo": [ints], "O2_exo": [ints]}``
    Universal: no SMARTS, no ligand-name knowledge.
    """
    n = len(syms)
    P = np.asarray(P, dtype=float)
    metal_pos = np.asarray(metal_pos, dtype=float)
    if n == 0:
        return []
    donors = [int(d) for d in donor_idxs if 0 <= int(d) < n]
    # Only oxygen donors are interesting for oxalate
    O_donors = [d for d in donors if syms[d] == "O"]
    if len(O_donors) < 2:
        return []
    # M-O coordination cut (covalent-sum × 1.35 if metal_sym known)
    if metal_sym is not None:
        try:
            md_ideal = float(_bd._ideal_bond(metal_sym, "O"))
            mo_cut = 1.35 * md_ideal
        except Exception:
            mo_cut = 2.6
    else:
        mo_cut = 2.6
    # All donor Os must be within mo_cut of the metal
    for d in O_donors:
        if float(np.linalg.norm(P[d] - metal_pos)) > mo_cut:
            return []
    # Candidate carboxylate carbons: C bonded to ≥2 O AND to one donor-O
    # First map each donor O to its carbon neighbour
    o_to_c = {}
    for o in O_donors:
        Cs = [j for j in range(n) if syms[j] == "C" and _bonded(syms, P, o, j)]
        if len(Cs) != 1:
            return []
        o_to_c[o] = Cs[0]
    results: List[dict] = []
    seen: set = set()
    for i, oi in enumerate(O_donors):
        for oj in O_donors[i + 1:]:
            ci = o_to_c[oi]
            cj = o_to_c[oj]
            if ci == cj:                # both donor Os on the same C → carboxylate κ²,
                continue                # not the M-O-C-C-O-M oxalate ring.
            if not _bonded(syms, P, ci, cj):
                continue
            # ci must carry ≥2 Os, same for cj (so each carbon has an exocyclic O)
            Os_i = [j for j in range(n) if syms[j] == "O" and _bonded(syms, P, ci, j)]
            Os_j = [j for j in range(n) if syms[j] == "O" and _bonded(syms, P, cj, j)]
            if len(Os_i) < 2 or len(Os_j) < 2:
                continue
            key = tuple(sorted([oi, ci, cj, oj]))
            if key in seen:
                continue
            seen.add(key)
            results.append({
                "O1": int(oi), "C1": int(ci),
                "C2": int(cj), "O2": int(oj),
                "O1_exo": [int(x) for x in Os_i if x != oi],
                "O2_exo": [int(x) for x in Os_j if x != oj],
            })
    return results


# ---------------------------------------------------------------------------
# Geometric helpers
# ---------------------------------------------------------------------------

def _best_fit_plane(pts: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """Return (centroid, unit_normal) of the LS plane through ``pts``."""
    c = pts.mean(axis=0)
    _, _, Vt = np.linalg.svd(pts - c, full_matrices=False)
    n = Vt[-1]
    nn = np.linalg.norm(n)
    if nn < 1e-12:
        return c, np.array([0.0, 0.0, 1.0])
    return c, n / nn


def _project_onto_plane(p: np.ndarray, c: np.ndarray, n: np.ndarray) -> np.ndarray:
    """Project point ``p`` onto plane (centroid ``c``, unit normal ``n``)."""
    return p - n * float(np.dot(p - c, n))


def _max_plane_dev(pts: np.ndarray) -> float:
    """Maximum |signed distance| from the LS plane through ``pts``."""
    if pts.shape[0] < 3:
        return 0.0
    c, n = _best_fit_plane(pts)
    return float(np.max(np.abs((pts - c) @ n)))


def _solve_planar_5ring(P_m: np.ndarray, P_O1: np.ndarray, P_C1: np.ndarray,
                        P_C2: np.ndarray, P_O2: np.ndarray,
                        md_O: float,
                        bite_deg: float = _TARGET_BITE_DEG,
                        cc: float = _TARGET_CC,
                        co: float = _TARGET_CO_ENDO,
                        ) -> Tuple[np.ndarray, np.ndarray, np.ndarray,
                                   np.ndarray, np.ndarray]:
    """Compute idealised in-plane 2D positions for the 5 ring atoms, then
    re-embed them into 3D using the *current* ring orientation.

    The ring topology constrains all 5 positions once we fix:
      - M at origin,
      - axis bisecting O-M-O along +x,
      - M-O distance = ``md_O``,
      - O-M-O bite angle,
      - C-C and C-O endo bond lengths.

    We solve for C1 and C2 in the plane by triangulation from M and Oi.

    Returns
    -------
    (M, O1, C1, C2, O2) — five 3D points, in the SAME orientation as the
    input (plane-normal preserved).  The mid-O-O bisector is mapped back to
    the input mid-O-O bisector so the placed donors stay near their assigned
    vertices; bite is contracted/expanded to the ideal.
    """
    # Current geometry: compute ring plane and an in-plane frame.
    pts = np.array([P_m, P_O1, P_C1, P_C2, P_O2], dtype=float)
    cen, normal = _best_fit_plane(pts)
    # In-plane axes anchored at metal:
    # x_axis = bisector(M->O1, M->O2) projected into the plane
    v1 = P_O1 - P_m
    v2 = P_O2 - P_m
    # Project bisector into the plane
    bis = v1 + v2
    bis_in = bis - normal * float(np.dot(bis, normal))
    bn = np.linalg.norm(bis_in)
    if bn < 1e-9:
        # Degenerate (O1 and O2 collinear through M opposite) — pick an axis
        # orthogonal to v1 inside the plane.
        v1_in = v1 - normal * float(np.dot(v1, normal))
        if np.linalg.norm(v1_in) < 1e-9:
            return P_m, P_O1, P_C1, P_C2, P_O2
        x_axis = v1_in / np.linalg.norm(v1_in)
    else:
        x_axis = bis_in / bn
    # y_axis = normal × x  (right-handed in-plane frame)
    y_axis = np.cross(normal, x_axis)
    ny = np.linalg.norm(y_axis)
    if ny < 1e-9:
        return P_m, P_O1, P_C1, P_C2, P_O2
    y_axis = y_axis / ny

    # Determine which side (sign in y) each donor sits on, so we preserve
    # the orientation of the chelate (mirror-invariant solution).
    sign_O1 = 1.0 if np.dot(v1, y_axis) >= 0 else -1.0
    sign_O2 = -sign_O1   # O2 on the opposite side

    half_bite = np.radians(bite_deg) * 0.5
    cos_h, sin_h = np.cos(half_bite), np.sin(half_bite)
    # 2D positions (x along bisector, y perpendicular within ring plane)
    O1_2d = np.array([md_O * cos_h, sign_O1 * md_O * sin_h])
    O2_2d = np.array([md_O * cos_h, sign_O2 * md_O * sin_h])
    # C1: triangulate from O1 (distance co) and "constraint" that |C1-C2| = cc.
    # Use symmetry: C1 and C2 share the same x-coordinate (mirror across y=0),
    # and (C1.y - C2.y) = sign_O1 * cc.
    # Distance C1-O1 = co constrains:
    #   (C.x - O1.x)^2 + (C1.y - O1.y)^2 = co^2
    # with C1.y = sign_O1 * cc/2, C2.y = -sign_O1 * cc/2.
    C1y = sign_O1 * cc * 0.5
    C2y = -C1y
    dy = C1y - O1_2d[1]   # = sign_O1 * (cc/2 - md_O * sin_h)
    dx2 = co * co - dy * dy
    if dx2 <= 1e-9:
        # Bite too narrow for these bond lengths — give up, return input.
        return P_m, P_O1, P_C1, P_C2, P_O2
    Cx = O1_2d[0] + float(np.sqrt(dx2))  # carbon is further from metal than O on +x axis
    C1_2d = np.array([Cx, C1y])
    C2_2d = np.array([Cx, C2y])

    # Reconstruct 3D: anchor metal back at its input position, x_axis & y_axis
    # define the plane (passing through P_m).
    M3 = P_m.copy()
    O1_3 = P_m + O1_2d[0] * x_axis + O1_2d[1] * y_axis
    O2_3 = P_m + O2_2d[0] * x_axis + O2_2d[1] * y_axis
    C1_3 = P_m + C1_2d[0] * x_axis + C1_2d[1] * y_axis
    C2_3 = P_m + C2_2d[0] * x_axis + C2_2d[1] * y_axis
    return M3, O1_3, C1_3, C2_3, O2_3


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def apply(syms: Sequence[str], P: np.ndarray,
          donor_idxs: Sequence[int],
          metal_pos: Optional[np.ndarray] = None,
          metal_sym: Optional[str] = None,
          plane_trigger_A: float = _PLANE_DEV_TRIGGER_A,
          ) -> Tuple[np.ndarray, List[dict]]:
    """Project oxalate 5-rings onto a planar ideal geometry.

    Parameters
    ----------
    syms, P
        Atom symbols and 3D coords of the (oriented + placed) chelate
        LIGAND block (metal is NOT a row in ``P`` — its position is
        passed separately via ``metal_pos``).
    donor_idxs
        Indices of the assigned donors (must include the two oxalate O
        atoms for the fix to fire).
    metal_pos
        Metal position in the same frame as ``P``.  Default
        ``np.zeros(3)`` (metal-frame: matches the assemble_complex
        convention where the chelate is in metal-at-origin frame).
    metal_sym
        Optional metal element symbol; used only to compute the M-O
        coordination cutoff.  Default: a generous fixed 2.6 Å gate.
    plane_trigger_A
        Only act when the existing ring is non-planar BY MORE THAN this
        threshold AND the bite differs from the ideal by > 4° (skip
        already-good rings to stay byte-identical for the easy cases).
        Default 0.05 Å.

    Returns
    -------
    (P_out, applied)
        ``P_out`` = updated coords (copy on every call).
        ``applied`` = list of dicts describing each ring touched
        (``{"ring": {...}, "plane_dev_before": float, "bite_before": float,
        "plane_dev_after": float, "bite_after": float}``) — empty when the
        env-gate is off or no oxalate fragment qualifies.

    Determinism: identical inputs → identical outputs (no RNG).  Byte-
    identical to HEAD when ``DELFIN_FFFREE_OXALATE_5RING_PLANAR`` is unset.
    """
    P = np.asarray(P, dtype=float)
    if not flag_active() or P.size == 0:
        return P.copy(), []
    if metal_pos is None:
        M_pos = np.zeros(3, dtype=float)
    else:
        M_pos = np.asarray(metal_pos, dtype=float)
    rings = detect_oxalate_5rings(syms, P, donor_idxs, M_pos, metal_sym=metal_sym)
    if not rings:
        return P.copy(), []
    Q = P.copy()
    applied: List[dict] = []
    for ring in rings:
        o1, c1, c2, o2 = ring["O1"], ring["C1"], ring["C2"], ring["O2"]
        ring_pts_before = np.array([M_pos, Q[o1], Q[c1], Q[c2], Q[o2]])
        dev_before = _max_plane_dev(ring_pts_before)
        bite_before = _bite(Q[o1] - M_pos, Q[o2] - M_pos)
        if dev_before < plane_trigger_A and abs(bite_before - _TARGET_BITE_DEG) < 4.0:
            # Already planar AND bite is close-to-ideal: nothing to do, stay
            # byte-identical.
            continue
        md_O = 0.5 * (float(np.linalg.norm(Q[o1] - M_pos))
                      + float(np.linalg.norm(Q[o2] - M_pos)))
        if md_O < 1e-6:
            continue
        # Capture old C positions so exocyclic Os move WITH their parent C
        # (preserves C=O direction modulo plane projection).
        old_C1 = Q[c1].copy()
        old_C2 = Q[c2].copy()
        new_M, new_O1, new_C1, new_C2, new_O2 = _solve_planar_5ring(
            M_pos, Q[o1], Q[c1], Q[c2], Q[o2], md_O=md_O,
        )
        # Compute the new ring plane to push the exocyclic Os onto.
        new_ring_pts = np.array([new_M, new_O1, new_C1, new_C2, new_O2])
        cen, normal = _best_fit_plane(new_ring_pts)
        # Apply ring updates.
        Q[o1] = new_O1
        Q[c1] = new_C1
        Q[c2] = new_C2
        Q[o2] = new_O2
        # Exocyclic oxygens: shift by the parent-C delta, then project onto
        # the ring plane so the C=O π system stays coplanar with the ring.
        # C-O bond length is preserved.
        delta_C1 = new_C1 - old_C1
        delta_C2 = new_C2 - old_C2
        for exo in ring["O1_exo"]:
            shifted = Q[exo] + delta_C1
            Q[exo] = _project_exocyclic_O(shifted, Q[c1], cen, normal)
        for exo in ring["O2_exo"]:
            shifted = Q[exo] + delta_C2
            Q[exo] = _project_exocyclic_O(shifted, Q[c2], cen, normal)
        ring_pts_after = np.array([M_pos, Q[o1], Q[c1], Q[c2], Q[o2]])
        dev_after = _max_plane_dev(ring_pts_after)
        bite_after = _bite(Q[o1] - M_pos, Q[o2] - M_pos)
        applied.append({
            "ring": ring,
            "plane_dev_before": float(dev_before),
            "plane_dev_after": float(dev_after),
            "bite_before": float(bite_before),
            "bite_after": float(bite_after),
        })
    return Q, applied


def _bite(u: np.ndarray, v: np.ndarray) -> float:
    """Angle O-M-O (deg)."""
    nu = float(np.linalg.norm(u))
    nv = float(np.linalg.norm(v))
    if nu < 1e-9 or nv < 1e-9:
        return 0.0
    cos = float(np.dot(u, v) / (nu * nv))
    cos = max(-1.0, min(1.0, cos))
    return float(np.degrees(np.arccos(cos)))


def _project_exocyclic_O(p_after_shift: np.ndarray, p_parent_C: np.ndarray,
                         plane_centroid: np.ndarray,
                         plane_normal: np.ndarray) -> np.ndarray:
    """Place an exocyclic O so that:
      - C-O bond length is preserved (same as before, modulo the C shift),
      - O lies in the ring plane,
      - O sits on the side of C OPPOSITE to the ring-side O (so the C=O π
        system stays planar with the ring, classic carboxylate geometry).

    ``p_after_shift`` is the exocyclic O position translated by the same
    delta as its parent C (so the C-O bond vector is preserved from the
    free embed).  We then project into the ring plane to remove any
    out-of-plane component, and rescale so |C-O| matches the input length.
    """
    v = p_after_shift - p_parent_C
    L = float(np.linalg.norm(v))
    if L < 1e-9:
        return p_after_shift
    # Strip the out-of-plane component
    v_proj = v - plane_normal * float(np.dot(v, plane_normal))
    nv = float(np.linalg.norm(v_proj))
    if nv < 1e-9:
        # The original C-O direction was perpendicular to the new plane;
        # pick any in-plane direction that points away from the metal.
        # Fall back to leaving the atom alone.
        return p_after_shift
    return p_parent_C + v_proj / nv * L
