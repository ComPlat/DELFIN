"""Iter-24 (2026-05-20) — post-UFF aromatic-ring planarity enforcement.

Re-grounding forensic (project_iter24_pi_planarity_forensic) showed that in
hapto / multi_hapto structures, TRUE aromatic rings (M_coord chelate rings
excluded) pucker badly: 72-75 % of rings violate the 0.10 Å out-of-plane
tolerance, mean OOP 0.34-0.38 Å.  e6761e4's soft-donor port (Iter-19) only
covered sigma / multi_sigma so it never helped the hapto aromatic rings.

This module flattens puckered aromatic rings onto their SVD best-fit plane,
**preserving the ring centroid** (so the metal-ring distance / M-D invariant
is untouched) and dragging ring-attached H rigidly afterwards.  It operates
on XYZ text only (no RDKit dependency for the geometry step) so it is robust
to any mol/atom-order drift — aromatic rings are identified geometrically:

  * 5/6-membered ring of C/N/O/S (no metal, no H)
  * mean intra-ring heavy-bond length < ``_AROMATIC_BOND_MAX`` Å
    (aromatic/conjugated ~1.34-1.43; excludes saturated cyclohexane ~1.54)
  * current max OOP in (``_OOP_TOL``, ``_PUCKER_CAP``) — puckered but still
    recognisably planar-intended (avoids touching genuinely 3-D cages)

Per-frame rollback if the flatten breaks any pre-existing M-D bond by more
than the Iter-15 hard invariant (_MD_INVARIANT_TOL).  Bit-exact to input
when no qualifying ring is found.

Reuses the proven machinery from ``delfin._pi_h_projector``.
"""
from __future__ import annotations

from typing import List, Tuple

import numpy as np

from delfin._pi_h_projector import (
    _parse_xyz,
    _format_xyz,
    _build_geometric_adjacency,
    _snapshot_md_bonds,
    _md_invariant_violated,
    _is_metal_sym,
    project_ring_h_atoms,
    _MD_INVARIANT_TOL,
)

_AROMATIC_ELIGIBLE = {"C", "N", "O", "S"}
_OOP_TOL: float = 0.10          # below this the ring is already flat
_PUCKER_CAP: float = 0.80       # above this it is not a planar-intended ring
_AROMATIC_BOND_MAX: float = 1.46  # mean intra-ring heavy-bond length gate (Å)
_M_COORD_DIST: float = 2.60       # ring atom within this of a metal → coordinated (skip)


def _detect_aromatic_rings(
    syms: List[str],
    pts: np.ndarray,
    nbrs: List[List[int]],
) -> List[Tuple[int, ...]]:
    """Geometric 5/6-membered C/N/O/S rings with aromatic-range mean bond
    length.  Saturated rings (mean bond ~1.54) and metal-containing rings
    are excluded.  Returns canonical sorted ring tuples."""
    n = len(syms)
    heavy_nbrs: List[List[int]] = [
        [j for j in nbrs[i]
         if syms[j] in _AROMATIC_ELIGIBLE and syms[i] in _AROMATIC_ELIGIBLE
         and not _is_metal_sym(syms[j])]
        for i in range(n)
    ]
    rings_canon: set = set()
    for start in range(n):
        if syms[start] not in _AROMATIC_ELIGIBLE:
            continue
        stack = [(start, [start])]
        while stack:
            cur, path = stack.pop()
            if len(path) > 6:
                continue
            for nx in heavy_nbrs[cur]:
                if nx == path[0] and len(path) >= 5:
                    rings_canon.add(tuple(sorted(path)))
                    continue
                if nx in path:
                    continue
                if len(path) < 6:
                    stack.append((nx, path + [nx]))

    out: List[Tuple[int, ...]] = []
    for ring in rings_canon:
        ring = tuple(ring)
        # mean intra-ring heavy-bond length (consecutive bonded pairs only)
        rset = set(ring)
        bond_lens: List[float] = []
        for i in ring:
            for j in heavy_nbrs[i]:
                if j in rset and j > i:
                    bond_lens.append(float(np.linalg.norm(pts[i] - pts[j])))
        if not bond_lens:
            continue
        if (sum(bond_lens) / len(bond_lens)) >= _AROMATIC_BOND_MAX:
            continue  # saturated ring — leave its (correct) pucker alone
        out.append(ring)
    return out


def _flatten_ring(pts: np.ndarray, ring: Tuple[int, ...]) -> float:
    """Project ring heavy atoms onto their SVD best-fit plane through the
    centroid (centroid preserved → M-ring distance invariant).  Mutates
    ``pts`` in place.  Returns the pre-flatten max OOP (0.0 if not eligible)."""
    ring_pts = np.array([pts[i] for i in ring])
    centroid = ring_pts.mean(axis=0)
    centered = ring_pts - centroid
    try:
        _, _, vh = np.linalg.svd(centered, full_matrices=False)
    except np.linalg.LinAlgError:
        return 0.0
    normal_raw = vh[-1]
    nn = float(np.linalg.norm(normal_raw))
    if nn < 1e-9:
        return 0.0
    normal = normal_raw / nn
    oop = float(np.max(np.abs(centered @ normal)))
    if oop <= _OOP_TOL or oop > _PUCKER_CAP:
        return 0.0
    for idx in ring:
        v = pts[idx] - centroid
        pts[idx] = pts[idx] - float(np.dot(v, normal)) * normal
    return oop


def correct_results(
    mol,
    results: List[Tuple[str, str]],
) -> List[Tuple[str, str]]:
    """Flatten puckered aromatic rings in every frame of ``results``.

    ``mol`` is accepted for signature parity with the other baustein
    correctors but is NOT required (geometry-only).  Per-frame M-D
    invariant rollback; bit-exact when no qualifying ring is present.
    """
    if not results:
        return results
    out: List[Tuple[str, str]] = []
    for item in results:
        try:
            xyz, label = item
        except Exception:
            out.append(item)
            continue
        try:
            syms, pts, lines = _parse_xyz(xyz)
        except Exception:
            out.append((xyz, label))
            continue
        if len(syms) < 5:
            out.append((xyz, label))
            continue
        try:
            nbrs = _build_geometric_adjacency(syms, pts)
            rings = _detect_aromatic_rings(syms, pts, nbrs)
        except Exception:
            out.append((xyz, label))
            continue
        if not rings:
            out.append((xyz, label))
            continue
        # Identify metals; only flatten PENDANT aromatic rings (no ring atom
        # within M-D bonding distance).  η-coordinated rings (Cp/arene) cannot
        # be flattened atom-wise without changing individual M-C distances and
        # tripping the M-D invariant (centroid-preservation only protects
        # M-centroid, not per-bond) — they need a centroid-based treatment,
        # deferred.  Pendant rings carry no M-D constraint → safe + cover the
        # majority of violations (forensic: pendant 3873 vs coordinated 1599).
        metals = [i for i, s in enumerate(syms) if _is_metal_sym(s)]

        def _is_coordinated(ring: Tuple[int, ...]) -> bool:
            for ai in ring:
                for mi in metals:
                    if float(np.linalg.norm(pts[ai] - pts[mi])) < _M_COORD_DIST:
                        return True
            return False

        pre_md = _snapshot_md_bonds(pts, syms)
        work = pts.copy()
        moved = False
        for ring in rings:
            if _is_coordinated(ring):
                continue
            if _flatten_ring(work, ring) > 0.0:
                moved = True
        if not moved:
            out.append((xyz, label))
            continue
        # Drag ring-attached H onto the now-flat ring planes (reuse projector).
        try:
            project_ring_h_atoms(syms, work)
        except Exception:
            pass
        # Hard invariant: never break a pre-existing M-D bond.
        if _md_invariant_violated(pre_md, work, tol=_MD_INVARIANT_TOL):
            out.append((xyz, label))  # rollback whole frame
            continue
        try:
            out.append((_format_xyz(lines, syms, work), label))
        except Exception:
            out.append((xyz, label))
    return out
