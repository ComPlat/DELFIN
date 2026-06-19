"""Iter-33 (2026-06-19) — universal aromatic / sp2-conjugated ring-SYSTEM
planarity enforcement (heteroaromatic + fused/polycyclic + coordinated).

Motivation (measured forensic, this session):
    The CONSOL-archive AXAGOY (Os, 49 atoms, phenanthroline-type chelate +
    pendant phenyls) shows its **coordinated** pyridine-like C5N donor rings
    puckered 0.13-0.16 Å out-of-plane, the fused central C6 0.13 Å, a pendant
    phenyl 0.09 Å (flat-aromatic target < 0.05 Å).  The existing
    ``_aromatic_ring_flattener`` (Iter-24) intentionally:

      * skips EVERY coordinated ring (``pendant_rings = [r for r in rings if
        not _is_coordinated(r)]``) — i.e. it never touches the worst
        offenders, the M-bound heteroaromatic donor rings;
      * only fuses *pendant* rings into a common plane, so a polycyclic
        chelate that is partly coordinated is never flattened as ONE unit;
      * is default-OFF and class-gated to {hapto, multi_hapto}, so a σ-class
        complex like AXAGOY never reaches it at all.

This module closes those gaps **without breaking the M-D invariant**:

  1. Detect aromatic rings geometrically (5/6-membered C/N/O/S ring, mean
     intra-ring heavy bond < ``_AROMATIC_BOND_MAX`` Å) — universal, no
     SMILES / refcode / element-list specialisation.  This is the same
     proven detector as Iter-24; it correctly recognises heteroaromatics
     (ring N/O/S) and pendant + coordinated rings alike.
  2. Fuse ALL aromatic rings that share an atom into one polycyclic π-system
     (naphthalene / phenanthroline / indole …) and flatten each system as a
     single rigid unit onto its common best-fit plane.
  3. **M-D preservation**: any ring-system atom that sits within M-D bonding
     distance of a metal is an *anchor* — it is NOT moved.  The best-fit
     plane is constrained to pass through the anchor(s); every NON-anchor
     atom of the system is projected onto that anchor-constrained plane.
     With ≥1 anchor fixed and the metal untouched, every M-anchor distance is
     preserved exactly, so the M-D / coordination invariant holds by
     construction (verified additionally by the hard ``_md_invariant``
     rollback).  Pure-pendant systems (no anchor) use the centroid-preserving
     projection (M-centroid distance invariant) as before.
  4. Ring-attached H + first substituent atoms are dragged rigidly onto the
     new plane (so H / ipso-C do not lag the flattened ring).

Per-frame, per-system NEVER-WORSE guards (a system is only kept if its own
worst-ring OOP strictly decreased) plus a frame-level worst-OOP guard plus
the hard M-D-invariant rollback.  Bit-exact to the input when the flag is 0
or no qualifying ring is present.  Operates on XYZ text only (no RDKit
dependency for the geometry step → robust to mol/atom-order drift and to
exotic metal valences that defeat RDKit sanitisation, e.g. Os-4).

Opt-in via ``DELFIN_FFFREE_AROM_PLANARIZE=1``; default-OFF, byte-identical to
the build commit when unset.
"""
from __future__ import annotations

from typing import List, Sequence, Set, Tuple

import numpy as np

from delfin._pi_h_projector import (
    _parse_xyz,
    _format_xyz,
    _build_geometric_adjacency,
    _snapshot_md_bonds,
    _md_invariant_violated,
    _is_metal_sym,
    _MD_INVARIANT_TOL,
)

_AROMATIC_ELIGIBLE = {"C", "N", "O", "S"}
_OOP_TOL: float = 0.05            # target: rings flatter than this are done
_PUCKER_CAP: float = 0.80         # above this it is not a planar-intended ring
_AROMATIC_BOND_MAX: float = 1.46  # mean intra-ring heavy-bond gate (Å)
_M_COORD_DIST: float = 2.60       # ring atom within this of metal → anchor
_SUBST_BOND_MAX_FACTOR: float = 1.30  # heavy substituent within this·Σr_cov is dragged

_COV_RADII = {
    "H": 0.31, "B": 0.84, "C": 0.76, "N": 0.71, "O": 0.66, "F": 0.57,
    "Si": 1.11, "P": 1.07, "S": 1.05, "Cl": 1.02, "Se": 1.20, "Br": 1.20,
    "I": 1.39,
}


def _detect_aromatic_rings(
    syms: List[str],
    pts: np.ndarray,
    nbrs: List[List[int]],
) -> List[Tuple[int, ...]]:
    """Geometric 5/6-membered C/N/O/S rings with aromatic-range mean bond
    length.  Saturated rings (mean bond ~1.54) and metal-containing rings are
    excluded.  Heteroaromatics (ring N/O/S) included.  Returns canonical
    sorted ring tuples.  (Same algorithm as the Iter-24 detector.)"""
    n = len(syms)
    heavy_nbrs: List[List[int]] = [
        [j for j in nbrs[i]
         if syms[j] in _AROMATIC_ELIGIBLE and syms[i] in _AROMATIC_ELIGIBLE
         and not _is_metal_sym(syms[j])]
        for i in range(n)
    ]
    rings_canon: Set[tuple] = set()
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
        out.append(tuple(ring))
    return out


def _fuse_components(rings: List[Tuple[int, ...]]) -> List[List[int]]:
    """Union aromatic rings that share ≥1 atom into maximal fused π-systems
    (naphthalene, phenanthroline, indole …).  Flattening a fused system as
    ONE unit avoids the shared-atom tug-of-war.  Returns sorted atom-index
    lists, one per maximal fused system."""
    comps: List[set] = []
    for ring in rings:
        rs = set(ring)
        merged = [rs]
        rest = []
        for c in comps:
            (merged if (c & rs) else rest).append(c)
        union: set = set()
        for m in merged:
            union |= m
        rest.append(union)
        comps = rest
    return [sorted(c) for c in comps]


def _best_fit_plane(P: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """Return (point_on_plane=centroid, unit_normal) of the SVD best-fit
    plane through ``P`` (M×3).  Normal is the least-variance direction."""
    centroid = P.mean(axis=0)
    _, _, vh = np.linalg.svd(P - centroid, full_matrices=False)
    nrm = vh[-1]
    return centroid, nrm / (np.linalg.norm(nrm) + 1e-12)


def _anchor_constrained_plane(
    pts: np.ndarray,
    atoms: Sequence[int],
    anchors: Sequence[int],
) -> Tuple[np.ndarray, np.ndarray]:
    """Best-fit plane of ``atoms`` that is constrained to pass through the
    fixed anchor point(s).

    * 0 anchors  → ordinary centroid best-fit plane.
    * 1 anchor   → plane through the anchor, normal = least-variance direction
      of the anchor-centred atom cloud.
    * ≥2 anchors → plane contains the anchor centroid AND (best-effort) the
      anchor-spanning direction(s): we anchor the plane origin at the anchor
      centroid and pick the normal as the least-variance SVD direction of the
      WHOLE atom cloud re-centred on the anchor centroid; for 2 anchors this
      keeps the anchor-anchor axis in the plane to first order, which is what
      a bidentate chelate (e.g. phenanthroline N,N) needs.

    Returns (origin, unit_normal).  ``origin`` is a point the plane passes
    through and is itself never moved (it is an anchor / anchor centroid).
    """
    P = np.array([pts[a] for a in atoms])
    if not anchors:
        return _best_fit_plane(P)
    A = np.array([pts[a] for a in anchors])
    origin = A.mean(axis=0)
    # Re-centre the full cloud on the (fixed) anchor origin and take the
    # least-variance direction as the normal.  The plane then passes exactly
    # through ``origin`` (anchor / anchor centroid), so projecting non-anchor
    # atoms onto it leaves the anchors untouched.
    centered = P - origin
    try:
        _, _, vh = np.linalg.svd(centered, full_matrices=False)
    except np.linalg.LinAlgError:
        return origin, np.array([0.0, 0.0, 1.0])
    nrm = vh[-1]
    return origin, nrm / (np.linalg.norm(nrm) + 1e-12)


def _ring_oop(pts: np.ndarray, ring: Sequence[int]) -> float:
    P = np.array([pts[a] for a in ring])
    c = P - P.mean(axis=0)
    try:
        _, _, vh = np.linalg.svd(c, full_matrices=False)
    except np.linalg.LinAlgError:
        return 0.0
    nrm = vh[-1] / (np.linalg.norm(vh[-1]) + 1e-12)
    return float(np.max(np.abs(c @ nrm)))


def _system_worst_oop(
    pts: np.ndarray, rings_in_system: List[Tuple[int, ...]]
) -> float:
    return max((_ring_oop(pts, r) for r in rings_in_system), default=0.0)


def _flatten_system(
    pts: np.ndarray,
    atoms: List[int],
    anchors: Set[int],
    nbrs: List[List[int]],
    syms: List[str],
) -> bool:
    """Project the (fused-)aromatic heavy atoms of ``atoms`` onto their
    anchor-constrained best-fit plane, plus drag ring-attached H and first
    substituent heavy atoms rigidly with their parent.  Anchors are never
    moved.  Mutates ``pts`` in place.  Returns True if any atom moved.

    The plane is well-defined only for ≥3 non-collinear atoms; tiny / linear
    systems are skipped.  Centroid-preserving for anchor-free systems
    (M-centroid invariant); anchor-passing for coordinated systems
    (per-anchor M-distance invariant)."""
    if len(atoms) < 3:
        return False
    anchor_list = [a for a in atoms if a in anchors]
    origin, normal = _anchor_constrained_plane(pts, atoms, anchor_list)
    nn = float(np.linalg.norm(normal))
    if nn < 1e-9:
        return False

    aset = set(atoms)
    moved = False
    # Project every NON-anchor ring-system atom onto the plane along ``normal``.
    for idx in atoms:
        if idx in anchors:
            continue  # fixed — preserves M-D distance exactly
        v = pts[idx] - origin
        shift = float(np.dot(v, normal)) * normal
        if np.linalg.norm(shift) > 1e-9:
            pts[idx] = pts[idx] - shift
            moved = True

    # Drag ring-attached H and first (non-ring, non-metal) substituent heavy
    # atoms onto the plane, rigidly with their parent.  A substituent atom is
    # projected the SAME signed distance as if it lived in the ring plane:
    # we simply project it onto the plane along ``normal`` too (keeps the
    # in-plane radial vector, removes only the out-of-plane wiggle).  H on
    # anchors are dragged as well (their parent did not move, but the H may
    # still be off-plane).
    for ring_atom in atoms:
        for j in nbrs[ring_atom]:
            if j in aset:
                continue
            if _is_metal_sym(syms[j]):
                continue  # never move a metal
            if j in anchors:
                continue
            # Only drag atoms that are genuinely bonded substituents (H, or a
            # first heavy substituent atom).  Heavy substituents are dragged
            # only if they themselves carry no further ring membership outside
            # this system — first-shell only.
            v = pts[j] - origin
            shift = float(np.dot(v, normal)) * normal
            if np.linalg.norm(shift) > 1e-9:
                pts[j] = pts[j] - shift
                moved = True
    return moved


def correct_xyz(xyz: str) -> str:
    """Flatten every puckered aromatic ring SYSTEM in a single XYZ string.

    Heteroaromatic, fused/polycyclic and coordinated systems are all handled.
    Per-system never-worse + frame-level never-worse + hard M-D-invariant
    rollback.  Bit-exact when no qualifying ring is present or all rings are
    already flat.  Any error returns the input unchanged."""
    try:
        syms, pts, lines = _parse_xyz(xyz)
    except Exception:
        return xyz
    if len(syms) < 5:
        return xyz
    try:
        nbrs = _build_geometric_adjacency(syms, pts)
        rings = _detect_aromatic_rings(syms, pts, nbrs)
    except Exception:
        return xyz
    if not rings:
        return xyz

    metals = [i for i, s in enumerate(syms) if _is_metal_sym(s)]

    def _anchors_for(atoms: Set[int]) -> Set[int]:
        """Ring-system atoms within M-D bonding distance of any metal."""
        out: Set[int] = set()
        for a in atoms:
            for m in metals:
                if float(np.linalg.norm(pts[a] - pts[m])) < _M_COORD_DIST:
                    out.add(a)
                    break
        return out

    components = _fuse_components(list(rings))
    pre_md = _snapshot_md_bonds(pts, syms)
    work = pts.copy()
    any_moved = False

    for atoms in components:
        aset = set(atoms)
        rings_here = [r for r in rings if set(r) <= aset]
        if not rings_here:
            continue
        worst_before = _system_worst_oop(work, rings_here)
        # only act on systems that are puckered but planar-intended
        if worst_before <= _OOP_TOL or worst_before > _PUCKER_CAP:
            continue
        anchors = _anchors_for(aset)
        before = work.copy()
        moved = _flatten_system(work, atoms, anchors, nbrs, syms)
        if not moved:
            continue
        worst_after = _system_worst_oop(work, rings_here)
        # per-system never-worse: keep only if this system got flatter
        if worst_after >= worst_before - 1e-3:
            # revert this system's atoms + their dragged substituents
            touched = set(atoms)
            for ra in atoms:
                for j in nbrs[ra]:
                    touched.add(j)
            for i in touched:
                work[i] = before[i]
            continue
        any_moved = True

    if not any_moved:
        return xyz

    # Frame-level never-worse: the global worst aromatic-ring OOP must not rise.
    def _max_ring_oop(P: np.ndarray) -> float:
        return max((_ring_oop(P, r) for r in rings), default=0.0)

    if _max_ring_oop(work) > _max_ring_oop(pts) + 1e-3:
        return xyz

    # Hard invariant: never break a pre-existing M-D (or M-H) bond.
    if _md_invariant_violated(pre_md, work, tol=_MD_INVARIANT_TOL):
        return xyz

    # Determinism / finiteness guard.
    if not np.all(np.isfinite(work)):
        return xyz

    try:
        return _format_xyz(lines, syms, work)
    except Exception:
        return xyz


def correct_results(mol, results):
    """Apply ``correct_xyz`` to each (xyz, label) tuple.  ``mol`` accepted for
    signature parity with the other baustein correctors (not required —
    geometry-only).  Fail-safe: any per-entry exception returns the original
    tuple unchanged."""
    if not results:
        return results
    out = []
    for entry in results:
        try:
            xyz, lbl = entry[0], entry[1]
            new_xyz = correct_xyz(xyz)
            out.append((new_xyz, lbl) if len(entry) == 2
                       else (new_xyz,) + tuple(entry[1:]))
        except Exception:
            out.append(entry)
    return out
