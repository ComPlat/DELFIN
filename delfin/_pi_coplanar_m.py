"""Iter-34 (2026-06-19) — coordinated planar π-donor co-planar-M orienter.

ROOT-CAUSE (eye-flagged ABIZIW).  A coordinated, internally-planar aromatic /
conjugated π-donor that binds the metal through an **in-plane sp2 σ lone pair**
(a pyridine-type ring N, an amidinate N, …) must lie so that its ring
mean-plane CONTAINS the metal: the donor→M vector lies IN the ring plane and
the donor's in-plane lone-pair direction points straight at M.

``DELFIN_FFFREE_AROM_PLANARIZE`` (Iter-33) makes such a ring *internally* flat
but does NOT rotate the coordinated ring plane THROUGH the metal — so the ring
sits tilted, the metal out of its plane (ABIZIW's coordinated C5N donor rings:
M out-of-plane 0.9-1.4 Å, donor→M to ring-plane angle 22-47°).  This module
closes that gap.

WHAT IT DOES (universal, geometry/graph-only — no SMILES / refcode / element
specialisation):

  1. Detect every internally-planar 5/6-membered aromatic/conjugated ring
     (C/N/O/S, aromatic-range mean bond length) — same proven detector as
     Iter-24/33.
  2. For each such ring, decide from GEOMETRY whether it is an **in-plane σ
     donor** (one ring atom bonded to M, the donor's outward in-plane radial /
     lone-pair points roughly toward M, M near the ring plane direction) or a
     **hapto / η π-face donor** (M sits over the ring face — its
     foot-of-perpendicular lands near the centroid, several ring atoms roughly
     equidistant from M).  HAPTO RINGS ARE EXCLUDED (they bind perpendicular;
     their plane must stay ⟂ to M→centroid, untouched).
  3. For an in-plane σ donor, rotate the RIGID ring (all ring-system atoms +
     ring-attached H + first non-ring substituents) about the FIXED donor atom
     so the donor's body-fixed in-plane lone-pair direction aligns with the
     (fixed) donor→M direction.  Donor and metal are NEVER moved, so the M–D
     bond length is preserved exactly.  Aligning the in-plane lone-pair with
     donor→M automatically brings M into the ring mean-plane (M out-of-plane
     → 0, angle → 0).

NEVER-WORSE / SAFETY (deterministic, never raises, default-OFF byte-id):
  * Donor + metal frozen — rotation is about the donor point only, M–D intact
    by construction (hard ``_md_invariant`` rollback as belt-and-braces).
  * Per-ring never-worse on M-out-of-plane: keep the rotation only if M's
    out-of-plane distance strictly decreased.
  * Inter-ligand clash never-worse: keep the rotation only if the minimum
    distance between any moved atom and any non-moved atom did not get worse
    (below a hard floor).  Otherwise skip this ring's rotation.
  * Frame-level + finiteness guards; any exception returns the input unchanged.

Opt-in via ``DELFIN_FFFREE_PI_COPLANAR_M=1``; default-OFF, byte-identical to
the build commit when unset.  Spec = this docstring.
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
    _COV_RADII,
    _MD_INVARIANT_TOL,
)
from delfin._arom_planarize import _detect_aromatic_rings, _fuse_components

# Chemically-coordinating donor elements (in-plane σ lone-pair donors): an
# atom of one of these elements that lies within M-bond distance is treated as a
# GENUINE M-D donor whose bond must be preserved.  A ring CARBON that is merely
# within M-bond distance because the ring is tilted (an artifact, e.g. a biaryl
# linker carbon dragged toward M) is NOT a genuine donor — it moves rigidly with
# its ring and its spurious M-contact is not protected.  (Genuine organometallic
# M-C — carbene / carbanion / CO — is still protected via the closest-contact
# rule in ``_genuine_md_atoms``: the carbon nearest M within its fragment.)
_DONOR_ELEMENTS = {"N", "O", "P", "S", "As", "Se", "F", "Cl", "Br", "I", "B"}
# A ring atom is treated as bonded to the metal if within this · Σr_cov.
_M_BOND_FACTOR: float = 1.30
# An in-plane σ donor must have M within this distance to count as coordinated.
_M_COORD_DIST: float = 3.20
# A planar-ring's max out-of-plane must be below this to be a real flat π-ring.
_RING_FLAT_TOL: float = 0.30
# Hapto discriminator: if M's foot-of-perpendicular onto the ring plane lands
# closer than this (Å) to the ring centroid, the metal sits over the ring face
# (η π-face donor) → EXCLUDE.  An in-plane σ donor's foot lands near the donor,
# i.e. ~1 ring radius (≈1.3-1.4 Å) from the centroid, far above this floor.
_HAPTO_FOOT_FRAC: float = 0.55
# Hapto discriminator: number of ring atoms within M-bonding distance.  An
# in-plane σ donor coordinates through ONE atom; an η ring coordinates through
# several.  A tilted (defective) σ ring may show 2-3 spuriously-close atoms, so
# this is combined with the foot-of-perpendicular test (face-on is the real
# discriminator), not used alone.
_HAPTO_MIN_COORD: int = 3
# Already-good gate: if M is already within this of the ring plane, leave it.
_ALREADY_GOOD: float = 0.05
# Hard inter-ligand clash floor: a moved atom may not come closer than this ·
# Σr_cov to any non-moved heavy/H atom (else the rotation is rejected).
_CLASH_FACTOR: float = 0.75
# Contacts at or above this ·Σr_cov ratio are comfortably non-clashing — a
# rotation that keeps the min contact above this is allowed even if the ratio
# dips slightly (a sub-vdW shuffle between far-apart atoms is not a clash).
# Only contacts that fall BELOW this AND get tighter trip the never-worse gate.
_CLASH_SAFE: float = 1.00
_SUBST_BOND_MAX_FACTOR: float = 1.30


def _ring_plane(pts: np.ndarray, ring: Sequence[int]) -> Tuple[np.ndarray, np.ndarray, float]:
    """Return (centroid, unit_normal, max_out_of_plane) of the SVD best-fit
    plane through the ring atoms."""
    P = np.array([pts[a] for a in ring])
    c = P.mean(axis=0)
    cen = P - c
    _, _, vh = np.linalg.svd(cen, full_matrices=False)
    nrm = vh[-1]
    nn = float(np.linalg.norm(nrm))
    if nn < 1e-12:
        return c, np.array([0.0, 0.0, 1.0]), 0.0
    nrm = nrm / nn
    oop = float(np.max(np.abs(cen @ nrm)))
    return c, nrm, oop


def _rotation_aligning(u: np.ndarray, v: np.ndarray) -> np.ndarray:
    """Return the 3×3 rotation matrix that rotates unit vector ``u`` onto unit
    vector ``v`` (shortest-arc / Rodrigues).  Deterministic; handles the
    parallel and antiparallel degenerate cases."""
    u = u / (np.linalg.norm(u) + 1e-15)
    v = v / (np.linalg.norm(v) + 1e-15)
    c = float(np.dot(u, v))
    if c > 1.0 - 1e-12:
        return np.eye(3)
    axis = np.cross(u, v)
    s = float(np.linalg.norm(axis))
    if s < 1e-12:
        # antiparallel: rotate 180° about any axis ⟂ u (deterministic pick)
        # choose the most-orthogonal cardinal axis
        cand = np.eye(3)[int(np.argmin(np.abs(u)))]
        axis = np.cross(u, cand)
        axis = axis / (np.linalg.norm(axis) + 1e-15)
        # 180° rotation about axis
        K = np.array([[0, -axis[2], axis[1]],
                      [axis[2], 0, -axis[0]],
                      [-axis[1], axis[0], 0]])
        return np.eye(3) + 2.0 * (K @ K)
    axis = axis / s
    K = np.array([[0, -axis[2], axis[1]],
                  [axis[2], 0, -axis[0]],
                  [-axis[1], axis[0], 0]])
    return np.eye(3) + s * K + (1.0 - c) * (K @ K)


def _in_plane_lone_pair(
    pts: np.ndarray,
    donor: int,
    ring_nbrs: Sequence[int],
    centroid: np.ndarray,
    normal: np.ndarray,
) -> np.ndarray:
    """In-plane sp2 lone-pair direction at the donor: outward bisector of the
    donor→ring-neighbour bonds, projected into the ring plane and normalised.
    Falls back to the outward radial (centroid→donor) if the bisector is
    degenerate."""
    bis = np.zeros(3)
    for j in ring_nbrs:
        v = pts[j] - pts[donor]
        nv = float(np.linalg.norm(v))
        if nv > 1e-9:
            bis += v / nv
    lp = -bis  # away from the ring neighbours = outward sp2 lone pair
    lp = lp - float(np.dot(lp, normal)) * normal
    nlp = float(np.linalg.norm(lp))
    if nlp < 1e-6:
        rad = pts[donor] - centroid
        rad = rad - float(np.dot(rad, normal)) * normal
        nrad = float(np.linalg.norm(rad))
        if nrad < 1e-6:
            return np.zeros(3)
        return rad / nrad
    return lp / nlp


def _linked_systems(
    rings: List[Tuple[int, ...]],
    nbrs: List[List[int]],
    syms: List[str],
) -> List[Set[int]]:
    """Group aromatic rings into maximal RIGID systems: rings are merged if they
    share an atom (fused) OR are joined by a single inter-ring heavy bond (a
    biaryl / linker bond — a conjugated planar polypyridyl like terpyridine is
    one rigid coplanar unit).  A per-ring rotation would tear such a linker, so
    the whole linked system must rotate together (or not at all).  Returns one
    atom-set per maximal system."""
    ring_sets = [set(r) for r in rings]
    n = len(ring_sets)
    parent = list(range(n))

    def find(x: int) -> int:
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    def union(a: int, b: int) -> None:
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[rb] = ra

    for i in range(n):
        for j in range(i + 1, n):
            if ring_sets[i] & ring_sets[j]:
                union(i, j)
                continue
            # biaryl / linker bond: an atom of ring i bonded to an atom of ring j
            linked = False
            for a in ring_sets[i]:
                for b in nbrs[a]:
                    if b in ring_sets[j]:
                        linked = True
                        break
                if linked:
                    break
            if linked:
                union(i, j)
    groups: dict = {}
    for i in range(n):
        groups.setdefault(find(i), set()).update(ring_sets[i])
    return list(groups.values())


def _genuine_md_atoms(
    pts: np.ndarray,
    syms: List[str],
    metals: Sequence[int],
    nbrs: List[List[int]],
) -> Set[int]:
    """Donor atoms whose M-D bond must be preserved.

    A genuine donor is an atom within M-bond distance of a metal that is EITHER
    a heteroatom/halide donor element (``_DONOR_ELEMENTS``) OR a carbon that is
    the closest M-contact within its own bonded ligand fragment (a genuine
    organometallic M-C: carbene / carbanion / CO).  A ring carbon that is
    M-close only because its ring is tilted — while a heteroatom of the SAME
    ligand sits even closer to M — is excluded: it is an artifact and must be
    free to move rigidly with its ring so the ring can be oriented through M.

    Distinguishing genuine M-C from a tilt artifact is done per metal: among the
    carbons within M-bond distance of a metal, a carbon is a genuine M-C donor
    only if NO heteroatom donor of the SAME connected (heavy) fragment is
    closer to that metal.  (For a real σ-aryl/carbene/CO the carbon IS the
    closest contact; for a tilted N-heterocycle the ring N is closer.)"""
    md: Set[int] = set()
    n = len(syms)
    # connected-fragment id over the non-metal heavy+H graph (so we can ask
    # "is there a closer heteroatom donor in the SAME ligand fragment?")
    frag = [-1] * n
    fid = 0
    for s in range(n):
        if frag[s] != -1 or _is_metal_sym(syms[s]):
            continue
        stack = [s]
        frag[s] = fid
        while stack:
            c = stack.pop()
            for j in nbrs[c]:
                if _is_metal_sym(syms[j]) or frag[j] != -1:
                    continue
                frag[j] = fid
                stack.append(j)
        fid += 1
    for m in metals:
        mr = _COV_RADII.get(syms[m], 1.5)
        # all heavy atoms within M-bond distance, with distances
        contacts = []
        for j in range(n):
            if j == m or _is_metal_sym(syms[j]) or syms[j] == "H":
                continue
            d = float(np.linalg.norm(pts[j] - pts[m]))
            if d <= _M_BOND_FACTOR * (mr + _COV_RADII.get(syms[j], 1.5)):
                contacts.append((d, j))
        for d, j in contacts:
            if syms[j] in _DONOR_ELEMENTS:
                md.add(j)
                continue
            if syms[j] != "C":
                # other element within range (rare) — protect conservatively
                md.add(j)
                continue
            # carbon: genuine M-C only if no heteroatom of the SAME fragment is
            # closer to THIS metal (else this C is a tilt artifact).
            closer_hetero = any(
                syms[k] in _DONOR_ELEMENTS
                and frag[k] == frag[j]
                and dk < d
                for dk, k in contacts
            )
            if not closer_hetero:
                md.add(j)
    return md


def _collect_rigid_atoms(
    ring_atoms: Set[int],
    nbrs: List[List[int]],
    syms: List[str],
    pts: np.ndarray,
    donor: int,
    all_ring_atoms: Set[int],
    frozen: Set[int],
) -> Set[int]:
    """Atoms that move rigidly with the ring: the (fused) ring-system atoms
    plus their first ring-attached H / non-ring substituent atoms.  The donor,
    every metal and every other frozen (M-D-protected) atom are EXCLUDED so
    they are never moved.  Substituents that belong to ANOTHER aromatic ring
    system are excluded (they are oriented by their own pass)."""
    rigid: Set[int] = set(ring_atoms)
    for ra in ring_atoms:
        for j in nbrs[ra]:
            if j in ring_atoms:
                continue
            if _is_metal_sym(syms[j]):
                continue
            if j == donor:
                continue
            if j in frozen:
                continue
            if j in all_ring_atoms:
                continue
            # only attach FIRST-shell substituents (bonded within covalent range)
            ri = _COV_RADII.get(syms[ra], 1.5)
            rj = _COV_RADII.get(syms[j], 1.5)
            d = float(np.linalg.norm(pts[ra] - pts[j]))
            if d <= _SUBST_BOND_MAX_FACTOR * (ri + rj):
                rigid.add(j)
    rigid.discard(donor)
    # never move a frozen / metal atom even if it slipped in
    rigid = {a for a in rigid if a not in frozen and not _is_metal_sym(syms[a])}
    return rigid


def _min_pair_dist(
    pts: np.ndarray,
    syms: List[str],
    moved: Set[int],
    others: Sequence[int],
    nbrs: List[List[int]],
) -> float:
    """Minimum (distance / Σr_cov) ratio between any moved atom and any non-moved
    atom that is NOT bonded to it (bonded pairs excluded so a covalent contact
    is not mistaken for a clash).  Lower = closer contact."""
    best = float("inf")
    moved_l = sorted(moved)
    for i in moved_l:
        ri = _COV_RADII.get(syms[i], 1.5)
        nbr_i = nbrs[i]
        for j in others:
            if j in moved:
                continue
            if j in nbr_i:
                continue
            d = float(np.linalg.norm(pts[i] - pts[j]))
            ratio = d / (ri + _COV_RADII.get(syms[j], 1.5) + 1e-12)
            if ratio < best:
                best = ratio
    return best


def _boundary_bond_broken(
    pre: np.ndarray,
    post: np.ndarray,
    moved: Set[int],
    nbrs: List[List[int]],
    tol: float = _MD_INVARIANT_TOL,
) -> bool:
    """Return True if any covalent bond that crosses the moved/non-moved
    boundary changed length by more than ``tol``.  This catches a rigid ring
    rotation that would tear an inter-ring / linker bond of a FUSED or LINKED
    polydentate chelate (e.g. a terpyridine biaryl C-C): only the donor pivot
    is shared, so rotating one ring drags its linker carbon away from the
    (stationary) carbon of the adjacent ring.  Such rotations are rejected
    (never-worse) — a rigid planar polydentate chelate's M-out-of-plane is fixed
    by its bite geometry, not freely adjustable by per-ring rotation."""
    for i in moved:
        for j in nbrs[i]:
            if j in moved:
                continue  # both endpoints moved rigidly — length preserved
            pre_d = float(np.linalg.norm(pre[i] - pre[j]))
            post_d = float(np.linalg.norm(post[i] - post[j]))
            if abs(post_d - pre_d) > tol:
                return True
    return False


def correct_xyz(xyz: str) -> str:
    """Orient each coordinated, internally-planar aromatic/conjugated π-DONOR
    SYSTEM so its mean-plane contains the metal: rotate the RIGID system about
    its single coordinating donor so the donor's in-plane sp2 lone-pair points
    at M (M brought into the system plane).  Donor + metal are FROZEN, so the
    M–D bond length is preserved exactly.

    SCOPE (universal, geometry/graph-only):
      * The rigid body is the maximal FUSED-AND-LINKED aromatic system (rings
        sharing an atom OR joined by a biaryl/linker bond — a terpyridine is one
        rigid coplanar unit), plus ring-attached H / first substituents.
      * Acted on only when the system has EXACTLY ONE genuine M-D donor (a
        heteroatom/halide donor, or a genuine organometallic carbon).  A single
        rigid rotation about that one frozen donor can bring M into the plane
        without disturbing any other M-D bond.  A multi-donor planar chelate
        (e.g. a chelating terpyridine) is PINNED by its several frozen M-D bonds
        — its M-out-of-plane is fixed by the bite/placement, not adjustable by a
        rigid rotation — so it is left untouched (any rotation would move the
        other donors and is rejected by the hard M-D invariant).
      * HAPTO / η π-face donors (M over the ring face — foot-of-perpendicular
        near the centroid) are EXCLUDED: they bind perpendicular and must stay.

    NEVER-WORSE / SAFETY: per-system never-worse on M-out-of-plane; linker-bond
    preservation; inter-ligand clash floor + already-tight never-worse; hard M-D
    invariant rollback; finiteness.  Bit-exact when no qualifying system is
    present.  Any error → input unchanged."""
    try:
        syms, pts, lines = _parse_xyz(xyz)
    except Exception:
        return xyz
    if len(syms) < 5:
        return xyz
    metals = [i for i in range(len(syms)) if _is_metal_sym(syms[i])]
    if not metals:
        return xyz
    try:
        nbrs = _build_geometric_adjacency(syms, pts)
        rings = _detect_aromatic_rings(syms, pts, nbrs)
    except Exception:
        return xyz
    if not rings:
        return xyz

    pre_md = _snapshot_md_bonds(pts, syms)
    # GENUINE donors (heteroatom/halide, or genuine organometallic C) — used as
    # anchors/pivots.  A ring carbon merely M-close because its ring is tilted is
    # NOT genuine and is free to move rigidly with its ring.
    genuine = _genuine_md_atoms(pts, syms, metals, nbrs)
    all_ring_atoms: Set[int] = {a for r in rings for a in r}
    # Maximal fused-AND-linked rigid aromatic systems (terpyridine = one unit).
    systems = _linked_systems(list(rings), nbrs, syms)

    work = pts.copy()
    any_moved = False

    for sys_atoms in systems:
        # system must be internally flat (a real planar π-system)
        sys_list = sorted(sys_atoms)
        _, sys_nrm0, sys_oop = _ring_plane(work, sys_list)
        if sys_oop > _RING_FLAT_TOL:
            continue

        for m in metals:
            mr = _COV_RADII.get(syms[m], 1.5)
            # genuine donors of THIS system to THIS metal
            sys_donors = [
                a for a in sys_atoms
                if a in genuine
                and float(np.linalg.norm(work[a] - work[m]))
                <= _M_BOND_FACTOR * (mr + _COV_RADII.get(syms[a], 1.5))
            ]
            if len(sys_donors) != 1:
                # 0 donors → not coordinated to this metal; ≥2 → pinned chelate.
                continue
            donor = sys_donors[0]
            dM = work[m] - work[donor]
            ndM = float(np.linalg.norm(dM))
            if ndM < 1e-6 or ndM > _M_COORD_DIST:
                continue
            dM_u = dM / ndM

            cen, nrm, _ = _ring_plane(work, sys_list)
            m_oop = abs(float(np.dot(work[m] - cen, nrm)))

            # ---- hapto / η π-face discriminator (face-on => EXCLUDE) --------
            foot = work[m] - float(np.dot(work[m] - cen, nrm)) * nrm
            foot_to_cen = float(np.linalg.norm(foot - cen))
            sys_radius = float(np.mean(
                [np.linalg.norm(work[a] - cen) for a in sys_atoms]
            ))
            n_close = sum(
                1 for a in sys_atoms
                if float(np.linalg.norm(work[a] - work[m]))
                <= _M_BOND_FACTOR * (mr + _COV_RADII.get(syms[a], 1.5))
            )
            face_on = (
                foot_to_cen < _HAPTO_FOOT_FRAC * (sys_radius + 1e-12)
                and n_close >= _HAPTO_MIN_COORD
            )
            if face_on:
                continue  # η π-face donor — binds perpendicular, leave it

            if m_oop <= _ALREADY_GOOD:
                continue  # already co-planar

            ring_nbrs = [
                j for j in nbrs[donor] if j in sys_atoms and j != donor
            ]
            if len(ring_nbrs) < 1:
                continue
            lp = _in_plane_lone_pair(work, donor, ring_nbrs, cen, nrm)
            if float(np.linalg.norm(lp)) < 1e-6:
                continue

            rigid = _collect_rigid_atoms(
                sys_atoms, nbrs, syms, work, donor, all_ring_atoms, genuine,
            )
            if not rigid:
                continue

            R = _rotation_aligning(lp, dM_u)
            pivot = work[donor].copy()
            trial = work.copy()
            for a in rigid:
                trial[a] = pivot + R @ (work[a] - pivot)

            # never tear a bond crossing the moved/non-moved boundary (a linker
            # to a non-moved fragment, or a genuine donor held fixed).
            if _boundary_bond_broken(work, trial, rigid, nbrs):
                continue

            # per-system never-worse: M out-of-plane must strictly improve
            cen2, nrm2, _ = _ring_plane(trial, sys_list)
            m_oop_after = abs(float(np.dot(trial[m] - cen2, nrm2)))
            if m_oop_after >= m_oop - 1e-3:
                continue

            # inter-ligand clash gate: reject only if the rotation drives a
            # contact below the hard floor, OR makes an ALREADY-tight contact
            # (below the comfortably-safe ratio) strictly tighter.
            others = [i for i in range(len(syms)) if i not in rigid]
            after_min = _min_pair_dist(trial, syms, rigid, others, nbrs)
            if after_min < _CLASH_FACTOR:
                continue  # hard clash floor
            if after_min < _CLASH_SAFE:
                before_min = _min_pair_dist(work, syms, rigid, others, nbrs)
                if after_min < before_min - 1e-3:
                    continue  # worsened an already-tight contact

            # hard M-D invariant (every genuine M-D / M-H bond preserved)
            if _md_invariant_violated(pre_md, trial, tol=_MD_INVARIANT_TOL):
                continue
            if not np.all(np.isfinite(trial)):
                continue

            work = trial
            any_moved = True

    if not any_moved:
        return xyz
    if not np.all(np.isfinite(work)):
        return xyz
    if _md_invariant_violated(pre_md, work, tol=_MD_INVARIANT_TOL):
        return xyz
    try:
        return _format_xyz(lines, syms, work)
    except Exception:
        return xyz


def correct_results(mol, results):
    """Apply ``correct_xyz`` to each (xyz, label) tuple.  ``mol`` accepted for
    signature parity with the other baustein correctors (geometry-only).
    Fail-safe: any per-entry exception returns the original tuple unchanged."""
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
