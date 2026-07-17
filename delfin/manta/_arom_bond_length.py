"""delfin.manta._arom_bond_length — resonance-aware aromatic bond-length corrector.

MANTA seats ligand skeletons from distance-geometry / covalent priors, which for
AROMATIC rings leaves the ring bonds (C–C / C–N / …) drifting toward the SINGLE-
bond covalent length and/or ALTERNATING (Kekulé-localised: ~1.34 / ~1.50) rather
than sitting at the DELOCALISED, mesomerism-equalised length.  Measured against
the eye's organic-bond-length axis this is the single largest systematic organic-
geometry defect on the champion pool (aromatic ring bonds too long / bond-length
alternation too high).

Physical target — bond-length equalisation from resonance (Entartung):
    * benzene: all six C–C EQUAL at ~1.39 Å (full delocalisation), NOT 1.34/1.50.
    * heteroaromatic / fused π: every ring bond at its own delocalised length.
The delocalised bond order of a fully aromatic ring is the Hückel/Pauling benzene
value n ≈ 5/3, i.e. an interpolation fraction f = 2/3 between the single- and
double-bond covalent radii.  The per-bond target is therefore the Pyykkö
single↔double covalent-radius interpolation at f = 2/3:

    L(A–B) = r_f(A) + r_f(B),   r_f = r_single − f·(r_single − r_double)

giving C–C 1.393, C–N 1.333, C–O 1.289, C–S 1.671, … from FIRST PRINCIPLES
(Pyykkö radii) for ANY element pair — a Pauling-type L(order), no per-pattern
table, no SMILES specialisation.

Mechanism (per frame, geometry-only → robust to mol/XYZ atom-order drift, the
same reason _arom_planarize / _bond_decollapse are geometry-only — a mol-index
version detects 0 bonds because the emitted XYZ order need not match the mol):
  1. Perceive aromatic ring systems geometrically — planar 5/6-rings of C/N/O/S
     with an aromatic-band mean bond — reusing the _arom_planarize detector, so
     exactly the rings that pass the aromatic gate are touched; fused rings are
     unioned into one π-system.
  2. Reshape ONLY the ring atoms toward the per-bond delocalised target with a
     Position-Based-Dynamics distance solver, projected INTO the ring's best-fit
     plane, with any metal-coordinated ring atom ANCHORED (M–D invariant hard-
     preserved).  Ring angles move only TOWARD the regular-polygon ideal
     (equalisation) and planarity is held by the in-plane projection — LENGTHS
     ONLY, never a term that distorts angles or planarity.
  3. Every NON-ring atom (substituents, ring-H, the metal side) is RIGIDLY
     translated by the displacement of the ring atom it hangs off, so all
     substituent internal geometry (bonds / angles / planarity) is preserved
     byte-exactly; only the aromatic ring reshapes.
  4. Per-frame NEVER-WORSE rollback: the reshaped frame is kept only if the ring
     bond-length deviation from target STRICTLY drops AND no M–D bond breaks AND
     no new heavy–heavy clash appears; otherwise the ORIGINAL frame text is
     returned unchanged (byte-identical).  Already-equalised rings (max bond
     error below the activation tol) are left byte-identical.

Opt-in via ``DELFIN_FFFREE_AROM_BOND_LENGTH=1``; default-OFF → the dispatch never
calls this module, so the champion build is byte-identical.
"""
from __future__ import annotations

from typing import Dict, List, Sequence, Set, Tuple

import numpy as np

from delfin.manta._pi_h_projector import (
    _parse_xyz,
    _format_xyz,
    _build_geometric_adjacency,
    _snapshot_md_bonds,
    _md_invariant_violated,
    _is_metal_sym,
)
from delfin.manta._arom_planarize import (
    _detect_aromatic_rings,
    _fuse_components,
    _COV_RADII,
)
from delfin.manta.polyhedra import PYYKKO_SINGLE, PYYKKO_DOUBLE

# Hückel/Pauling delocalised bond order of a fully aromatic ring ≈ 5/3 → the
# single↔double covalent-radius interpolation fraction is 2/3 (benzene C–C 1.393).
_AROM_FRACTION: float = 2.0 / 3.0
# A ring atom within this of a metal is anchored (its M–D bond is the invariant);
# matches _arom_planarize._M_COORD_DIST.
_M_COORD_DIST: float = 2.60
_MAX_PASSES: int = 200            # PBD sweeps (small rings converge in < 30)
_CONVERGE_TOL: float = 1e-4       # Å — stop when the worst ring-bond error is below
_STEP_CLAMP: float = 0.15         # Å — max per-atom move per sweep (stability)
_ACTIVATE_TOL: float = 0.025      # Å — a ring system is only reshaped if its worst
                                  # ring-bond error exceeds this (already-equalised
                                  # rings stay byte-identical)
_CLASH_FACTOR: float = 0.80       # heavy–heavy overlap (clash) if d < factor·Σr_cov
_IMPROVE_EPS: float = 1e-4        # min ring-bond-deviation drop to accept a frame


# ── first-principles delocalised target length ────────────────────────────────
def _arom_radius(sym: str) -> float:
    """Delocalised (f = 2/3 single↔double) Pyykkö covalent radius for an aromatic
    ring atom.  Falls back to the single-bond covalent radius when the element is
    outside the Pyykkö double table → target = single sum → no shortening (safe)."""
    rs = PYYKKO_SINGLE.get(sym)
    if rs is None:
        return _COV_RADII.get(sym, 0.75)
    rd = PYYKKO_DOUBLE.get(sym)
    if rd is None:
        return rs
    return rs - _AROM_FRACTION * (rs - rd)


def _target_length(a: str, b: str) -> float:
    return _arom_radius(a) + _arom_radius(b)


# ── graph helpers ─────────────────────────────────────────────────────────────
def _ring_bonds_of_systems(
    systems: List[List[int]], nbrs: List[List[int]], syms: List[str]
) -> Tuple[List[Tuple[int, int]], Dict[int, int]]:
    """Return (ring_bond_pairs, atom→system-id).  A ring bond is a heavy–heavy
    edge whose two atoms lie in the SAME fused aromatic system (an inter-system
    biaryl edge is NOT a ring bond and is left untouched)."""
    sys_of: Dict[int, int] = {}
    for sid, atoms in enumerate(systems):
        for a in atoms:
            sys_of[a] = sid
    seen: Set[Tuple[int, int]] = set()
    bonds: List[Tuple[int, int]] = []
    for a in sys_of:
        for b in nbrs[a]:
            if b <= a or syms[b] == "H" or _is_metal_sym(syms[b]):
                continue
            if sys_of.get(b) == sys_of.get(a):
                key = (a, b)
                if key not in seen:
                    seen.add(key)
                    bonds.append(key)
    return bonds, sys_of


def _collect_substituent_owners(
    ring_atoms: Set[int], metal_set: Set[int], nbrs: List[List[int]]
) -> Dict[int, int]:
    """Assign every NON-ring atom to the ring atom it hangs off (BFS over the
    heavy+H graph, blocked by ring atoms AND metals).  Blocking metals is
    essential: without it a coordinated donor substituent could be reached from a
    DIFFERENT ligand's ring via the metal and be dragged off its M–D bond.
    Deterministic: ring atoms are visited in ascending index, first owner wins.
    Atoms behind a metal (other ligands, the metal itself) are left unowned →
    untouched."""
    blocked = ring_atoms | metal_set
    owner: Dict[int, int] = {}
    for a in sorted(ring_atoms):
        stack = [nb for nb in nbrs[a] if nb not in blocked and nb not in owner]
        while stack:
            x = stack.pop()
            if x in blocked or x in owner:
                continue
            owner[x] = a
            for nb in nbrs[x]:
                if nb not in blocked and nb not in owner:
                    stack.append(nb)
    return owner


# ── plane projection ──────────────────────────────────────────────────────────
def _project_to_plane(
    work: np.ndarray, atoms: Sequence[int], frozen: Set[int]
) -> None:
    """Project the non-frozen ``atoms`` onto their best-fit plane in place.  The
    plane is anchored through the frozen (coordinated) atoms so those never move
    (M–D invariant); with no frozen atom the centroid best-fit plane is used."""
    idx = list(atoms)
    if len(idx) < 3:
        return
    anchors = [a for a in idx if a in frozen]
    P = work[idx]
    origin = work[anchors].mean(axis=0) if anchors else P.mean(axis=0)
    try:
        _, _, vh = np.linalg.svd(P - origin, full_matrices=False)
    except np.linalg.LinAlgError:
        return
    normal = vh[-1]
    nn = float(np.linalg.norm(normal))
    if nn < 1e-9:
        return
    normal = normal / nn
    for a in idx:
        if a in frozen:
            continue
        work[a] = work[a] - float((work[a] - origin) @ normal) * normal


# ── PBD ring-bond equalisation ────────────────────────────────────────────────
def _relax_ring_bonds(
    pts: np.ndarray,
    syms: List[str],
    systems: List[List[int]],
    active_atoms: Set[int],
    ring_bonds: List[Tuple[int, int]],
    target: Dict[Tuple[int, int], float],
    frozen: Set[int],
) -> np.ndarray:
    """Return a copy of ``pts`` with the ACTIVE ring atoms reshaped toward the
    per-bond delocalised target (Jacobi PBD, in-plane, frozen anchors fixed)."""
    work = pts.copy()
    act_bonds = [
        (i, j) for (i, j) in ring_bonds
        if (i in active_atoms or j in active_atoms)
    ]
    for _ in range(_MAX_PASSES):
        disp = np.zeros_like(work)
        cnt = np.zeros(len(syms))
        worst = 0.0
        for (i, j) in act_bonds:
            v = work[j] - work[i]
            d = float(np.linalg.norm(v))
            if d < 1e-6:
                continue
            err = target[(i, j)] - d          # +ve → too short; −ve → too long
            worst = max(worst, abs(err))
            u = v / d
            fi, fj = i in frozen, j in frozen
            if fi and fj:
                continue
            if fi:
                disp[j] += err * u; cnt[j] += 1
            elif fj:
                disp[i] -= err * u; cnt[i] += 1
            else:
                disp[i] -= 0.5 * err * u; cnt[i] += 1
                disp[j] += 0.5 * err * u; cnt[j] += 1
        moved = False
        for a in active_atoms:
            if a in frozen or cnt[a] == 0:
                continue
            dv = disp[a] / cnt[a]
            m = float(np.linalg.norm(dv))
            if m > _STEP_CLAMP:
                dv = dv / m * _STEP_CLAMP
            if m > 1e-9:
                work[a] += dv
                moved = True
        for atoms in systems:
            if any(a in active_atoms for a in atoms):
                _project_to_plane(work, atoms, frozen)
        if worst < _CONVERGE_TOL or not moved:
            break
    return work


# ── clash proxy (never-worse guard) ───────────────────────────────────────────
def _clash_count(pts: np.ndarray, syms: List[str]) -> int:
    n = len(syms)
    c = 0
    for i in range(n):
        if syms[i] == "H" or _is_metal_sym(syms[i]):
            continue
        ri = _COV_RADII.get(syms[i], 0.75)
        for j in range(i + 1, n):
            if syms[j] == "H" or _is_metal_sym(syms[j]):
                continue
            d = float(np.linalg.norm(pts[i] - pts[j]))
            if d < _CLASH_FACTOR * (ri + _COV_RADII.get(syms[j], 0.75)):
                c += 1
    return c


def _ring_bond_dev(
    pts: np.ndarray, bonds: Sequence[Tuple[int, int]],
    target: Dict[Tuple[int, int], float],
) -> float:
    return float(sum(
        abs(float(np.linalg.norm(pts[i] - pts[j])) - target[(i, j)])
        for (i, j) in bonds
    ))


# ── public API ────────────────────────────────────────────────────────────────
def correct_xyz(xyz: str) -> str:
    """Equalise aromatic ring bond lengths to their delocalised targets.  Returns
    the ORIGINAL text unchanged (byte-identical) when there is no aromatic ring,
    nothing is off-target, or the reshape does not strictly improve the ring
    bond-length deviation without breaking an M–D bond / adding a clash."""
    syms, pts, lines = _parse_xyz(xyz)
    n = len(syms)
    if n < 5 or not any(s in ("C", "N", "O", "S") for s in syms):
        return xyz

    nbrs = _build_geometric_adjacency(syms, pts)
    rings = _detect_aromatic_rings(syms, pts, nbrs)
    if not rings:
        return xyz
    systems = _fuse_components(rings)
    ring_atoms: Set[int] = set().union(*[set(s) for s in systems])
    ring_bonds, _sys_of = _ring_bonds_of_systems(systems, nbrs, syms)
    if not ring_bonds:
        return xyz

    target = {(i, j): _target_length(syms[i], syms[j]) for (i, j) in ring_bonds}

    # activation: only reshape systems whose worst ring-bond error exceeds the tol.
    sys_worst: Dict[int, float] = {}
    for sid, atoms in enumerate(systems):
        aset = set(atoms)
        errs = [
            abs(float(np.linalg.norm(pts[i] - pts[j])) - target[(i, j)])
            for (i, j) in ring_bonds if i in aset and j in aset
        ]
        sys_worst[sid] = max(errs) if errs else 0.0
    active_atoms: Set[int] = set()
    for sid, atoms in enumerate(systems):
        if sys_worst[sid] > _ACTIVATE_TOL:
            active_atoms.update(atoms)
    if not active_atoms:
        return xyz                                   # every ring already equalised

    metals = [i for i in range(n) if _is_metal_sym(syms[i])]
    metal_set = set(metals)
    coordinated = {
        a for a in range(n)
        if a not in metal_set and any(nb in metal_set for nb in nbrs[a])
    }
    # ownership map (non-ring atom → the ring atom it hangs off), computed up front
    # so the freeze rule can see the WHOLE substituent subtree of each ring atom.
    owner = _collect_substituent_owners(ring_atoms, metal_set, nbrs)

    # freeze ring atoms whose position is pinned by the coordination sphere so the
    # M–D invariant is preserved BY CONSTRUCTION — else the rigid substituent drag
    # would pull a coordinated donor off the metal:
    #   (1) a ring atom itself near/bonded to a metal (pyridine-type ring-N donor);
    #   (2) a ring atom DIRECTLY bonded to a coordinated donor (catecholate /
    #       phenolate O — its ipso carbon);
    #   (3) a ring atom whose substituent SUBTREE contains a coordinated donor
    #       (salen/salophen CH=N–M imine: the coordinated N sits two bonds off the
    #       ring, so freezing only direct neighbours would let the ring drag it and
    #       drift Cu–N).
    frozen: Set[int] = set()
    for r in ring_atoms:
        if any(float(np.linalg.norm(pts[r] - pts[m])) < _M_COORD_DIST
               for m in metals):
            frozen.add(r)
            continue
        if any(nb not in ring_atoms and nb in coordinated for nb in nbrs[r]):
            frozen.add(r)
    for x, a in owner.items():                       # rule (3): deep-subtree donor
        if x in coordinated:
            frozen.add(a)

    # NEVER-WORSE SCOPE (2026-07-17, arom_v1 A/B): a ring SYSTEM that touches the
    # coordination sphere (ANY frozen atom) is left BYTE-IDENTICAL — do not reshape
    # it at all.  Anchoring only the donor still lets the rest of a coordinated ring
    # reshape + re-planarise, which SILENTLY distorts the coordination geometry: the
    # M–D-bond guard checks only donor distances and the clash proxy excludes metals,
    # so neither sees it — but the eye does (round-trip lost / polyhedron distorted /
    # graph-geometry regressed on HEGCEC, ZEYMUL, EFEFAZ, SEJFOF, VIBTAC, ERIBEM,
    # RETFON, LUMTAP).  The mesomeric equalisation is applied ONLY to FREE organic
    # aromatic systems, where the rigid substituent drag cannot perturb an M–D
    # relationship.  Coordinated aromatics need their correct delocalised lengths
    # seated at CONSTRUCTION time (root fix), which is deferred.
    _act_old = active_atoms
    active_atoms = set()
    for atoms in systems:
        aset = set(atoms)
        if (aset & _act_old) and not (aset & frozen):
            active_atoms |= aset
    if not active_atoms:
        return xyz                                   # only coordinated rings were off-target

    work = _relax_ring_bonds(
        pts, syms, systems, active_atoms, ring_bonds, target, frozen
    )

    # assemble the candidate frame: ring atoms take their reshaped position; every
    # substituent subtree is rigidly translated by its owner ring atom's shift.
    new_pts = pts.copy()
    for a in ring_atoms:
        new_pts[a] = work[a]
    for x, a in owner.items():
        new_pts[x] = pts[x] + (work[a] - pts[a])

    # never-worse rollback (ring-bond deviation must strictly drop; M–D intact;
    # no new heavy clash).
    act_bonds = [
        (i, j) for (i, j) in ring_bonds
        if i in active_atoms or j in active_atoms
    ]
    if _ring_bond_dev(new_pts, act_bonds, target) >= \
            _ring_bond_dev(pts, act_bonds, target) - _IMPROVE_EPS:
        return xyz
    if _md_invariant_violated(_snapshot_md_bonds(pts, syms), new_pts):
        return xyz
    if _clash_count(new_pts, syms) > _clash_count(pts, syms):
        return xyz
    return _format_xyz(lines, syms, new_pts)


def correct_results(mol, results):
    """Apply :func:`correct_xyz` to each ``(xyz, label)`` frame.  ``mol`` is
    accepted for dispatch-signature parity but NOT used — perception is geometric
    so the corrector is robust to mol/XYZ atom-order drift."""
    if not results:
        return results
    out = []
    for item in results:
        try:
            xyz, label = item
        except Exception:
            out.append(item)
            continue
        try:
            out.append((correct_xyz(xyz), label))
        except Exception:
            out.append((xyz, label))
    return out
