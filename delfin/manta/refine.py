"""delfin.manta.refine — geometry refiner (NO force field).

Deterministic coordinate descent that directly minimises a structural-defect
count (close-contact clashes / hydrogen over-coordination / collapsed bonds),
moving only the ligand periphery with the metal + donor atoms FROZEN.  Each pass
is accepted only if the defect count decreases, otherwise the step is halved and
rolled back.

This makes the rigid geometric placement clean without a force field on the
metal: the placement already gives correct connectivity and metal-donor
distances; this step removes the residual close contacts that rigid placement
leaves behind.
"""
from __future__ import annotations
from typing import List, Set, Tuple
import os
import numpy as np
import delfin.manta._bond_decollapse as bd
from delfin.manta import cod_ideals as _CODI

# Iter 29 (construction-driver): refine heavy-heavy bonds toward COD-empirical lengths
# (real crystal p50) instead of the generic covalent-sum _ideal_bond, which is
# systematically too long for aromatic/conjugated ligands (C-C 1.52 vs COD 1.40) and
# was the dominant F3_bond / fragfit gap of FF-free vs UFF.  Env-gated, default OFF
# (byte-identical when unset).  CCDC-ready: cod_ideals is swappable for a CSD source.
_USE_COD_BONDS = os.environ.get("DELFIN_FFFREE_COD_BONDS", "0") == "1"

# Iter 30 (User 2026-05-28, RUJSIY eye-validation): when a heavy atom (C/N/O of a ring)
# is moved by the refiner, the H atoms bonded to it MUST move along (rigid-H drag),
# else the X–H bond stretches from 1.08 Å to ~1.3 Å as the heavy drifts away — user
# observed "in den frames bewegen sich nur die C und N atome im ring aber die H bleiben
# an ihrer stelle und dadurch entstehen unrealitische strukturen".  Env-gated default
# OFF (byte-identical when unset); the propagation copies the heavy's accumulated
# displacement onto its bonded H atoms before the step is applied.
_RIGID_H_DRAG = os.environ.get("DELFIN_FFFREE_RIGID_H_DRAG", "0") == "1"

# Iter 31 (aromatic bond-length ROOT seat, User 2026-07-17): seat aromatic ring
# bonds at their delocalised CCDC targets DURING this coordinate descent, so
# metal-COORDINATED aromatic rings are built at the right length (coherent
# relaxation) instead of being dragged post-hoc.  The frozen metal+donor atoms
# (``fixed_idx``) keep the coordination invariant while the FREE ring atoms move
# and the ring<->substituent junction angles co-adapt (no rigid substituent
# translation).  Env-gated default OFF: when unset ``arom_targets`` is never
# built and every ``_violations`` call gets ``None`` -> bit-for-bit identical to
# the current behaviour.  This is the metal-path lever the post-hoc corrector
# (_arom_bond_length, which now skips all coordinated rings) cannot reach.
_AROM_SEAT = os.environ.get("DELFIN_FFFREE_AROM_SEAT", "0") == "1"
_AROM_SEAT_TOL = 0.02      # Å — a ring bond is seated once |d − target| ≤ this
_AROM_SEAT_GAIN = 0.4      # per-pass correction fraction (matches the distort gain)

# Frame-level coordination/validity GUARANTEE for the seat (never-worse by
# construction): after seating, the frame is KEPT only if it stays at-least-as-
# valid as the no-seat refine — MANTA's own topology_preserved (M–D + bond graph +
# donor amine-H) AND no new heavy–heavy clash AND coordination polyhedron preserved
# AND ring planarity not worse; otherwise the seat is reverted BYTE-IDENTICALLY to
# the no-seat result for that frame.  This is what lets a COORDINATED aromatic ring
# seat only when the coordination survives, else be left byte-identical (the LUMTAP
# v1→v2 lesson: donor–M–donor angles are frozen-invariant, so the ONLY signal of a
# coordinated-ring reshape distortion is the clash/planarity of the reshaped ring).
_AROM_PLANARITY_TOL = 0.05     # Å — a seated ring's max out-of-plane may not worsen by more
_POLY_MD_TOL = 0.05            # Å — donor–metal distance must be preserved within this
_POLY_ANGLE_TOL = 3.0          # deg — donor–metal–donor angle must be preserved within this
# A ring SYSTEM that touches the coordination sphere is NEVER seated (pre-excluded):
# reshaping a coordinated ring silently distorts the coordination geometry the eye
# perceives, and the donor–M–donor angles are frozen-invariant so no cheap in-refine
# check can catch it (LUMTAP).  This hard scope is the GUARANTEE; the frame-level
# guard above is the additional belt for the FREE rings that DO seat.  "Touches" =
# a ring atom is within _M_COORD_DIST of a metal, is a frozen donor, or is bonded to
# either — matches _arom_bond_length / _arom_planarize coordinated-ring scope.
_M_COORD_DIST = 2.6            # Å — atom within this of a metal counts as coordinated


def _bond_ideal(a, b):
    if _USE_COD_BONDS:
        v = _CODI.cod_ideal_bond(a, b)
        if v is not None:
            return v
    return bd._ideal_bond(a, b)

_VDW = {"H": 1.20, "C": 1.70, "N": 1.55, "O": 1.52, "F": 1.47, "P": 1.80, "S": 1.80,
        "Cl": 1.75, "Br": 1.85, "I": 1.98, "B": 1.92, "Si": 2.10, "Se": 1.90,
        "As": 1.85, "Te": 2.06}
_VDW_D = 1.7
_CLASH = 0.70          # non-bonded pair < this * (vdW_i+vdW_j) = clash
_COLLAPSE = 0.70       # bond < this * ideal = collapsed


def _vdw(s):
    return _VDW.get(s, _VDW_D)


def _bonds_adj(syms, P):
    b = bd._geometric_bonds(syms, P)
    bonded = {(min(i, j), max(i, j)) for i, j in b}
    adj: List[List[int]] = [[] for _ in syms]
    for i, j in b:
        adj[i].append(j); adj[j].append(i)
    return b, bonded, adj


def _violations(syms, P, bonded, adj, arom_targets=None):
    """Return (loss, moves) where moves = list of (atom_idx, displacement) that
    would relieve a violation.  loss = #clashes + #collapses + #h_overcoord.

    ``arom_targets`` (None unless DELFIN_FFFREE_AROM_SEAT is on) maps each
    perceived aromatic ring bond ``(i, j)`` → its delocalised target length; when
    provided, an additional continuous term seats those bonds (see below).  When
    None, this function is bit-for-bit identical to the pre-seat behaviour."""
    n = len(syms)
    loss = 0
    moves: List[Tuple[int, np.ndarray]] = []
    # clashes (non-bonded, incl. H)
    for i in range(n):
        for j in range(i + 1, n):
            if (i, j) in bonded:
                continue
            d = float(np.linalg.norm(P[i] - P[j]))
            if d < 1e-6:
                continue
            tgt = _CLASH * (_vdw(syms[i]) + _vdw(syms[j]))
            if d < tgt:
                loss += 1
                u = (P[i] - P[j]) / d
                push = 0.5 * (tgt - d)
                moves.append((i, +push * u)); moves.append((j, -push * u))
    # collapsed bonds -> stretch to ideal; distorted heavy-heavy bonds -> toward ideal
    for i, j in bonded:
        d = float(np.linalg.norm(P[i] - P[j]))
        if d <= 1e-6:
            continue
        ideal = _bond_ideal(syms[i], syms[j])
        if d < _COLLAPSE * ideal:
            loss += 1
            u = (P[j] - P[i]) / d
            stretch = 0.5 * (ideal - d)
            moves.append((j, +stretch * u)); moves.append((i, -stretch * u))
        elif syms[i] != "H" and syms[j] != "H":
            dev = (d - ideal) / ideal
            if dev < -0.25 or dev > 0.085:        # bond_distort band
                loss += 1
                u = (P[j] - P[i]) / d
                corr = 0.4 * (ideal - d)           # +stretch if short, -compress if long
                moves.append((j, +corr * u)); moves.append((i, -corr * u))
    # H over-coordination -> push H off the non-parent heavy
    for h in range(n):
        if syms[h] != "H":
            continue
        heavies = [k for k in range(n) if k != h and syms[k] != "H"
                   and not bd._is_metal(syms[k])
                   and np.linalg.norm(P[h] - P[k]) < 1.3 * bd._ideal_bond(syms[h], syms[k])]
        if len(heavies) >= 2:
            loss += 1
            parent = min(heavies, key=lambda k: np.linalg.norm(P[h] - P[k]))
            for k in heavies:
                if k == parent:
                    continue
                d = float(np.linalg.norm(P[h] - P[k]))
                if d > 1e-6:
                    moves.append((h, 0.4 * (P[h] - P[k]) / d))
    # Aromatic ring-bond seating (DELFIN_FFFREE_AROM_SEAT).  Pull each perceived
    # aromatic ring bond toward its delocalised target.  Moves on frozen
    # metal/donor atoms are discarded by the caller, so a coordinated ring bond
    # is seated by moving ONLY its free atom while the donor stays put (M–D
    # invariant preserved by construction).  The penalty is the CONTINUOUS |err|
    # in Å — always far smaller than a unit clash/collapse/H-overcoord defect, so
    # it never overrides coordination integrity, and any partial seating strictly
    # lowers the loss so the accept-if-better gate settles it monotonically.
    if arom_targets:
        for (i, j), tgt in arom_targets.items():
            if (i, j) not in bonded:
                continue
            d = float(np.linalg.norm(P[i] - P[j]))
            if d <= 1e-6:
                continue
            err = tgt - d                       # +ve → too short; −ve → too long
            if abs(err) > _AROM_SEAT_TOL:
                loss += abs(err)
                u = (P[j] - P[i]) / d
                corr = _AROM_SEAT_GAIN * err
                moves.append((j, +corr * u)); moves.append((i, -corr * u))
    return loss, moves


def _precompute_arom_targets(syms, P, fixed=None):
    """Perceive the SEATABLE aromatic ring bonds ONCE (geometry-only, robust to
    mol/XYZ atom-order drift — the same reason the other MANTA correctors are
    geometric).

    Return ``(targets, systems)`` where ``targets`` is
    ``{(i, j): delocalised_target_length}`` for every intra-system aromatic ring
    bond of a FREE ring system that has a vendored target (``None`` when nothing is
    seatable), and ``systems`` is the list of the corresponding free fused ring-atom
    groups (consumed by the frame-level planarity guard).  Both empty/None → the
    caller treats the seat as inactive (byte-identical to flag-off).

    Ring systems that TOUCH the coordination sphere are EXCLUDED (never seated):
    any ring atom within ``_M_COORD_DIST`` of a metal, in the frozen ``fixed``
    set (metal+donors), or bonded to such an atom.  Reshaping a coordinated ring
    distorts the coordination geometry the eye perceives while the donor–M–donor
    angles stay frozen-invariant (LUMTAP) — so coordinated rings are left byte-
    identical by construction and only free rings seat.

    Imports are lazy so the module's import cost/order is unchanged when the flag
    is off, and so no import cycle is introduced on the always-on path."""
    try:
        from delfin.manta._pi_h_projector import (
            _build_geometric_adjacency, _is_metal_sym,
        )
        from delfin.manta._arom_planarize import (
            _detect_aromatic_rings, _fuse_components,
        )
        from delfin.manta._arom_bond_length import _ring_bonds_of_systems
        from delfin.fffree.aromatic_bond_targets import aromatic_target
    except Exception:
        return None, []
    try:
        syms = list(syms)
        Pn = np.asarray(P, float)
        nbrs = _build_geometric_adjacency(syms, Pn)
        rings = _detect_aromatic_rings(syms, Pn, nbrs)
        if not rings:
            return None, []
        systems = _fuse_components(rings)
        ring_bonds, _sys_of = _ring_bonds_of_systems(systems, nbrs, syms)
        # coordination sphere = frozen metal+donors (fixed) ∪ any atom within
        # _M_COORD_DIST of a metal (self-contained; catches a donor NOT in `fixed`).
        coord = set(fixed) if fixed else set()
        metals = [i for i in range(len(syms)) if _is_metal_sym(syms[i])]
        for m in metals:
            for a in range(len(syms)):
                if a != m and float(np.linalg.norm(Pn[a] - Pn[m])) < _M_COORD_DIST:
                    coord.add(a)
        # keep only FREE ring systems (no ring atom is, or is bonded to, a
        # coordination-sphere atom) — coordinated rings are left byte-identical.
        seat_systems = []
        for atoms in systems:
            touches = any(
                (a in coord) or any(nb in coord for nb in nbrs[a]) for a in atoms
            )
            if not touches:
                seat_systems.append(atoms)
        seat_atoms = set().union(*[set(s) for s in seat_systems]) if seat_systems else set()
        out = {}
        for (i, j) in ring_bonds:
            if i not in seat_atoms or j not in seat_atoms:
                continue
            t = aromatic_target(syms[i], syms[j])
            if t is not None:
                out[(min(i, j), max(i, j))] = t
        return (out or None), seat_systems
    except Exception:
        return None, []


def _ring_oop(P, atoms) -> float:
    """Max out-of-plane distance (Å) of the ring atoms from their SVD best-fit
    plane — a planarity proxy for a seated ring."""
    idx = list(atoms)
    if len(idx) < 3:
        return 0.0
    pts = np.asarray(P, float)[idx]
    c = pts.mean(axis=0)
    try:
        _, _, vh = np.linalg.svd(pts - c, full_matrices=False)
    except np.linalg.LinAlgError:
        return 0.0
    nrm = vh[-1] / (float(np.linalg.norm(vh[-1])) + 1e-12)
    return float(np.max(np.abs((pts - c) @ nrm)))


def _heavy_clash_count(syms, P) -> int:
    """Heavy–heavy NON-bonded pairs closer than the refine clash floor (metals and
    H excluded).  Includes 1–3 ring pairs, so a ring CONTRACTION that pulls a meta
    atom pair into the clash band — the LUMTAP coordinated-ring failure mode — is
    counted."""
    P = np.asarray(P, float)
    b = bd._geometric_bonds(syms, P)
    bonded = {(min(i, j), max(i, j)) for i, j in b}
    n = len(syms)
    c = 0
    for i in range(n):
        if syms[i] == "H" or bd._is_metal(syms[i]):
            continue
        for j in range(i + 1, n):
            if syms[j] == "H" or bd._is_metal(syms[j]):
                continue
            if (i, j) in bonded:
                continue
            d = float(np.linalg.norm(P[i] - P[j]))
            if d < _CLASH * (_vdw(syms[i]) + _vdw(syms[j])):
                c += 1
    return c


def _polyhedron_preserved(syms, P0, P1, fixed) -> bool:
    """True iff every donor–metal distance and donor–metal–donor angle is preserved
    between ``P0`` and ``P1`` within tolerance.  Donors are frozen in the descent so
    this is normally exact; it is a belt that also catches a donor that was NOT in
    ``fixed`` (detection miss) and got dragged."""
    P0 = np.asarray(P0, float); P1 = np.asarray(P1, float)
    metals = [i for i in fixed if bd._is_metal(syms[i])]
    if not metals:
        metals = [i for i in range(len(syms)) if bd._is_metal(syms[i])]
    if not metals:
        return True
    donors = sorted(a for a in fixed if not bd._is_metal(syms[a]))
    for m in metals:
        for a in donors:
            if abs(float(np.linalg.norm(P1[a] - P1[m]))
                   - float(np.linalg.norm(P0[a] - P0[m]))) > _POLY_MD_TOL:
                return False
        near = [a for a in donors if float(np.linalg.norm(P0[a] - P0[m])) < 3.0]
        for x in range(len(near)):
            for y in range(x + 1, len(near)):
                a, b2 = near[x], near[y]

                def _ang(P, _a=a, _b=b2):
                    va = P[_a] - P[m]; vb = P[_b] - P[m]
                    cc = float(np.dot(va, vb) / (
                        np.linalg.norm(va) * np.linalg.norm(vb) + 1e-12))
                    return float(np.degrees(np.arccos(max(-1.0, min(1.0, cc)))))
                if abs(_ang(P1) - _ang(P0)) > _POLY_ANGLE_TOL:
                    return False
    return True


def _seat_frame_ok(syms, P_noseat, P_seat, systems, fixed) -> bool:
    """Frame-level never-worse GUARANTEE for the aromatic seat: accept the seated
    frame ONLY if it is at-least-as-valid as the no-seat refine by EVERY check
    below; otherwise the caller reverts to the no-seat frame byte-identically.
    Conservative by design (any doubt → reject):

      1. MANTA's OWN topology/coordination predicate (``topology_preserved``):
         M–D distances + M–D edge set + covalent bond multiset + donor amine-H
         realism all unchanged.
      2. No NEW heavy–heavy clash (a coordinated-ring contraction pulling a ring
         pair into the clash band — the LUMTAP failure — is rejected).
      3. Coordination polyhedron (donor–M distances + donor–M–donor angles) kept.
      4. No seated ring loses planarity.
    """
    try:
        from delfin.manta._topology_hash import topology_preserved
        if not topology_preserved(list(syms), P_noseat, P_seat).passed:
            return False
    except Exception:
        return False                       # cannot verify → conservative reject
    if _heavy_clash_count(syms, P_seat) > _heavy_clash_count(syms, P_noseat):
        return False
    if not _polyhedron_preserved(syms, P_noseat, P_seat, fixed):
        return False
    for atoms in systems:
        if _ring_oop(P_seat, atoms) > _ring_oop(P_noseat, atoms) + _AROM_PLANARITY_TOL:
            return False
    return True


def _refine_core(syms, P, fixed, arom_targets, max_passes: int, damp: float):
    """The coordinate-descent refiner proper.  ``arom_targets`` None → bit-for-bit
    the pre-seat refiner (so the flag-off path is byte-identical)."""
    P = np.array(P, float).copy()
    n = len(syms)
    _, bonded, adj = _bonds_adj(syms, P)
    best_loss, _ = _violations(syms, P, bonded, adj, arom_targets)
    if best_loss == 0:
        return P
    for _ in range(max_passes):
        _, bonded, adj = _bonds_adj(syms, P)
        loss, moves = _violations(syms, P, bonded, adj, arom_targets)
        if loss == 0:
            break
        disp = np.zeros((n, 3))
        cnt = np.zeros(n)
        for idx, d in moves:
            if idx in fixed:
                continue
            disp[idx] += d; cnt[idx] += 1
        for i in range(n):
            if cnt[i] > 0:
                disp[i] /= cnt[i]
        # Iter 30: rigid-H drag — propagate heavy-atom displacement to bonded H atoms
        # so X-H bonds preserve length instead of stretching as the heavy drifts.
        # H atoms must not already have an own displacement (else they get double-pushed
        # — overwrite only when cnt[h]==0).
        if _RIGID_H_DRAG:
            # bonded H map: for each H, identify its closest heavy (parent) using the
            # current bonded graph that produced this pass (`adj`)
            for hi in range(n):
                if syms[hi] != "H" or hi in fixed or cnt[hi] > 0:
                    continue
                parents = [k for k in adj[hi] if syms[k] != "H" and not bd._is_metal(syms[k])]
                if not parents:
                    continue
                par = parents[0] if len(parents) == 1 else min(
                    parents, key=lambda k: float(np.linalg.norm(P[hi] - P[k])))
                if cnt[par] > 0:
                    disp[hi] = disp[par]      # rigid drag: H rides with its heavy parent
        trial = P + damp * disp
        _, tb, ta = _bonds_adj(syms, trial)
        tloss, _ = _violations(syms, trial, tb, ta, arom_targets)
        if tloss < best_loss:          # accept-if-better gate
            P = trial; best_loss = tloss
        else:
            damp *= 0.6                # smaller step
            if damp < 0.05:
                break
    return P


def refine(syms, P, fixed_idx: Set[int], max_passes: int = 80, damp: float = 0.6):
    """Coordinate descent minimising the defect loss; metal+donors frozen.
    Deterministic.  Returns refined P (or the original if no improvement).

    When DELFIN_FFFREE_AROM_SEAT is on, aromatic ring bonds are additionally seated
    toward their delocalised targets DURING the descent — but the seated frame is
    KEPT only if it stays at-least-as-valid as the no-seat refine (``_seat_frame_ok``);
    otherwise it reverts BYTE-IDENTICALLY to the no-seat result for that frame.  So a
    coordinated aromatic ring is seated only when the coordination survives, else it
    is left byte-identical (LUMTAP).  When off, this is bit-for-bit the original
    refiner (the no-seat core with ``arom_targets`` None)."""
    fixed = set(fixed_idx)
    if not _AROM_SEAT:
        return _refine_core(syms, P, fixed, None, max_passes, damp)
    P0 = np.array(P, float).copy()
    targets, systems = _precompute_arom_targets(syms, P0, fixed)
    if not targets:
        return _refine_core(syms, P0, fixed, None, max_passes, damp)
    P_noseat = _refine_core(syms, P0, fixed, None, max_passes, damp)
    P_seat = _refine_core(syms, P0, fixed, targets, max_passes, damp)
    if _seat_frame_ok(syms, P_noseat, P_seat, systems, fixed):
        return P_seat
    return P_noseat
