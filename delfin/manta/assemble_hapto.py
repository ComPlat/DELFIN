"""Hapto assembly of the FF-free constructor: eta ring placement, piano-stool leg tilt, ring spins, slip modes, puckers and the hapto ensemble.

Moved verbatim from delfin/manta/assemble_complex.py (MANTA split, 2026-10); every
statement is the original text, only the imports below are new.
"""
from __future__ import annotations

import numpy as np
import os
from rdkit import Chem

from delfin.manta import metal_sphere_builder as MSB
from delfin.manta.assemble_chelate import (
    _place_chelate_block,
)
from delfin.manta.assemble_donor_plane import (
    _canonical_arm_order,
)
from delfin.manta.assemble_fold_fp import (
    _fold_fp,
    _fold_fp_enabled,
    _fold_mc_arms,
    _fold_rings_with_mc,
    _fold_same,
)
from delfin.manta.assemble_ligand_confs import (
    _clash_count,
    _ligand_confs_from_mol,
    _refine_guarded,
)
from delfin.manta.assemble_ligand_embed import (
    _kabsch_rot,
    _rot_align,
)
from delfin.manta.assemble_orient import (
    _orient_chelate_to_vertices,
    _vsepr_reconstruct,
)


def _eta_centroid_distance(metal, eta_n):
    """Crystallographic metal→ring-centroid distance (Å) for an η-face.  Pulls the
    open-source literature-averaged table from smiles_converter (a hard-coded table
    of published averages; reads no proprietary database file) and falls back to its
    geometric estimator.  Deterministic; always finite."""
    try:
        from delfin.smiles_converter import _target_mc_dist
        d = float(_target_mc_dist(metal, int(eta_n)))
        if np.isfinite(d) and d > 0.5:
            return d
    except Exception:
        pass
    # last-resort: covalent M-C sum minus a ring-radius slip (always finite)
    from delfin.manta import polyhedra as PLY
    rmc = PLY.COV.get(metal, 1.5) + PLY.COV.get("C", 0.76)
    if eta_n >= 3:
        rr = 1.40 / (2.0 * np.sin(np.pi / max(eta_n, 3)))
        v = rmc * rmc - rr * rr
        return float(np.sqrt(v)) if v > 0.25 else 0.80 * rmc
    return 0.85 * rmc


def _rot_about_axis(axis, theta):
    """Rodrigues rotation matrix: rotate by ``theta`` (rad) about unit ``axis``."""
    axis = np.asarray(axis, float)
    n = np.linalg.norm(axis)
    if n < 1e-12:
        return np.eye(3)
    axis = axis / n
    c, s = float(np.cos(theta)), float(np.sin(theta))
    x, y, z = axis
    K = np.array([[0.0, -z, y], [z, 0.0, -x], [-y, x, 0.0]])
    return np.eye(3) * c + s * K + (1.0 - c) * np.outer(axis, axis)


def _piano_leg_tilt_enabled() -> bool:
    """Piano-stool σ-leg TILT correction (default OFF -> byte-identical).

    A half-sandwich / piano-stool (ONE η-face + a small set of σ-legs) is NOT a
    regular tetrahedron: the η-ring is a fat 3-electron-pair FACE donor, and the legs
    splay DOWN, away from the ring, so the real (ring-centroid)-M-L angle is ~120-130
    deg (not the ~90-110 deg a generic CN4/CN5 polyhedron places the legs at).  When
    ON, the legs are rigidly tilted to that piano-stool angle; OFF, the block is
    skipped and the build is byte-identical to the historic polyhedron placement."""
    return os.environ.get("DELFIN_FFFREE_PIANO_LEG_TILT", "0") == "1"


def _piano_leg_target_deg() -> float:
    """Target (η-ring-centroid)-M-(σ-leg) tilt angle in degrees for the piano-stool
    leg-tilt correction.  Default 125 deg -- the canonical crystallographic value for
    half-sandwich piano-stools: (arene)Cr(CO)3 centroid-Cr-C ~125.7 deg, CpMn(CO)3
    centroid-Mn-C ~125 deg, CpFe(CO)2X ~121-126 deg.  Override via
    DELFIN_FFFREE_PIANO_LEG_DEG; the legs splay DOWN (away from the ring), so values
    >90 deg point the legs to the FAR hemisphere from the centroid axis."""
    try:
        v = float(os.environ.get("DELFIN_FFFREE_PIANO_LEG_DEG", "125.0"))
    except (TypeError, ValueError):
        return 125.0
    # clamp to a sane piano-stool window so a bad env value can never flip the tripod
    return min(160.0, max(95.0, v))


def _tilt_piano_legs(P, prim_V, sigma_blocks, target_deg):
    """Rigidly tilt each σ-leg block to the piano-stool (ring-centroid)-M-L angle.

    ``P`` is the assembled coordinate array with the metal at the origin (row 0).
    ``prim_V`` is the UNIT M->(primary η-ring centroid) direction.  ``sigma_blocks``
    is the list of (start, end) atom-row half-open ranges of the σ-co-ligand blocks
    (the legs).  ``target_deg`` is the desired centroid-M-L angle.

    For each leg, the representative donor (the block atom CLOSEST to the metal) gives
    the leg's current M->L direction.  The whole leg block is rotated RIGIDLY about an
    axis through the metal (the origin) that is PERPENDICULAR to the plane spanned by
    ``prim_V`` and the leg direction, by exactly the angle needed to bring the
    centroid-M-L angle to ``target_deg``.  Because the rotation passes through the
    metal it preserves every atom's metal distance (M-L INVARIANT) and the leg's
    internal geometry; it only swings the leg DOWN, away from the η-ring.  Each leg
    keeps its own azimuth about ``prim_V`` (the rotation axis is radial-tangential),
    so a 3-fold-symmetric tripod stays 3-fold symmetric.  The η-ring is never touched.

    Universal over leg count (1..n) and ring type; geometry-only, deterministic,
    never raises (returns ``P`` unchanged on any degeneracy).  Returns the modified
    array (in-place edits applied to ``P``)."""
    prim_V = np.asarray(prim_V, float)
    nV = np.linalg.norm(prim_V)
    if nV < 1e-9 or not sigma_blocks:
        return P
    prim_V = prim_V / nV
    tgt = np.radians(float(target_deg))
    M = P[0]
    for (s, e) in sigma_blocks:
        if e <= s:
            continue
        blk = P[s:e]
        # representative leg donor = block atom nearest the metal (the σ-donor).
        rel = blk - M
        dists = np.linalg.norm(rel, axis=1)
        di = int(np.argmin(dists))
        L = rel[di]
        nL = np.linalg.norm(L)
        if nL < 1e-9:
            continue
        Lu = L / nL
        cur = float(np.arccos(np.clip(np.dot(prim_V, Lu), -1.0, 1.0)))
        delta = tgt - cur
        if abs(delta) < 1e-6:
            continue
        # rotation axis = prim_V x leg-direction (perpendicular to the swing plane).
        # Sign of `delta` then moves the leg toward/away from the centroid axis; a
        # positive delta (target > current) opens the angle -> legs splay DOWN.
        ax = np.cross(prim_V, Lu)
        nax = np.linalg.norm(ax)
        if nax < 1e-8:
            # leg is (anti)parallel to the centroid axis -> tilt direction undefined;
            # use any axis perpendicular to prim_V (deterministic) so it still splays.
            ref = np.array([1.0, 0.0, 0.0]) if abs(prim_V[0]) < 0.9 else np.array([0.0, 1.0, 0.0])
            ax = np.cross(prim_V, ref)
            nax = np.linalg.norm(ax)
            if nax < 1e-8:
                continue
        ax = ax / nax
        R = _rot_about_axis(ax, delta)
        # rotate the WHOLE leg block about the metal at the origin (no recentre).
        P[s:e] = blk @ R.T
    return P


def _place_eta_ring(metal, lsyms, lP, eta_idxs, Vunit, mc_dist,
                    spin=0.0, pucker=0.0):
    """Rigidly seat an η-face so its ring CENTROID sits on the polyhedron vertex
    direction ``Vunit`` at distance ``mc_dist`` from the metal (origin), with the
    ring plane PERPENDICULAR to the M→centroid axis.  The whole ligand (ring +
    substituents) is moved as one rigid body, so substituents are dragged with the
    ring and the η-carbons sit at the ring radius AROUND the centroid (never on the
    metal).  Returns the transformed coords, or None on degeneracy.

    Geometry-only, deterministic.  The in-plane spin is fixed canonically (first
    ring atom placed at azimuth 0) so two builds of the same input are identical.

    ``spin`` (rad): additional rigid rotation of the whole ligand ABOUT the
    M→centroid axis — generates the η-ring ROTAMER orientations (rotates substituent
    positions around the ring while the ring itself maps onto itself up to symmetry).
    ``pucker`` (Å): out-of-plane displacement amplitude of alternating ring atoms
    (Cremer-Pople-style envelope/twist mode) applied to NON-aromatic η-rings before
    seating, so diene / allyl / cyclohexadienyl rings emit puckered conformers.  Both
    default to 0.0 -> byte-identical to the canonical single rigid build."""
    lP = np.asarray(lP, float).copy()
    ring = [int(i) for i in eta_idxs]
    if len(ring) < 2:
        return None
    R = lP[ring]
    cen = R.mean(axis=0)                                   # ring centroid (ligand frame)
    X = R - cen
    # ring-plane normal via SVD (smallest singular direction); robust for any n>=2.
    try:
        _, _, Vt = np.linalg.svd(X)
    except Exception:
        return None
    if len(ring) == 2:
        # η2 (alkene): "plane normal" is ambiguous; use the C=C bond direction as the
        # in-plane axis and pick a normal perpendicular to it (deterministic).
        nrm = Vt[2] if Vt.shape[0] > 2 else np.cross(Vt[0], np.array([0.0, 0.0, 1.0]))
    else:
        nrm = Vt[2]
    nn = np.linalg.norm(nrm)
    if nn < 1e-8:
        return None
    nrm = nrm / nn
    # Cremer-Pople-style pucker: displace ring atoms along the ring NORMAL with an
    # alternating sign (B/T-like envelope), in the LIGAND frame, before seating.  Only
    # the ring atoms move; substituents stay rigid relative to the (unpuckered) ring
    # carbon they hang off — this is a deterministic conformer, not an MMFF relax.
    if abs(pucker) > 1e-9 and len(ring) >= 4:
        for k, ai in enumerate(ring):
            lP[ai] = lP[ai] + nrm * (float(pucker) * (1.0 if (k % 2 == 0) else -1.0))
        # recompute centroid after the symmetric in-out displacement (≈ unchanged)
        cen = lP[ring].mean(axis=0)
    Vunit = np.asarray(Vunit, float)
    Vunit = Vunit / np.linalg.norm(Vunit)
    target_centroid = Vunit * float(mc_dist)
    # Rotate so the ring-plane normal aligns with the M→centroid axis: the metal
    # then sits along the ring normal (face-on coordination), exactly as in a real
    # piano-stool / metallocene.  Sign chosen so the ring is on the FAR side of the
    # metal (centroid points away from origin along +Vunit).
    Rrot = _rot_align(nrm, Vunit)
    Q = (lP - cen) @ Rrot.T + target_centroid
    if abs(spin) > 1e-9:
        # rotate the seated ligand about the M→centroid axis through the centroid:
        # spins substituents to a distinct rotamer; the bare ring maps onto itself.
        Spin = _rot_about_axis(Vunit, float(spin))
        Q = (Q - target_centroid) @ Spin.T + target_centroid
    if not np.all(np.isfinite(Q)):
        return None
    return Q


def assemble_hapto(metal, geometry, d, variant=None):
    """Build a hapto complex on the FF-free path: η-faces are placed as RIGID rings
    on their centroid vertices; σ-donors are oriented onto the remaining vertices
    (homoleptic/heteroleptic monodentate).  ``d`` is a rigid-hapto decompose dict
    (has per-ligand 'is_eta'/'eta_local_idxs'/'eta_n').  v1: chelating σ-arms are
    placed by the existing rigid chelate block.  Returns (syms, P, donors) where
    ``donors`` are the global indices fffree treats as the coordination shell (one
    representative atom per η-face + every σ-donor), or None on failure.

    Universal, geometry-only, deterministic (fixed seeds, canonical in-plane spin).

    Returns the tuple (syms, P, donors, exempt_pairs) — exempt_pairs are the global
    index pairs of GENUINE multiple bonds (C≡O carbonyls, C≡N nitriles, C=C, …) whose
    short length is chemically correct, so the self-gate does not mistake them for a
    collapse.

    ``variant`` (optional dict, default None = byte-identical canonical build): a
    per-build perturbation spec used by ``assemble_hapto_ensemble`` to enumerate the
    rigid ENSEMBLE.  Recognised keys::

        variant["eta"][li] = {"spin": rad, "pucker": Å, "eta_idxs": [...], "eta_n": k}

    where ``spin`` rotates the seated η-ligand about the M→centroid axis (η-ring
    ROTAMER), ``pucker`` applies a Cremer-Pople envelope to non-aromatic rings, and
    ``eta_idxs``/``eta_n`` (optional) override the ring-atom set + hapticity for a
    deterministic ring-SLIP (η6→η4→η2 / η5→η3→η1) isomer.  All variation is geometric
    + deterministic; an absent/None variant reproduces the historic single build."""
    variant = variant or {}
    eta_var = variant.get("eta", {})
    # Hebel A (DELFIN_FFFREE_HAPTO_AXIS_ROT): rigid rotation (rad) of every NON-
    # primary-η atom block (the CO-tripod / co-ligands / a 2nd ring) about the
    # PRIMARY η-face M→centroid axis.  The metal sits at the origin, so any rotation
    # about an axis through the origin preserves every atom's distance to the metal
    # (M-centroid + all M-D distances INVARIANT); only the ring-vs-rest azimuthal
    # CLOCK changes — the η-ring internal geometry is untouched.  0.0 = byte-identical
    # to the historic build (the rotation is skipped entirely).
    _axis_rot = float(variant.get("axis_rot", 0.0) or 0.0)
    ref = MSB._ref_vectors(geometry)
    n_vert = len(ref)
    ligands = d["ligands"]
    # how many vertices each ligand occupies (η-face = 1, σ = denticity)
    occ = [(1 if lg.get("is_eta") else lg["denticity"]) for lg in ligands]
    if sum(occ) != n_vert:
        return None
    # Deterministic vertex assignment: η-faces first (lowest vertex indices), then
    # σ-donors in ligand order.  Simple + reproducible; isomer enumeration is a
    # later increment (v1 emits ONE faithful build per hapto complex).
    order = sorted(range(len(ligands)), key=lambda i: (0 if ligands[i].get("is_eta") else 1, i))
    out_syms = [metal]
    placed = [np.zeros(3)]
    placed_syms = [metal]
    donors = []
    relax_frags = []
    exempt_pairs = []        # global index pairs that are genuine double/triple bonds
    pos = 1
    vi = 0
    fixed = {0}
    # Hebel-A bookkeeping: per-η-face (atom-block-start, atom-block-end, centroid-unit-
    # direction); and the contiguous atom blocks of every NON-η ligand (CO/σ co-ligands).
    eta_blocks = []          # [(start, end, Vunit), ...] in emission (= placement) order
    nonprimary_eta_blocks = []   # (start, end) of η faces 2..n  -> rotated with the rest
    sigma_blocks = []        # (start, end) of σ ligands

    # LENGTH-GATED HAPTO EXEMPTION (DELFIN_FFFREE_HAPTO_EXEMPT_LENGTH, default OFF -> byte-identical).
    # The hapto path collects every heavy-heavy bond of order >= 1.5 into a LIST, and
    # converter_backend._is_exempt reads the list form as an UNCONDITIONAL pass ("if key in _ex:
    # return True") -- so once a bond is multiple/aromatic, ANY length survives the collapse
    # self-gate, however crushed.  Measured: AJUWUY ships a C-C at 0.495 A through this door.
    # The dict form is already implemented on the reading side (_ex_len + DELFIN_FFFREE_MULTIBOND_TOL
    # 0.15) and gates on `d >= ideal - tol`, which still passes every genuine short multiple bond
    # (C=O 1.13, C#N 1.16, aromatic C~C 1.39) while catching the collapses.  Class: 908 systems,
    # 86.6% of them no_valid -- it tracks the W 69% / Re 53% / Rh 48% / Mo 44% failure rates.
    # Pairs with no known ideal map to 0.0 = unconditional, i.e. exactly today's behaviour.
    _hapto_len_gate = os.environ.get("DELFIN_FFFREE_HAPTO_EXEMPT_LENGTH", "0") == "1"
    if _hapto_len_gate:
        exempt_pairs = {}

    def _collect_exempt(frag_mol, lig_offset):
        """Local heavy-atom double/triple bonds -> global index pairs (the AddHs
        ligand block starts at lig_offset+1 in the assembled coords)."""
        for b in frag_mol.GetBonds():
            bt = b.GetBondTypeAsDouble()
            if bt >= 1.5:                       # aromatic(1.5)/double(2)/triple(3)
                a1, a2 = b.GetBeginAtomIdx(), b.GetEndAtomIdx()
                at1 = frag_mol.GetAtomWithIdx(a1)
                at2 = frag_mol.GetAtomWithIdx(a2)
                if at1.GetAtomicNum() > 1 and at2.GetAtomicNum() > 1:
                    g1 = lig_offset + 1 + a1
                    g2 = lig_offset + 1 + a2
                    key = (min(g1, g2), max(g1, g2))
                    if _hapto_len_gate:
                        from delfin.manta.converter_backend import multibond_ideal as _mbi
                        exempt_pairs[key] = _mbi(at1.GetSymbol(), at2.GetSymbol(), bt) or 0.0
                    else:
                        exempt_pairs.append(key)
    for li in order:
        lg = ligands[li]
        lig_offset = pos - 1
        if lg.get("is_eta"):
            confs = _ligand_confs_from_mol(lg["mol"])
            if confs is None:
                return None
            lsyms, coords_list, lmol = confs
            Vunit = ref[vi] / np.linalg.norm(ref[vi])
            # per-build variation (ensemble): η-ring spin/pucker + optional ring-slip.
            _ev = eta_var.get(li, {})
            _spin = float(_ev.get("spin", 0.0))
            _pucker = float(_ev.get("pucker", 0.0))
            _eta_idxs = _ev.get("eta_idxs", lg["eta_local_idxs"])
            _eta_n = int(_ev.get("eta_n", lg["eta_n"]))
            mc = _eta_centroid_distance(metal, _eta_n)
            best_Q, best_clash = None, 1e18
            for lP in coords_list:
                Q = _place_eta_ring(metal, lsyms, lP, _eta_idxs, Vunit, mc,
                                    spin=_spin, pucker=_pucker)
                if Q is None:
                    continue
                cl = _clash_count(Q, np.array(placed), lsyms, placed_syms)
                if cl < best_clash:
                    best_clash, best_Q = cl, Q
                if cl == 0:
                    break
            if best_Q is None:
                return None
            out_syms += lsyms
            for row in best_Q:
                placed.append(row)
            placed_syms += lsyms
            # representative donor = ring atom closest to the metal (for the self-gate).
            # Uses the ACTIVE ring set (_eta_idxs) so ring-slip isomers freeze + report
            # the η-carbons they actually coordinate, not the full pre-slip face.
            ring_globals = [lig_offset + 1 + int(j) for j in _eta_idxs]
            rep = min(ring_globals, key=lambda g: float(np.linalg.norm(placed[g])))
            donors.append(rep)
            fixed.update(lig_offset + 1 + int(j) for j in _eta_idxs)
            relax_frags.append((Chem.AddHs(lg["mol"]), lig_offset))
            _collect_exempt(lg["mol"], lig_offset)
            _blk = (lig_offset + 1, lig_offset + 1 + len(lsyms))
            eta_blocks.append((_blk[0], _blk[1], Vunit.copy()))
            if len(eta_blocks) > 1:               # secondary η face -> rotates with rest
                nonprimary_eta_blocks.append(_blk)
            pos += len(lsyms)
            vi += 1
        elif lg["denticity"] == 1:
            confs = _ligand_confs_from_mol(lg["mol"])
            if confs is None:
                return None
            lsyms, coords_list, lmol = confs
            di = lg["donor_local_idxs"][0]
            Vunit = ref[vi] / np.linalg.norm(ref[vi])
            md = MSB.md_distance(metal, lsyms[di],
                                 atom=lmol.GetAtomWithIdx(di), mol=lmol)
            best_Q, best_clash = None, 1e18
            for lP in coords_list:
                if len(lsyms) == 1:
                    Q = (Vunit * md).reshape(1, 3)
                else:
                    lPv, lp = _vsepr_reconstruct(lsyms, lP, lmol, di)
                    Q = (lPv - lPv[di]) @ _rot_align(lp, -Vunit).T + Vunit * md
                cl = _clash_count(Q, np.array(placed), lsyms, placed_syms)
                if cl < best_clash:
                    best_clash, best_Q = cl, Q
                if cl == 0:
                    break
            if best_Q is None:
                return None
            out_syms += lsyms
            for row in best_Q:
                placed.append(row)
            placed_syms += lsyms
            donors.append(lig_offset + 1 + di)
            fixed.add(lig_offset + 1 + di)
            relax_frags.append((Chem.AddHs(lg["mol"]), lig_offset))
            _collect_exempt(lg["mol"], lig_offset)
            sigma_blocks.append((lig_offset + 1, lig_offset + 1 + len(lsyms)))
            pos += len(lsyms)
            vi += 1
        else:
            # chelating σ-arm: rigid-fit onto the next `denticity` vertices.
            dent = lg["denticity"]
            confs = _ligand_confs_from_mol(lg["mol"])
            if confs is None:
                return None
            lsyms, coords_list, lmol = confs
            dons_d = _canonical_arm_order(lg, dent)
            verts = list(range(vi, vi + dent))
            targets = [ref[verts[i]] / np.linalg.norm(ref[verts[i]])
                       * MSB.md_distance(metal, lsyms[dons_d[i]],
                                         atom=lmol.GetAtomWithIdx(dons_d[i]), mol=lmol)
                       for i in range(dent)]
            best_Q, best_clash = None, 1e18
            for lP in coords_list:
                if dent == 2:
                    Q = _place_chelate_block(metal, lsyms, lP, dons_d[0], dons_d[1],
                                             targets[0], targets[1], mol=lg["mol"])
                else:
                    _rigid_seat = (os.environ.get("DELFIN_FFFREE_LIGAND_RIGID", "0") == "1"
                                   and dent >= 3)
                    Q = _orient_chelate_to_vertices(lP, dons_d, targets, asym=True,
                                                    rigid=_rigid_seat, lsyms=lsyms)
                if Q is None:
                    continue
                cl = _clash_count(Q, np.array(placed), lsyms, placed_syms)
                if cl < best_clash:
                    best_clash, best_Q = cl, Q
                if cl == 0:
                    break
            if best_Q is None:
                return None
            out_syms += lsyms
            for row in best_Q:
                placed.append(row)
            placed_syms += lsyms
            for di in dons_d:
                donors.append(lig_offset + 1 + di)
                fixed.add(lig_offset + 1 + di)
            relax_frags.append((Chem.AddHs(lg["mol"]), lig_offset))
            _collect_exempt(lg["mol"], lig_offset)
            sigma_blocks.append((lig_offset + 1, lig_offset + 1 + len(lsyms)))
            pos += len(lsyms)
            vi += dent
    P = np.vstack([np.zeros((1, 3))] + [np.array(placed[1:], float)])
    # --- Hebel A: rigid η-axis rotation of the non-primary-η atoms --------------
    # Rotate the CO-tripod / co-ligands / a secondary ring as ONE rigid body about
    # the PRIMARY η-face M→centroid axis (through the metal at the origin).  This is
    # the missing relative DOF for piano-stool / sandwich systems: the η-ring sits
    # right but the tripod/2nd-ring CLOCK is not yet at the crystal rotamer.  Because
    # the axis passes through the metal, every atom's metal distance is preserved
    # (M-centroid + M-D INVARIANT); only the azimuth changes.  Strictly additive: the
    # canonical build (axis_rot==0) is byte-identical (this block is skipped).
    if abs(_axis_rot) > 1e-9 and eta_blocks:
        prim_V = np.asarray(eta_blocks[0][2], float)
        nV = np.linalg.norm(prim_V)
        if nV > 1e-9:
            prim_V = prim_V / nV
            Rax = _rot_about_axis(prim_V, _axis_rot)
            rot_blocks = list(sigma_blocks) + list(nonprimary_eta_blocks)
            for (s, e) in rot_blocks:
                P[s:e] = P[s:e] @ Rax.T            # metal at origin -> no recentre
            if not np.all(np.isfinite(P)):
                return None
    # --- Piano-stool σ-leg TILT (DELFIN_FFFREE_PIANO_LEG_TILT, default OFF) -------
    # A half-sandwich (EXACTLY one η-face + σ-legs, no second ring) is not a regular
    # polyhedron: the η-ring is a fat FACE donor and the legs splay DOWN, so the real
    # (ring-centroid)-M-L angle is ~125 deg, not the ~110 deg (T-4) / ~90 deg (SP-4)
    # the generic CNn polyhedron places them at.  Rigidly tilt each σ-leg block about
    # the metal so that angle is correct; the η-ring + every M-L distance are invariant.
    # Skipped (byte-identical) when the flag is off, when there is not exactly one
    # η-face (full sandwich / metallocene untouched), or when there are no σ-legs.
    if (_piano_leg_tilt_enabled() and len(eta_blocks) == 1 and sigma_blocks):
        try:
            prim_V = np.asarray(eta_blocks[0][2], float)
            P = _tilt_piano_legs(P, prim_V, sigma_blocks, _piano_leg_target_deg())
        except Exception:
            pass
        if not np.all(np.isfinite(P)):
            return None
    # FF-free geometric clash-relief (η-ring + σ-donors all frozen so the rigid
    # ring + constructed coordination are preserved; periphery relaxes only).
    P = _refine_guarded(out_syms, P, fixed)
    if not np.all(np.isfinite(P)):
        return None
    return out_syms, P, sorted(set(donors)), exempt_pairs


# ---------------------------------------------------------------------------
# Rigid-hapto ENSEMBLE (deterministic, RMSD-deduplicated) — used by the FF-free
# converter backend when DELFIN_FFFREE_RIGID_HAPTO=1.  Every member is a fully
# rigid build (no ring collapse); the variation is purely geometric + enumerated
# deterministically, so two runs of the same SMILES emit byte-identical frames.
# ---------------------------------------------------------------------------

# Symmetry-distinct in-plane spin counts per hapticity for a BARE ring (the ring
# carbons map onto themselves under the cyclic group C_n; a SUBSTITUTED ring breaks
# that symmetry, so spinning to n_fold orientations samples the substituent rotamers
# without emitting bare-ring duplicates — the RMSD dedup removes any that coincide).
_ETA_NFOLD = {2: 2, 3: 3, 4: 4, 5: 5, 6: 6}


def _ring_spins(eta_n):
    """Deterministic, symmetry-reduced list of η-ring spin angles (rad), starting
    at the canonical 0.  For an η-n face the carbons recur every 2π/n; we sample one
    full inter-carbon sector in n_fold steps (so a substituted ring visits each
    distinct substituent-vs-coligand clock position exactly once).  Duplicate frames
    (symmetric rings) are removed later by the RMSD dedup."""
    nf = _ETA_NFOLD.get(int(eta_n), max(2, int(eta_n)))
    sector = 2.0 * np.pi / float(nf)
    # sample the sector at nf sub-steps -> the substituent sweeps a full inter-carbon
    # gap; canonical 0 first so member[0] is the historic single build.
    return [sector * (k / float(nf)) for k in range(nf)]


def _axis_rot_angles(n_rot, eta_n):
    """Deterministic η-axis rotation angles (rad) for Hebel A: rotate the non-η
    block (CO-tripod / co-ligands / 2nd ring) about the primary η M→centroid axis.

    The canonical 0 is NOT included here (the base ensemble already carries it as
    axis_rot==0); these are the ADDED frames.  We sweep the full 0–2π in ``n_rot``
    equal steps and drop the 0 step, so the extra frames are a uniform azimuthal
    sweep of the ring-vs-rest clock.  Bare-ring / symmetric duplicates are removed
    downstream by the RMSD dedup; the sweep is intentionally NOT symmetry-reduced
    here because the relevant symmetry is the PRODUCT of the ring C_n and the tripod
    C_m (system-dependent), and over-reducing would skip the crystal rotamer."""
    n_rot = max(2, int(n_rot))
    step = 2.0 * np.pi / float(n_rot)
    return [step * k for k in range(1, n_rot)]


def _slip_modes(lg):
    """Deterministic ring-SLIP isomers for one η-ligand, gated by valence sanity.
    Returns a list of (eta_idxs, eta_n) the ring may bind through, lowest-disruption
    first and ALWAYS including the full face.  An η-n face can ring-slip to a
    CONTIGUOUS sub-arc of the ring (η6→η4→η2 for arenes/dienes, η5→η3→η1 for Cp/allyl)
    — the canonical odd/even hapticity ladder.  We only emit a slip when the sub-arc
    atoms are mutually contiguous through ring bonds (graph-checked), so no impossible
    mode is invented; substituents stay rigidly attached (whole ligand still rigid)."""
    full = list(lg["eta_local_idxs"])
    n = len(full)
    modes = [(full, n)]
    if n < 3:
        return modes
    mol = lg["mol"]
    # ring-bond adjacency among the η carbons (graph-only)
    fs = set(full)
    adj = {i: [] for i in full}
    for b in mol.GetBonds():
        a1, a2 = b.GetBeginAtomIdx(), b.GetEndAtomIdx()
        if a1 in fs and a2 in fs:
            adj[a1].append(a2)
            adj[a2].append(a1)
    # canonical ring traversal order (follow adjacency from the lowest index)
    seq = [full[0]]
    seen = {full[0]}
    while len(seq) < n:
        nxt = next((x for x in adj.get(seq[-1], []) if x not in seen), None)
        if nxt is None:
            return modes                      # not a simple ring path -> no slip
        seq.append(nxt)
        seen.add(nxt)
    # hapticity ladder: drop 2 carbons at a time (one from each end of the arc) so the
    # bound arc stays contiguous + centred; preserve odd/even parity (η6→η4→η2,
    # η5→η3→η1).  Stop at η1 (single σ-C, no longer a face — handled as edge case).
    k = n - 2
    while k >= 1:
        drop = (n - k) // 2
        arc = seq[drop: drop + k] if k > 1 else [seq[n // 2]]
        if len(arc) == k and k >= 1:
            modes.append((sorted(arc), k))
        k -= 2
    return modes


def _cp_pucker_amps(eta_n, aromatic):
    """Cremer-Pople-style out-of-plane pucker amplitudes (Å) for a NON-aromatic
    η-ring (diene / allyl / cyclohexadienyl).  Aromatic faces (Cp, arene) are planar
    -> no pucker (amp 0 only).  Returns [0.0] for planar faces (single conformer) and
    a small symmetric set for puckerable rings."""
    if aromatic or int(eta_n) < 4:
        return [0.0]
    return [0.0, 0.15, -0.15]


def _eta_is_aromatic(lg):
    """True if every η-ring atom is aromatic (planar face: Cp / arene)."""
    mol = lg["mol"]
    try:
        return all(mol.GetAtomWithIdx(int(j)).GetIsAromatic()
                   for j in lg["eta_local_idxs"])
    except Exception:
        return True


def _rmsd_aligned(A, B):
    """Heavy-rigid coordinate RMSD between two SAME-ORDERING coordinate sets after a
    proper-rotation Kabsch superposition (both already metal-centred at origin).
    Deterministic; used only to dedup ensemble members of ONE complex."""
    A = np.asarray(A, float); B = np.asarray(B, float)
    if A.shape != B.shape or A.shape[0] < 3:
        return 1e9
    Ac = A - A.mean(0); Bc = B - B.mean(0)
    R = _kabsch_rot(Ac, Bc)
    D = Ac @ R.T - Bc
    return float(np.sqrt((D * D).sum() / len(A)))


def _hapto_fold_rings(d, syms):
    """Global ring index lists for the eta builds -- the offset arithmetic of
    ``assemble_hapto`` RETRACED and then CHECKED against reality.

    ``assemble_hapto`` returns only ``(syms, P, donors, exempt)``; the offsets
    stay local there.  They are therefore recomputed exactly as they arise there:
      * emission order is NOT the list order of ``d["ligands"]``,
        but eta ligands first, then the rest, each ascending (line 3621).
      * block size is ``Chem.AddHs(lg["mol"]).GetNumAtoms()``, NOT
        ``lg["mol"].GetNumAtoms()`` -- the ring H exist in ``lg["mol"]`` only as a
        NumExplicitHs property and become atoms only through AddHs.
      * index 0 is the metal, so the first ligand starts at 1.

    ⚠ RECOMPUTING IS AN ASSUMPTION, SO IT IS CHECKED.  If the total number
      of atoms does not match, there is NO verdict (``None``) and the predicate falls
      back to the pure RMSD -- exactly the case the kekulize detour
      (DELFIN_FFFREE_KEKULIZE_SPLIT) can produce when it neutralises an aromatic N+
      and thereby changes the H count.  ``_fold_rings_from_blocks``
      afterwards still checks every single ring atom against its element symbol."""
    try:
        ligs = d["ligands"]
        order = sorted(range(len(ligs)),
                       key=lambda i: (0 if ligs[i].get("is_eta") else 1, i))
        blocks = []
        off = 1                                    # index 0 is the metal
        for i in order:
            m = Chem.AddHs(ligs[i]["mol"])
            # third entry = the chelate arms for the metallacycle fingerprint
            # (``None`` for eta and for denticity 1 -> no ring).  MC switch OFF
            # -> `_fold_rings_with_mc` never reads it.
            blocks.append((off, m, _fold_mc_arms(ligs[i])))
            off += m.GetNumAtoms()
        if off != len(syms):
            return None                            # offset model does not fit: no verdict
        return _fold_rings_with_mc(blocks, syms, 0)
    except Exception:
        return None


def _dedup_builds(builds, rmsd_tol=0.25, fold_rings=None):
    """RMSD-deduplicate a list of (syms, P, donors, exempt) builds (same complex, so
    identical atom ordering).  Keeps the FIRST occurrence (emission order = canonical
    build first), dropping any later build within ``rmsd_tol`` Å of a kept one.
    Deterministic.  Distinct-by-atom-count builds (ring-slip changes nothing in the
    atom list, so counts always match) are compared directly.

    ⚠ WHY A FOLD FINGERPRINT IS NEEDED HERE.  The variant list that runs into
      this function contains ``_cp_pucker_amps``: the Cremer-Pople fold
      of the eta face with an amplitude of only +/-0.15 A, compared over ALL
      atoms of the complex.  That lies safely below ``rmsd_tol`` = 0.25 -- the
      eta fold axis was completely eaten at this site, and both
      callers (RIGID_HAPTO, HAPTO_AXIS_ROT) are champion switches.
      ``fold_rings=None`` (default) -> predicate byte-identical to before."""
    kept = []
    fps = []                                       # fingerprint per kept build
    for b in builds:
        P = b[1]
        dup = False
        fp = None                                  # lazy: only on RMSD proximity
        for ki, kb in enumerate(kept):
            if kb[1].shape == P.shape and _rmsd_aligned(kb[1], P) < rmsd_tol:
                if fold_rings is None:
                    dup = True                     # switch OFF -> old predicate
                    break
                if fp is None:
                    fp = _fold_fp(P, fold_rings)
                if fps[ki] is None:
                    fps[ki] = _fold_fp(kb[1], fold_rings)
                if _fold_same(fps[ki], fp):
                    dup = True
                    break
        if not dup:
            kept.append(b)
            fps.append(None)
    return kept


def _hapto_base_variants(d):
    """Deterministic base variant grid for the rigid-hapto ensemble: η-ring rotamers
    (symmetry-reduced spin), valence-gated η/σ ring-slip isomers, Cremer-Pople pucker
    (non-aromatic faces).  Returns ([variant, ...], eta_lis) with the canonical
    all-base variant first, or (None, None) if the decompose dict carries no η-face.
    Graph-only, deterministic; factored out so Hebel A can reuse the same base grid."""
    ligands = d["ligands"]
    eta_lis = [i for i, lg in enumerate(ligands) if lg.get("is_eta")]
    if not eta_lis:
        return None, None
    per_eta_choices = {}
    for li in eta_lis:
        lg = ligands[li]
        aromatic = _eta_is_aromatic(lg)
        slips = _slip_modes(lg)                       # [(eta_idxs, eta_n), ...]
        choices = []
        for (eidx, en) in slips:
            spins = _ring_spins(en)
            for sp in spins:
                for pk in _cp_pucker_amps(en, aromatic):
                    choices.append({"spin": sp, "pucker": pk,
                                    "eta_idxs": eidx, "eta_n": en})
        per_eta_choices[li] = choices
    # Vary ONE η-ligand at a time off the canonical base; canonical all-base leads.
    base_choice = {li: per_eta_choices[li][0] for li in eta_lis}
    variants = [{"eta": dict(base_choice)}]           # member 0 = canonical
    for li in eta_lis:
        for ch in per_eta_choices[li][1:]:
            ev = dict(base_choice)
            ev[li] = ch
            variants.append({"eta": ev})
    return variants, eta_lis


def assemble_hapto_ensemble(metal, geometry, d, max_builds=30):
    """Deterministic RIGID-hapto ENSEMBLE: enumerate η-ring rotamers (symmetry-
    reduced spin), η/σ ring-slip isomers (valence-gated), and Cremer-Pople ring
    pucker (non-aromatic faces), assembling each as a fully rigid build (no collapse),
    then RMSD-deduplicate.  Returns a list of (syms, P, donors, exempt_pairs) builds
    (canonical build first), or None on total failure.

    Universal, graph-only, deterministic (fixed seeds; canonical ordering); the
    canonical build (member 0) is byte-identical to ``assemble_hapto(...)``."""
    variants, eta_lis = _hapto_base_variants(d)
    if variants is None:
        return None

    builds = []
    for var in variants:
        if len(builds) >= max_builds * 3:             # generous cap before dedup
            break
        try:
            b = assemble_hapto(metal, geometry, d, variant=var)
        except Exception:
            b = None
        if b is None:
            continue
        builds.append(b)
    if not builds:
        return None
    # Default OFF -> `fold_rings` stays None -> dedup byte-identical.
    builds = _dedup_builds(
        builds,
        fold_rings=(_hapto_fold_rings(d, builds[0][0]) if _fold_fp_enabled() else None))
    return builds[:max_builds] if builds else None


def assemble_hapto_axis_rotants(metal, geometry, d, n_axis=8, max_builds=60):
    """Hebel A (DELFIN_FFFREE_HAPTO_AXIS_ROT) — STRICTLY ADDITIVE η-axis rotamers.

    For every base variant (same grid as ``assemble_hapto_ensemble``) emit the rigid
    builds that additionally rotate the non-primary-η block (CO-tripod / co-ligands /
    a 2nd ring) about the PRIMARY η M→centroid axis by a uniform 0–2π sweep (the
    canonical 0 step is EXCLUDED — those are exactly the base ensemble builds, which
    the caller emits separately).  The rotation axis runs through the metal at the
    origin, so M-centroid + every M-D distance is INVARIANT and the η-ring internal
    geometry is untouched; only the ring-vs-rest azimuthal clock changes.

    Returns a list of (syms, P, donors, exempt_pairs) builds (NEW orientations only),
    or [] if no η-face / no rotatable block.  RMSD-deduplicated.  Deterministic."""
    variants, eta_lis = _hapto_base_variants(d)
    if variants is None:
        return []
    prim_n = int(d["ligands"][eta_lis[0]].get("eta_n", 6))
    angles = _axis_rot_angles(n_axis, prim_n)
    builds = []
    for var in variants:
        if len(builds) >= max_builds * 3:
            break
        for a in angles:
            v = {"eta": dict(var.get("eta", {})), "axis_rot": float(a)}
            try:
                b = assemble_hapto(metal, geometry, d, variant=v)
            except Exception:
                b = None
            if b is not None:
                builds.append(b)
            if len(builds) >= max_builds * 3:
                break
    if not builds:
        return []
    # Default OFF -> `fold_rings` stays None -> dedup byte-identical.
    builds = _dedup_builds(
        builds,
        fold_rings=(_hapto_fold_rings(d, builds[0][0]) if _fold_fp_enabled() else None))
    return builds[:max_builds]
