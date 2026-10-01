#!/usr/bin/env python3
"""assemble_complex.py — assemble a full 3D TMC from a metal + ligand(s) +
geometry, by orienting each ligand so its donor sits on the placed polyhedron
vertex with its lone pair pointing at the metal.  Metal-FF-free: the sphere is
geometric (metal_sphere_builder); ligands are rigidly oriented; (constrained MMFF
relax of the organic periphery with the core frozen = next step).

This ties the pieces together: enumeration (which donors / which vertex) +
sphere placement (where) + this (orient + merge) -> a DFT-startable structure.

Prototype: homoleptic monodentate [M(L)n].  Deterministic.
"""
from __future__ import annotations
import math
import os
import itertools
from typing import List, Tuple
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem
from delfin.manta import metal_sphere_builder as MSB
import delfin.manta._bond_decollapse as _bd

from delfin.manta.assemble_donor_plane import (
    _BETA_BAND,
    _COMBO_MATERIALISE_MAX,
    _DP_RESID_SLACK,
    _DP_STEPS,
    _H_CONTACT_FLOOR,
    _TRILAT_RESCUE,
    _beta_score,
    _bite_aware_targets,
    _canonical_arm_order,
    _collapsed_heavy_bonds_strict,
    _donor_follow_weights,
    _donor_plane_beta,
    _donor_plane_relax,
    _ffree_flag,
    _has_collapsed_heavy_bonds,
    _hydrogens_riding_on,
    _ranked_combos,
    _riders_that_may_move,
    _rigid_planar_central_arm,
    _trilat_targets_on,
    _trilaterate_donor_targets,
    trilat_rescue_enabled,
    trilaterate_rescue,
)  # noqa: F401

from delfin.manta.assemble_ligand_embed import (
    SEED,
    _POLY_FIDELITY_BIN,
    _bite_free_bounds,
    _coplanar_metal_centered_conformer,
    _embed_metallacycle,
    _kabsch_rot,
    _ligand_3d,
    _planar_mer_cn5_enabled,
    _planar_polydentate_place_enabled,
    _rigid_ligand_cavity_conformer,
    _ring_bounds_enabled,
    _rot_align,
    _tighten_ring_bounds,
)  # noqa: F401

from delfin.manta.assemble_orient import (
    _CHALCOGENS,
    _DONOR_BEND_DEG,
    _PNICTOGENS,
    _axis_rot,
    _diatomic_donor_partner,
    _diatomic_orient_enabled,
    _donor_and_lp,
    _donor_bend_angle,
    _donor_c_angle,
    _orient_chelate_to_vertices,
    _orient_diatomic_block,
    _straighten_sp_chain,
    _subtree,
    _vsepr_reconstruct,
)  # noqa: F401

from delfin.manta.assemble_chelate import (
    _BITE_MAX_RING,
    _chelate_ring_size,
    _lone_pair_dir,
    _lp_aligned_angle,
    _lp_orient_seated_bidentate,
    _measured_bite,
    _min_interligand,
    _place_chelate_block,
    assemble_chelate,
    assemble_monodentate,
    assemble_multichelate,
    bite_close,
    clash_relief,
)  # noqa: F401

from delfin.manta.assemble_ligand_confs import (
    _CONF_CACHE,
    _CONF_CACHE_MAX,
    _clash_count,
    _joint_declash_frame,
    _ligand_3d_from_mol,
    _ligand_block_bonds,
    _ligand_confs_from_mol,
    _refine_guarded,
    _relax_confs_ffree,
    _sphere_flex_frame,
    _symmetrize_degenerate,
    _torsion_relax_frame,
)  # noqa: F401

from delfin.manta.assemble_fold_fp import (
    _FOLD_FP_MAXRING,
    _FOLD_FP_MCMAX,
    _FOLD_FP_QMIN,
    _FOLD_FP_RINGMAX,
    _complex_rmsd,
    _fold_fp,
    _fold_fp_enabled,
    _fold_fp_mc_enabled,
    _fold_mc_arms,
    _fold_mc_rings,
    _fold_rings_from_blocks,
    _fold_rings_with_mc,
    _fold_same,
    _lig_path,
)  # noqa: F401


def assemble_heteroleptic_from_mols(metal: str, geometry: str, vertex_specs,
                                    refine: bool = True):
    """vertex_specs[i] = (frag_mol, donor_local_idx).  Takes ligand MOLS directly
    (preserves donor index).  refine=True runs the geometry refiner (deterministic
    coordinate descent, metal+donors frozen) to remove rigid-placement clashes."""
    ref = MSB._ref_vectors(geometry)
    if len(vertex_specs) != len(ref):
        raise ValueError("vertex_specs count != vertices")
    out_syms = [metal]; blocks = [np.zeros((1, 3))]
    placed_P = [np.zeros(3)]; placed_syms = [metal]
    fixed = {0}                       # metal frozen
    pos = 1
    block_specs = []                  # (offset, lmol, donor_local) for the torsion relax
    for i, (frag, di) in enumerate(vertex_specs):
        Vunit = ref[i] / np.linalg.norm(ref[i])
        confs = _ligand_confs_from_mol(frag)
        if confs is None:
            return None
        lsyms, coords_list, lmol = confs
        md = MSB.md_distance(metal, lsyms[di],
                             atom=lmol.GetAtomWithIdx(di), mol=lmol)
        vertex = Vunit * md
        # UNIVERSAL: pick the conformer whose placement clashes least with the
        # metal + already-placed ligands (defect-count-driven conformer selection).
        best_Q, best_clash = None, 1e18
        for lP in coords_list:
            if len(lsyms) == 1:
                Q = vertex.reshape(1, 3)
            else:
                lPv, lp = _vsepr_reconstruct(lsyms, lP, lmol, di)
                R = _rot_align(lp, -Vunit)
                Q = (lPv - lPv[di]) @ R.T + vertex
                # diatomic-orientation guard (env-gated, default OFF -> byte-identical):
                # keep the SMILES-bonded donor atom (not its partner) at the vertex.
                if _diatomic_orient_enabled():
                    _dp = _diatomic_donor_partner(
                        {"mol": frag, "donor_local_idxs": [di], "denticity": 1})
                    if _dp is not None:
                        Q = _orient_diatomic_block(
                            Q, lsyms, _dp[0], _dp[1], np.zeros(3), vertex)
            cl = _clash_count(Q, np.array(placed_P), lsyms, placed_syms)
            # THE FREE DEGREE OF FREEDOM, CALL SITE 2 of 2 (purely monodentate path).
            # Here the azimuth about M-D is the only free DOF, and on this path it is
            # sampled NOWHERE -- the existing axis rotation (CN2_SPINS)
            # sits in the ensemble sibling assemble_heteroleptic_ensemble, not
            # here.  Default OFF -> None -> byte-identical.
            _fd = _free_dof_reseat(Q, lsyms, [di], np.array(placed_P), placed_syms, cl)
            if _fd is not None:
                Q, cl = _fd
            if cl < best_clash:
                best_clash, best_Q = cl, Q
            if cl == 0:
                break
        block_specs.append((pos, lmol, di))
        out_syms += lsyms; blocks.append(best_Q)
        for row in best_Q:
            placed_P.append(row)
        placed_syms += lsyms
        if len(lsyms) == 1:
            fixed.add(pos); pos += 1
        else:
            fixed.add(pos + di); pos += len(lsyms)
    P = np.vstack(blocks)
    if refine:
        # ===== THE SEATING HANDS THE RELAXATION ITS ASSERTION =====================
        # Measured 2026-08-18 on 30921 systems: net 988 systems flow from the
        # octahedron into the trigonal prism (McNemar X2 = 860.8), the builder produces
        # 2.99 times too many prisms, while every other shape stays between 0.86 and
        # 1.29.  But they are NOT prisms -- the CShM mass lies unimodally at 8
        # to 12 instead of at 16.7, i.e. a HALF Bailar twist.  And it is not a
        # selection: poly_match is false in 1061 of 1061 cases, although the eye reads
        # the best frame over the WHOLE manifold -- in the whole manifold there is
        # no octahedron.  The signal is the interlocking, monotone from 2.49 % at zero
        # chelate rings to 16.12 % at five; metal and d-count are flat.
        #
        # ⇒ The seating sets the polyhedron correctly, and the chelate pull twists it
        # out afterwards.  That is exactly what the assertion protocol is built for: the
        # construction says WHAT it claims, and the relaxation must not break it.
        #
        # The existing countermeasure in the source (smiles_converter.py:38074, the
        # UFF angle targets on the opposing donor pairs) hangs via :27451
        # on apply_uff and lies behind the FF-free return -- it NEVER ran on this
        # path.  This here is its FF-free counterpart, and it does not correct,
        # it FORBIDS: if the relaxation breaks the assertion, the frame BEFORE the
        # relaxation stands.  Never-worse by construction, no target value, no threshold
        # that could be fine-tuned to a pool.
        #
        # ⚠ WHY THIS IS NOT A REPAIRER.  The module census of 18.08. measured
        # that 20 of 21 repairers do nothing anyway and the one remaining
        # (the unconditional post-relaxation) is LOAD-BEARING -- without it every axis
        # gets worse (uffoffE).  So the answer is not "no optimisation",
        # but "no BLIND optimisation".  This block takes nothing away from the
        # relaxation; it only gives it what it did not know before.
        #
        # DELFIN_FFFREE_ASSERT_ENFORCE (default 0 -> byte-identical).
        P = _refine_guarded(out_syms, P, fixed)
        # #308 whole-complex torsion-space clash relax (env-gated, default-OFF
        # byte-id): when rigid M-D-axis selection is not enough and ligand-internal
        # rotation is needed, jointly optimise all rotatable single bonds of the
        # assembled complex (metal + donors = `fixed` frozen).  Torsion-only -> bond
        # lengths/angles preserved exactly; never-worse.  No-op when flag unset.
        P = _torsion_relax_frame(out_syms, P, fixed, block_specs)
        # JOINT inter-ligand declash (env-gated, default-OFF byte-id): whole-ligand
        # M-D-axis rotations + internal torsions jointly minimising the GLOBAL
        # inter-ligand heavy-heavy clash (the self-gate blocker for class-B), core
        # frozen.  Runs after #308; no-op when DELFIN_FFFREE_JOINT_DECLASH unset.
        P = _joint_declash_frame(out_syms, P, fixed, block_specs, geom=geometry)
        # SOFT coordination-sphere flex (env-gated, default-OFF byte-id): let the
        # donors breathe a hard-bounded amount to open the residual MILD inter-ligand
        # heavy-heavy clashes that rotation/frozen-donor refine cannot (clash forensik
        # 2026-06-29).  Restored toward the ideal vertex+M-D; never-worse on clash.
        P = _sphere_flex_frame(out_syms, P, fixed, block_specs)
    return out_syms, P


def assemble_heteroleptic_ensemble(metal: str, geometry: str, vertex_specs,
                                   n_frames: int = 6, rmsd_dedup: float = 0.5,
                                   per_lig_confs: int = 6, refine: bool = True):
    """Ensemble variant of ``assemble_heteroleptic_from_mols``: instead of keeping
    only the single clash-minimal conformer per ligand, emit a small RMSD-deduped
    ENSEMBLE of full-complex frames that vary ligand internal conformation +
    substituent/co-ligand orientation while keeping the ideal coordination core
    (vertex directions + M-D distances) FIXED.

    This reuses the SAME FF-free Layer-3 machinery as the single path —
    ``_ligand_confs_from_mol`` (deterministic ETKDG conformer pool + MMFF),
    ``_vsepr_reconstruct``/``_rot_align`` orientation, ``_clash_count`` scoring, and
    the geometric ``refine`` — only it RETAINS several distinct low-clash conformer
    *combinations* across ligands rather than collapsing to one.  Best-of-ensemble
    crystal-recall then gets a comparable conformer spray to the legacy multi-frame
    path.  Deterministic (fixed SEED, single-thread), graph-only, never non-finite.

    Returns a list of (syms, P) frames (>=1) or None on failure.  The first frame is
    byte-identical to ``assemble_heteroleptic_from_mols`` (the clash-minimal pick).
    """
    ref = MSB._ref_vectors(geometry)
    if len(vertex_specs) != len(ref):
        raise ValueError("vertex_specs count != vertices")
    # Per-vertex candidate placements: for each ligand build its conformer pool and
    # orient every conformer onto the vertex; keep the distinct low-clash candidates
    # (clash vs the metal only -> placement-order-independent + diverse).  An
    # all-monatomic ligand set has a single rigid placement (-> 1 frame).
    out_syms = [metal]
    fixed = {0}
    per_vertex_cands = []     # per vertex: list of (Q, internal_clash) sorted best-first
    vertex_lsyms = []
    block_specs = []          # (offset, lmol, donor_local) for the #308 torsion relax
    pos = 1
    metal_P = np.zeros((1, 3))
    metal_sym = [metal]
    for i, (frag, di) in enumerate(vertex_specs):
        Vunit = ref[i] / np.linalg.norm(ref[i])
        confs = _ligand_confs_from_mol(frag, k=max(per_lig_confs, 1))
        if confs is None:
            return None
        lsyms, coords_list, lmol = confs
        block_specs.append((pos, lmol, di))
        md = MSB.md_distance(metal, lsyms[di],
                             atom=lmol.GetAtomWithIdx(di), mol=lmol)
        vertex = Vunit * md
        cand = []                                   # (Q, clash_vs_metal)
        seen_local = []
        for lP in coords_list:
            if len(lsyms) == 1:
                Q = vertex.reshape(1, 3)
            else:
                lPv, lp = _vsepr_reconstruct(lsyms, lP, lmol, di)
                R = _rot_align(lp, -Vunit)
                Q = (lPv - lPv[di]) @ R.T + vertex
                # diatomic-orientation guard (env-gated, default OFF -> byte-identical)
                if _diatomic_orient_enabled():
                    _dp = _diatomic_donor_partner(
                        {"mol": frag, "donor_local_idxs": [di], "denticity": 1})
                    if _dp is not None:
                        Q = _orient_diatomic_block(
                            Q, lsyms, _dp[0], _dp[1], np.zeros(3), vertex)
            if not np.all(np.isfinite(Q)):
                continue
            cl = _clash_count(Q, metal_P, lsyms, metal_sym)
            # dedup conformers of THIS ligand by intra-ligand RMSD (identity corr.)
            dup = False
            for Qs in seen_local:
                if Qs.shape == Q.shape:
                    c0 = Q - Q.mean(axis=0); c1 = Qs - Qs.mean(axis=0)
                    if float(np.sqrt(((c0 - c1) ** 2).sum(axis=1).mean())) < rmsd_dedup:
                        dup = True
                        break
            if dup:
                continue
            seen_local.append(Q)
            cand.append((Q, cl))
        if not cand:
            return None
        # Axial-spin rotamers (substituent/co-ligand orientation): for a near-rigid
        # donor (e.g. =S/=P thiourea/phosphine CN2 ligands) the dominant pose DOF vs
        # the crystal is the ROTATION of the whole ligand about the M-donor axis, which
        # the ETKDG internal-conformer pool does not sample.  Spin each candidate block
        # about its M-D axis (= Vunit through the donor) at deterministic increments,
        # keeping distinct low-clash variants.  Pure geometry (the linear core/M-D stay
        # fixed -> only orientation changes).  Opt-out via DELFIN_FFFREE_CN2_NOSPIN=1.
        if len(lsyms) > 1 and os.environ.get("DELFIN_FFFREE_CN2_NOSPIN", "0") != "1":
            n_spin = int(os.environ.get("DELFIN_FFFREE_CN2_SPINS", "6"))
            spun = []
            for Q0, _cl0 in list(cand):
                for s in range(1, max(n_spin, 1)):
                    ang = 2.0 * np.pi * s / max(n_spin, 1)
                    Rs = _axis_rot(Vunit, ang)
                    Qs = (Q0 - vertex) @ Rs.T + vertex   # spin about the donor vertex
                    if not np.all(np.isfinite(Qs)):
                        continue
                    cls = _clash_count(Qs, metal_P, lsyms, metal_sym)
                    dup = False
                    for Qe, _ in (cand + spun):
                        if Qe.shape == Qs.shape:
                            c0 = Qs - Qs.mean(axis=0); c1 = Qe - Qe.mean(axis=0)
                            if float(np.sqrt(((c0 - c1) ** 2).sum(axis=1).mean())) < rmsd_dedup:
                                dup = True
                                break
                    if not dup:
                        spun.append((Qs, cls))
            cand += spun
        cand.sort(key=lambda t: t[1])               # low clash-vs-metal first
        per_vertex_cands.append(cand)
        vertex_lsyms.append(lsyms)
        out_syms += lsyms
        if len(lsyms) == 1:
            fixed.add(pos); pos += 1
        else:
            fixed.add(pos + di); pos += len(lsyms)

    # Enumerate full-complex frames over the Cartesian product of per-vertex
    # conformer choices, greedily ordered (sum of candidate ranks -> the best
    # combinations first), score each by total inter-ligand clash, RMSD-dedup at the
    # complex level, keep up to n_frames.  Capped product keeps it bounded+fast.
    import itertools as _it
    rank_lists = [list(range(len(c))) for c in per_vertex_cands]
    _nprod = 1
    for _rl in rank_lists:
        _nprod *= max(1, len(_rl))
        if _nprod > _COMBO_MATERIALISE_MAX:
            break
    if _nprod > _COMBO_MATERIALISE_MAX:      # never materialise a 10^8 product to use 64
        combos = _ranked_combos(rank_lists, 256)
    else:
        combos = list(_it.product(*rank_lists))
    # order combos by total rank (frame 0 = all best = the single-path pick), then
    # lexicographically -> deterministic.
    combos.sort(key=lambda cb: (sum(cb), cb))

    # ---- inter-ligand-clash-aware selection (env-gated, default OFF = byte-id) ----
    # ROOT FIX (#306 inter-ligand): the per-vertex candidate scoring above clashes vs
    # the METAL ONLY, and this combo loop orders by metal-clash RANK SUM and NEVER
    # scores the assembled frame for inter-ligand overlap.  For homoleptic bulky
    # ligands every vertex then independently picks the same metal-optimal pose ->
    # systematic inter-ligand clash, identical in every frame (e.g. Zr(C6H5)6: min
    # inter-aryl H-H = 0.73 A in all frames).  The axial-spin candidates that COULD
    # interleave the ligands already exist in per_vertex_cands -- they are just never
    # selected for.  When DELFIN_FFFREE_INTERLIG_RANK=1: (a) lead frame 0 with a
    # GREEDY sequential placement that picks, at each vertex, the candidate minimising
    # _clash_count vs the already-placed ligands (mirrors assemble_heteroleptic_from_
    # mols), and (b) emit the remaining frames in ASCENDING total inter-ligand clash.
    # RIGID only: we re-select among existing conformer/spin candidates; the core
    # (vertex directions + M-D distances) stays fixed.  Deterministic, geometry-only.
    _interlig = os.environ.get("DELFIN_FFFREE_INTERLIG_RANK", "0") == "1"
    MAX_EVAL = 64                                    # bound the build work

    def _interlig_clash(cb):
        """Total inter-ligand _clash_count for combo cb: sum over each ligand block
        of its clash vs the union of PREVIOUSLY-placed ligand blocks (metal at origin
        excluded; same H-inclusive measure used everywhere).  Rigid-block geometry."""
        placed = []                       # accumulated already-placed ligand atoms
        placed_syms = []
        total = 0
        for vi, ci in enumerate(cb):
            Q = per_vertex_cands[vi][ci][0]
            lsyms = vertex_lsyms[vi]
            if placed:
                total += _clash_count(Q, np.array(placed), lsyms, placed_syms)
            for row in Q:
                placed.append(row)
            placed_syms += lsyms
        return total

    if _interlig:
        MAX_EVAL = 128                               # modest, capped bump when ON
        # (a) greedy inter-ligand-minimal lead combo: fixed deterministic vertex order,
        #     each vertex takes the candidate minimising clash vs already-placed blocks.
        greedy = []
        gplaced = []; gplaced_syms = []
        for vi in range(len(per_vertex_cands)):
            cands_vi = per_vertex_cands[vi]
            lsyms = vertex_lsyms[vi]
            best_ci, best_cl = 0, None
            for ci, (Q, _mc) in enumerate(cands_vi):
                if not gplaced:
                    cl = 0
                else:
                    cl = _clash_count(Q, np.array(gplaced), lsyms, gplaced_syms)
                # strictly-less keeps the FIRST (lowest-index = lowest metal-clash)
                # candidate on ties -> deterministic, and =all-best when no clash.
                if best_cl is None or cl < best_cl:
                    best_cl, best_ci = cl, ci
                if cl == 0:
                    break
            greedy.append(best_ci)
            Qsel = cands_vi[best_ci][0]
            for row in Qsel:
                gplaced.append(row)
            gplaced_syms += lsyms
        greedy = tuple(greedy)
        # rank the bounded combo window by total inter-ligand clash (ascending), then
        # by the original metal-rank ordering as a stable deterministic tiebreak; put
        # the greedy lead combo first (dedup so it is not double-evaluated).
        window = combos[:MAX_EVAL]
        if greedy not in window:
            window = [greedy] + window[:MAX_EVAL - 1]
        scored = [(_interlig_clash(cb), sum(cb), cb) for cb in window]
        # greedy combo forced to the front via a sentinel key; rest by (clash, rank).
        scored.sort(key=lambda t: (0 if t[2] == greedy else 1, t[0], t[1], t[2]))
        eval_order = [t[2] for t in scored]
    else:
        eval_order = combos[:MAX_EVAL]

    frames = []                                      # (syms, P) kept (deduped)
    # FOLD FINGERPRINT (see block at `_fold_fp_enabled`).  Default OFF ->
    # `_fold_rings` stays None -> the predicate below is literally the old one.
    # `block_specs` here already carries (global offset, lmol, donor-local).
    #
    # ⛔ THE METALLACYCLE FINGERPRINT (DELFIN_FFFREE_DEDUP_FOLD_FP_MC) IS
    #    DELIBERATELY NOT WIRED HERE, and the reason is structural, not caution:
    #    the metal IS at index 0 (`out_syms = [metal]` :3262), so it is not
    #    missing.  What is missing is the SECOND donor.  `vertex_specs` is a sequence
    #    of `(frag, di)` with EXACTLY ONE donor index per vertex (:3270/:3276), and
    #    both callers build it from `lig_ref[lab] = (lg["mol"],
    #    lg["donor_local_idx"])` -- singular -- in the MONODENTATE branch of
    #    `converter_backend` (:2615, calls :2727/:3050).  Chelates do not reach
    #    this branch, they go to `assemble_from_config` beforehand.
    #    A chelate ring needs two donors FROM THE SAME block; two vertices
    #    are two separate blocks without a bond between them.
    #    ⇒ A wiring here could never fire.  It would be exactly the mistake
    #      this campaign has found five times in one day since 19.08.:
    #      a mechanism that is built in and whose reach is zero.
    #      If `vertex_specs` ever becomes multidentate, `_fold_rings_with_mc` belongs
    #      here -- not before.
    _fold_rings = (_fold_rings_from_blocks([(o, m) for (o, m, _dl) in block_specs],
                                           out_syms)
                   if _fold_fp_enabled() else None)
    _fold_kept = []                                  # fingerprint per kept frame
    for cb in eval_order:
        blocks = [np.zeros((1, 3))]
        placed_P = [np.zeros(3)]; placed_syms = [metal]
        ok = True
        for vi, ci in enumerate(cb):
            Q = per_vertex_cands[vi][ci][0]
            blocks.append(Q)
            for row in Q:
                placed_P.append(row)
            placed_syms += vertex_lsyms[vi]
        P = np.vstack(blocks)
        if not np.all(np.isfinite(P)):
            continue
        if refine:
            P = _refine_guarded(out_syms, P, fixed)
            # #308 whole-complex torsion-space clash relax (env-gated, default-OFF
            # byte-id); torsion-only, never-worse, metal+donors (`fixed`) frozen.
            P = _torsion_relax_frame(out_syms, P, fixed, block_specs)
            # JOINT inter-ligand declash (env-gated, default-OFF byte-id): global
            # inter-ligand heavy-heavy minimisation, core frozen.  After #308.
            P = _joint_declash_frame(out_syms, P, fixed, block_specs, geom=geometry)
            # SOFT coordination-sphere flex (env-gated, default-OFF byte-id): donors
            # breathe a bounded amount to open residual mild inter-ligand clashes.
            P = _sphere_flex_frame(out_syms, P, fixed, block_specs)
        if not np.all(np.isfinite(P)):
            continue
        # complex-level RMSD dedup vs already-kept frames
        dup = False
        _fp = None                                   # lazy: only on RMSD proximity
        for _ki, (_, Pk) in enumerate(frames):
            if Pk.shape == P.shape and _complex_rmsd(out_syms, P, Pk) < rmsd_dedup:
                if _fold_rings is None:
                    dup = True                       # switch OFF -> old predicate
                    break
                if _fp is None:
                    _fp = _fold_fp(P, _fold_rings)
                if _fold_kept[_ki] is None:
                    _fold_kept[_ki] = _fold_fp(Pk, _fold_rings)
                if _fold_same(_fold_kept[_ki], _fp):
                    dup = True
                    break
        if dup:
            continue
        frames.append((list(out_syms), P))
        _fold_kept.append(None)                      # same length as `frames`
        if len(frames) >= n_frames:
            break
    if not frames:
        return None
    return frames


def _constrained_uff_relax(ligmol, fixed_idx, max_its=500):
    """METAL-FREE constrained UFF relax (the _hapto_rigid_v2 pattern): the ligands-
    only mol (no metal, no M-D bonds) is UFF-minimized with the donor atoms FIXED
    at their placed vertex positions and inter-fragment vdW ON, so ligand
    peripheries relax + inter-ligand clashes/H-overcoord resolve WITHOUT UFF having
    to type the metal (which it can't reliably do for 4d/5d)."""
    try:
        ff = AllChem.UFFGetMoleculeForceField(ligmol, ignoreInterfragInteractions=False)
        if ff is None:
            return False
        ff.Initialize()                       # MUST precede AddFixedPoint, else the fixed
        for i in fixed_idx:                   # points are cleared -> donors drift (M-D break)
            ff.AddFixedPoint(int(i))
        ff.Minimize(maxIts=max_its)
        return True
    except Exception:
        return False


def build_and_relax(metal: str, geometry: str, vertex_specs, relax: bool = True):
    """Heteroleptic monodentate assembly + metal-free constrained relax (donors
    pinned at vertices).  vertex_specs[i] = (frag_mol, donor_local_idx).
    Returns (syms, P) = [metal] + relaxed ligand atoms.  Deterministic."""
    ref = MSB._ref_vectors(geometry)
    if len(vertex_specs) != len(ref):
        raise ValueError("vertex_specs count != vertices")
    lig = Chem.RWMol()              # ligands-only mol (no metal) for the FF
    conf_xyz = []
    fixed = []
    for i, (frag, di) in enumerate(vertex_specs):
        Vunit = ref[i] / np.linalg.norm(ref[i])
        emb = _ligand_3d_from_mol(frag)
        if emb is None:
            return None
        lsyms, lP, lmol = emb
        md = MSB.md_distance(metal, lsyms[di],
                             atom=lmol.GetAtomWithIdx(di), mol=lmol)
        vertex = Vunit * md
        if len(lsyms) == 1:
            Q = vertex.reshape(1, 3)
        else:
            lp = _donor_and_lp(lsyms, lP, lmol, di)
            R = _rot_align(lp, -Vunit)
            Q = (lP - lP[di]) @ R.T + vertex
        offset = lig.GetNumAtoms()
        for a in lmol.GetAtoms():
            lig.AddAtom(Chem.Atom(a.GetAtomicNum()))
        for b in lmol.GetBonds():
            lig.AddBond(b.GetBeginAtomIdx() + offset, b.GetEndAtomIdx() + offset,
                        b.GetBondType())
        fixed.append(offset + di)          # donor pinned at its vertex
        for row in Q:
            conf_xyz.append(np.asarray(row, float))
    if lig.GetNumAtoms() == 0:
        return None
    conf = Chem.Conformer(lig.GetNumAtoms())
    for k, xyz in enumerate(conf_xyz):
        conf.SetAtomPosition(k, [float(xyz[0]), float(xyz[1]), float(xyz[2])])
    lig.AddConformer(conf, assignId=True)
    if relax:
        try:
            Chem.SanitizeMol(lig, catchErrors=True)
        except Exception:
            pass
        _constrained_uff_relax(lig, fixed)
    LP = lig.GetConformer().GetPositions()
    syms = [metal] + [a.GetSymbol() for a in lig.GetAtoms()]
    P = np.vstack([np.zeros((1, 3)), np.array(LP, float)])
    return syms, P


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


# ===== GLOBAL DONOR SEATING: the distance nobody in the construction is looking at ===========
#
# WHAT IS MISSING.  Every seating quantity this builder computes is computed for ONE ligand at
# a time: _orient_chelate_to_vertices takes ONE ligand's coordinates, fits ITS donors onto ITS
# assigned vertices and resets ITS M-D radii.  The bite is an intra-ligand distance, r(M-D) is
# a radius, beta is a local plane.  Not one term in the whole seating has a distance between
# TWO DIFFERENT ligands in it.  The conformer PICK sees the neighbour (_clash_count vs
# `placed`), but it is greedy and sequential -- ligand 1 is chosen while ligands 2 and 3 do not
# exist yet -- and once a conformer is picked, nothing can move it relative to its neighbours.
#
# MEASURED, three times independently on 2026-08-02:
#   KEJCUZ dissected: C4...O28 and O15...C39 both sit at 1.85 A = 57 % of the vdW sum, and BOTH
#     pairs are INTER-LIGAND (donor of one chelate against backbone carbon of another).  The
#     per-ligand bite is blind to them by construction.
#   donorplane3 (187 systems): rotating a seated ligand about its donor BEST-FIT line breaks
#     the topology from THREE donors on -- the line passes through none of them.
#   donorplane4 (187 systems): even at TWO donors the backbone swing cost 12 systems
#     (ccdc_arrangement_lost / ccdc_backbone_lost) when the objective was beta.
#
# The lesson of the last two is NOT "never move a seated ligand".  It is "do not move it for a
# reason unrelated to what is broken": beta fired on EVERY system, so it paid the backbone cost
# everywhere and collected a benefit only sometimes.  This fires ONLY where an inter-ligand
# contact is actually below the floor -- a frame with no such contact is returned untouched and
# unchanged, which is most frames.
#
# WHAT IS SOLVED FOR, and why the two things we already get right cannot be traded away.
# The unknowns are the positions of ALL donors of ALL ligands at once.  The two constraints the
# seating already satisfies are enforced by the PARAMETERISATION, not as penalty terms that an
# optimiser could sell:
#     r(M-D)      is invariant under any rotation about the METAL (the metal is the origin in
#                 this frame), so every M-D band stays exactly where md_distance put it;
#     the BITE    -- and with it every intra-ligand distance -- is invariant under ANY rigid
#                 motion of the whole ligand body.  d_DD does not change, so by the ring-closure
#                 law cos(bite) = (r1^2 + r2^2 - d_DD^2)/(2 r1 r2) the bite ANGLE cannot change
#                 either.  Nothing is re-solved; it is carried.
# What is left free is precisely the quantity that was never in any objective: how the ligands
# sit relative to EACH OTHER.  Both invariants are re-CHECKED numerically per move (_gd_move_ok)
# rather than asserted -- one instrument that nobody holds against a second says nothing.
#
# WHY THIS IS NOT joint_declash.  delfin/manta/joint_declash.py already minimises an
# inter-ligand heavy-heavy objective, but it runs with the metal AND ALL DONORS FROZEN, and
# torsion_relax.identify_dofs drops any rotation whose moving half contains a frozen atom other
# than the pivot.  For a CHELATE the M-D1 whole-body spin moves D2, which is frozen -> the DOF
# is dropped, and the D1-D2 line is not a bond so it is not a torsion axis either.  A rigid
# chelate (acac, oxalate, bipy, phen) therefore has NO degree of freedom there at all -- exactly
# the KEJCUZ case.  This pass exists because the donors have to be allowed to move, together.
_GD_CLASH_F = 0.75      # the same clash factor the self-gate and joint_declash use
_GD_H_W = 0.05          # H contacts are a tie-breaker only; the gate rejects on heavy-heavy
_GD_RESID = 0.05        # A, furthest a donor may travel from the vertex it was seated on
# 2026-08-02: 0.25 A was MEASURED and lost.  n=53 affected, valid 25->27, cap_gained=4 but
# cap_LOST=2 (KEGMEP, KIQNUT) and good_regr=1 -> topology_floor=False.  0.25 A at r=2.1 A is
# ~7 deg, which is enough for a donor to leave the polyhedron vertex it was enumerated onto,
# and the topology floor is exactly the term that notices.  0.05 A is ~1.4 deg: the bounded
# stage survives only as a nudge, and what is left is essentially the axis on which the
# donors PROVABLY do not move.  Set DELFIN_FFFREE_GD_RESID to sweep it; 0 disables the
# bounded stage outright (free axis only).
_GD_FREE_STEPS = 24     # 15 deg grid on the zero-cost axis; the objective is smooth
_GD_BOUND_STEPS = 6     # steps each way inside the residual cap
_GD_PASSES = 4          # coordinate-descent sweeps over the ligands
_GD_EPS = 1e-6          # numerical slack on "did not move" -- a rigid rotation about an axis
                        # THROUGH a donor leaves it put only to float precision, so resid=0
                        # must still admit that, or the free axis rejects itself.


def _global_donor_seat_enabled() -> bool:
    """THE one place DELFIN_FFFREE_GLOBAL_DONORS is read (default OFF -> byte-identical)."""
    return os.environ.get("DELFIN_FFFREE_GLOBAL_DONORS", "0") == "1"


def _gd_resid() -> float:
    """THE one place DELFIN_FFFREE_GD_RESID is read.  Only ever consulted from inside
    _global_donor_seat, i.e. only when DELFIN_FFFREE_GLOBAL_DONORS is on -- with the flag off
    this is dead code and the frame is byte-identical whatever the variable says."""
    try:
        return max(0.0, float(os.environ["DELFIN_FFFREE_GD_RESID"]))
    except Exception:
        return _GD_RESID


def _gd_loss(X, mh, ml, fl):
    """(loss, worst inter-ligand heavy contact).  loss = sum of squared shortfalls below the
    vdW floor over inter-ligand pairs, heavy-heavy dominant.  Same shape as the self-gate's
    own measure, so a move that lowers this is a move toward passing the gate."""
    D = np.linalg.norm(X[:, None, :] - X[None, :, :], axis=2)
    over = np.clip(fl - D, 0.0, None)
    L = float((over[mh] ** 2).sum() + _GD_H_W * (over[ml] ** 2).sum())
    return L, (float(D[mh].min()) if mh.any() else float("inf"))


def _gd_move_ok(T, X0, blocks, resid):
    """Re-measure the two invariants instead of trusting the parameterisation.

    r(M-D): the metal is the origin, so |x_d| must be unchanged to numerical precision.
    The BITE: every donor-donor distance INSIDE a ligand must be unchanged likewise.
    Plus the one quantity this pass is allowed to spend -- how far a donor has drifted from
    the vertex the enumeration seated it on -- capped at ``resid`` against the ORIGINAL
    frame (not the previous step), so repeated sweeps cannot accumulate a walk-away.
    ``resid`` is passed in rather than read here: this runs once per candidate angle per
    donor, and an environment lookup in that loop would be the most expensive line in it."""
    cap = resid + _GD_EPS
    for _st, _ln, dn in blocks:
        for a in range(len(dn)):
            da = dn[a]
            if abs(float(np.linalg.norm(T[da])) - float(np.linalg.norm(X0[da]))) > 1e-6:
                return False                                   # M-D band broken
            if float(np.linalg.norm(T[da] - X0[da])) > cap:
                return False                                   # drifted off its vertex
            for b in range(a + 1, len(dn)):
                db = dn[b]
                if abs(float(np.linalg.norm(T[da] - T[db]))
                       - float(np.linalg.norm(X0[da] - X0[db]))) > 1e-6:
                    return False                               # bite broken
    return True


def _global_donor_seat(syms, P, blocks):
    """Move every ligand body -- donors included -- against every OTHER ligand, jointly.

    ``blocks``: one ``(start, n_atoms, [donor global indices])`` per placed ligand; the metal
    is index 0 and never moves.  Returns the improved frame, or None when there is nothing
    below the floor / nothing improved (the caller then keeps the seated frame verbatim).

    Two kinds of rigid motion, tried in this order because the first one is FREE:

      1) rotation about the ligand's OWN donor axis -- the M-D bond for a monodentate, the
         line THROUGH both donors for a bidentate.  The donors lie ON the axis, so they do
         not move at all: r(M-D), the bite AND the vertex alignment are all untouched, and
         only the backbone swings.  Cost to everything already right: exactly zero.
         NOT offered from three donors on: there the "axis" would be a best-fit line through
         none of them, which is the donorplane3 failure verbatim.

      2) rotation of the whole ligand about the METAL, capped so no donor leaves its vertex
         by more than _GD_RESID.  This is the step where the donors genuinely re-seat, and it
         is the only one available to a rigid tri-/tetradentate -- and it is also the step
         that lost the 0.25 A measurement, because a donor that leaves its vertex is exactly
         what the topology floor is watching for.  At resid 0 it is not offered at all and
         only (1) remains.

    Coordinate descent, fixed ligand order, fixed angular grid, no RNG, accept-only-if-better,
    with a never-worse floor on the WORST inter-ligand heavy contact so the sum objective
    cannot buy three mild reliefs by crushing one comfortable pair.
    """
    try:
        X0 = np.array(P, float)
    except Exception:
        return None
    n = len(syms)
    if len(blocks) < 2 or X0.shape != (n, 3) or not np.all(np.isfinite(X0)):
        return None
    try:
        from delfin.manta.refine import _vdw
    except Exception:
        return None
    lig = np.full(n, -1, dtype=int)
    for bi, (st, ln, _dn) in enumerate(blocks):
        if st < 1 or st + ln > n:
            return None                                  # bookkeeping mismatch -> do nothing
        lig[st:st + ln] = bi
    vdw = np.array([_vdw(s) for s in syms], float)
    fl = _GD_CLASH_F * (vdw[:, None] + vdw[None, :])
    isH = np.array([s == "H" for s in syms])
    inter = ((lig[:, None] != lig[None, :]) & (lig[:, None] >= 0) & (lig[None, :] >= 0)
             & np.triu(np.ones((n, n), bool), 1))
    mh = inter & ~isH[:, None] & ~isH[None, :]
    ml = inter & (isH[:, None] | isH[None, :])
    if not mh.any() and not ml.any():
        return None
    L0, h0 = _gd_loss(X0, mh, ml, fl)
    if L0 <= 1e-12:
        return None            # nothing inter-ligand under the floor -> the frame is returned
    hfloor = h0 - 1e-6         # ... unchanged, which is what makes this cheap on clean frames
    resid = _gd_resid()        # read ONCE, outside every loop
    Xc = X0.copy()
    bestL = L0
    for _p in range(_GD_PASSES):
        improved = False
        for st, ln, dn in blocks:
            if ln < 1 or not dn:
                continue
            sl = slice(st, st + ln)
            cands = []                         # (origin, axis, angles) in try-order
            # 1) THE FREE AXIS (see above): donors stay exactly put.  Needs a BODY to swing
            #    -- for a monatomic ligand (Cl-, the atom IS the donor) it moves nothing, so
            #    such a ligand only gets the bounded stage below.
            if len(dn) == 1 and ln >= 2:
                cands.append((np.zeros(3), np.asarray(Xc[dn[0]], float),
                              [2.0 * math.pi * k / _GD_FREE_STEPS
                               for k in range(1, _GD_FREE_STEPS)]))
            elif len(dn) == 2:
                cands.append((np.asarray(Xc[dn[0]], float),
                              np.asarray(Xc[dn[1]], float) - np.asarray(Xc[dn[0]], float),
                              [2.0 * math.pi * k / _GD_FREE_STEPS
                               for k in range(1, _GD_FREE_STEPS)]))
            # 2) THE BOUNDED AXES through the metal: the donors move, together, by at most
            #    the chord _GD_RESID allows.  The chord of a given rotation grows with the
            #    radius, so it is the LONGEST M-D in this ligand that decides the angle for
            #    the whole body -- taking the shortest would let the outer donors overrun the
            #    cap (_gd_move_ok would then throw those candidates away, silently).
            rmax = max(float(np.linalg.norm(Xc[d])) for d in dn)
            if rmax > 1e-6 and resid > _GD_EPS:
                tmax = 2.0 * math.asin(min(1.0, resid / (2.0 * rmax)))
                angs = [s * m for s in
                        [tmax * k / _GD_BOUND_STEPS for k in range(1, _GD_BOUND_STEPS + 1)]
                        for m in (1.0, -1.0)]
                for e in (np.array([1.0, 0.0, 0.0]), np.array([0.0, 1.0, 0.0]),
                          np.array([0.0, 0.0, 1.0])):
                    cands.append((np.zeros(3), e, angs))
            for org, ax, angs in cands:
                na = float(np.linalg.norm(ax))
                if na < 1e-6:
                    continue
                a_ = np.asarray(ax, float) / na
                locL, locX = bestL, None
                for th in angs:
                    R = _axis_rot(a_, th)
                    T = Xc.copy()
                    T[sl] = (Xc[sl] - org) @ R.T + org
                    if not np.all(np.isfinite(T)):
                        continue
                    Lt, ht = _gd_loss(T, mh, ml, fl)
                    if Lt < locL - 1e-9 and ht >= hfloor and _gd_move_ok(T, X0, blocks, resid):
                        locL, locX = Lt, T
                if locX is not None:
                    Xc, bestL, improved = locX, locL, True
        if not improved:
            break
    if bestL >= L0 - 1e-9 or not np.all(np.isfinite(Xc)):
        return None
    return Xc


# ===== THE HALF BAILAR TWIST, TURNED BACK FF-FREE ===========================
# Measured 18.08.2026 on 30921 systems: net +988 systems flow from the octahedron
# into the trigonal prism (McNemar X2 = 860.8 on 1 df).  TPR-6 is built 2.99 times
# as often as it really occurs, while EVERY other shape lies between 0.86 and 1.29
# -- the largest single defect of the polyhedron axis, 3.7 times the mass of the
# second-largest pair.
#
# THREE MEASUREMENTS SAY WHAT IT IS NOT:
#   * It is NOT a ligand-field effect.  Metal flat, d-count flat; the signal is
#     solely the INTERLOCKING by chelate rings, monotone from 2.49 % at zero rings
#     to 16.12 % at five (6.5-fold).  So geometry, not chemistry.
#   * They are NOT real prisms.  CShM(OC-6) lies at a median of 11.03 instead of
#     16.7, as an ideal TPR would have -- a HALF Bailar twist, stopped
#     halfway.
#   * It is NOT a selection question.  poly_match is false in 1061 of 1061 cases,
#     although the eye reads poly_build as the minimum over ALL realistic frames.
#     In the whole manifold there is no octahedron -- so none is built.
#
# ⛔ WHY THE EXISTING REPAIR IS NOT ENOUGH.  It exists twice, both in
# smiles_converter.py: DELFIN_FFFREE_CN6_OH_ADD (:27452) begins with `apply_uff and`,
# and DELFIN_FFFREE_CN6_OH_ANGLES (:38153) sits in
# _build_coordination_constraints_from_xyz, i.e. in the UFF constraint machinery.  Both
# give UFF octahedral angle targets (90/180 degrees) so that UFF relaxes the twist out.
# On the FF-free path no UFF runs ⇒ reach zero.  They cannot be
# wired in; they must arise anew on the construction side, and that is this block.
#
# ⚠ WHAT IS DELIBERATELY CARRIED OVER HERE, AND WHY.  The three safeguards of the UFF
# version are not trimmings, they are the reason it is isomer-safe:
#   1) The three trans pairs come from the FRAME ITSELF (greedy: each donor with its
#      most opposite one), NOT from an enumerator permutation.  The
#      PERM variant was measured on 2026-07-14 and COLLAPSED isomers
#      (VOYWUD lost all-trans + trans-OH) -- it forces the wrong trans set on some
#      arrangements.  What the frame already has stays: fac stays fac,
#      cis stays cis.  It is a twist correction, not an arrangement change.
#   2) Only if ALL three pairs are clearly trans (min > 120 degrees).  A valid
#      TPR/OC frame sits at 140-180; an ambiguous one does not -> skipped, so that
#      nothing collapses.
#   3) Only at CN 6 and _PREFERRED_CN6_GEOMETRY.get(metal, 'OH') == 'OH'.
#
# ⚠ AND WHAT IS DIFFERENT -- that is the reason why this can land.  By the cost law
# measured on 18.08. at three points, ordering/selection costs about 0,
# isometry +0.98 pp, a rigid rotation with NEW conformation +6.57 pp and a
# re-embedding +11.9 pp.  What ADDS, without inventing new geometry, lands.
# That is why here the WHOLE ligand arm is rotated RIGIDLY about the metal and not the
# single donor atom shifted: a rotation about the metal leaves r(M-D) exact and
# the bite exact, it invents no conformation.  Shifting single donors
# would tear bonds -- exactly the re-embedding that was measured to be the most expensive.
#
# A chelate with a 78-degree bite CANNOT give a perfect octahedron, and it is not
# meant to: the Kabsch fit puts its donors as close to the ideal directions as
# its bite allows, and the bite wins.  The goal is not CShM 0, the goal
# is "no more half twist".
_OC6_TRANS_MIN = 120.0   # degrees; below it the pair is not unambiguously trans -> abort
_OC6_AXIS_MIN = 0.20     # smallest singular value of the three axes: below it they are
                         # almost coplanar and the orthonormalisation would be guessed
_OC6_RIGID_TOL = 1e-6    # Angstrom; re-measurement of r(M-D) and bite AFTER the rotation

# ⚠ ONE EXIT, ONE LINE -- and `call` right at the front.  The function has ten ways
# to come back with None, and every single one means something different: "not my case",
# "ambiguous, hands off", "rotated, but it achieved nothing".  Without the denominator
# they would all be the same zero in the report -- the mistake that was made here on
# 10.08. and once more on 14.08.  Only recorded if the corrector is called at all,
# i.e. only behind oc6_twist=True.  Read by _self_test_oc6_twist and by the
# sibling self-test in converter_backend.
_OC6_SEAT_CENSUS = dict.fromkeys(
    ("call", "not_oc6", "metal_pref", "shape", "book", "not_cn6", "zero_md",
     "ambiguous", "coplanar_axes", "rot_bad", "rigid_broken", "cshm_flat", "ok"), 0)
# (CShM before, CShM after, smallest trans angle) per call -- the raw numbers from
# which to read WHETHER the seating is twisted at all.  A counter alone could not
# say that: "not improved" means either "already correct" or "too
# bad to rescue", and those are opposite findings.  Only behind oc6_twist.
_OC6_SEAT_CSHM = []
_OC6_CSHM_KEEP = 4096   # cap; the counters above stay complete, only the
                        # raw-value list stops growing at some point


def _oc6_twist_seat_enabled() -> bool:
    """THE one read site of DELFIN_FFFREE_OC6_TWIST_SEAT (default 0 -> byte-identical).

    It does NOT decide about the primary frame.  ``assemble_from_config`` executes the
    correction exclusively on the keyword ``oc6_twist=True``, which
    is False by default -- so the built primary frame is byte-identical, no
    matter what is in the environment.  This switch only says whether the caller
    additionally builds a SIBLING FRAME."""
    return os.environ.get("DELFIN_FFFREE_OC6_TWIST_SEAT", "0") == "1"


def _oc6_trans_pairs(u, donors):
    """The three trans pairs from the frame itself: each donor with its most
    opposite one, greedy in fixed index order (deterministic, no RNG).

    ``u``: dict donor index -> unit vector from the metal.  Returns
    ``(pairs, smallest_trans_angle_deg)`` or ``(None, 0.0)``.  Word for word the
    pairing of the UFF version in smiles_converter.py:38160 -- not out of convenience,
    but because EXACTLY this pairing is the isomer-safe one (see head note)."""
    rem = list(donors)
    pairs = []
    min_trans = 180.0
    while len(rem) >= 2:
        a = rem[0]
        b = min(rem[1:], key=lambda x: float(np.dot(u[a], u[x])))
        c = max(-1.0, min(1.0, float(np.dot(u[a], u[b]))))
        min_trans = min(min_trans, math.degrees(math.acos(c)))
        pairs.append((a, b))
        rem.remove(a)
        rem.remove(b)
    if len(pairs) != 3:
        return None, 0.0
    return pairs, min_trans


def _oc6_ideal_axes(u, pairs):
    """The three measured trans axes, pulled onto the NEAREST orthonormal triad
    (polar decomposition, ``A = U S Vt`` -> ``U Vt``).

    WHY POLAR DECOMPOSITION AND NOT GRAM-SCHMIDT: Gram-Schmidt is order-dependent
    -- the first axis would stay untouched, the third would carry the whole error.  The
    polar decomposition minimises the sum of squares over all three simultaneously and
    is thereby independent of which pair was found first.  That matters,
    because the greedy pairing above has an index order, but the chemistry does not.

    The handedness is NOT corrected.  What is sought are three mutually perpendicular
    unit vectors; {±e1, ±e2, ±e3} is the same octahedron whether the triad is right-
    or left-handed.  A det correction would be no protection here, but an
    additional, unnecessary rotation.  Returns ``(E, smallest_singular_value)``."""
    A = []
    for a, b in pairs:
        ax = u[a] - u[b]
        na = float(np.linalg.norm(ax))
        if na < 1e-9:
            return None, 0.0
        A.append(ax / na)
    A = np.asarray(A, float)
    try:
        U, S, Vt = np.linalg.svd(A)
    except Exception:
        return None, 0.0
    if not np.all(np.isfinite(U)) or not np.all(np.isfinite(Vt)):
        return None, 0.0
    return U @ Vt, float(S[-1])


def _oc6_twist_seat(syms, P, blocks, metal, geometry):
    """Rotate the half twist out -- rigidly, ligand by ligand, about the metal.

    ``blocks``: per seated ligand one ``(start, n_atoms, [global donor indices])``,
    the same bookkeeping ``_global_donor_seat`` uses; atom 0 is the metal.
    Returns: the corrected frame, or ``None`` if one of the safeguards
    triggers OR the correction does not measurably reduce the twist.  The caller
    then keeps the seated frame verbatim.

    Procedure:
      1) CN 6, OC-6 requested, metal prefers OH -- otherwise nothing.
      2) trans pairs from the frame, all three clearly trans (> 120 degrees).
      3) orthonormalise the three axes -> the octahedron NEAREST to the frame.
         Not the lab-axes octahedron: the nearest is the one requiring the least
         movement, and movement is exactly what is paid for under the cost
         law.
      4) per ligand ONE rigid rotation about the metal that puts its donors onto
         their target directions in the Kabsch sense.  A monodentate ligand gets
         the minimal rotation (Rodrigues), from two donors on the Kabsch fit -- which
         is NOT degenerate for two points, because the metal at the origin is held
         along and the covariance thus has rank 2, whose null direction is unique.
      5) RE-MEASURE instead of trust: r(M-D) and every intra-ligand donor-donor
         distance must be unchanged to 1e-6, and CShM(OC-6) must have STRICTLY
         fallen.  Both are measurements on the result, not thresholds tuned to a
         pool -- a rotation that does not reduce the twist
         is rejected, instead of being dressed up."""
    def _no(reason):
        _OC6_SEAT_CENSUS[reason] += 1
        return None

    _OC6_SEAT_CENSUS["call"] += 1
    if not str(geometry).startswith("OC-6"):
        return _no("not_oc6")              # a REQUESTED TPR-6 stays a TPR-6
    # Metal preference.  Deferred import as in _finish_config_frame; if the
    # module fails, the default 'OH' applies -- the same the table itself gives.
    try:
        from delfin.smiles_converter import _PREFERRED_CN6_GEOMETRY as _PCN6
        if _PCN6.get(str(metal), 'OH') != 'OH':
            return _no("metal_pref")
    except Exception:
        pass
    try:
        X0 = np.asarray(P, float)
    except Exception:
        return _no("shape")
    n = len(syms)
    if X0.shape != (n, 3) or not np.all(np.isfinite(X0)) or not blocks:
        return _no("shape")
    donors = []
    for st, ln, dn in blocks:
        if st < 1 or st + ln > n:
            return _no("book")             # bookkeeping does not fit -> do nothing
        donors += [int(x) for x in dn]
    if len(donors) != 6 or len(set(donors)) != 6:
        return _no("not_cn6")              # CN 6, and every donor exactly once
    M = X0[0].copy()
    u = {}
    r = {}
    for d in sorted(donors):
        v = X0[d] - M
        nv = float(np.linalg.norm(v))
        if nv < 1e-6:
            return _no("zero_md")
        u[d] = v / nv
        r[d] = nv
    pairs, min_trans = _oc6_trans_pairs(u, sorted(donors))
    if pairs is None or min_trans <= _OC6_TRANS_MIN:
        return _no("ambiguous")            # ambiguous -> skip, nothing collapses
    E, smin = _oc6_ideal_axes(u, pairs)
    if E is None or smin < _OC6_AXIS_MIN:
        return _no("coplanar_axes")        # almost coplanar axes -> the triad would be guessed
    tgt = {}
    for i, (a, b) in enumerate(pairs):
        e = np.asarray(E[i], float)
        ne = float(np.linalg.norm(e))
        if ne < 1e-9:
            return _no("coplanar_axes")
        e = e / ne
        if float(np.dot(u[a], e)) < 0.0:
            e = -e                         # the axis points towards a, not away from a
        tgt[a] = e * r[a]
        tgt[b] = -e * r[b]
    Xc = X0.copy()
    for st, ln, dn in blocks:
        dn = [int(x) for x in dn]
        if ln < 1 or not dn:
            continue
        obs = np.asarray([X0[d] - M for d in dn], float)
        tar = np.asarray([tgt[d] for d in dn], float)
        if len(dn) == 1:
            R = _rot_align(obs[0], tar[0])
        else:
            R = _kabsch_rot(obs, tar)
        if R is None or not np.all(np.isfinite(R)):
            return _no("rot_bad")
        Xc[st:st + ln] = (X0[st:st + ln] - M) @ R.T + M
    if not np.all(np.isfinite(Xc)):
        return _no("rot_bad")
    # 5a) RE-MEASURE the two invariants.  A rotation about the metal holds them
    #     mathematically; they are measured anyway, because a degenerate Kabsch matrix
    #     could silently smuggle in a mirroring exactly here.
    for st, ln, dn in blocks:
        dn = [int(x) for x in dn]
        for i, da in enumerate(dn):
            if abs(float(np.linalg.norm(Xc[da] - M)) - r[da]) > _OC6_RIGID_TOL:
                return _no("rigid_broken")             # r(M-D) broken
            for db in dn[i + 1:]:
                if abs(float(np.linalg.norm(Xc[da] - Xc[db]))
                       - float(np.linalg.norm(X0[da] - X0[db]))) > _OC6_RIGID_TOL:
                    return _no("rigid_broken")         # bite broken
    # 5b) and the one number that matters.  If it does not fall, the block has nothing
    #     to offer and returns the frame unchanged (the caller keeps it).
    try:
        from delfin.manta import polyhedra as _PH
        before = _PH.cshm([X0[d] - M for d in sorted(donors)], "OC-6 octahedron")
        after = _PH.cshm([Xc[d] - M for d in sorted(donors)], "OC-6 octahedron")
        if len(_OC6_SEAT_CSHM) < _OC6_CSHM_KEEP:      # capped: a 30-hour run
            _OC6_SEAT_CSHM.append(                    # must not drag a list along
                (float(before), float(after), float(min_trans)))
    except Exception:
        return _no("cshm_flat")
    if not (after < before - 1e-9):
        return _no("cshm_flat")
    _OC6_SEAT_CENSUS["ok"] += 1
    return Xc


# ===== THE LIGAND'S LAST DEGREE OF FREEDOM ===================================
# THE FINDING THAT FORCES THIS BLOCK (18.08., 492 clean remaining systems).
# `org_bond` -- the largest defect mass of the organic geometry -- does NOT hang on
# the bond class (within a molecule the localisation is massive, but identical between
# hits and failures, max |delta| 0.14).  Two quantities genuinely enrich:
# `worst_n` (number of simultaneously bent organic bonds) 3.83x, and the
# COORDINATION NUMBER 1.62x (CN>=5 51.3 % against CN<=4 31.7 %).  And CN acts THROUGH
# `worst_n`: mean `worst_n` rises 3.46 -> 7.07 (CN 2..6) at practically
# constant molecule size, with a CORRECTLY built polyhedron 3.6 against 7.1.
# ⇒ At equal ligand size a CN-6 centre bends twice as many
# organic bonds as a CN-4 centre.  The polyhedron is enforced, the ligand pays.
#
# THE BUILDER'S DOF BALANCE (re-measured, see _self_test_ligand_dof):
#   * A ligand block is seated RIGIDLY -- 6 rigid-body DOF.
#   * Monodentate: `_rot_align(lp, -Vunit)` fixes 5 of them (3 translation via
#     `+ Vunit*md`, 2 direction via the Rodrigues rotation).  The SIXTH -- the
#     azimuth about the M-D axis -- is chemically completely undetermined and is
#     nonetheless nailed down, namely to the random value that the MINIMAL rotation of
#     `_rot_align` happens to deliver.  It is NEVER sampled on the configuration path.
#   * Bidentate: the Kabsch fit in `_orient_chelate_to_vertices` puts both donors
#     onto their vertices; exactly ONE rotation remains -- about the donor-donor axis.
#     The fit cannot see it (both donors lie ON the axis).  Today it is
#     used only by `_lp_orient_seated_bidentate` to optimise a BOND angle
#     (default OFF, measured negative as a seating) -- for the PACKING nobody
#     has ever used it.
#   * Tridentate and higher: three non-collinear donors fix the rigid body
#     COMPLETELY.  ZERO free rotation.  The register is right: 0 DOF.
#
# WHAT THIS AXIS PRESERVES -- and that is the reason it is built.  Both
# cases rotate about a line that CONTAINS EVERY donor of this block.  The donors
# are thereby POINTWISE fixed: M-D distance, M-D DIRECTION, bite, vertex angles and CShM
# do not change by one bit.  It is an isometry of the ligand block on an
# already built conformation -- no new embedding, no new conformation.
# By the cost law measured at four points (ordering/selection ~0, isometry
# +0.98 pp, rigid rotation with NEW conformation +6.57 pp, re-embedding +11.9 pp)
# that is the cheapest class that can change anything about the packing at all.
#
# ⚠ DISTINCTION FROM DELFIN_FFFREE_LIGAND_SWING (joint_declash.py:305).  The swing is
# a DIFFERENT movement at a DIFFERENT place: it rotates the whole ligand about the
# axis M -> donor centroid, AFTER assembly, in the declash pass, and in doing so THE
# DONORS MOVE -- that is why it needs a 3-degree cap and a CShM budget, and
# that is why it can damage the polyhedron at all (at 8 degrees it did).  This
# axis here does not move the donors, therefore needs no angle cap and may sample the
# full circle; and it acts at the SELECTION, where the candidate is still being
# decided, not afterwards on the finished frame.  Extending the swing would mean
# rebuilding it into a movement it is not -- and in a file that does not belong to this
# territory.  Hence a second axis, not an extension.
_LIGAND_DOF_AXIS_TOL = 1.0e-6      # collinearity of the donors / axis proximity (Angstrom)


def _ligand_dof_seat_enabled() -> bool:
    """THE one place where DELFIN_FFFREE_LIGAND_DOF_SEAT is read
    (default 0 -> byte-identical)."""
    return os.environ.get("DELFIN_FFFREE_LIGAND_DOF_SEAT", "0") == "1"


def _ligand_dof_seat_steps() -> int:
    """Number of sampled angles on the full circle (default 6, like the already
    existing CN2 axis rotation DELFIN_FFFREE_CN2_SPINS).  Only read when the
    switch above is on."""
    try:
        n = int(os.environ.get("DELFIN_FFFREE_LIGAND_DOF_SEAT_N", "6"))
    except Exception:
        n = 6
    return max(2, min(n, 36))


def _free_rigid_axis(Q, donor_locals, metal_pos):
    """The ONLY line about which this seated ligand block may still be rotated
    rigidly without changing the coordination sphere by even one bit -- or ``None``.

    Condition: the axis must contain EVERY donor of the block, then all donors are
    pointwise fixed and thus M-D distance, M-D direction, bite and CShM exactly preserved.
      * 1 donor  -> the line M--D (among all lines through the donor the only one that
                    additionally preserves EVERY M-X distance of the ligand, because the
                    metal then lies on the axis itself -- the ligand cannot swing into
                    the metal).
      * >=2 donors -> the line through the donors, but ONLY if they are collinear.
      * otherwise -> ``None``.  Three non-collinear donors fix the rigid body
                    completely; rotating anything here would mean touching the polyhedron.
    Returns ``(origin, unit_axis)``; ``origin`` is a donor point, so that this donor
    stays BIT-EXACT under the rotation.
    """
    try:
        d = [np.asarray(Q[int(i)], float) for i in donor_locals]
    except Exception:
        return None
    if not d or any(not np.all(np.isfinite(x)) for x in d):
        return None
    if len(d) == 1:
        a = d[0] - np.asarray(metal_pos, float)
        n = float(np.linalg.norm(a))
        if n < _LIGAND_DOF_AXIS_TOL:
            return None
        return (d[0], a / n)
    a = d[-1] - d[0]
    n = float(np.linalg.norm(a))
    if n < _LIGAND_DOF_AXIS_TOL:
        return None
    a = a / n
    for p in d[1:-1]:                       # Kollinearitaet ALLER Donoren
        w = p - d[0]
        if float(np.linalg.norm(w - float(np.dot(w, a)) * a)) > _LIGAND_DOF_AXIS_TOL:
            return None
    return (d[0], a)


def _free_dof_reseat(Q, lsyms, donor_locals, existing, existing_syms, base_clash,
                     metal_pos=None):
    """Samples the one free degree of freedom of this ligand block and returns
    ``(Q_rotated, clash)`` -- or ``None`` if nothing is STRICTLY better.

    NEVER-WORSE BY CONSTRUCTION, at three places:
      1) Switch off -> immediately ``None``, the caller sees nothing.  Byte-identical.
      2) The same quantity is measured that the selection reads anyway as its first
         key (`_clash_count` against the already seated atoms) -- no
         new criterion that could win against the old one.
      3) Only a STRICTLY smaller clash is taken; on a tie the
         historic pose stays.  And only a block that ALREADY collides
         (`base_clash > 0`) is touched at all -- a collision-free seating stays
         untouched, there is nothing to gain in the packing there.
    Deterministic: fixed angle list, smallest index wins on a tie.
    """
    if not _ligand_dof_seat_enabled():
        return None
    if base_clash <= 0 or len(Q) < 2:
        return None
    ax = _free_rigid_axis(Q, donor_locals,
                          np.zeros(3) if metal_pos is None else metal_pos)
    if ax is None:
        return None
    o, u = ax
    W = np.asarray(Q, float) - o
    perp = W - np.outer(W @ u, u)
    if float(np.max(np.linalg.norm(perp, axis=1))) < _LIGAND_DOF_AXIS_TOL:
        return None                     # everything lies ON the axis -> rotation = identity
    n = _ligand_dof_seat_steps()
    best = None
    for k in range(1, n):
        try:
            R = _axis_rot(u, 2.0 * np.pi * float(k) / float(n))
        except Exception:
            continue
        Qk = W @ R.T + o
        if not np.all(np.isfinite(Qk)):
            continue
        ck = _clash_count(Qk, existing, lsyms, existing_syms)
        if ck < base_clash and (best is None or ck < best[1]):
            best = (Qk, ck)
    if best is None:
        return None
    # THE ASSERTION IS RE-MEASURED, NOT CLAIMED.  Mathematically a rotation about an
    # axis through all donors holds every donor point fixed -- it is measured
    # anyway, because exactly here a degenerate axis could silently do something
    # else, and the price would be the coordination sphere.  If the assertion breaks,
    # the rotation is REJECTED (the caller keeps its historic pose).
    _mp = np.zeros(3) if metal_pos is None else np.asarray(metal_pos, float)
    for _d in donor_locals:
        _d = int(_d)
        if float(np.linalg.norm(best[0][_d] - Q[_d])) > 1.0e-9:
            return None                                   # donor has moved
        if abs(float(np.linalg.norm(best[0][_d] - _mp))
               - float(np.linalg.norm(np.asarray(Q[_d], float) - _mp))) > 1.0e-9:
            return None                                   # r(M-D) broken
    _tp = os.environ.get("DELFIN_LIGAND_DOF_TRACE", "")
    if _tp and _tp != "0":
        # ⚠️ TRACE, BECAUSE "byte-identical" HAS SEVERAL CAUSES HERE: the block does
        # not run, there is no free axis, nothing collided, or the rotation found
        # no better pose.  Without this line one could not say which applies.
        try:
            with open(_tp, "a") as _fh:
                _fh.write("[LIGDOF] ndon=%d nat=%d clash %d -> %d\n"
                          % (len(list(donor_locals)), len(Q), int(base_clash),
                             int(best[1])))
        except Exception:
            pass
    return best


def assemble_from_config(metal, geometry, config, ligands, refine=True,
                         n_frames=1, per_lig_confs=6, rmsd_dedup=0.5,
                         planar_bite=None, planar_coplanar=None, prefer_beta=False,
                         lp_orient=False, oc6_twist=False):
    """Build a 3D complex from a chelate-isomer config (vertex -> (ligand_idx,
    arm_idx)) and the decomposed ligand list.  Chelating ligands are Kabsch-fit
    onto their two assigned vertices; monodentate ligands are oriented onto their
    vertex.  Multi-conformer selection per ligand + constrained refine.

    ``n_frames``: with the default ``1`` the single clash-minimal frame is returned
    (byte-identical to the historic behaviour) as ``(syms, P, donors)``.  With
    ``n_frames > 1`` an RMSD-deduped ENSEMBLE of full-complex frames is returned as
    a list ``[(syms, P, donors), ...]`` (canonical clash-minimal frame first) that
    varies ligand internal / chelate-ring conformation while keeping the
    constructed coordination (vertex directions + M-D distances) fixed -- the same
    best-of-ensemble lever proven for CN2 / rigid-hapto, generalised to the chelate
    σ sub-path (Task A.1).  Returns ``None`` on failure (in either mode).
    Deterministic (fixed SEED, single-thread); every emitted frame is internally
    refined; ETKDG / metallacycle conformer pool already samples chelate-ring
    pucker (Cremer-Pople) across conformers, so distinct puckers fall out naturally.
    """
    ref = MSB._ref_vectors(geometry)
    # group config by ligand instance: lig_idx -> [(vertex, arm), ...]
    by_lig = {}
    for v, (li, arm) in config.items():
        by_lig.setdefault(li, []).append((v, arm))
    out_syms = [metal]; placed = [np.zeros(3)]; placed_syms = [metal]
    fixed = {0}; pos = 1
    relax_frags = []          # (AddHs(lg.mol), ligands-only offset) for the internal relax
    lig_blocks = []           # (start, n_atoms, [donor global idxs]) per placed ligand --
                              # which atoms form ONE rigid body, for the global donor seating
    # ENSEMBLE refactor (Task A.1): collect per-ligand candidate placements
    # (Q, clash_vs_metal) deduped by intra-ligand RMSD instead of only the single
    # clash-minimal pick, then enumerate full-complex combinations.  n_frames==1
    # keeps only the all-best combo -> byte-identical to the historic single path.
    ensemble = int(n_frames) > 1
    _k_confs = max(int(per_lig_confs), 1) if ensemble else 6
    per_lig_cands = []        # per ligand: [(Q, clash_vs_metal), ...] sorted best-first
    per_lig_syms = []         # per ligand: lsyms (for combo assembly)
    metal_P = np.zeros((1, 3)); metal_sym = [metal]
    # ===== THE MOST RIGID LIGAND FIRST =======================================
    # Until 18.08.2026 this loop ran in the insertion order of `config`,
    # i.e. in VERTEX order -- a monodentate could be seated before a tetradentate.
    # That is the reverse of what the task demands:
    #
    #   * An already seated monodentate counts for EVERY later candidate in
    #     `_clash_count`.  The tetradentate then arrives with almost no freedom --
    #     its bite fixes four vertices -- and the candidate that wins the
    #     selection is the one that TWISTS AWAY.
    #   * Conversely the rigid ligand sits first on its ideal vertices, and the
    #     monodentates have full rotational freedom to dodge it.  A monodentate
    #     can almost always dodge, a chelate ring never.
    #
    # That is the classic rule "most constrained variable first",
    # and it fits the measured signal exactly: the error rate of the
    # octahedron-to-prism confusion is MONOTONE in the number of chelate rings --
    # 2.49 % at zero, 8.38 % at three, 16.12 % at five (30921 systems, net +988
    # systems, McNemar X2 = 860.8).  The more rigid ligands compete for the same
    # vertices, the more often the polyhedron twists.
    #
    # ⚠ WHY THIS IS THE CHEAPEST CLASS OF ALL: it only changes the
    # ORDER in which already existing candidates are evaluated.  No
    # new geometry, not even a new selection quantity.  By the law measured today at
    # three points (isometry +0.98 pp, rigid rotation +6.57 pp,
    # re-embedding +11.9 pp) that lies even below isometry.
    #
    # ⚠ DETERMINISM: the secondary key is the ligand index, not chance.
    # Same denticity -> same order as before.
    #
    # DELFIN_FFFREE_SEAT_RIGID_FIRST (default 0 -> byte-identical).
    _lig_order = list(by_lig.items())
    if os.environ.get("DELFIN_FFFREE_SEAT_RIGID_FIRST", "0") == "1":
        def _dent_of(_li):
            try:
                return int(ligands[_li].get("denticity") or len(by_lig[_li]))
            except Exception:
                return len(by_lig[_li])
        _before = [kv[0] for kv in _lig_order]
        _lig_order.sort(key=lambda kv: (-_dent_of(kv[0]), kv[0]))
        # ⚠️ TRACE, BECAUSE "byte-identical" HAS THREE CAUSES HERE: the block does
        # not run, the order was already right, or it changes and the
        # result converges anyway.  The first smoke test delivered 19 of 19
        # identical -- without this line one could not say which of them applies,
        # and "achieves nothing" would be a claim instead of a measurement.
        # The case "was already right" is not a failure, but the
        # answer: then the builder already does it implicitly.
        _tp = os.environ.get("DELFIN_SEAT_ORDER_TRACE", "")
        if _tp and _tp != "0":
            _after = [kv[0] for kv in _lig_order]
            try:
                with open(_tp, "a") as _fh:
                    _fh.write("[SEATORD] nlig=%d dents=%s changed=%d\n"
                              % (len(_before), [_dent_of(i) for i in _before],
                                 1 if _after != _before else 0))
            except Exception:
                pass
    for li, va in _lig_order:
        lg = ligands[li]
        dons = lg["donor_local_idxs"]
        lig_offset = pos - 1                       # start index in the ligands-only frame
        # Chelates: embed the metallacycle (ligand + dative-bonded metal) so the
        # ring forms with correct geometry (backbone clears the metal) instead of
        # rigid-fitting a non-chelating free-ligand conformer.  Fall back to the
        # free-ligand path if the metallacycle embed fails.
        ring_confs = None
        lmol = None
        if lg["denticity"] >= 2:
            # ideal donor positions at the assigned polyhedron vertices (metal at origin)
            # -> the constrained metallacycle embed pins the chelate bite to match them.
            _dent = lg["denticity"]
            # CONFIG-FAITHFUL: arm index i -> the i-th CANONICAL (element-sorted)
            # donor, so the enumerator's arm->vertex assignment is realised the
            # same way for every instance of a ligand type (see _canonical_arm_order).
            _dons_d = _canonical_arm_order(lg, _dent)
            _vts = [v for v, arm in sorted(va, key=lambda x: x[1])]
            # element of the i-th canonical arm (for the ideal M-D target distance)
            _delems = [lg["mol"].GetAtomWithIdx(int(di)).GetSymbol() for di in _dons_d]
            try:
                _dtp = [ref[_vts[i]] / np.linalg.norm(ref[_vts[i]])
                        * MSB.md_distance(metal, _delems[i],
                                          atom=lg["mol"].GetAtomWithIdx(int(_dons_d[i])),
                                          mol=lg["mol"]) for i in range(_dent)]
            except Exception:
                _dtp = None
            # CHELATE-BACKBONE hardening is SCOPED to the newly-admitted LARGE backbones
            # (per-arm heavy > the historic cap of 8) so every existing native chelate
            # (per-arm <= 8: en, acac, bipy, terpy, ...) keeps the EXACT historic embed
            # -> byte-identical to OFF.  Only an over-cap arm (BIQCOV/ABEZAJ-class, only
            # reachable at all because the flag lifted the decompose cap) gets the soft
            # donor-donor windows + more conformers.  Flag-gated; OFF => harden=False
            # everywhere => byte-identical.
            _nheavy_lg = sum(1 for a in lg["mol"].GetAtoms() if a.GetAtomicNum() > 1)
            _harden = (os.environ.get("DELFIN_FFFREE_CHELATE_BACKBONE", "0") == "1"
                       and _nheavy_lg / max(_dent, 1) > 8.0)
            # RIGID PLANAR tridentate on CN5 (TBP-5 / SPY-5): force the DG bite
            # constraint so the embed carries the meridional bite (the rigid-Kabsch
            # orient cannot open a folded ~110deg bite onto the TBP/SPY meridian).
            # Scoped to CN5 + rigid_planar + dent 3.  planar_bite overrides the env
            # gate (the never-worse caller builds BOTH ways); when None the env flag
            # DELFIN_FFFREE_PLANAR_MER_CN5 decides (else byte-identical: force_bite False).
            _eligible_rp_cn5 = bool(
                lg.get("rigid_planar") and _dent == 3
                and geometry in ("TBP-5 trigonal bipyramid", "SPY-5 square pyramid"))
            if planar_bite is None:
                _force_bite = bool(_eligible_rp_cn5 and _planar_mer_cn5_enabled())
            else:
                _force_bite = bool(_eligible_rp_cn5 and planar_bite)
            ring_confs = _embed_metallacycle(lg["mol"], _dons_d, metal,
                                             donor_target_pos=_dtp, harden=_harden,
                                             force_bite=_force_bite)
            # RIGID PLANAR polydentate (terpy / pincer, dent >= 3): the metallacycle
            # embed FOLDS the rigid backbone so the metal lifts ~1.0 A OUT of the
            # donors' own plane.  When enabled, build a metal-COPLANAR conformer
            # (the metal solved INTO the rigid donor plane) so the in-plane meridional
            # pose is realised at placement time.  planar_coplanar overrides the env
            # gate (the never-worse caller builds BOTH ways); when None the env flag
            # DELFIN_FFFREE_PLANAR_POLYDENTATE_PLACE decides (else byte-identical:
            # the coplanar conformer is never built).  Scoped to rigid_planar dent>=3;
            # falls back to the folded metallacycle embed if the coplanar solve fails.
            _eligible_rp_cop = bool(lg.get("rigid_planar") and _dent >= 3
                                    and _dtp is not None)
            _use_cop = (bool(_eligible_rp_cop and _planar_polydentate_place_enabled())
                        if planar_coplanar is None
                        else bool(_eligible_rp_cop and planar_coplanar))
            if _use_cop:
                _mds = [float(np.linalg.norm(p)) for p in _dtp]
                _cop = _coplanar_metal_centered_conformer(
                    lg["mol"], _dons_d, metal, _mds)
                if _cop is not None:
                    ring_confs = _cop
            elif (_dent >= 2 and _dtp is not None
                  and os.environ.get("DELFIN_FFFREE_PI_RIGID_PLACE", "0") == "1"):
                # PHASE 2 — FLEXIBLE multi-arm polydentate (NOT rigid_planar, e.g. a
                # tripodal with pyridyl arms on an sp3 backbone = 91 % of the metal-
                # out-of-plane defects).  Choose the backbone conformer that lets the
                # metal sit COPLANAR with EVERY aromatic arm's ring (per-ring solve;
                # _multiarm_coplanar.best_conformer).  Self-gating: returns None when
                # the donors are not in aromatic rings -> falls back to the metallacycle
                # embed.  Default-OFF / byte-identical when the flag is unset.
                try:
                    from delfin.manta._multiarm_coplanar import embed_coplanar_multiarm as _ma_emb
                    _ma = _ma_emb(
                        lg["mol"], _dons_d, metal,
                        [float(np.linalg.norm(p)) for p in _dtp],
                        donor_target_pos=_dtp)
                    if _ma is not None:
                        # ADD the coplanar-M conformer to the ensemble (preserve
                        # completeness — do NOT replace the existing conformers;
                        # V16 ADD-not-REPLACE lesson).  _ma is only returned when
                        # genuinely coplanar (accept-if-better gate in the module),
                        # so no regression: ranking/clean-gate pick the best frame.
                        if (ring_confs is not None and ring_confs[0] == _ma[0]):
                            ring_confs = (ring_confs[0],
                                          list(ring_confs[1]) + list(_ma[1]))
                        else:
                            ring_confs = _ma
                except Exception:
                    pass
            # UNIVERSAL LIGAND-INVARIANT cavity seating (DELFIN_FFFREE_RIGID_LIGAND_SEAT,
            # default OFF -> byte-identical).  For ANY rigid polydentate (dent>=3:
            # κ3 pincer / κ4 macrocycle / κ6 cage) replace the metallacycle/coplanar
            # conformer with the FREE-ligand conformer whose donor cavity best admits the
            # metal at the ideal M-D radii (3-D distance solve; NO backbone-flatness
            # requirement -> works where _coplanar_metal_centered_conformer bails).  Paired
            # with the rigid-body orient below (_rigid_seat forced True for this ligand) the
            # ligand internal geometry is preserved EXACTLY and the polyhedron is EMERGENT.
            # Overrides _use_cop / PI_RIGID_PLACE when set (the more general construction).
            # Self-gating: None on solve failure -> keeps the metallacycle embed (never-worse).
            _rls = (_dent >= 3 and _dtp is not None
                    and os.environ.get("DELFIN_FFFREE_RIGID_LIGAND_SEAT", "0") == "1")
            if _rls:
                _mds = [float(np.linalg.norm(p)) for p in _dtp]
                _cav = _rigid_ligand_cavity_conformer(lg["mol"], _dons_d, metal, _mds)
                if _cav is not None:
                    ring_confs = _cav
        if ring_confs is not None:
            lsyms, coords_list = ring_confs
        else:
            confs = _ligand_confs_from_mol(lg["mol"], k=_k_confs)
            if confs is None:
                return None
            lsyms, coords_list, lmol = confs
        # Build every valid placement for this ligand across the conformer pool.
        # Single mode: pick the clash-minimal vs already-placed atoms (historic,
        # byte-identical).  Ensemble mode: keep distinct low-clash-vs-metal
        # candidates (placement-order-independent + diverse), deduped intra-ligand.
        #
        # COLLAPSE-AWARE SELECTION (DELFIN_FFFREE_COLLAPSE_AWARE_SELECT, default OFF ->
        # byte-identical).  The pick below ranks conformers by CLASH alone, and is blind to
        # the one criterion the self-gate will later reject the whole complex for: a bonded
        # heavy-heavy pair below 0.82 x the covalent sum.  Measured 2026-08-02: 171 of the
        # 215 CHELATE_EMPTY systems die on exactly that, and the mode of the surviving
        # frames' tightest bond sits in the first bin ABOVE the floor -- the gate cuts
        # through the middle of the distribution, it does not trim a tail.
        #
        # _collapsed_heavy_bonds_strict is the gate's OWN predicate, same 0.82, same
        # bonded-pair rule, already unconditional -- it is merely scoped to rigid-planar
        # tridentates further down.  This widens the KNOWLEDGE of it to every chelate
        # without widening the REJECTION: a collapsed conformer is only DEPRIORITISED, never
        # discarded.  If every conformer collapses, the same one wins as today, so the set
        # of buildable systems cannot shrink -- never-worse by construction, not by measurement.
        _csel = os.environ.get("DELFIN_FFFREE_COLLAPSE_AWARE_SELECT", "0") == "1"
        # BETA-AWARE SELECTION (DELFIN_FFFREE_BETA_AWARE_SELECT, default OFF -> byte-id).
        #
        # pyramidal_sp2 -- the metal sitting OUT of its donor's own substituent plane -- is the
        # largest single defect class we carry: 18-19 % of frames against 0.98 % in 509 clean
        # crystals, a 19-fold over-representation.  And it is a PATH difference, not physics:
        # the monodentate path aligns the metal onto the donor's lone-pair axis by construction
        # (_vsepr_reconstruct + _rot_align) and reaches 0.9 deg, BETTER than the crystals; the
        # chelate path never looks at the donor plane at all and reaches 9.3 / 15.9 deg.
        #
        # WHY THIS IS A SELECTION AND NOT A CONSTRAINT.  beta is SECOND ORDER in distance:
        # with d^2 = r^2 + a^2 + 2*r*a*cos(beta)*cos(half) the first derivative vanishes at
        # beta=0, so d changes only as beta^2.  A 0.10 A bound still admits ~37 deg, a 0.05 A
        # bound ~26 deg, and our worst 15.9 deg costs all of 0.018 A.  Holding beta to 5 deg
        # would need a 0.002 A tolerance, which triangle smoothing will not survive.  NO set of
        # distances to atoms IN the donor plane can pin beta to first order -- the distance
        # from a point to coplanar anchors is stationary against perpendicular displacement.
        # (The in-plane azimuth theta IS first order, 0.014 A/deg; that one belongs in bounds.)
        #
        # WHY IT REACHES THE FRAME.  The per-donor radial rescale below scales a donor along
        # M->D, a line that lies IN the plane, so plane(D',X1,X2) is the same plane and still
        # contains M: beta is preserved.  The Kabsch fit is rigid, likewise.  So whatever beta
        # the chosen conformer has is the beta that ships -- the lever belongs here and nowhere
        # after.
        #
        # beta decides ONLY among conformers of EQUAL clash.  It cannot displace the criterion
        # that decides today, so it cannot introduce a clash regression at all.
        # ``prefer_beta`` lets a CALLER ask for the beta-optimal pick without touching the
        # environment.  THE POINT IS THAT IT IS ADDITIVE: the caller builds the frame twice
        # and appends the second as a SIBLING, so the primary is untouched by construction.
        #
        # That shape is not a preference, it is what the record says works.  Every flag that
        # LANDED in this project adds something -- D8_SQ_ADD "a PURELY ADDITIVE sibling ...
        # the PRIMARY frame is untouched", CN6_OH_ADD "adds the OC-6 isomer as a PURELY
        # ADDITIVE sibling", STEREOCENTER_ENUM "additively builds every buildable fold",
        # CN4_BOTH "native-additive".  Everything measured on 2026-08-02 that CHOSE instead
        # of ADDING died: beta-as-replacement cost 3-4 capabilities, collapse-as-replacement
        # died on all three pools, and five reach levers that replaced legacy's frame lost
        # 78/52/52/15/3.  The one lever still alive that night -- the ring pucker -- is the
        # only one that appends and leaves the primary byte-identical.
        _bsel = bool(prefer_beta) or os.environ.get(
            "DELFIN_FFFREE_BETA_AWARE_SELECT", "0") == "1"
        best_Q, best_clash = None, 1e18
        best_coll = True                            # a collapsed pick loses to a clean one
        best_pol = 10 ** 9                          # ... then the polyhedron fidelity (18.08.)
        best_beta = float("inf")                    # ... and among equals, the flatter donor
        # Default OFF -> _pol is 0 for EVERY candidate, so the tuple compare falls
        # back exactly onto the historic order: byte-identical.
        _polysel = os.environ.get("DELFIN_FFFREE_POLY_FIDELITY_SEAT", "0") == "1"
        cands = []                                  # (Q, clash_vs_metal) for ensemble
        seen_local = []                             # intra-ligand RMSD dedup
        for lP in coords_list:
            dent_targets = None       # only the polydentate branch sets target vertices
            if lg["denticity"] == 1:
                v = va[0][0]
                Vunit = ref[v] / np.linalg.norm(ref[v])
                md = MSB.md_distance(metal, lsyms[dons[0]],
                                     atom=lg["mol"].GetAtomWithIdx(int(dons[0])),
                                     mol=lg["mol"])
                if len(lsyms) == 1:
                    Q = (Vunit * md).reshape(1, 3)
                else:
                    lPv, lp = _vsepr_reconstruct(lsyms, lP, lmol, dons[0])
                    Q = (lPv - lPv[dons[0]]) @ _rot_align(lp, -Vunit).T + Vunit * md
                    # Post-placement diatomic-orientation guard (env-gated, default OFF
                    # => byte-identical): for a linear diatomic donor (M-C#O carbonyl /
                    # M-C#N cyanide / M-N=O nitrosyl) make sure the SMILES-bonded DONOR
                    # atom — not its partner — sits at the vertex.  A shared conformer
                    # cache (keyed by canonical SMILES) can hand back a block whose atom
                    # ordering differs from this ligand's graph, flipping the unit to an
                    # M-O-C isocarbonyl; the guard reflects it donor-first.  Connectivity-
                    # derived donor element, no per-case logic.
                    if _diatomic_orient_enabled():
                        _dp = _diatomic_donor_partner(lg)
                        if _dp is not None:
                            Q = _orient_diatomic_block(
                                Q, lsyms, _dp[0], _dp[1], np.zeros(3), Vunit * md)
            else:
                # chelate (bi-/tri-dentate): the d assigned vertices host the d donors;
                # arm index i (config order) seats on verts[i] in FIXED correspondence
                # (config-faithful), using the canonical element-sorted arm order so the
                # enumeration bijects onto the distinct stereoisomers.
                dent = lg["denticity"]
                dons_d = _canonical_arm_order(lg, dent)
                verts = [v for v, arm in sorted(va, key=lambda x: x[1])]
                targets = [ref[verts[i]] / np.linalg.norm(ref[verts[i]])
                           * MSB.md_distance(metal, lsyms[dons_d[i]],
                                             atom=lg["mol"].GetAtomWithIdx(int(dons_d[i])),
                                             mol=lg["mol"]) for i in range(dent)]
                # The target vertices of this chelate -- from here on the selection knows
                # WHERE the donors actually belong, and can score fidelity to that.
                dent_targets = targets
                # config-faithful seating is only meaningful for ASYMMETRIC chelates
                # (distinct donor elements); for symmetric chelates every arm seating
                # is the same isomer, so keep the best-fit (lower ring strain).
                _de = [lsyms[di] for di in dons_d]
                _asym = len(set(_de)) > 1
                # RIGID PLANAR tridentate: the CENTRAL donor MUST seat on the central
                # meridian vertex (arm 1 -> targets[1]); the best-fit permutation would
                # scramble central/outer (a flat terpy is element-symmetric all-N), so
                # force config-faithful seating to realise the meridional arrangement.
                if lg.get("rigid_planar") and dent == 3:
                    _asym = True
                _rigid_planar = bool(lg.get("rigid_planar")) and dent == 3
                Q = None
                # κ4 chain: when the chelate's metallacycle embed failed (ring_confs
                # is None, e.g. an unkekulizable aromatic-N⁺ scorpionate cap) the
                # tridentate was NEVER seated -> skipped -> whole complex to legacy.
                # Seat the FREE-ligand conformer too (Kabsch-fit the rigid donor
                # triangle onto the assigned vertices).  Gated default-OFF byte-id.
                _orient_free = (ring_confs is None
                                and os.environ.get("DELFIN_FFFREE_KEKULIZE_SPLIT", "0") == "1")
                if ring_confs is not None or _orient_free:
                    # rigid-body orient (no per-donor rescale) when LIGAND_RIGID OR the
                    # universal cavity seating (RIGID_LIGAND_SEAT) is active for a dent>=3
                    # rigid polydentate -> ligand internal geometry preserved, polyhedron
                    # emergent.  Default OFF on both -> byte-identical (per-donor rescale).
                    _rigid_seat = (dent >= 3 and (
                        os.environ.get("DELFIN_FFFREE_LIGAND_RIGID", "0") == "1"
                        or os.environ.get("DELFIN_FFFREE_RIGID_LIGAND_SEAT", "0") == "1"))
                    Q = _orient_chelate_to_vertices(lP, dons_d, targets, asym=_asym,
                                                    rigid=_rigid_seat, lsyms=lsyms)
                    # THE LEFTOVER ROTATION, SET BY THE LONE PAIRS INSTEAD OF BY THE SVD.
                    # Seating two donors onto two vertices leaves exactly one rotational
                    # freedom -- about the donor-donor axis -- and the Kabsch fit cannot see
                    # it, because both donors lie ON that axis and do not move.  See
                    # _lp_orient_seated_bidentate for why that one angle IS beta, and for the
                    # 995-system measurement behind the law.  Donors, M-D distances, bite and
                    # arm-to-vertex correspondence are all untouched; only the backbone turns.
                    #
                    # AS A SEATING THIS IS MEASURED NEGATIVE, AS A SIBLING IT IS OPEN.
                    # lpseat2 (187 systems, the rotation applied unconditionally): reach 49,
                    # but cap_LOST 2, valid 28->27, mean +0.181 and ten red terms.  The turned
                    # backbones do collide, the self-gate drops those colorings, and the whole
                    # system hands over to legacy.  So the source's own claim at
                    # _lp_aligned_angle -- "the plane wins, without a steric escape" -- is
                    # refuted AT THE LIVE SEATING: a bonding argument does not survive contact
                    # with a co-ligand that is also real.
                    # A guarded version is no answer either: with the collapse guard the same
                    # run measured affected=0, i.e. exactly the "changed almost nothing" the
                    # docstring already recorded once.  Guarded it does nothing, unguarded it
                    # costs capability.
                    # Hence ``lp_orient``: the CALLER asks for the lp-oriented frame and
                    # appends it as a SIBLING, leaving the primary untouched.  Then the plane
                    # can win in a frame of its own without any config ever being dropped,
                    # which is the one shape that has ever landed here.
                    if (Q is not None and dent == 2
                            and (lp_orient
                                 or os.environ.get("DELFIN_FFFREE_LP_SEAT", "0") == "1")):
                        _Ql = _lp_orient_seated_bidentate(Q, lg.get("mol"),
                                                          dons_d[0], dons_d[1], lsyms=lsyms)
                        if _Ql is not None:
                            Q = _Ql
                    # Collapse-rigid-fallback (DELFIN_FFFREE_CHELATE_RIGID_FALLBACK,
                    # default OFF -> byte-id).  The per-donor RADIAL rescale moves each
                    # donor independently onto its exact ideal radius, which on many
                    # chelates distorts the donor-backbone bonds into a collapse near
                    # the metal (the #1 build-coverage gap: 15% selfgate-collapse, 77%
                    # within 3A of the metal).  When the rescaled seating collapsed,
                    # retry the RIGID-body Kabsch fit (no rescale): ligand geometry is
                    # preserved exactly, the M-D distances come out emergent (cavity-
                    # dictated, refined downstream), and the backbone does NOT collapse
                    # -> the config is recovered instead of skipped to legacy.
                    if (Q is not None and not _rigid_seat
                            and os.environ.get("DELFIN_FFFREE_CHELATE_RIGID_FALLBACK", "0") == "1"
                            and _has_collapsed_heavy_bonds(lsyms, Q)):
                        _Qr = _orient_chelate_to_vertices(lP, dons_d, targets,
                                                          asym=_asym, rigid=True, lsyms=lsyms)
                        if (_Qr is not None and np.all(np.isfinite(_Qr))
                                and not _has_collapsed_heavy_bonds(lsyms, _Qr)):
                            Q = _Qr
                    # iter-32e (YILNUF oxalate-collapse class): the DG metallacycle
                    # embed sometimes produces a backbone with collapsed C-C / C-O bonds.
                    # _orient_chelate_to_vertices preserves that (rigid Kabsch+rescale),
                    # so the resulting Q poisons _build_is_clean and the whole complex
                    # falls back to legacy.  Reject such Q here → the bidentate rigid
                    # fallback below picks up; for tridentate we skip the config (the
                    # build-gate would reject anyway).  Env-gated default OFF.
                    if Q is not None and _has_collapsed_heavy_bonds(lsyms, Q):
                        Q = None
                    # RIGID PLANAR tridentate (DELFIN_FFFREE_PLANAR_MER): the DG
                    # metallacycle embed yields chelating conformers (outer-outer ~150-
                    # 158deg) but a few carry a collapsed donor-backbone bond.  Reject
                    # such conformers UNCONDITIONALLY (not env-gated) so a CLEAN one is
                    # selected from the pool -> the meridional placement survives the
                    # self-gate instead of bailing the whole complex to legacy.
                    if Q is not None and _rigid_planar and _collapsed_heavy_bonds_strict(lsyms, Q):
                        Q = None
                if Q is None and dent == 2:         # embed/orient failed -> rigid fallback (bidentate only)
                    Q = _place_chelate_block(metal, lsyms, lP, dons_d[0], dons_d[1],
                                             targets[0], targets[1], mol=lg["mol"])
                if Q is None:                       # tridentate embed failure -> skip this config
                    continue
            if not np.all(np.isfinite(Q)):
                continue
            cl = _clash_count(Q, np.array(placed), lsyms, placed_syms)
            # THE FREE DEGREE OF FREEDOM, CALL SITE 1 of 2 (the main road: this
            # path builds every chelate configuration).  It acts ONLY if this block
            # already collides with something already seated -- exactly the case in
            # which the register says "the ligand pays".  Default OFF -> `_fd` is always
            # None, `Q` and `cl` stay untouched: byte-identical.
            _fd = _free_dof_reseat(Q, lsyms, dons, np.array(placed), placed_syms, cl)
            if _fd is not None:
                Q, cl = _fd
            # (collapsed, clash) beats (clash) alone: a clean conformer outranks a colliding
            # one, and among equals the historic clash order is untouched.  With the flag OFF
            # _coll is False for every candidate, so the tuple compare degenerates to the
            # historic `cl < best_clash` exactly -- byte-identical.
            _coll = bool(_csel and _collapsed_heavy_bonds_strict(lsyms, Q))
            _bta = _beta_score(lsyms, Q, dons) if _bsel else 0.0
            # ===== POLYHEDRON FIDELITY AS A SELECTION CRITERION ================
            # THE FINDING THAT FORCES THIS TERM (18.08., 30921 systems): net
            # +988 systems flow from the octahedron into the trigonal prism, McNemar
            # X2 = 860.8; the builder produces 2.99 times too many prisms, while every
            # other shape stays between 0.86 and 1.29.  They are not prisms --
            # CShM median 11.03 instead of 16.7, i.e. a HALF Bailar twist.  And the
            # signal is the INTERLOCKING, monotone from 2.49 % at zero chelate rings to
            # 16.12 % at five.  Metal and d-count are flat.
            #
            # ⚠ WHY HERE AND NOT IN THE RELAXATION -- that is DECIDED today,
            # not suspected.  The assertion protocol was hung on all four
            # refine call sites and caught the twist six times (109 calls,
            # all six trans).  Nevertheless 19 of 19 archive files stayed
            # byte-identical.  The discriminator: the rollback additionally shifted by 100
            # Angstrom -- a frame that every gate rejects and every
            # byte comparison would see.  The archive stayed identical.  ⇒ The frames whose
            # assertion breaks are rejected anyway; the SURVIVORS are
            # already twisted when they arrive at refine.  The relaxation is thereby
            # excluded as the cause -- the chelate pull acts HERE, at
            # placement.
            #
            # The builder puts the donors onto `targets` (ideal vertices times
            # md_distance), but a rigid chelate ring cannot reach them all:
            # _orient_chelate_to_vertices lays it on as well as possible, and the
            # rest is deviation.  Until now NOTHING decided about that -- the selection
            # read (collapse, clash, beta).  A conformer that twists the polyhedron by 30
            # degrees won against a faithful one as soon as it had one clash less
            # or even just came earlier.
            #
            # ⚠ WHY THIS IS CHEAP: pure SELECTION among already built
            # candidates -- no new geometry arises.  By the law measured
            # today (isometry +0.98 pp, rigid rotation +6.57 pp,
            # re-embedding +11.9 pp) that is the cheapest class of all.
            #
            # ⚠ WHY BINNED: `_pol` is a floating-point number.  Unbinned it would decide
            # practically every comparison and `_bta` would never get its turn again --
            # that would be a silent switch-off of a landed term.  The band of
            # 0.05 Angstrom lets it speak only when the difference is large
            # enough to mean something chemically.
            #
            # DELFIN_FFFREE_POLY_FIDELITY_SEAT (default 0 -> byte-identical).
            _pol = 0
            if _polysel and dent_targets is not None:
                try:
                    _dev = [float(np.linalg.norm(Q[dons_d[i]] - dent_targets[i]))
                            for i in range(len(dent_targets))]
                    _pol = int(round((sum(d * d for d in _dev) / len(_dev)) ** 0.5
                                     / _POLY_FIDELITY_BIN))
                except Exception:
                    _pol = 0
            if (_coll, cl, _pol, _bta) < (best_coll, best_clash, best_pol, best_beta):
                best_coll, best_clash, best_pol, best_beta, best_Q = (
                    _coll, cl, _pol, _bta, Q)
            if ensemble:
                # dedup conformers of THIS ligand by intra-ligand RMSD (identity corr.)
                dup = False
                for Qs in seen_local:
                    if Qs.shape == Q.shape:
                        c0 = Q - Q.mean(axis=0); c1 = Qs - Qs.mean(axis=0)
                        if float(np.sqrt(((c0 - c1) ** 2).sum(axis=1).mean())) < rmsd_dedup:
                            dup = True
                            break
                if not dup:
                    seen_local.append(Q)
                    cands.append((Q, _clash_count(Q, metal_P, lsyms, metal_sym)))
            elif cl == 0 and not _coll and not _bsel:
                # the early exit must not fire on a conformer that is clash-free but carries
                # a bond the gate will kill the whole complex for -- that is the exact trade
                # this flag exists to stop.  Flag OFF -> _coll is False -> historic break.
                # It must also not fire under BETA_AWARE_SELECT at all: taking the FIRST
                # clash-free conformer is precisely what makes the beta tie-break unable to
                # see the others.  At most per_lig_confs candidates, so the cost is bounded.
                break
        if best_Q is None:                  # no conformer could be placed -> bail to legacy
            return None
        if ensemble:
            if not cands:
                return None
            cands.sort(key=lambda t: t[1])          # low clash-vs-metal first
            per_lig_cands.append(cands)
            per_lig_syms.append(lsyms)
        out_syms += lsyms
        placed_syms += lsyms
        relax_frags.append((Chem.AddHs(lg["mol"]), lig_offset))   # bonds for internal relax
        # in single mode the placed coords come from best_Q; the ensemble assembles
        # placed coords per combo below (placed[] is only the single-mode buffer)
        if not ensemble:
            for row in best_Q:
                placed.append(row)
        # freeze the donor atoms at their vertices
        for d in dons:
            fixed.add(pos + d)
        # ... and record which atoms belong to this one ligand.  Pure bookkeeping: nothing
        # reads it unless DELFIN_FFFREE_GLOBAL_DONORS is on, so the frame is unaffected.
        lig_blocks.append((pos, len(lsyms), sorted(pos + int(d) for d in dons)))
        pos += len(lsyms)
    donors = sorted(fixed - {0})              # global indices of the constructed donor atoms
    if not ensemble:
        P = np.vstack([np.zeros((1, 3))] + [np.array(placed[1:], float)])
        # OC-6 TWIST CORRECTION IN THE SEATING (only on ``oc6_twist=True``; the
        # keyword is False by default, so the primary frame is byte-identical,
        # independent of any environment variable).  It stands BEFORE the global
        # donor seating, because the drift permitted there (0.05 A) is measured against
        # the frame it finds: first set the polyhedron right, then de-distort --
        # the other way round the de-distortion would have to give up its own work again.
        # And BEFORE _finish_config_frame, because the relaxation there nails down the donors.
        if oc6_twist:
            _t = _oc6_twist_seat(out_syms, P, lig_blocks, metal, geometry)
            if _t is not None:
                P = _t
        # GLOBAL DONOR SEATING (default OFF -> byte-identical).  Every ligand up to here was
        # seated ALONE; this is the first and only point where all of them exist at once, so
        # it is the first point where an inter-ligand distance can even be written down.  It
        # runs BEFORE _finish_config_frame on purpose: the relax there pins the donors, so
        # whatever this leaves them at is what the coordination is finished around.
        if _global_donor_seat_enabled():
            _g = _global_donor_seat(out_syms, P, lig_blocks)
            if _g is not None:
                P = _g
        P = _finish_config_frame(out_syms, P, fixed, relax_frags, refine, geom=geometry)
        return out_syms, P, donors

    # ENSEMBLE assembly (Task A.1): enumerate the Cartesian product of per-ligand
    # conformer choices, greedily ordered (sum of candidate ranks -> the best
    # combinations first = frame 0 == the single-path pick), refine each, RMSD-dedup
    # at the complex level, keep up to n_frames.  Capped product keeps it bounded.
    import itertools as _it
    rank_lists = [list(range(len(c))) for c in per_lig_cands]
    _nprod = 1
    for _rl in rank_lists:
        _nprod *= max(1, len(_rl))
        if _nprod > _COMBO_MATERIALISE_MAX:
            break
    if _nprod > _COMBO_MATERIALISE_MAX:      # never materialise a 10^8 product to use 64
        combos = _ranked_combos(rank_lists, 256)
    else:
        combos = list(_it.product(*rank_lists))
    combos.sort(key=lambda cb: (sum(cb), cb))      # deterministic; frame 0 = all-best
    MAX_EVAL = 64
    frames = []                                    # (syms, P) kept (deduped)
    # FOLD FINGERPRINT (see block at `_fold_fp_enabled`).  Default OFF ->
    # `_fold_rings` stays None -> the predicate below is literally the old one.
    # `relax_frags` carries (AddHs(lg.mol), LIGANDS-ONLY offset); globally the
    # block starts at `lig_offset + 1`, because index 0 is the metal -- the same convention
    # `_collect_exempt` and `_finish_config_frame` already use.
    # THIRD ENTRY = the chelate arms for the metallacycle (switch, see above).  The
    # ligand for it is in `ligands[_lig_order[i][0]]`: the loop above appends per
    # iteration EXACTLY ONE `relax_frags` entry (:5550, no `continue`
    # before it), so index i is the same.  ⚠ That is an assumption about the
    # loop, therefore it is CHECKED -- if the length does not match, there are no
    # chelate arms and the fingerprint is exactly the organic one from before.
    _fold_rings = None
    if _fold_fp_enabled():
        _fb = [(int(o) + 1, m, None) for (m, o) in relax_frags]
        if len(relax_frags) == len(_lig_order):
            _fb = [(int(o) + 1, m, _fold_mc_arms(ligands[_lig_order[i][0]]))
                   for i, (m, o) in enumerate(relax_frags)]
        _fold_rings = _fold_rings_with_mc(_fb, out_syms, 0)
    _fold_kept = []                                # fingerprint per kept frame
    for cb in combos[:MAX_EVAL]:
        blocks = [np.zeros((1, 3))]
        ok = True
        for vi, ci in enumerate(cb):
            Q = per_lig_cands[vi][ci][0]
            blocks.append(Q)
        Pc = np.vstack(blocks)
        if not np.all(np.isfinite(Pc)):
            continue
        # the same twist correction for the ensemble combinations: the ligand spans
        # are the same across all combinations (same conformer atom counts), so lig_blocks
        # applies unchanged.  Again only on the keyword.
        if oc6_twist:
            _t = _oc6_twist_seat(out_syms, Pc, lig_blocks, metal, geometry)
            if _t is not None:
                Pc = _t
        # same global seating for the ensemble combos; the per-ligand spans are identical
        # across combos (same conformer atom counts), so lig_blocks applies unchanged.
        if _global_donor_seat_enabled():
            _g = _global_donor_seat(out_syms, Pc, lig_blocks)
            if _g is not None:
                Pc = _g
        Pc = _finish_config_frame(out_syms, Pc, fixed, relax_frags, refine, geom=geometry)
        if not np.all(np.isfinite(Pc)):
            continue
        dup = False
        _fp = None                                 # lazy: only on RMSD proximity
        for _ki, (_, Pk) in enumerate(frames):
            if Pk.shape == Pc.shape and _complex_rmsd(out_syms, Pc, Pk) < rmsd_dedup:
                if _fold_rings is None:
                    dup = True                     # switch OFF -> old predicate
                    break
                if _fp is None:
                    _fp = _fold_fp(Pc, _fold_rings)
                if _fold_kept[_ki] is None:
                    _fold_kept[_ki] = _fold_fp(Pk, _fold_rings)
                if _fold_same(_fold_kept[_ki], _fp):
                    dup = True
                    break
        if dup:
            continue
        frames.append((list(out_syms), Pc))
        _fold_kept.append(None)                    # same length as `frames`
        if len(frames) >= int(n_frames):
            break
    if not frames:
        return None
    return [(syms, P, donors) for syms, P in frames]


def _finish_config_frame(out_syms, P, fixed, relax_frags, refine=True, geom=None):
    """Shared finishing tail for a single assembled chelate-config frame: the
    FUNDAMENTAL internal relaxation (division-of-labor doctrine) — build the
    ligands-only mol (NO metal) at the placed coords and UFF-relax it with the
    DONORS FIXED + inter-fragment vdW ON.  All organic internals (bond lengths,
    angles, funcgroup/aromatic planarity, H geometry) recover to MM-ideal and
    inter-ligand clashes resolve, while the constructed coordination is preserved
    (donors pinned -> M-D invariant; metal never enters the force field).  This
    replaces the weak defect-count refine(), which left ETKDG-rough internals.
    INTERNAL finishing.  THREE-WAY RACE design:
      Track 3 (FF-FREE, DEFAULT): pure construction -- GEOMETRIC correctors only
        (sp2-flatten projection + defect-count clash-relief).  No force field.
      Track 2 (ligand-FF, opt-in DELFIN_FFFREE_LIGANDFF=1): additionally a constrained
        UFF relax of the organic internals (metal + donors frozen) before the geometric
        correctors.  Track 1 (UFF) is the legacy path (DELFIN_FFFREE_BUILDER=0).
    The complex mol (metal + ligand bonds) is the sp2-flatten template (and the optional
    UFF mol); metal + donors stay frozen so the constructed coordination is preserved.
    Returns the (possibly relaxed) coordinate array; never raises."""
    if not refine:
        return P
    try:
        cm = Chem.RWMol()
        cm.AddAtom(Chem.Atom(out_syms[0]))
        for frag, _off in relax_frags:
            base = cm.GetNumAtoms()
            for a in frag.GetAtoms():
                cm.AddAtom(Chem.Atom(a.GetAtomicNum()))
            for b in frag.GetBonds():
                cm.AddBond(b.GetBeginAtomIdx() + base, b.GetEndAtomIdx() + base, b.GetBondType())
        donor_globals = sorted(fixed - {0})
        for dg in donor_globals:
            cm.AddBond(0, dg, Chem.BondType.SINGLE)
        if cm.GetNumAtoms() == len(out_syms):
            conf = Chem.Conformer(cm.GetNumAtoms())
            for kk in range(cm.GetNumAtoms()):
                x = P[kk]; conf.SetAtomPosition(kk, [float(x[0]), float(x[1]), float(x[2])])
            cm.AddConformer(conf, assignId=True)
            try:
                Chem.SanitizeMol(cm, catchErrors=True)
            except Exception:
                pass
            if os.environ.get("DELFIN_FFFREE_LIGANDFF", "0") == "1":   # Track 2: FF relax
                if _constrained_uff_relax(cm, [0] + donor_globals):
                    Pn = np.array(cm.GetConformer().GetPositions(), float)
                    if Pn.shape == P.shape:
                        P = Pn
                        c = cm.GetConformer()
                        for kk in range(len(P)):
                            c.SetAtomPosition(kk, [float(P[kk][0]), float(P[kk][1]), float(P[kk][2])])
            # FF-FREE geometric planarity corrector (BOTH tracks): project sp2 atoms onto
            # their neighbour plane (aromatic/amide/funcgroup planarity).  Pure geometry.
            try:
                from delfin.smiles_converter import _flatten_sp2_atoms_xyz
                xyz_str = "\n".join(f"{s} {float(p[0]):.6f} {float(p[1]):.6f} {float(p[2]):.6f}"
                                    for s, p in zip(out_syms, P))
                flat = _flatten_sp2_atoms_xyz(xyz_str, cm)
                if flat:
                    newP = np.array([[float(x) for x in ln.split()[1:4]]
                                     for ln in flat.splitlines() if ln.strip()], float)
                    if newP.shape == P.shape:
                        P = newP
            except Exception:
                pass
    except Exception:
        pass
    # FF-free geometric clash-relief (both tracks)
    P = _refine_guarded(out_syms, P, fixed)
    # #308 whole-complex torsion-space clash relax (env-gated, default-OFF byte-id):
    # joint multi-axis torsion of all rotatable single bonds, metal + donors (`fixed`)
    # frozen.  Chelate ring + M-D arms are ring bonds -> kept rigid by construction;
    # only the σ-arms/substituents rotate.  Torsion-only, never-worse.  No-op when unset.
    try:
        from delfin.manta import torsion_relax as _TR
        bp = None
        try:
            # true connectivity: ligand-internal bonds (lig_offset shifts the AddHs
            # mol's local indices to global = lig_offset+1, atom 0 = metal) + M-D bonds
            # to every constructed donor (so the relaxer keeps chelate metallacycles
            # rigid via ring detection and never mis-perceives a crowded contact).
            pairs = []
            for frag, lig_off in relax_frags:
                base = lig_off + 1                 # global index of the frag's atom 0
                for b in frag.GetBonds():
                    pairs.append((base + b.GetBeginAtomIdx(), base + b.GetEndAtomIdx()))
            for dg in sorted(fixed - {0}):
                pairs.append((0, dg))
            bp = pairs or None
        except Exception:
            bp = None
        P = np.asarray(_TR.relax_if_enabled(out_syms, P, fixed, bond_pairs=bp),
                       dtype=float)
        # JOINT inter-ligand declash (env-gated, default-OFF byte-id): global
        # inter-ligand heavy-heavy minimisation, core frozen.  After #308.  Reuses
        # the same true-connectivity bond list.  No-op when the flag is unset.
        from delfin.manta import joint_declash as _JD
        P = np.asarray(_JD.declash_if_enabled(out_syms, P, fixed, geom=geom, bond_pairs=bp),
                       dtype=float)
        # SOFT coordination-sphere radial flex (env-gated, default-OFF byte-id):
        # let crowded monodentate ligands translate radially outward a bounded
        # amount to open residual mild inter-ligand clashes (real-crystal 0.85*vdw
        # target).  This is the assemble_from_config path (the main metal-complex
        # builder); never-worse on clash.  No-op when DELFIN_FFFREE_SPHERE_FLEX unset.
        from delfin.manta import sphere_flex as _SF
        P = np.asarray(_SF.flex_if_enabled(out_syms, P, fixed, bond_pairs=bp),
                       dtype=float)
    except Exception:
        pass
    return P


def assemble_heteroleptic(metal: str, geometry: str, vertex_specs):
    """vertex_specs[i] = (ligand_smiles, donor_idx) for polyhedron vertex i.
    Heteroleptic monodentate assembly (different ligand per vertex) — the basis
    for building enumerated coordination isomers (e.g. cis/trans MA4B2)."""
    ref = MSB._ref_vectors(geometry)
    if len(vertex_specs) != len(ref):
        raise ValueError("vertex_specs count != vertices")
    out_syms = [metal]; blocks = [np.zeros((1, 3))]
    for i, (smi, di) in enumerate(vertex_specs):
        Vunit = ref[i] / np.linalg.norm(ref[i])
        lsyms, lP, lmol = _ligand_3d(smi)
        md = MSB.md_distance(metal, lsyms[di],
                             atom=lmol.GetAtomWithIdx(di), mol=lmol)
        vertex = Vunit * md
        if len(lsyms) == 1:                       # monatomic ligand (e.g. Cl-)
            out_syms += lsyms; blocks.append(vertex.reshape(1, 3)); continue
        lp = _donor_and_lp(lsyms, lP, lmol, di)
        R = _rot_align(lp, -Vunit)
        Q = (lP - lP[di]) @ R.T + vertex
        out_syms += lsyms; blocks.append(Q)
    return out_syms, np.vstack(blocks)


def generate_complex_conformers(metal, ligand_smiles, donor_idx, geometry,
                                 max_lig_conf=6):
    """Wire L3 into the complex: enumerate the LIGAND's conformers (conformer_enum),
    build the homoleptic complex from each, dedup. Demonstrates conformers-per-isomer
    (the L3 dimension of generate-gate-floor)."""
    from delfin.manta import conformer_enum as CE
    _, confs = CE.enumerate_conformers(ligand_smiles)
    ref = MSB._ref_vectors(geometry); n = len(ref)
    out = []
    for e, m in confs[:max_lig_conf]:
        lsyms = [a.GetSymbol() for a in m.GetAtoms()]
        lP = m.GetConformer().GetPositions()
        md = MSB.md_distance(metal, lsyms[donor_idx],
                             atom=m.GetAtomWithIdx(donor_idx), mol=m)
        lp = _donor_and_lp(lsyms, lP, m, donor_idx)
        syms = [metal]; blocks = [np.zeros((1, 3))]
        for i in range(n):
            Vunit = ref[i] / np.linalg.norm(ref[i])
            R = _rot_align(lp, -Vunit)
            Q = (lP - lP[donor_idx]) @ R.T + Vunit * md
            syms += lsyms; blocks.append(Q)
        out.append((round(e, 2), syms, np.vstack(blocks)))
    return out


# ---------------------------------------------------------------------------
# Self-test — does the DERIVED bite actually preserve the ligand's own d_DD?
# ---------------------------------------------------------------------------

def _self_test_bite_law() -> None:
    """The claim is not "the angle is nicer" -- it is that d_DD is PRESERVED.

    mdAB measured what happens when it is not: correcting the M-D length while the vertex
    angle stays at the ideal 90 deg forces d_DD = |u1*r1 - u2*r2| to a value the ligand
    does not have, and the ligand absorbs the difference (valid 155 -> 146).  So the test
    that matters compares the donor-donor distance the ligand HAS against the one the
    seating ASKS FOR, with the law off and on.  A lever that leaves that gap unchanged
    would be dead code that looks alive.
    """
    cases = [
        ("NCCN",              [0, 3], "ethylenediamine, sp3 backbone"),
        ("c1ccc(-c2ccccn2)nc1", None, "2,2'-bipyridine, aromatic backbone"),
        ("CC(=O)CC(C)=O",     None,   "acetylacetone, 6-ring"),
    ]
    # WHERE THE MISMATCH ACTUALLY SHOWS.  _place_chelate_block fits the ligand RIGIDLY, so
    # d_DD is preserved either way -- the first version of this test measured that and
    # learned nothing.  A rigid body cannot satisfy two vertex directions AND two M-D
    # lengths when the ligand's own bite disagrees with the vertex spacing, so the
    # compromise lands in the M-D DISTANCES: the donors end up nearer or further than the
    # length that was requested.  That is also the sharper reading of mdAB's failure --
    # setting r is pointless if the rigid fit then moves the donor off it again.
    print("derived bite -- does the seat DELIVER the M-D length it asked for?")
    print(f"{'ligand':>34} | {'d_DD':>6} | {'M-D err off':>11} | {'M-D err on':>10} | "
          f"{'angle off':>9} {'angle on':>8} | OK")
    metal_sym = "Ni"
    for smi, dons, label in cases:
        try:
            lsyms, lP, lmol = _ligand_3d(smi)
        except Exception as exc:
            print(f"{label:>34} | SKIP ({type(exc).__name__})")
            continue
        if dons is None:                      # first two N/O heavy donors in the SMILES
            dons = [a.GetIdx() for a in lmol.GetAtoms()
                    if a.GetSymbol() in ("N", "O")][:2]
        if len(dons) != 2:
            print(f"{label:>34} | SKIP (no donor pair)")
            continue
        d1, d2 = dons
        dd_lig = float(np.linalg.norm(lP[d1] - lP[d2]))
        _want1 = MSB.md_distance(metal_sym, lsyms[d1],
                                 atom=lmol.GetAtomWithIdx(d1), mol=lmol)
        _want2 = MSB.md_distance(metal_sym, lsyms[d2],
                                 atom=lmol.GetAtomWithIdx(d2), mol=lmol)
        errs, angs = [], []
        for flag in ("0", "1"):
            prev = os.environ.get("DELFIN_FFREE_BITE_LAW")
            os.environ["DELFIN_FFREE_BITE_LAW"] = flag
            try:
                # OC-6 vertices 0/1 are TRANS ([1,0,0] / [-1,0,0]); a chelate needs a CIS
                # edge, so 0/2 ([1,0,0] / [0,1,0], 90 deg apart).
                syms, P = assemble_multichelate(metal_sym, "OC-6 octahedron",
                                                [(smi, [d1, d2], [0, 2])])
                # metal sits at row 0; the ligand block follows in ligand order
                v1, v2 = P[1 + d1] - P[0], P[1 + d2] - P[0]
                g1, g2 = float(np.linalg.norm(v1)), float(np.linalg.norm(v2))
                errs.append(max(abs(g1 - _want1), abs(g2 - _want2)))
                angs.append(math.degrees(math.acos(max(-1.0, min(1.0,
                            float(np.dot(v1, v2)) / (g1 * g2))))))
            except Exception as exc:
                errs.append(float("nan"))
                angs.append(float("nan"))
                print(f"    ({label}: {type(exc).__name__}: {exc})")
            finally:
                if prev is None:
                    os.environ.pop("DELFIN_FFREE_BITE_LAW", None)
                else:
                    os.environ["DELFIN_FFREE_BITE_LAW"] = prev
        # NEVER-WORSE is the criterion, not "always changes something".  Three outcomes are
        # all correct: the law repairs a real gap (en), it reproduces what an already-exact
        # rigid fit does (bipy -- a confirmation, not a no-op), or it DECLINES because the
        # triangle does not close (acac, whose freely embedded conformer is the open form,
        # not the chelating one -- d_DD must come from the chelating basin).
        _c = ((_want1 ** 2 + _want2 ** 2 - dd_lig ** 2) / (2.0 * _want1 * _want2)
              if _want1 > 0 and _want2 > 0 else 2.0)
        note = ("declined (cos=%.2f)" % _c if not -1.0 <= _c <= 1.0
                else ("repairs" if errs[0] - errs[1] > 1.0e-6 else "reproduces"))
        ok = errs[1] <= errs[0] + 1.0e-9
        print(f"{label:>34} | {dd_lig:6.3f} | {errs[0]:11.4f} | {errs[1]:10.4f} | "
              f"{angs[0]:9.2f} {angs[1]:8.2f} | {str(ok):>5} {note}")


def _run_self_tests() -> None:
    _self_test_bite_law()
    _self_test_lp_orient()
    _self_test_trilateration()
    _self_test_oc6_twist()
    _self_test_ligand_dof()


def _self_test_ligand_dof() -> None:
    """HOW MUCH FREEDOM DOES THE LIGAND STILL HAVE -- and does using it cost the sphere?

    Four questions, and each may say NO:
      1) CENSUS: how many free rigid-body rotations remain per denticity?
         Expectation from reading the code: monodentate 1, bidentate 1, tridentate 0.
      2) PRESERVATION: under this rotation do donor position, M-D distance, bite and
         the ENTIRE internal ligand geometry stay exact?  Measured in Angstrom, not
         claimed.  (Isometry: only rounding residues may remain.)
      3) DEFAULT OFF: does `_free_dof_reseat` without the switch guarantee None?
      4) REACH ON: does the rotation resolve overlap on a REAL, truly colliding
         ligand pair -- or does it find nothing?  Either is a finding.
    """
    print("\nfreier Ligand-Freiheitsgrad -- Zensus, Erhaltung, Vorgabe, Reichweite")
    metal = "Fe"
    ref = MSB._ref_vectors("OC-6 octahedron")
    _ANG = 137.0                       # odd angle: no symmetry coincidence

    def _seat_mono(smi, want=("N", "P", "S", "O")):
        lsyms, lP, lmol = _ligand_3d(smi)
        di = next(a.GetIdx() for a in lmol.GetAtoms() if a.GetSymbol() in want)
        V = ref[0] / np.linalg.norm(ref[0])
        md = MSB.md_distance(metal, lsyms[di], atom=lmol.GetAtomWithIdx(di), mol=lmol)
        lPv, lp = _vsepr_reconstruct(lsyms, lP, lmol, di)
        return lsyms, (lPv - lPv[di]) @ _rot_align(lp, -V).T + V * md, [di]

    def _seat_chelate(smi, dons):
        lmol = Chem.AddHs(Chem.MolFromSmiles(smi))
        tp = [ref[v] / np.linalg.norm(ref[v])
              * MSB.md_distance(metal, lmol.GetAtomWithIdx(int(d)).GetSymbol(),
                                atom=lmol.GetAtomWithIdx(int(d)), mol=lmol)
              for v, d in zip((0, 2, 4)[:len(dons)], dons)]
        rc = _embed_metallacycle(lmol, list(dons), metal, donor_target_pos=tp)
        if rc is None:
            return None
        lsyms, confs = rc
        Q = _orient_chelate_to_vertices(confs[0], list(dons), tp, asym=False)
        return (lsyms, Q, list(dons)) if Q is not None else None

    cases = []
    try:
        cases.append(("pyridin (monodentat)",) + _seat_mono("c1ccncc1"))
    except Exception as exc:
        print("   pyridin: SKIP (%s)" % type(exc).__name__)
    for label, smi, dons in (("ethylendiamin (bidentat)", "NCCN", [0, 3]),
                             ("diethylentriamin (tridentat)", "NCCNCCN", [0, 3, 6])):
        try:
            s = _seat_chelate(smi, dons)
        except Exception as exc:
            s = None
            print("   %s: SKIP (%s)" % (label, type(exc).__name__))
        if s is not None:
            cases.append((label,) + s)

    print(f"{'Ligand':>30} | {'#Don':>4} | {'freie Achse':>11} | {'d(Donor)':>9} | "
          f"{'d(M-D)':>8} | {'d(Biss)':>8} | {'d(intern)':>9}")
    M = np.zeros(3)
    for label, lsyms, Q, dons in cases:
        ax = _free_rigid_axis(Q, dons, M)
        if ax is None:
            print(f"{label:>30} | {len(dons):4d} | {'KEINE (0 DOF)':>11} | "
                  f"{'--':>9} | {'--':>8} | {'--':>8} | {'--':>9}")
            continue
        o, u = ax
        Qr = (np.asarray(Q, float) - o) @ _axis_rot(u, math.radians(_ANG)).T + o
        d_don = max(float(np.linalg.norm(Qr[d] - Q[d])) for d in dons)
        d_md = max(abs(float(np.linalg.norm(Qr[d] - M))
                       - float(np.linalg.norm(Q[d] - M))) for d in dons)
        d_bite = 0.0
        for i, a in enumerate(dons):
            for b in dons[i + 1:]:
                d_bite = max(d_bite, abs(float(np.linalg.norm(Qr[a] - Qr[b]))
                                         - float(np.linalg.norm(Q[a] - Q[b]))))
        D0 = np.linalg.norm(Q[:, None, :] - Q[None, :, :], axis=-1)
        D1 = np.linalg.norm(Qr[:, None, :] - Qr[None, :, :], axis=-1)
        d_int = float(np.max(np.abs(D1 - D0)))
        print(f"{label:>30} | {len(dons):4d} | {'1 (%.0f Grad)' % _ANG:>11} | "
              f"{d_don:9.2e} | {d_md:8.2e} | {d_bite:8.2e} | {d_int:9.2e}")

    # --- 3) DEFAULT OFF ------------------------------------------------------
    _prev = os.environ.pop("DELFIN_FFFREE_LIGAND_DOF_SEAT", None)
    try:
        _off = None
        if cases:
            _l, _s, _Q, _d = cases[0]
            _off = _free_dof_reseat(_Q, _s, _d, np.zeros((1, 3)), [metal], 99)
        print("   Vorgabe AUS -> _free_dof_reseat liefert %s  (%s)"
              % (_off, "OK" if _off is None else "FEHLER"))
    finally:
        if _prev is not None:
            os.environ["DELFIN_FFFREE_LIGAND_DOF_SEAT"] = _prev

    # --- 4) REACH: a REAL colliding ligand pair -------------------------------
    # Two bulky phosphines on ADJACENT octahedron vertices.  Both are seated exactly
    # as assemble_from_config does it (VSEPR + _rot_align), the second
    # sees the first in `_clash_count` -- and afterwards has only the azimuth left.
    try:
        smi = "P(C(C)(C)C)(C(C)(C)C)C(C)(C)C"
        lsyms, lP, lmol = _ligand_3d(smi)
        di = next(a.GetIdx() for a in lmol.GetAtoms() if a.GetSymbol() == "P")
        md = MSB.md_distance(metal, "P", atom=lmol.GetAtomWithIdx(di), mol=lmol)
        lPv, lp = _vsepr_reconstruct(lsyms, lP, lmol, di)
        blocks = []
        for v in (0, 2):
            V = ref[v] / np.linalg.norm(ref[v])
            blocks.append((lPv - lPv[di]) @ _rot_align(lp, -V).T + V * md)
        ex = np.vstack([np.zeros((1, 3)), blocks[0]])
        ex_s = [metal] + list(lsyms)
        c0 = _clash_count(blocks[1], ex, lsyms, ex_s)
        os.environ["DELFIN_FFFREE_LIGAND_DOF_SEAT"] = "1"
        try:
            res = _free_dof_reseat(blocks[1], lsyms, [di], ex, ex_s, c0)
        finally:
            if _prev is None:
                os.environ.pop("DELFIN_FFFREE_LIGAND_DOF_SEAT", None)
            else:
                os.environ["DELFIN_FFFREE_LIGAND_DOF_SEAT"] = _prev
        if c0 <= 0:
            print("   Reichweite: das Paar kollidiert gar nicht (clash 0) -> die Achse"
                  " wird hier nicht gebraucht; das ist ein Befund, kein Erfolg.")
        elif res is None:
            print("   Reichweite: clash %d, aber KEINE Drehung ist besser -> der"
                  " Freiheitsgrad traegt an diesem Paar nicht." % c0)
        else:
            Qn, c1 = res
            _dmd = abs(float(np.linalg.norm(Qn[di])) - float(np.linalg.norm(blocks[1][di])))
            print("   Reichweite: P(tBu)3 / P(tBu)3 cis  clash %d -> %d, "
                  "M-D-Drift %.2e A" % (c0, c1, _dmd))
    except Exception as exc:
        print("   Reichweite: SKIP (%s: %s)" % (type(exc).__name__, exc))


def _self_test_oc6_twist() -> None:
    """Does the corrector really rotate the half twist out -- and in doing so leave
    everything standing that it must leave standing?

    Measured on SYNTHETIC frames whose truth is known, not on a
    pool: an octahedron twisted about the C3 axis by a known angle.
    phi = 0 is the octahedron, phi = 60 degrees the ideal trigonal prism, and the
    measured median CShM 11.03 lies in between -- the half twist this is about.

    Five questions, and each of them can say NO:
      1) does CShM(OC-6) fall -- and to ~0, not just a bit?
      2) does r(M-D) stay exact (the rotation is about the metal)?
      3) does the bite stay exact (the arm rotates RIGIDLY)?
      4) does it leave a REQUESTED TPR-6 alone?  (otherwise it destroys an isomer)
      5) does it leave an ambiguous frame alone (min-trans <= 120 degrees)?
    And the sixth, which is not a question but the default: with the switch OFF
    nothing happens at all.
    """
    from delfin.manta import polyhedra as _PH
    print("\nOC-6 Twist-Korrektor -- der halbe Bailar-Twist, FF-frei zurueckgedreht")
    print(f"  Schalter DELFIN_FFFREE_OC6_TWIST_SEAT gelesen als: "
          f"{_oc6_twist_seat_enabled()}  (Vorgabe muss False sein)")

    def _twisted_oct(phi_deg, rad=2.10):
        """An octahedron twisted about the C3 axis [1,1,1]: the upper triangular face
        by +phi/2, the lower by -phi/2.  phi=0 -> OC-6, phi=60 -> ideal TPR-6."""
        V = _PH.ref_vectors("OC-6 octahedron")
        c3 = np.array([1.0, 1.0, 1.0]) / math.sqrt(3.0)
        top, bot = [], []
        for v in V:                       # the two triangles perpendicular to the C3 axis
            (top if float(np.dot(v, c3)) > 0 else bot).append(v)
        Rt = _axis_rot(c3, math.radians(+phi_deg / 2.0))
        Rb = _axis_rot(c3, math.radians(-phi_deg / 2.0))
        return ([np.asarray(v) @ Rt.T * rad for v in top]
                + [np.asarray(v) @ Rb.T * rad for v in bot])

    def _mk_monodentate_frame(dirs):
        """Metal + six monodentate arms (donor + one backbone atom pointing outwards), so
        that every block really has a BODY that the rotation must carry along."""
        syms = ["Fe"]
        P = [np.zeros(3)]
        blocks = []
        for v in dirs:
            st = len(P)
            syms.append("N")
            P.append(np.asarray(v, float))
            syms.append("C")
            P.append(np.asarray(v, float) * 1.65)          # radially outward
            blocks.append((st, 2, [st]))
        return syms, np.asarray(P, float), blocks

    print(f"\n{'Frame':>26} | {'min-trans':>9} | {'CShM vor':>9} | {'CShM nach':>9} | "
          f"{'dM-D':>8} | Urteil")
    ok_all = True
    for phi in (0.0, 15.0, 30.0, 45.0, 60.0):
        dirs = _twisted_oct(phi)
        syms, P, blocks = _mk_monodentate_frame(dirs)
        u = {b[2][0]: P[b[2][0]] / float(np.linalg.norm(P[b[2][0]])) for b in blocks}
        _pairs, _mt = _oc6_trans_pairs(u, sorted(u))
        before = _PH.cshm([P[b[2][0]] for b in blocks], "OC-6 octahedron")
        Xc = _oc6_twist_seat(syms, P, blocks, "Fe", "OC-6 octahedron")
        if Xc is None:
            verdict = ("uebersprungen (min-trans <= 120)" if _mt <= _OC6_TRANS_MIN
                       else "uebersprungen")
            print(f"{'OC-6 twist %4.1f' % phi:>26} | {_mt:8.1f}d | {before:9.3f} | "
                  f"{'--':>9} | {'--':>8} | {verdict}")
            # phi=0 is ALREADY the octahedron -> CShM cannot fall -> None is correct
            if phi == 0.0 and before < 1e-6:
                continue
            if _mt > _OC6_TRANS_MIN and before > 1.0:
                ok_all = False              # should have acted
            continue
        after = _PH.cshm([Xc[b[2][0]] for b in blocks], "OC-6 octahedron")
        dmd = max(abs(float(np.linalg.norm(Xc[b[2][0]]))
                      - float(np.linalg.norm(P[b[2][0]]))) for b in blocks)
        # the arm must have COME ALONG: the backbone atom keeps its distance to the donor
        darm = max(abs(float(np.linalg.norm(Xc[b[0]] - Xc[b[0] + 1]))
                       - float(np.linalg.norm(P[b[0]] - P[b[0] + 1]))) for b in blocks)
        good = (after < before - 1e-9) and dmd < 1e-6 and darm < 1e-6
        ok_all = ok_all and good
        print(f"{'OC-6 twist %4.1f' % phi:>26} | {_mt:8.1f}d | {before:9.3f} | "
              f"{after:9.3f} | {dmd:8.1e} | {'OK' if good else 'FEHLER'}"
              f"  (Arm {darm:.1e})")

    # 4) a REQUESTED TPR-6 must stay untouched -- otherwise an isomer is lost
    dirs = [np.asarray(v, float) * 2.10 for v in _PH.ref_vectors("TPR-6 trigonal prism")]
    syms, P, blocks = _mk_monodentate_frame(dirs)
    tpr = _oc6_twist_seat(syms, P, blocks, "Fe", "TPR-6 trigonal prism")
    print(f"{'TPR-6 angefordert':>26} | {'--':>9} | {'--':>9} | {'--':>9} | {'--':>8} | "
          f"{'OK (nicht angefasst)' if tpr is None else 'FEHLER: Isomer zerstoert'}")
    ok_all = ok_all and (tpr is None)

    # 5) ambiguous: a frame whose best pairing stays below 120 degrees (all six
    #    donors crowded into ONE hemisphere) -- the corrector must not guess
    amb = []
    for k in range(6):
        a = 2.0 * math.pi * k / 6.0
        v = np.array([math.cos(a) * 0.80, math.sin(a) * 0.80, 0.60])
        amb.append(v / float(np.linalg.norm(v)) * 2.10)
    syms, P, blocks = _mk_monodentate_frame(amb)
    u = {b[2][0]: P[b[2][0]] / float(np.linalg.norm(P[b[2][0]])) for b in blocks}
    _p2, _mt2 = _oc6_trans_pairs(u, sorted(u))
    ambr = _oc6_twist_seat(syms, P, blocks, "Fe", "OC-6 octahedron")
    print(f"{'mehrdeutig (Halbkugel)':>26} | {_mt2:8.1f}d | {'--':>9} | {'--':>9} | "
          f"{'--':>8} | {'OK (uebersprungen)' if ambr is None else 'FEHLER: geraten'}")
    ok_all = ok_all and (ambr is None)

    # 6) a CHELATE: two donors on ONE rigid body.  It may only rotate RIGIDLY
    #    -- the bite is ligand geometry, not polyhedron geometry.
    dirs = _twisted_oct(30.0)
    syms = ["Fe"]; P = [np.zeros(3)]; blocks = []
    for i in range(0, 6, 2):
        st = len(P)
        syms += ["N", "C", "N"]
        P += [np.asarray(dirs[i], float),
              (np.asarray(dirs[i], float) + np.asarray(dirs[i + 1], float)) * 0.62,
              np.asarray(dirs[i + 1], float)]
        blocks.append((st, 3, [st, st + 2]))
    P = np.asarray(P, float)
    bite0 = [float(np.linalg.norm(P[b[2][0]] - P[b[2][1]])) for b in blocks]
    before = _PH.cshm([P[d] for b in blocks for d in b[2]], "OC-6 octahedron")
    Xc = _oc6_twist_seat(syms, P, blocks, "Fe", "OC-6 octahedron")
    if Xc is None:
        print(f"{'3 Chelate, Biss fest':>26} | {'--':>9} | {before:9.3f} | {'--':>9} | "
              f"{'--':>8} | uebersprungen")
        ok_all = False
    else:
        after = _PH.cshm([Xc[d] for b in blocks for d in b[2]], "OC-6 octahedron")
        bite1 = [float(np.linalg.norm(Xc[b[2][0]] - Xc[b[2][1]])) for b in blocks]
        dbite = max(abs(x - y) for x, y in zip(bite0, bite1))
        dmd = max(abs(float(np.linalg.norm(Xc[d])) - float(np.linalg.norm(P[d])))
                  for b in blocks for d in b[2])
        good = (after < before - 1e-9) and dbite < 1e-6 and dmd < 1e-6
        ok_all = ok_all and good
        print(f"{'3 Chelate, Biss fest':>26} | {'--':>9} | {before:9.3f} | {after:9.3f} | "
              f"{dmd:8.1e} | {'OK' if good else 'FEHLER'}  (Biss {dbite:.1e})")
    print(f"\n  Gesamturteil: "
          f"{'ALLE BESTANDEN' if ok_all else 'MINDESTENS EINER FEHLGESCHLAGEN'}")


def _self_test_lp_orient() -> None:
    """How far do the donor lone pairs point PAST the metal?

    That angle IS the flapping.  A conjugated donor's lone pair lies in its pi plane, so a
    correctly seated planar chelate has the metal IN that plane and the angle is zero;
    there is no such thing as a tilted aromatic chelate.  The old orientation was chosen by
    a 36-step sweep maximising the minimum backbone-to-metal distance -- collision
    avoidance where the constraint is bonding -- so this measures what that costs.
    """
    cases = [
        ("c1ccc(-c2ccccn2)nc1",        "2,2'-bipyridine (planar, aromatic)"),
        ("c1cnc2c(c1)ccc1cccnc12",     "1,10-phenanthroline (fused, rigid)"),
        ("NCCN",                       "ethylenediamine (sp3)"),
        ("CC(=O)[O-]",                 "acetate (carboxylate O donors)"),
    ]
    print("\nlone-pair orientation -- degrees the lone pairs miss the metal by")
    print("  'converge' = angle between the two lone pairs measured on the FREE ligand.")
    print("  A chelating conformer aims both at one point, so the two directions must")
    print("  CONVERGE (obtuse, near 180-bite).  If they diverge, no rotation of the whole")
    print("  ligand can aim them at a metal -- the conformer itself is the wrong one, and")
    print("  the seat can only flap it into place.")
    print(f"{'ligand':>36} | {'conv':>5} | {'--':>6} {'lp':>6} {'bite':>6} "
          f"{'both':>6} | {'best':>7} | OK")
    for smi, label in cases:
        try:
            lsyms, lP, lmol = _ligand_3d(smi)
        except Exception as exc:
            print(f"{label:>36} | SKIP ({type(exc).__name__})")
            continue
        dons = [a.GetIdx() for a in lmol.GetAtoms() if a.GetSymbol() in ("N", "O")][:2]
        if len(dons) != 2:
            print(f"{label:>36} | SKIP (no donor pair)")
            continue
        d1, d2 = dons
        _l1, _l2 = _lone_pair_dir(lP, lmol, d1), _lone_pair_dir(lP, lmol, d2)
        conv = (math.degrees(math.acos(max(-1.0, min(1.0, float(np.dot(_l1, _l2))))))
                if (_l1 is not None and _l2 is not None) else float("nan"))
        # Four configurations, because the sweep was only ONE candidate root.  A rigid
        # ligand like phenanthroline cannot have a wrong conformer, so if its lone pairs
        # still miss, the cause is upstream: the donors are pinned at an angle and a length
        # the ligand does not aim at.  That is precisely what the derived bite removes, so
        # it has to be in the comparison.
        miss = []
        for flag, bite in (("0", "0"), ("1", "0"), ("0", "1"), ("1", "1")):
            prev = os.environ.get("DELFIN_FFREE_LP_ORIENT")
            prevb = os.environ.get("DELFIN_FFREE_BITE_LAW")
            os.environ["DELFIN_FFREE_LP_ORIENT"] = flag
            os.environ["DELFIN_FFREE_BITE_LAW"] = bite
            try:
                syms, P = assemble_multichelate("Ni", "OC-6 octahedron",
                                                [(smi, [d1, d2], [0, 2])])
                worst = 0.0
                for d in (d1, d2):
                    lp = _lone_pair_dir(P[1:], lmol, d)
                    if lp is None:
                        continue
                    t = P[0] - P[1 + d]
                    nt = float(np.linalg.norm(t))
                    if nt < 1.0e-9:
                        continue
                    c = max(-1.0, min(1.0, float(np.dot(lp, t / nt))))
                    worst = max(worst, math.degrees(math.acos(c)))
                miss.append(worst)
            except Exception as exc:
                miss.append(float("nan"))
                print(f"    ({label}: {type(exc).__name__}: {exc})")
            finally:
                for _k, _v in (("DELFIN_FFREE_LP_ORIENT", prev),
                               ("DELFIN_FFREE_BITE_LAW", prevb)):
                    if _v is None:
                        os.environ.pop(_k, None)
                    else:
                        os.environ[_k] = _v
        best = min(m for m in miss if m == m)
        ok = best <= miss[0] + 1.0e-6             # never-worse is the criterion
        print(f"{label:>36} | {conv:4.0f}d | {miss[0]:6.1f}d {miss[1]:6.1f}d "
              f"{miss[2]:6.1f}d {miss[3]:6.1f}d | {miss[0]-best:+6.1f}d | {ok}")




def _self_test_trilateration() -> None:
    """Does trilateration really hold BOTH sides, where the two old branches each drop one?

    The seating documents the trade-off as unavoidable: per-donor radial keeps the M-D
    lengths and splays the ring, rigid keeps the ring and lets M-D drift.  So the test
    reports exactly those two errors for all three placements.  A lever that merely moves
    the error from one column to the other would be a rename, not a fix.
    """
    from delfin.manta import polyhedra as _PH
    cases = [
        ("tridentate, 2.1 A radii", "OC-6 octahedron", [0, 2, 4], 2.10),
        ("tridentate, 2.4 A radii", "OC-6 octahedron", [0, 2, 4], 2.40),
        ("tetradentate",            "OC-6 octahedron", [0, 2, 4, 1], 2.10),
    ]
    print("\ntrilateration -- both constraints at once?  (errors in Angstrom)")
    print(f"{'case':>24} | {'placement':>14} | {'M-D err':>8} | {'D-D err':>8}")
    for label, geom, verts, rad in cases:
        ref = _PH.ref_vectors(geom)
        # a synthetic rigid donor set: an equilateral-ish arm spacing the ideal vertices
        # cannot satisfy, which is the situation a real polydentate is in
        L = np.array([[0.0, 0.0, 0.0], [2.55, 0.0, 0.0], [1.27, 2.30, 0.0],
                      [1.27, 0.77, 2.10]])[:len(verts)]
        T = [ref[v] / np.linalg.norm(ref[v]) * rad for v in verts]
        D = {(i, j): float(np.linalg.norm(L[i] - L[j]))
             for i in range(len(verts)) for j in range(i + 1, len(verts))}

        def _err(X):
            md = max(abs(float(np.linalg.norm(X[i])) - rad) for i in range(len(verts)))
            dd = max(abs(float(np.linalg.norm(X[i] - X[j])) - D[(i, j)])
                     for (i, j) in D)
            return md, dd

        # (a) ideal vertices, per-donor radial: radii exact, ligand distances ignored
        a_md, a_dd = _err(np.array(T))
        # (b) rigid body onto the ideal vertices: ligand exact, radii drift
        dmu, tmu = L.mean(0), np.array(T).mean(0)
        Rk = _kabsch_rot(L - dmu, np.array(T) - tmu)
        Xb = (L - dmu) @ Rk.T + tmu
        b_md, b_dd = _err(Xb)
        # (c) trilateration, then the same rigid fit onto the moved targets
        tri = _trilaterate_donor_targets(L, list(range(len(verts))), T)
        if tri is None:
            print(f"{label:>24} | trilateration returned None")
            continue
        tgt = np.array(tri)
        Rk2 = _kabsch_rot(L - dmu, tgt - tgt.mean(0))
        Xc = (L - dmu) @ Rk2.T + tgt.mean(0)
        c_md, c_dd = _err(Xc)
        for nm, (md, dd) in (("radial(ideal)", (a_md, a_dd)),
                             ("rigid(ideal)", (b_md, b_dd)),
                             ("TRILATERATED", (c_md, c_dd))):
            print(f"{label:>24} | {nm:>14} | {md:8.3f} | {dd:8.3f}")


_REAL_FRAME_CASES = [
    ("cisplatin CN4",        "N[Pt](N)(Cl)Cl"),
    ("CoCl3(NH3)3 CN6",      "[NH3][Co]([NH3])([NH3])([Cl])([Cl])[Cl]"),
    ("Fe(en)3 CN6, 3 Chel",  "[Fe]123([N]CC[N]1)([N]CC[N]2)[N]CC[N]3"),
    ("CoCl(en)2 CN6, 2 Chel", "[Cl][Co]12([NH3])([N]CC[N]1)[N]CC[N]2"),
    ("Pd(PtBu3)Cl2",         "[Cl][Pd]([Cl])[P](C(C)(C)C)(C(C)(C)C)C(C)(C)C"),
    ("Ni(Cl)2(NMe3) CN4",    "[Ni](Cl)(Cl)N(C)(C)C"),
    ("CoCl2(en)2 CN6",       "[NH2]CC[NH2][Co]1([NH2]CC[NH2]1)([Cl])[Cl]"),
    ("Fe-Citrat, polydentat",
     "O=C1C(CC(O)=O)(CC(O)=O)O[Fe]23(OC(C(CC(O)=O)(CC(O)=O)O2)=O)"
     "(OC(CC(O)=O)(CC(O)=O)C(O3)=O)O1"),
    ("Pt(en)(NH2CH2CH2NH2)", "NCCN[Pt]1NCCN1"),
]


def _real_frame_hashes() -> None:
    """REAL frames, not synthetic ones: builds every system of the set via the
    normal converter path and prints per system a SHA1 over ALL emitted
    XYZ strings.  Two uses:

      * BYTE IDENTITY.  The same call against the unchanged state (switch
        not set) must deliver character for character the same hashes.
      * THE FREEZE CENSUS.  Per system: how many atoms does the builder nail down
        (metal + donors, see `fixed` in assemble_from_config) and how many
        stay movable?  That is the number the DOF question is about.
    """
    import hashlib
    from delfin.manta import converter_backend as CB
    from delfin.manta import decompose as DEC
    print("%-24s | %4s | %5s | %-8s | %s"
          % ("System", "#Iso", "CN", "fix/Atome", "sha1(alle XYZ)"))
    for label, smi in _REAL_FRAME_CASES:
        try:
            r = CB._fffree_isomers(smi)
        except Exception as exc:
            print("%-24s | FEHLER %s: %s" % (label, type(exc).__name__, exc))
            continue
        if not r:
            print("%-24s | %4s | %5s | %-8s | %s" % (label, "-", "-", "-", "None (legacy)"))
            continue
        h = hashlib.sha1()
        for xyz, lab in r:
            h.update(lab.encode())
            h.update(xyz.encode())
        try:
            cn = int(DEC.decompose(smi).get("cn"))
        except Exception:
            cn = -1
        nat = 0
        for _ln in r[0][0].splitlines():
            _f = _ln.split()
            if len(_f) >= 4 and _f[0][:1].isalpha():
                try:
                    float(_f[1]); float(_f[2]); float(_f[3])
                except ValueError:
                    continue
                nat += 1
        print("%-24s | %4d | %5d | %8s | %s"
              % (label, len(r), cn, "%d/%d" % (cn + 1, nat), h.hexdigest()))


if __name__ == "__main__":       # pragma: no cover
    import sys as _sys
    if len(_sys.argv) > 1 and _sys.argv[1] == "realframes":
        _real_frame_hashes()
        raise SystemExit(0)
    _run_self_tests()
