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

from delfin.manta.assemble_hapto import (
    _ETA_NFOLD,
    _axis_rot_angles,
    _cp_pucker_amps,
    _dedup_builds,
    _eta_centroid_distance,
    _eta_is_aromatic,
    _hapto_base_variants,
    _hapto_fold_rings,
    _piano_leg_target_deg,
    _piano_leg_tilt_enabled,
    _place_eta_ring,
    _ring_spins,
    _rmsd_aligned,
    _rot_about_axis,
    _slip_modes,
    _tilt_piano_legs,
    assemble_hapto,
    assemble_hapto_axis_rotants,
    assemble_hapto_ensemble,
)  # noqa: F401

from delfin.manta.assemble_seat import (
    _GD_BOUND_STEPS,
    _GD_CLASH_F,
    _GD_EPS,
    _GD_FREE_STEPS,
    _GD_H_W,
    _GD_PASSES,
    _GD_RESID,
    _LIGAND_DOF_AXIS_TOL,
    _OC6_AXIS_MIN,
    _OC6_CSHM_KEEP,
    _OC6_RIGID_TOL,
    _OC6_SEAT_CENSUS,
    _OC6_SEAT_CSHM,
    _OC6_TRANS_MIN,
    _free_dof_reseat,
    _free_rigid_axis,
    _gd_loss,
    _gd_move_ok,
    _gd_resid,
    _global_donor_seat,
    _global_donor_seat_enabled,
    _ligand_dof_seat_enabled,
    _ligand_dof_seat_steps,
    _oc6_ideal_axes,
    _oc6_trans_pairs,
    _oc6_twist_seat,
    _oc6_twist_seat_enabled,
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
