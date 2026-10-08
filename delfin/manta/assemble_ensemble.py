"""Heteroleptic assembly from mols, the heteroleptic ensemble, the constrained UFF relax and build_and_relax of the FF-free constructor.

Moved verbatim from delfin/manta/assemble_complex.py (MANTA split, 2026-10); every
statement is the original text, only the imports below are new.
"""
from __future__ import annotations

import numpy as np
import os
from rdkit import Chem
from rdkit.Chem import AllChem

from delfin.manta import metal_sphere_builder as MSB
from delfin.manta.assemble_donor_plane import (
    _COMBO_MATERIALISE_MAX,
    _ranked_combos,
)
from delfin.manta.assemble_fold_fp import (
    _complex_rmsd,
    _fold_fp,
    _fold_fp_enabled,
    _fold_rings_from_blocks,
    _fold_same,
)
from delfin.manta.assemble_ligand_confs import (
    _clash_count,
    _joint_declash_frame,
    _ligand_3d_from_mol,
    _ligand_confs_from_mol,
    _refine_guarded,
    _sphere_flex_frame,
    _torsion_relax_frame,
)
from delfin.manta.assemble_ligand_embed import (
    _rot_align,
)
from delfin.manta.assemble_orient import (
    _axis_rot,
    _diatomic_donor_partner,
    _diatomic_orient_enabled,
    _donor_and_lp,
    _orient_diatomic_block,
    _vsepr_reconstruct,
)
from delfin.manta.assemble_seat import (
    _free_dof_reseat,
)


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
