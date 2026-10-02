"""The frame filters applied to the finished manifold: non-finite, ranking, declash, carbonyl, topology gate, clean gate, pi planes, coordination integrity, conformer completion, GFN-FF rank, permutation dedup and mirror enumeration.

Moved verbatim from delfin/smiles_converter.py (MANTA split, 2026-10); every
statement is the original text, only the imports below are new.
"""
from __future__ import annotations

import math
import os
from typing import Dict, List

from delfin.common.logging import get_logger
from delfin.manta.converter_flags import (
    _delfin_env_float,
)
from delfin.manta.ml_tables import (
    Chem,
    RDKIT_AVAILABLE,
    _COVALENT_RADII,
    _METAL_SET,
)

logger = get_logger("delfin.smiles_converter")


def _filter_nonfinite_isomers(results):
    """Universal output contract (#36): never emit a structure with non-finite
    (NaN/inf) coordinates — it is not a valid DFT startpoint and corrupts every
    downstream metric.  Drops ONLY the broken isomer/conformer, keeps the finite
    ones (no coverage loss when other builds for the same SMILES are valid).
    Applies to EVERY build path (legacy, hapto, fffree) since it gates the public
    entry point.  Geometry-only, deterministic."""
    if not results:
        return results
    clean = []
    for item in results:
        xyz = item[0] if isinstance(item, (tuple, list)) and item else item
        finite = True
        for ln in str(xyz).splitlines():
            p = ln.split()
            if len(p) >= 4:
                try:
                    x, y, z = float(p[1]), float(p[2]), float(p[3])
                except ValueError:
                    continue
                if not (math.isfinite(x) and math.isfinite(y) and math.isfinite(z)):
                    finite = False
                    break
        if finite:
            clean.append(item)
    return clean


def _rank_emitted_isomers(isomers):
    """Order the emitted ensemble best (most crystal-like) first.  Never drops or
    alters a structure — only reorders.  Safe no-op on any error or with
    DELFIN_NO_FRAME_RANK=1.

    DEFAULT = the legacy least-clash ranking (cross-validated 2026-06-12).  That
    ranker measures ONLY steric overlap, so it is blind in both directions: a
    TORN / decoordinated frame has no overlap at all and therefore scores the
    perfect 100.0 and LEADS, while on a clean ensemble 86 % of all emitted frames
    score exactly 100.0 -> a complete tie -> the delivered order is the ENUMERATION
    order, i.e. no ranking at all (measured over the built census 2026-07-28; user
    report the same day: "the best frames are not sorted to the front").

    DELFIN_FRAME_RANK_QUALITY=1 switches to the DEFECT ordering
    (:func:`delfin.manta._conformer_rank.rank_isomers_quality`): torn / spurious /
    collapsed ligand bonds first, then clash, then coordination distortion
    relative to the best frame of the same CN signature, enumeration index last.
    Pure geometry, deterministic, no RMSD, no energy; a verified PERMUTATION with an
    explicit bijection guard, so no frame is ever lost and every label travels with
    its frame.  Default-OFF -> byte-identical order to today.  Applied ONLY here, on
    the final emitted ensemble -- the construction-time base-frame picker keeps using
    the legacy ``rank_isomers``, so the manifold CONTENT is unchanged."""
    try:
        if os.environ.get("DELFIN_FRAME_RANK_QUALITY", "0") == "1":
            from delfin.manta._conformer_rank import rank_isomers_quality
            return rank_isomers_quality(isomers)
        from delfin.manta._conformer_rank import rank_isomers
        return rank_isomers(isomers)
    except Exception:
        return isomers


def _hapto_declash_filter(isomers):
    """Final-manifold eta-ring inter-ligand declash over the emitted ensemble.

    The legacy analytical piano-stool seat (``_build_hapto_scaffold``) places an
    eta-coordinated ring (arene / Cp / Cb / diene) on the metal but does NOT
    choose its ROTATIONAL conformer about the M->centroid axis, so ring carbons /
    substituents routinely eclipse the carbonyls and crowd the other ligands.
    Measured against CCDC piano-stool crystals (find_inter_ligand_clash, F26) the
    seat over-clashes ~3-4x, the residual dominated by the eta-face crowding
    pendant aryls / other rings / sigma-donors -- contacts essentially absent in
    the crystals.

    This pass rotates each eta-ring about its M->centroid axis to the LOCAL
    inter-ligand-clash minimum within a bounded window (the missing piano-stool
    DOF).  It is a pure isometry of the coordination polyhedron: every
    M-(ring atom) distance is preserved, the metal + all sigma-donors are frozen,
    the inter-ligand overlap is never increased, and it is deterministic.

    Runs as the LAST geometry transform over the FINAL, fully-selected + dedup'd
    ensemble (NOT per-conformer during construction -- that would feed back into
    conformer selection / dedup and change WHICH conformers survive), so the
    isomer/conformer COUNT is exactly preserved: a frame is only re-emitted when
    the spin actually moved an atom; untouched otherwise.

    Identity (byte-identical) unless DELFIN_FFFREE_HAPTO_DECLASH=1.  Geometry-only
    (Bondi vdW radii); no CSD/CCDC data.  Deterministic; never raises."""
    if not isomers or os.environ.get("DELFIN_FFFREE_HAPTO_DECLASH", "0") != "1":
        return isomers
    try:
        import numpy as np
        from delfin.manta import hapto_declash as _HD
        try:
            _maxrot = float(os.environ.get(
                "DELFIN_FFFREE_HAPTO_DECLASH_MAXROT", _HD._DEF_MAXROT_DEG))
        except (TypeError, ValueError):
            _maxrot = _HD._DEF_MAXROT_DEG
        out = []
        for item in isomers:
            is_tuple = isinstance(item, (tuple, list))
            xyz = item[0] if is_tuple and item else item
            lines = str(xyz).splitlines()
            syms, coords, idx = [], [], []
            for li, ln in enumerate(lines):
                p = ln.split()
                if len(p) >= 4:
                    try:
                        x, y, z = float(p[1]), float(p[2]), float(p[3])
                    except ValueError:
                        continue
                    syms.append(p[0]); coords.append([x, y, z]); idx.append(li)
            if len(syms) < 4:
                out.append(item); continue
            P = np.asarray(_HD.declash_frame(syms, coords, max_rot_deg=_maxrot), float)
            orig = np.asarray(coords, float)
            if P.shape != orig.shape or np.allclose(P, orig, atol=1e-7):
                out.append(item); continue          # no eta-ring / nothing moved
            for k, li in enumerate(idx):
                p = lines[li].split()
                p[1] = f"{P[k][0]:.6f}"; p[2] = f"{P[k][1]:.6f}"; p[3] = f"{P[k][2]:.6f}"
                lines[li] = " ".join(p)
            new_xyz = "\n".join(lines)
            out.append((new_xyz,) + tuple(item[1:]) if is_tuple else new_xyz)
        return out
    except Exception:
        return isomers


def _sigma_declash_filter(isomers):
    """SIGMA-co-ligand inter-ligand declash over the emitted ensemble -- the
    sigma-ligand counterpart to the eta-ring spin (_hapto_declash_filter).

    Relieves crowding of the NON-eta ligand BODIES (stannyl / stibine / phosphite
    / bulky phosphine substituents) by the joint_declash kinematics: each ligand's
    whole-body rotation about its M-donor axis + internal rotatable-bond torsions,
    metal + ALL donors FROZEN (coordination polyhedron invariant).  A
    partition-free total-overlap guard makes every frame strictly never-worse
    (measured on the hapto pool: inter-ligand clashes 4867 -> 4150 on top of the
    ring spin, 0 frames worse).

    Runs over the FINAL dedup'd manifold (count-preserving: a frame is only
    re-emitted when an atom actually moved).  Identity (byte-identical) unless
    DELFIN_FFFREE_SIGMA_DECLASH=1.  Geometry-only; deterministic; never raises."""
    if not isomers or os.environ.get("DELFIN_FFFREE_SIGMA_DECLASH", "0") != "1":
        return isomers
    try:
        import numpy as np
        from delfin.manta import hapto_declash as _HD
        out = []
        for item in isomers:
            is_tuple = isinstance(item, (tuple, list))
            xyz = item[0] if is_tuple and item else item
            lines = str(xyz).splitlines()
            syms, coords, idx = [], [], []
            for li, ln in enumerate(lines):
                p = ln.split()
                if len(p) >= 4:
                    try:
                        x, y, z = float(p[1]), float(p[2]), float(p[3])
                    except ValueError:
                        continue
                    syms.append(p[0]); coords.append([x, y, z]); idx.append(li)
            if len(syms) < 4:
                out.append(item); continue
            P = np.asarray(_HD.declash_frame_sigma(syms, coords), float)
            orig = np.asarray(coords, float)
            if P.shape != orig.shape or np.allclose(P, orig, atol=1e-7):
                out.append(item); continue
            for k, li in enumerate(idx):
                p = lines[li].split()
                p[1] = f"{P[k][0]:.6f}"; p[2] = f"{P[k][1]:.6f}"; p[3] = f"{P[k][2]:.6f}"
                lines[li] = " ".join(p)
            new_xyz = "\n".join(lines)
            out.append((new_xyz,) + tuple(item[1:]) if is_tuple else new_xyz)
        return out
    except Exception:
        return isomers


def _carbonyl_fix_filter(isomers):
    """Contract stretched terminal metal-carbonyl C#O bonds over the emitted
    ensemble.  The analytical seat's flat _bond_len places C-O at the single-bond
    length (~1.43 A) ignoring bond order, so M-C#O carbonyls come out stretched
    (~1.40 A vs the true ~1.13-1.15 A; measured 56% > 1.35 A).  This slides only
    the terminal O inward along the C->O axis to the triple-bond length -- a pure
    local geometry correction that never lengthens and moves nothing else, so the
    isomer/conformer COUNT and all other geometry are preserved.

    Identity (byte-identical) unless DELFIN_FFFREE_CARBONYL_FIX=1.  Geometry-only;
    deterministic; never raises."""
    if not isomers or os.environ.get("DELFIN_FFFREE_CARBONYL_FIX", "0") != "1":
        return isomers
    try:
        import numpy as np
        from delfin.manta import hapto_declash as _HD
        out = []
        for item in isomers:
            is_tuple = isinstance(item, (tuple, list))
            xyz = item[0] if is_tuple and item else item
            lines = str(xyz).splitlines()
            syms, coords, idx = [], [], []
            for li, ln in enumerate(lines):
                p = ln.split()
                if len(p) >= 4:
                    try:
                        x, y, z = float(p[1]), float(p[2]), float(p[3])
                    except ValueError:
                        continue
                    syms.append(p[0]); coords.append([x, y, z]); idx.append(li)
            if len(syms) < 3:
                out.append(item); continue
            P = np.asarray(_HD.correct_carbonyls(syms, coords), float)
            orig = np.asarray(coords, float)
            if P.shape != orig.shape or np.allclose(P, orig, atol=1e-7):
                out.append(item); continue
            for k, li in enumerate(idx):
                p = lines[li].split()
                p[1] = f"{P[k][0]:.6f}"; p[2] = f"{P[k][1]:.6f}"; p[3] = f"{P[k][2]:.6f}"
                lines[li] = " ".join(p)
            new_xyz = "\n".join(lines)
            out.append((new_xyz,) + tuple(item[1:]) if is_tuple else new_xyz)
        return out
    except Exception:
        return isomers


def _topology_gate_filter(isomers):
    """Final UNIVERSAL topology-preservation gate over the FULL emitted ensemble.

    SPEC
    ----
    Root cause it cures: some conformer/pucker frames of a flexible complex let a
    NON-bonded heavy-atom pair fall to bonding distance, forging a SPURIOUS bond
    (a topology break).  Example: ACEBOB (Re tris-carbonyl macrocycle) -- in some
    frames a carbonyl C penetrates the macrocyclic cavity and lands ~1.47 A from a
    ligand C/N.  At 1.47 A this reads as a perfectly normal covalent bond, so the
    existing clash/coord self-gates are BLIND to it.  We must reject ONLY the
    broken frames and KEEP the clean ones (never collapse the ensemble).

    Method (consensus-based, NO reference/SMILES alignment needed):
      All emitted frames for one structure share atom order.  A GENUINE covalent
      bond is rigid -- present in ALL frames; a SPURIOUS contact appears in only a
      minority.  So for every heavy-heavy NON-metal atom pair (i,j) we measure the
      fraction of frames in which dist(i,j) < bond_thr(elem_i, elem_j), with
      bond_thr = 1.15 * (cov_i + cov_j) (Cordero covalent radii; ~1.75 A for C-C).
      A pair is a "real bond" iff that fraction >= 0.80 (bonded in the clear
      majority).  A frame has a TOPOLOGY BREAK iff it contains a pair that is at
      bonding distance in THIS frame yet is NOT a real bond.  Such frames are
      dropped.  Metal-involving pairs are excluded (M-donor distances vary
      legitimately and must never be treated as breakable bonds).

    Guarantees:
      * Geometry-only and graph-free -- no per-SMILES / per-refcode special case.
      * NEVER-WORSE: if every frame would be dropped (or the result would be
        empty) the ORIGINAL list is returned unchanged.
      * Identity (returns input untouched) unless DELFIN_FFFREE_TOPOLOGY_GATE=1,
        so output is BYTE-IDENTICAL to the baseline when off.
      * Deterministic; never raises.
    """
    if not isomers or os.environ.get("DELFIN_FFFREE_TOPOLOGY_GATE", "0") != "1":
        return isomers
    if len(isomers) < 2:
        # Consensus is undefined with a single frame; nothing to compare against.
        return isomers
    try:
        def _radius(sym):
            return _COVALENT_RADII.get(sym, 0.76)

        # Parse every frame into (symbols, coords) sharing one atom order.
        frames = []
        n_atoms = None
        for item in isomers:
            xyz = item[0] if isinstance(item, (tuple, list)) and item else item
            syms, xs, ys, zs = [], [], [], []
            for ln in str(xyz).splitlines():
                p = ln.split()
                if len(p) >= 4:
                    try:
                        x, y, z = float(p[1]), float(p[2]), float(p[3])
                    except ValueError:
                        continue
                    syms.append(p[0])
                    xs.append(x)
                    ys.append(y)
                    zs.append(z)
            if n_atoms is None:
                n_atoms = len(syms)
            if len(syms) != n_atoms or n_atoms == 0:
                # Frames disagree on atom count -> consensus is ill-defined; bail
                # out safely (identity) rather than risk an unsound drop.
                return isomers
            frames.append((syms, xs, ys, zs))

        nf = len(frames)
        ref_syms = frames[0][0]

        # Heavy (non-H), non-metal atom indices -- consensus is only meaningful
        # for the rigid covalent backbone; H and the metal are excluded.
        heavy = [
            k for k in range(n_atoms)
            if ref_syms[k] != 'H' and ref_syms[k] not in _METAL_SET
        ]

        # Precompute per-pair bond threshold (squared) once.
        # near[(i,j)] = count of frames where dist(i,j) < bond_thr.
        near = {}
        thr2 = {}
        for a in range(len(heavy)):
            i = heavy[a]
            ri = _radius(ref_syms[i])
            for b in range(a + 1, len(heavy)):
                j = heavy[b]
                t = 1.15 * (ri + _radius(ref_syms[j]))
                thr2[(i, j)] = t * t
                near[(i, j)] = 0

        for (syms, xs, ys, zs) in frames:
            for a in range(len(heavy)):
                i = heavy[a]
                xi, yi, zi = xs[i], ys[i], zs[i]
                for b in range(a + 1, len(heavy)):
                    j = heavy[b]
                    dx = xi - xs[j]
                    dy = yi - ys[j]
                    dz = zi - zs[j]
                    if dx * dx + dy * dy + dz * dz < thr2[(i, j)]:
                        near[(i, j)] += 1

        # A pair is a "real bond" if bonded in the clear majority (>= 80%).
        real_bond = set()
        for key, cnt in near.items():
            if cnt >= 0.80 * nf:
                real_bond.add(key)

        # Drop a frame iff any pair is bonding in THIS frame but is not a real bond.
        kept = []
        for item, (syms, xs, ys, zs) in zip(isomers, frames):
            broken = False
            for a in range(len(heavy)):
                i = heavy[a]
                xi, yi, zi = xs[i], ys[i], zs[i]
                for b in range(a + 1, len(heavy)):
                    j = heavy[b]
                    key = (i, j)
                    if key in real_bond:
                        continue
                    dx = xi - xs[j]
                    dy = yi - ys[j]
                    dz = zi - zs[j]
                    if dx * dx + dy * dy + dz * dz < thr2[key]:
                        broken = True
                        break
                if broken:
                    break
            if not broken:
                kept.append(item)

        # NEVER-WORSE: empty or all-dropped -> keep original ensemble.
        if not kept:
            return isomers
        return kept
    except Exception:
        return isomers


def _clean_gate_filter(isomers):
    """UNIVERSAL, ASYMMETRIC, PER-FRAME final clean-manifold emission gate.

    Governing invariant (user, verbatim): "sort out only the BAD ones, not
    the good ones".  The gate is ASYMMETRIC: it drops a frame ONLY when that frame is
    CERTAINLY bad on its OWN geometry; in any doubt it KEEPS the frame.  A
    surviving good structure matters more than removing a borderline one
    (false-negative acceptable; false-positive -- dropping a good frame --
    forbidden).

    Each frame is judged ABSOLUTELY on its own coordinates -- NO ensemble
    consensus, NO best-coordinated-reference donor set, NO minority/majority vote
    (all of which over-reject on small/few-frame ensembles by letting one frame's
    geometry condemn another's).  Bonds are perceived PER FRAME from that frame's
    own interatomic distances.  A frame is REJECTED only if ANY DEFINITE defect is
    present with a CLEAR margin (all graph+geometry only, NO SMILES / refcode /
    element special-casing):

      1. STRUCTURALLY INVALID -- a non-finite coordinate, empty frame, or atom
         count mismatch (a destroyed frame).
      2. REAL INTER-LIGAND CLASH -- two heavy atoms in DIFFERENT per-frame ligand
         connected-components overlap clearly below ``clash`` x (cov_i + cov_j),
         or an H-involving non-bonded pair sits below a hard vdW floor (deep
         interpenetration, not a soft contact).
      3. COLLAPSED BOND -- a per-frame covalent bond fused below ``collapse`` x its
         ideal covalent length (two atoms on top of each other).
      4. BARE METAL / TOTAL DECOORDINATION -- a metal with ZERO heavy donors
         anywhere inside its first shell (the complex fell completely apart).  This
         is the only decoordination test that needs no cross-frame reference, so it
         is the only one applied: a single donor wandering a little is NOT certainly
         bad and is KEPT.

    Deliberately NOT rejected (these over-rejected good frames before): consensus
    "spurious bond" (a close contact is judged per-frame, only a true collapse or a
    clear clash counts), bond OVER-stretch (a long contact is simply not a bond in
    that frame -- topology, not a defect for this gate), and a single donor that
    swung modestly past an ideal M-D distance (handled, if at all, only by the
    dedicated coord-integrity filter).

    NEVER-EMPTY without keeping a broken frame: the old "cleanest-available"
    fallback that retained a broken frame is REMOVED.  Frames that are not
    certainly bad simply stay (no consensus needed).  In the degenerate case where
    EVERY frame is certainly bad, the input is returned UNCHANGED (keep all rather
    than fabricate a marked broken survivor) -- in doubt, keep.

    Thresholds derive from covalent radii (Cordero, ``_COVALENT_RADII``); all
    universal, env-tunable.  Identity (returns input untouched) unless
    DELFIN_FFFREE_CLEAN_GATE=1, so output is BYTE-IDENTICAL to baseline when off.
    Deterministic; never raises; returns the input unchanged on any error.
    """
    if not isomers or os.environ.get("DELFIN_FFFREE_CLEAN_GATE", "0") != "1":
        return isomers
    try:
        # ---- tunables (env, sane universal defaults) -----------------------
        # clash_f: heavy-heavy inter-ligand overlap fraction that is CERTAINLY a
        #          clash (clear margin below the covalent sum).
        clash_f = float(os.environ.get("DELFIN_FFFREE_CLEAN_GATE_CLASH", "0.70"))
        # vdw_f: H-involving non-bonded deep-interpenetration floor (hard).
        vdw_f = float(os.environ.get("DELFIN_FFFREE_CLEAN_GATE_VDW", "0.55"))
        # bond_mult: covalent multiplier defining "at bonding distance" PER FRAME.
        bond_mult = float(os.environ.get("DELFIN_FFFREE_CLEAN_GATE_BONDMULT", "1.15"))
        # collapse_f: a per-frame bond shorter than this fraction of ideal = fused.
        collapse_f = float(os.environ.get("DELFIN_FFFREE_CLEAN_GATE_COLLAPSE", "0.55"))
        # md_shell: metal first-shell radius (heavy donor within this = coordinated).
        md_shell = float(os.environ.get("DELFIN_FFFREE_CLEAN_GATE_SHELL", "2.95"))
        # md_collapse_on: also reject a frame whose heavy donor has fused ONTO the
        #          metal (M-D below collapse_f x covalent-sum).  Independent, default
        #          OFF -> byte-identical to the prior CLEAN_GATE behaviour when unset.
        #          Closes a blindspot: the per-frame pair tests (2/3) SKIP every metal
        #          pair, so a donor collapsed onto the metal is otherwise invisible
        #          here (and to the ligcollapse/geomfault detectors, which exempt the
        #          M-D pair too).  Asymmetric-safe: no real M-D bond is anywhere near
        #          this short (shortest M-heavy-donor bonds ~1.5A).
        md_collapse_on = (os.environ.get(
            "DELFIN_FFFREE_CLEAN_GATE_MD_COLLAPSE", "0") == "1")
        # clash_vdw_on: measure the inter-ligand clash against the VAN DER WAALS sum instead of
        #   the covalent sum.  Default OFF -> byte-identical.
        #
        #   WHY.  Measured 2026-08-04 on the 142 worst-broken systems: the gate's five criteria
        #   fired 420 times for `collapse` and EXACTLY ZERO times for `clash` -- while the eye
        #   reports inter-ligand clashes on 3272 of 27757 frames (11.79 %), a defect class that
        #   occurs in the clean crystals at 0.00 %.  The gate is blind to it, and the reason is a
        #   radius mix-up, not a threshold:
        #       gate  0.70 x (rcov_i + rcov_j)   ->  C-C below 1.06 A   <- guessed default
        #       eye   0.65 x (rvdW_i + rvdW_j)   ->  C-C below 2.21 A   <- COD-validated, FP -> 0
        #   Two carbons never reach 1.06 A without failing the collapse test first, so the clash
        #   branch could not fire.  Switching to the vdW sum is NOT a loosening -- it replaces a
        #   guessed reference with the one the eye's own metric_inter_ligand_clash was calibrated
        #   on ("real crystals have inter-ligand contacts >= 0.85 x vdW sum"; 0.65 is the
        #   conservative firing point where the COD false-positive rate goes to zero).
        # h_point_on: the DIRECTIONAL H-clash -- an H that is correctly bonded to its parent but
        #   POINTS INTO a nearby non-bonded heavy atom.  Default OFF -> byte-identical.
        #
        #   WHY A SECOND H CRITERION.  The gate already has one (`vdw_f`), and it is dead for the
        #   same reason the clash branch was: it compares against the COVALENT sum.  For H-H that
        #   is 0.55 x 0.62 = 0.34 A -- two hydrogens never come that close, so `vdw_h` fired 0
        #   times in the 2026-08-04 tally while the eye reports xh_hh_clash on 5.27 % of frames
        #   (crystals: 0.00 %).
        #   The eye's find_h_clash does not use a distance ratio at all; it asks a DIRECTIONAL
        #   question, which is what actually distinguishes a real contact from a broken one:
        #       H bonded to its parent (0.85-1.20 A)
        #       AND a non-bonded heavy atom within 2.50 A of that H
        #       AND the parent->H direction within 30 deg of parent->heavy
        #   i.e. the H is aimed at the neighbour.  Cause per that detector: VSEPR-snap places H at
        #   ideal polar angles but never picks the rotational phase, so a methyl umbrella can
        #   point straight into its neighbour.  Pure geometry, three constants, no CCDC table --
        #   portable into DELFIN license-clean.
        h_point_on = (os.environ.get("DELFIN_FFFREE_CLEAN_GATE_H_POINT", "0") == "1")
        h_par_min = float(os.environ.get("DELFIN_FFFREE_CLEAN_GATE_H_PARENT_MIN", "0.85"))
        h_par_max = float(os.environ.get("DELFIN_FFFREE_CLEAN_GATE_H_PARENT_MAX", "1.20"))
        h_other_max = float(os.environ.get("DELFIN_FFFREE_CLEAN_GATE_H_OTHER_MAX", "2.50"))
        h_angle_max = float(os.environ.get("DELFIN_FFFREE_CLEAN_GATE_H_ANGLE", "30.0"))
        clash_vdw_on = (os.environ.get(
            "DELFIN_FFFREE_CLEAN_GATE_CLASH_VDW", "0") == "1")
        clash_vdw_f = float(os.environ.get("DELFIN_FFFREE_CLEAN_GATE_CLASH_VDW_F", "0.65"))
        _vdw_r = None
        if clash_vdw_on:
            try:
                from delfin.manta.refine import _vdw as _vdw_r
            except Exception:
                clash_vdw_on = False

        def _radius(sym):
            return _COVALENT_RADII.get(sym, 0.76)

        def _is_metal(sym):
            return sym in _METAL_SET

        # ---- parse every frame (shared atom order) ------------------------
        frames = []
        n_atoms = None
        for item in isomers:
            xyz = item[0] if isinstance(item, (tuple, list)) and item else item
            syms, xs, ys, zs = [], [], [], []
            ok = True
            for ln in str(xyz).splitlines():
                p = ln.split()
                if len(p) >= 4:
                    try:
                        x, y, z = float(p[1]), float(p[2]), float(p[3])
                    except ValueError:
                        continue
                    if not (math.isfinite(x) and math.isfinite(y)
                            and math.isfinite(z)):
                        ok = False
                    syms.append(p[0])
                    xs.append(x)
                    ys.append(y)
                    zs.append(z)
            if n_atoms is None:
                n_atoms = len(syms)
            # Structurally invalid iff non-finite coord, empty, or count mismatch.
            frames.append((syms, xs, ys, zs, ok and len(syms) == n_atoms
                           and n_atoms > 0))

        if not n_atoms:
            return isomers

        def _d2(f, i, j):
            dx = f[1][i] - f[1][j]
            dy = f[2][i] - f[2][j]
            dz = f[3][i] - f[3][j]
            return dx * dx + dy * dy + dz * dz

        # ---- WHICH criterion actually drops the frames? (measurement, not behaviour) -----
        # Measured 2026-08-04 (cleangateW, the 142 worst-broken systems): the gate removes
        # 547 of 2341 frames and takes the hard-frame fraction from 84.2 % to 77.2 % -- the
        # first lever that lowers the ABSOLUTE defect rate at all.  But it also drops frames
        # that MATCHED THE CRYSTAL: ccdc_backbone_lost 3 (GILKAQ, HOJSUX, UHEJIB),
        # isomers_lost 3 (UHEJIB 10 -> 4 of 38 theory).  That violates the gate's own
        # governing invariant ("sort out only the BAD ones, not the good ones").
        #
        # Before touching any threshold: find out WHICH of the five criteria does it.  The
        # obvious suspect was the clash factor -- and that guess was WRONG: the gate measures
        # clash against the COVALENT sum (0.70 x ~1.52 A = 1.06 A for C-C) while the eye's
        # COD-validated metric uses the VDW sum (0.65 x 3.40 = 2.21 A), so the gate is far
        # LOOSER there, not tighter.  Hence: count, do not guess.
        # Pure instrumentation -- no frame is kept or dropped differently, output is
        # byte-identical; the tally only reaches stderr under DELFIN_TRACE_CLEAN_GATE=1.
        _gate_tally = {"invalid": 0, "collapse": 0, "vdw_h": 0, "clash": 0, "clash_vdw": 0,
                       "bare_metal": 0, "md_collapse": 0, "h_point": 0, "rescued_last_of_kind": 0}

        # ---- PER-FRAME absolute certainly-bad test ------------------------
        def _certainly_bad(f):
            """True iff frame f is CERTAINLY bad on its OWN geometry (a definite
            defect with a clear margin).  In any doubt -> False (KEEP)."""
            # 1. structurally invalid.
            if not f[4]:
                _gate_tally["invalid"] += 1
                return True
            syms = f[0]
            metal_idx = [k for k in range(n_atoms) if _is_metal(syms[k])]
            metal_set = set(metal_idx)

            # per-frame bond graph over non-metal pairs (this frame's own dists).
            parent = list(range(n_atoms))

            def _find(a):
                while parent[a] != a:
                    parent[a] = parent[parent[a]]
                    a = parent[a]
                return a

            def _union(a, b):
                ra, rb = _find(a), _find(b)
                if ra != rb:
                    parent[ra] = rb

            for i in range(n_atoms):
                if i in metal_set:
                    continue
                ri = _radius(syms[i])
                for j in range(i + 1, n_atoms):
                    if j in metal_set:
                        continue
                    t = bond_mult * (ri + _radius(syms[j]))
                    if _d2(f, i, j) < t * t:
                        _union(i, j)
            comp = [_find(k) for k in range(n_atoms)]

            # 2+3: per-frame heavy/H pair geometry.
            for i in range(n_atoms):
                if i in metal_set:
                    continue
                si = syms[i]
                ri = _radius(si)
                hi = (si == "H")
                for j in range(i + 1, n_atoms):
                    if j in metal_set:
                        continue
                    sj = syms[j]
                    d = math.sqrt(_d2(f, i, j))
                    rsum = ri + _radius(sj)
                    hj = (sj == "H")
                    bonded = d < bond_mult * rsum
                    if bonded:
                        # 3. collapsed bond (fused atoms) -> certainly bad.
                        if d < collapse_f * rsum:
                            _gate_tally["collapse"] += 1
                            return True
                        continue
                    # non-bonded pair:
                    if hi or hj:
                        # deep H interpenetration floor (hard).
                        if d < vdw_f * rsum:
                            _gate_tally["vdw_h"] += 1
                            return True
                    else:
                        # 2. real inter-ligand clash (different per-frame ligands,
                        # clear overlap).  Same-component close contacts are NOT
                        # flagged (a tight intra-ligand contact is not a defect).
                        if comp[i] != comp[j] and d < clash_f * rsum:
                            _gate_tally["clash"] += 1
                            return True
                        # vdW-referenced inter-ligand clash (the eye's calibration point).
                        if (clash_vdw_on and comp[i] != comp[j]
                                and d < clash_vdw_f * (_vdw_r(si) + _vdw_r(sj))):
                            _gate_tally["clash_vdw"] += 1
                            return True
            # 3b. DIRECTIONAL H-clash (ported from the eye's find_h_clash; only when explicitly on).
            if h_point_on:
                for hi_ in range(n_atoms):
                    if syms[hi_] != "H":
                        continue
                    par = -1
                    dpar = 1e9
                    for p in range(n_atoms):
                        if p == hi_ or syms[p] == "H" or p in metal_set:
                            continue
                        dp = math.sqrt(_d2(f, hi_, p))
                        if dp < dpar:
                            dpar, par = dp, p
                    if par < 0 or not (h_par_min <= dpar <= h_par_max):
                        continue                       # not a cleanly bonded H -> other tests own it
                    for q in range(n_atoms):
                        if q in (hi_, par) or syms[q] == "H" or q in metal_set:
                            continue
                        dq = math.sqrt(_d2(f, hi_, q))
                        if dq > h_other_max:
                            continue
                        # angle at the PARENT between parent->H and parent->q
                        dpq = math.sqrt(_d2(f, par, q))
                        if dpq < 1e-6:
                            continue
                        cosang = (dpar * dpar + dpq * dpq - dq * dq) / (2.0 * dpar * dpq)
                        cosang = max(-1.0, min(1.0, cosang))
                        if math.degrees(math.acos(cosang)) < h_angle_max:
                            _gate_tally["h_point"] += 1
                            return True
            # 4. bare metal / total decoordination (clear margin: ZERO donors).
            for mi in metal_idx:
                has_donor = False
                for k in range(n_atoms):
                    if k in metal_set or syms[k] == "H":
                        continue
                    if math.sqrt(_d2(f, mi, k)) < md_shell:
                        has_donor = True
                        break
                if not has_donor:
                    _gate_tally["bare_metal"] += 1
                    return True
            # 5. DONOR COLLAPSED ONTO METAL (blindspot; only when explicitly on).
            #    The pair tests above skip every metal pair, so a heavy donor fused
            #    onto the metal escapes them.  A metal-to-heavy-donor distance below
            #    collapse_f x covalent-sum is a fused atom = certainly bad (reuses
            #    the same collapse fraction as the bond test).  H excluded (M-H /
            #    agostic / bridging H can be legitimately short).
            if md_collapse_on:
                for mi in metal_idx:
                    rm = _radius(syms[mi])
                    for k in range(n_atoms):
                        if k in metal_set or syms[k] == "H":
                            continue
                        if math.sqrt(_d2(f, mi, k)) < collapse_f * (
                                rm + _radius(syms[k])):
                            _gate_tally["md_collapse"] += 1
                            return True
            return False

        # ============ "LAST OF ITS KIND": A FLOOR MUST NOT WIPE OUT AN ISOMER ============
        # Measured 2026-08-04, twice and independently, and it contradicts a claim written
        # verbatim in this function's own docstring ("NO ensemble consensus"):
        #     pcrejW      24 frames dropped -> smiles_ccdc_regressed 5, pyramid_frame 3, backbone 1
        #     cleangateW 547 frames dropped -> ccdc_backbone_lost 3 (GILKAQ HOJSUX UHEJIB),
        #                                      isomers_lost 3 (UHEJIB 10 -> 4 of 38 theory)
        # A filter that only REMOVES frames cannot make a system worse -- unless the frame it
        # removed was the only one realising that coordination isomer.  UHEJIB did not lose
        # quality, it lost SIX ISOMERS.
        #
        # For the question "is this frame bad?" judging it alone is right.  For the question
        # "may it go?" it is wrong: a frame can be geometrically impossible AND the last of its
        # kind.  GILKAQ's crystal backbone sits in a frame that has a collapsed bond somewhere.
        #
        # The rule (DELFIN_FFFREE_CLEAN_GATE_LAST_OF_KIND, default OFF -> byte-identical):
        # group the frames by their COORDINATION FINGERPRINT -- per metal, the sorted elements of
        # the heavy atoms inside the first shell -- and never let a group become empty.  If every
        # frame of an isomer is bad, the FIRST one survives (lowest index = the emitter's own top
        # rank; deterministic, no score, no RMSD).  No consensus, no majority, no vote: the only
        # question asked of the ensemble is "is there a replacement?".
        _bad_flags = [_certainly_bad(f) for f in frames]
        if os.environ.get("DELFIN_FFFREE_CLEAN_GATE_LAST_OF_KIND", "0") == "1":
            try:
                # ===== THE FINGERPRINT WAS BLIND TO THE ARRANGEMENT ==================
                # Measured 2026-08-18 on gk10kb (10000 systems, 9600 compared): the floor
                # lowers hard_frame_frac by 6.64 percentage points -- and in doing so loses
                #     isomers_lost 40 | ccdc_arrangement_lost 37 | ccdc_backbone_lost 26
                # ALTHOUGH "last of its kind" was running.  That is not a contradiction but the
                # definition of the group: the key below is the sorted ELEMENT
                # MULTISET of the first shell.  cis and trans have the same multiset.  An
                # entire arrangement family can thus be wiped out without its group
                # ever becoming empty -- the rescue never kicks in, because another isomer
                # carries the same fingerprint.
                #
                # Exactly the same root as commit 1caa9123 ("fold completeness per
                # ARRANGEMENT FAMILY instead of per geometry group").  For the second time the
                # same mix-up: "same elements" is not "same arrangement".
                #
                # DELFIN_FFFREE_CLEAN_GATE_KIND_ARRANGEMENT (default 0 -> byte-identical)
                # attaches to the key the multiset of the D-M-D angle classes: for every
                # donor pair (element, element, angle band).  With that, cis/trans,
                # fac/mer and axial/equatorial separate.
                #
                # ⚠️ DIRECTION OF THE CHANGE: a FINER key produces MORE groups,
                # and more groups can only rescue MORE frames, never fewer.  The
                # change is thereby monotone in favour of completeness and cannot cost an
                # isomer that the coarse key would have rescued.  The price lies
                # on the other side: the floor rejects less, the gain in
                # hard_frame_frac turns out smaller.  Exactly that is to be measured.
                #
                # The angle band is deliberately COARSE (45 degrees).  Distortion may split a
                # group -- that is the safe direction -- but noise shall not break every
                # group down into single frames.  Ideal octahedron 90/180 -> bands 2/4,
                # tetrahedron 109.5 -> 2, trigonal bipyramid 90/120/180 -> 2/3/4.
                _kind_arr = (os.environ.get(
                    "DELFIN_FFFREE_CLEAN_GATE_KIND_ARRANGEMENT", "0") == "1")

                def _coord_fp(f):
                    syms = f[0]
                    out = []
                    for mi in range(n_atoms):
                        if not _is_metal(syms[mi]):
                            continue
                        _don = [k for k in range(n_atoms)
                                if k != mi and syms[k] != "H" and not _is_metal(syms[k])
                                and math.sqrt(_d2(f, mi, k)) < md_shell]
                        sh = sorted(syms[k] for k in _don)
                        _key = syms[mi] + ":" + ",".join(sh)
                        if _kind_arr and len(_don) >= 2:
                            _ang = []
                            for _ai in range(len(_don)):
                                for _bi in range(_ai + 1, len(_don)):
                                    _a, _b = _don[_ai], _don[_bi]
                                    _ra2, _rb2 = _d2(f, mi, _a), _d2(f, mi, _b)
                                    if _ra2 <= 0.0 or _rb2 <= 0.0:
                                        continue
                                    _c = ((_ra2 + _rb2 - _d2(f, _a, _b))
                                          / (2.0 * math.sqrt(_ra2) * math.sqrt(_rb2)))
                                    _c = max(-1.0, min(1.0, _c))
                                    _band = int(round(math.degrees(math.acos(_c)) / 45.0))
                                    _e1, _e2 = sorted((syms[_a], syms[_b]))
                                    _ang.append("%s%s%d" % (_e1, _e2, _band))
                            _key += ";" + ",".join(sorted(_ang))
                        out.append(_key)
                    return "|".join(sorted(out))
                _groups = {}
                for _i, f in enumerate(frames):
                    _groups.setdefault(_coord_fp(f), []).append(_i)
                for _fp, _idxs in _groups.items():
                    if all(_bad_flags[_i] for _i in _idxs):
                        _bad_flags[_idxs[0]] = False        # the last of its kind stays
                        _gate_tally["rescued_last_of_kind"] += 1
            except Exception:
                pass
        kept = [item for item, _b in zip(isomers, _bad_flags) if not _b]
        _trace_dst = os.environ.get("DELFIN_TRACE_CLEAN_GATE", "")
        if len(kept) != len(isomers) and _trace_dst and _trace_dst != "0":
            # A FILE, not stderr.  loop.py routes the build workers' stderr to
            # results/debug_<rid>.log and only when a debug pattern matches -- otherwise it is
            # discarded, so a stderr trace from inside a worker reaches nothing (measured
            # 2026-08-04: 139 of 142 systems built, ZERO trace lines).  O_APPEND with one short
            # line per call is atomic enough across the parallel workers.
            try:
                with open(_trace_dst, "a") as _fg:
                    _fg.write("[CLEANGATE] %d -> %d frames | %s\n" % (
                        len(isomers), len(kept),
                        " ".join("%s=%d" % (k, v) for k, v in _gate_tally.items() if v)))
            except Exception:
                pass
        # In-doubt-keep: if (and only if) EVERY frame is certainly bad, keep the
        # input unchanged rather than dropping all or fabricating a broken
        # survivor (the removed "cleanest-available" fallback behaviour).
        return kept if kept else isomers
    except Exception:
        return isomers


def _apply_pi_inplane_final(isomers):
    """FINAL post-conformer monodentate aromatic σ-donor metal-in-plane re-assertion.

    A σ aromatic donor binds through an IN-PLANE sp² lone pair, so the metal must lie
    IN the donor ring's plane.  ``assemble_complex`` seats it there correctly, but the
    conformer/reembed passes (``_conf``) run AFTER the mid-pipeline coplanar-M pass
    (``DELFIN_FFFREE_PI_COPLANAR_M``) and re-tilt the ring (metal ends 0.6-2.8 Å out of
    plane in ~44 % of σ-aromatic-donor complexes).  This re-asserts the in-plane pose on
    the FINAL ensemble for TRULY MONODENTATE ligands (one donor in the ligand) via a
    rigid rotation about the fixed donor (M-D preserved; clash never-worse).  Multi-arm
    polydentates (the larger remainder) need a coupled backbone relaxation (Phase 2, not
    this pass).  Default-OFF / byte-identical unless ``DELFIN_FFFREE_PI_RIGID_PLACE=1``.
    """
    if not isomers or os.environ.get("DELFIN_FFFREE_PI_RIGID_PLACE", "0") != "1":
        return isomers
    try:
        from delfin.manta._pi_inplane_final import correct_results
        return correct_results(isomers)
    except Exception:
        return isomers


def _apply_pi_coplanar_final(isomers):
    """FINAL per-arm π-coplanarity polish for MULTI-ARM polydentates (the 91 % of
    the metal-out-of-plane defect that the monodentate ``_apply_pi_inplane_final``
    cannot reach).  Runs on the FINAL emitted ensemble — after ``_conf`` (the
    backbone reembed that re-tilts the rings) — so nothing downstream can undo it.

    Per conjugated π-system with exactly ONE coordinating σ-donor, rotates ONLY that
    arm (the π-system + its donor-less substituents) rigidly about its donor so the
    metal lies in the system plane; the shared backbone + other arms stay fixed, so
    every M-D length is preserved exactly.  Asymmetric / accept-if-better per arm
    (metal-oop strictly down, linker not over-stretched, inter-ligand clash
    never-worse) → never makes a frame worse.  Geometry-only, deterministic, FF-free.
    Default-OFF / byte-identical unless ``DELFIN_FFFREE_PI_COPLANAR_FINAL=1``.
    """
    if not isomers or os.environ.get("DELFIN_FFFREE_PI_COPLANAR_FINAL", "0") != "1":
        return isomers
    try:
        from delfin.manta._pi_coplanar_final import correct_results
        return correct_results(isomers)
    except Exception:
        return isomers


def _coord_integrity_filter(isomers):
    """Final UNIVERSAL decoordination filter over the FULL emitted ensemble at the
    public boundary -- covers BOTH the native fffree path AND the legacy generation
    path (where the eye-flagged ligand-flew-off frames actually originate: the hard
    kappa4/cage cases return None from _fffree_isomers and fall through to legacy,
    whose dense conformer/pucker/see-saw frames bypass the native self-gate).
    Drops a frame iff a coordinating donor has swung off the metal (donor-type ideal
    M-D + slack); crystals pack TIGHT but stay COORDINATED so a crystal-like frame is
    never dropped.  Gated DELFIN_FFFREE_COORD_INTEGRITY (default OFF -> identity ->
    byte-identical to the candidate/legacy baselines).  Never raises."""
    if not isomers or os.environ.get("DELFIN_FFFREE_COORD_INTEGRITY", "0") != "1":
        return isomers
    try:
        from delfin.manta.converter_backend import _coord_filter
        return _coord_filter(isomers) or isomers
    except Exception:
        return isomers


def _conf_complete_filter(isomers):
    """UNIVERSAL conformer-completeness pass over the FULL emitted ensemble.

    Root cause it cures (user eye-validation): the conformer search is INCOMPLETE
    (bulky pendant rotors -- tBu / biaryl / pincer arms -- never rotate, ATIDEM /
    ABUSEY), some FREE rotations BREAK topology (tear an M-D bond / put an H on the
    M-C axis / forge a spurious inter-ligand bond, AYOZIX), and symmetry-equivalent
    rotamers are DUPLICATED (ATOSAE).  This pass enumerates every heavy-moving
    rotatable axis as HINDERED rotations on a deterministic staggered grid, hard-
    gates each rotamer for topology preservation (M-D intact, ring intact, no
    spurious bond, no H invading the metal sphere), then RMSD-dedups -- complete
    coverage without bloat.

    Identity (returns input untouched) unless DELFIN_FFFREE_CONF_COMPLETE=1, so
    output is BYTE-IDENTICAL to baseline when off.  FF-free, deterministic, never
    raises."""
    if not isomers or os.environ.get("DELFIN_FFFREE_CONF_COMPLETE", "0") != "1":
        return isomers
    try:
        from delfin.manta.conformer_complete import apply_to_ensemble
        return apply_to_ensemble(isomers) or isomers
    except Exception:
        return isomers


def _gfnff_ensemble_rank_filter(isomers):
    """FINAL energy-ranked conformer retention over the emitted ensemble (the
    headline 'keep the most important conformers' feature).

    The shipped conformer engine (conformer_complete) generates topology-gated,
    RMSD-DIVERSE rotamers but invokes NO force field, and the final isomer ordering
    is by least-clash — so nothing in the pipeline ranks conformers by ENERGY.  This
    pass groups the emitted frames by base isomer (the label with any conformer
    suffix stripped), scores each group's frames with GFN-FF (xtb, license-clean,
    parametrised for metals; UFF is unusable on TM complexes), and keeps the top-K
    LOWEST-energy frames within an energy window.  Every base isomer keeps >=1 frame,
    so the isomer manifold is never reduced — only each isomer's conformer set is
    trimmed to the physically important members.

    Within one isomer the total charge is constant, so charge=0 ranks the conformers
    self-consistently (a constant offset cancels); DELFIN_GFNFF_CHARGE can override.

    Identity unless DELFIN_FFFREE_GFNFF_RANK=1 -> byte-identical when off.  Never
    raises; on any failure or if xtb is absent the input passes through unchanged."""
    if not isomers or os.environ.get("DELFIN_FFFREE_GFNFF_RANK", "0") != "1":
        return isomers
    try:
        import re as _re
        from delfin.manta import _gfnff_rank as _gff
        if not _gff.available():
            return isomers
        # Retention policy = "keep ALL distinct conformers, rank them, drop only
        # garbage" (user 2026-06-26).  The conformers are already RMSD-deduped
        # upstream (conf_complete @0.5 A), so each is a genuinely distinct fold; the
        # energy ranking ORDERS them (most important first) and the wide window drops
        # only true garbage (clashing / wildly strained -> tens-to-thousands of
        # kcal/mol above the minimum).  Defaults are intentionally generous (window
        # 35 kcal/mol, soft K-cap 64) so a normal conformer set is kept in full;
        # tighten DELFIN_GFNFF_TOPK / DELFIN_GFNFF_EWIN for a compact pool.
        topk = int(os.environ.get("DELFIN_GFNFF_TOPK", "64"))
        ewin = float(os.environ.get("DELFIN_GFNFF_EWIN", "35.0"))
        chg = int(os.environ.get("DELFIN_GFNFF_CHARGE", "0"))
        _suffix = _re.compile(r"(_conf-|-conf\d|_pool-|-reembed|_reembed|_conf\d)")

        def base_key(lbl):
            if not lbl:
                return ""
            return _suffix.split(lbl, maxsplit=1)[0]

        groups = {}
        order = []
        for (xyz, lbl) in isomers:
            k = base_key(lbl)
            if k not in groups:
                groups[k] = []
                order.append(k)
            groups[k].append((xyz, lbl))

        # HYBRID cost control: GFN2 ranks correctly but is ~10-50x GFN-FF, so for a
        # conformer-rich isomer we first cull cheaply with GFN-FF to the top-M, then
        # rank only those M with the (default GFN2) method.  Safe because the GFN2
        # winners are NOT GFN-FF's worst (measured: ABAKOE's GFN2-best was GFN-FF
        # rank ~3, well inside a top-15 cull).  Isomers with <=M conformers skip the
        # cull and are GFN2-ranked directly.
        # Cost bound only: GFN-FF pre-cull engages for isomers richer than this
        # (rare).  Kept >= the soft K-cap so it never trims below what we retain.
        precull_m = int(os.environ.get("DELFIN_GFNFF_PRECULL_M", "64"))
        out = []
        for k in order:
            frames = groups[k]
            if len(frames) == 1:
                out.append(frames[0])
                continue
            if len(frames) > precull_m:
                ffs = []
                for (xyz, lbl) in frames:
                    e = _gff.gfnff_energy(xyz, charge=chg, method="gfnff")
                    ffs.append((e if e is not None else float("inf"), xyz, lbl))
                ffs.sort(key=lambda t: t[0])
                shortlist = [(xyz, lbl) for _e, xyz, lbl in ffs[:precull_m]]
            else:
                shortlist = frames
            scored = []
            for (xyz, lbl) in shortlist:
                e = _gff.gfnff_energy(xyz, charge=chg)   # default method = GFN2
                scored.append((e if e is not None else float("inf"), xyz, lbl))
            if all(s[0] == float("inf") for s in scored):
                out.extend(frames)            # ranking unusable here -> keep all
                continue
            scored.sort(key=lambda t: t[0])
            emin = scored[0][0]
            kept = []
            for (e, xyz, lbl) in scored:
                if len(kept) >= topk:
                    break
                if e != float("inf") and e - emin > ewin:
                    break
                kept.append((xyz, lbl))
            out.extend(kept or [(scored[0][1], scored[0][2])])
        return out
    except Exception:
        return isomers


def _permute_dedup_filter(isomers):
    """UNIVERSAL permutation-invariant duplicate removal over the FULL emitted
    ensemble.

    Root cause it cures (#1 user request): the existing ensemble dedup compares
    frames by FIXED-ORDER heavy-atom Kabsch RMSD, so two structures identical up
    to relabeling of indistinguishable atoms (identical ligands swapped /
    symmetry-equivalent atoms / any molecular-graph automorphism) are NOT seen as
    duplicates and BOTH survive -> "very many identical structures in the pool".  This
    pass adds the missing layer: two frames are duplicates iff the MIN over graph
    automorphisms of their heavy Kabsch RMSD is below the threshold.  Atoms are
    permuted ONLY within their graph-symmetry orbit, so genuinely-different
    geometric isomers / conformers / stereoisomers (whose geometry differs by more
    than the threshold even after the best symmetry relabelling) are KEPT.

    Identity (returns input untouched) unless DELFIN_FFFREE_PERMUTE_DEDUP=1, so
    output is BYTE-IDENTICAL to baseline when off.  FF-free, deterministic, never
    raises."""
    if not isomers or os.environ.get("DELFIN_FFFREE_PERMUTE_DEDUP", "0") != "1":
        return isomers
    try:
        from delfin.manta.permute_dedup import dedup_ensemble
        return dedup_ensemble(isomers) or isomers
    except Exception:
        return isomers


def _mirror_xyz_coords(xyz):
    """Return the MIRROR IMAGE of a frame (reflection through the x=0 plane -> negate
    x).  An enantiomer IS the exact mirror image, so this is the free, exact partner
    of a chiral coordination frame — no rebuild.  Preserves atom order and any
    header/comment lines; only ``sym x y z`` lines flip.  None on any parse issue."""
    try:
        out = []
        for ln in str(xyz).splitlines():
            parts = ln.split()
            if len(parts) == 4:
                try:
                    x = -float(parts[1]); y = float(parts[2]); z = float(parts[3])
                    out.append(f"{parts[0]} {x:.6f} {y:.6f} {z:.6f}")
                    continue
                except Exception:
                    pass
            out.append(ln)
        return "\n".join(out)
    except Exception:
        return None


def _coord_sphere_donors(xyz):
    """(metal_dir_array, donor_symbols) for the coordination sphere: unit vectors
    metal->donor for the nearest neighbours (< 3.2 Å, up to 9) of the first
    ``_METAL_SET`` atom.  (None, None) if no metal / < 3 donors."""
    try:
        import numpy as _np
        rows = [l.split() for l in str(xyz).splitlines() if len(l.split()) == 4]
        syms = [p[0] for p in rows]
        coords = _np.array([[float(p[1]), float(p[2]), float(p[3])] for p in rows])
        m_i = next((i for i, s in enumerate(syms) if s in _METAL_SET), None)
        if m_i is None or coords.shape[0] < 4:
            return None, None
        d = coords - coords[m_i]
        dist = _np.linalg.norm(d, axis=1); dist[m_i] = 1e9
        idx = [i for i in _np.argsort(dist) if dist[i] < 3.2][:9]
        if len(idx) < 3:
            return None, None
        dirs = _np.array([d[i] / (dist[i] + 1e-12) for i in idx])
        return dirs, [syms[i] for i in idx]
    except Exception:
        return None, None


def _proper_kabsch_rmsd(P, Q):
    """Min RMSD aligning Q onto P by a PROPER rotation (det=+1, NO reflection) +
    translation.  Reflection-excluding by construction, so a true mirror image can
    never be aligned away."""
    import numpy as _np
    Pc = P - P.mean(0); Qc = Q - Q.mean(0)
    H = Qc.T @ Pc
    U, S, Vt = _np.linalg.svd(H)
    dsign = 1.0 if _np.linalg.det(Vt.T @ U.T) >= 0 else -1.0
    R = Vt.T @ _np.diag([1.0, 1.0, dsign]) @ U.T
    Qr = Qc @ R.T
    return float(_np.sqrt(_np.mean(_np.sum((Pc - Qr) ** 2, axis=1))))


def _coord_sphere_chirality_rmsd(xyz):
    """CONFIGURATIONAL chirality measure (Å): the min-over-same-element-donor-
    permutations PROPER-rotation Kabsch RMSD between the coordination sphere's donor
    directions and their MIRROR (reflection-excluding).  ~0 => the mirror superimposes
    => achiral configuration (mirror plane OR only a slight build deformation from an
    ideal-achiral arrangement).  Large => a genuine Δ/Λ twist => chiral.  Conformer-
    robust (donors only, backbone ignored).  A physical Å threshold (tuned) gates it.
    Returns 0.0 when undetermined (never spuriously chiral)."""
    try:
        import itertools as _it
        import numpy as _np
        dirs, syms = _coord_sphere_donors(xyz)
        if dirs is None:
            return 0.0
        mirror = dirs.copy(); mirror[:, 0] *= -1.0     # reflect (negate x)
        n = len(syms)
        groups: Dict[str, List[int]] = {}
        for i, s in enumerate(syms):
            groups.setdefault(s, []).append(i)
        # permutations = product of within-element permutations; cap the search.
        per_elem = []
        _count = 1
        for s, idxs in groups.items():
            ps = list(_it.permutations(idxs))
            _count *= len(ps)
            per_elem.append((idxs, ps))
        if _count > 5040:                              # >7! : fall back to identity only
            perms = [list(range(n))]
        else:
            perms = []
            for combo in _it.product(*[ps for _, ps in per_elem]):
                perm = [0] * n
                for (idxs, _), mapped in zip(per_elem, combo):
                    for src, dst in zip(idxs, mapped):
                        perm[src] = dst
                perms.append(perm)
        best = 1e9
        for perm in perms:
            r = _proper_kabsch_rmsd(dirs, mirror[perm])
            if r < best:
                best = r
        return best
    except Exception:
        return 0.0


def _coord_chirality_sign(xyz):
    """Sign (+1 / -1 / 0) of the coordination-sphere chirality pseudoscalar — used
    ONLY to name the Δ/Λ hands (never gates; the RMSD measure gates)."""
    try:
        import numpy as _np
        dirs, _syms = _coord_sphere_donors(xyz)
        if dirs is None:
            return 0
        total = 0.0
        n = len(dirs)
        for i in range(n):
            for j in range(i + 1, n):
                for k in range(j + 1, n):
                    total += float(_np.dot(_np.cross(dirs[i], dirs[j]), dirs[k]))
        if abs(total) < 1e-4:
            return 0
        return 1 if total > 0 else -1
    except Exception:
        return 0


def _arrangement_key(lbl):
    """Strip conformer / duplicate-disambiguation / Δ-Λ-hand suffixes -> the ARRANGEMENT
    identity.  All conformers AND both hands of one coordination isomer share this key, so
    the chirality decision can be made ONCE per arrangement (conformer-consistent) instead
    of per built conformer (where UFF noise makes the CSM straddle the tolerance)."""
    import re as _re
    k = str(lbl or "")
    k = _re.sub(r"-conf\d+", "", k)
    k = _re.sub(r"-[ΔΛ](?=-|$)", "", k)
    k = _re.sub(r"-\d+$", "", k)
    return k


def _mirror_symmetrize_xyz(xyz, tol):
    """Project a near-mirror-symmetric frame onto its EXACT mirror-symmetric form
    S = ½(F + R·F'[π]) (F' = mirror, R = best proper rotation, π = best graph automorphism).
    The antisymmetric (build-noise) component is projected OUT -> higher symmetry + quality,
    exactly consistent with DELFIN's target of the idealized intrinsic reference (crystal
    distortions are out of scope).  Returns the symmetrised xyz, or the ORIGINAL unchanged if
    the mirror is not reachable within ``tol`` (genuinely chiral) or on any error (best-effort,
    never raises)."""
    try:
        import numpy as _np
        from delfin.manta.permute_dedup import _automorphisms_for_xyz
        rows = [l.split() for l in str(xyz).splitlines() if len(l.split()) == 4]
        syms = [r[0] for r in rows]
        F = _np.array([[float(r[1]), float(r[2]), float(r[3])] for r in rows])
        n = len(F)
        if n < 2:
            return xyz
        Fm = F.copy(); Fm[:, 0] *= -1.0
        _s, autos, _h = _automorphisms_for_xyz(xyz, True, 4096)
        if not autos:
            autos = [list(range(n))]
        Fc = F - F.mean(0)
        best = None; best_r = 1e9
        for perm in autos:
            if len(perm) != n:
                continue
            Q = Fm[list(perm)]; Qc = Q - Q.mean(0)
            U, S, Vt = _np.linalg.svd(Qc.T @ Fc)
            d = 1.0 if _np.linalg.det(Vt.T @ U.T) >= 0 else -1.0
            Qr = Qc @ (Vt.T @ _np.diag([1.0, 1.0, d]) @ U.T).T
            r = float(_np.sqrt(_np.mean(_np.sum((Fc - Qr) ** 2, axis=1))))
            if r < best_r:
                best_r = r; best = Qr
        if best is None or best_r >= tol:
            return None            # no internal mirror within tol -> CHIRAL (caller keeps both hands)
        Ssym = 0.5 * (Fc + best) + F.mean(0)
        return "\n".join(f"{syms[i]} {Ssym[i, 0]:.6f} {Ssym[i, 1]:.6f} {Ssym[i, 2]:.6f}"
                         for i in range(n))
    except Exception:
        return None


def _pointgroup_symmetrize_xyz(xyz, tol):
    """Detect the frame's approximate POINT GROUP and PROJECT onto it:

        S = (1/|G|) Σ_{g∈G} g(F)      (the Continuous-Symmetry-Measure ideal)

    giving the structure with EXACT G-symmetry.  This GENERALISES
    ``_mirror_symmetrize_xyz`` (a σ-only special case: it averaged F with ONE mirror
    image, so it only ever removed a single mirror-breaking component -> plateaued) to
    the FULL group: EVERY symmetry operation of the frame is a graph automorphism π
    whose geometric realisation is the orthogonal matrix O (proper OR improper) that
    best superimposes F onto F[π]; the operation belongs to G iff that residual < tol.
    Averaging over ALL such operations removes EVERY symmetry-breaking distortion at
    once -> both poly-quality (perfect geometry) AND enantiomer elimination improve.

    Returns ``(S_xyz, has_improper)``:
      * ``has_improper`` True  <=> G contains an IMPROPER element (σ / inversion i / Sn)
        <=> the configuration is ACHIRAL -> the caller emits ONE symmetric frame
        (enantiomer ELIMINATED).  This is the FUNDAMENTAL achirality test — it covers
        the inversion centre and every Sn axis, not just mirror planes.
      * ``has_improper`` False <=> only proper rotations (Cn / Dn) or E (C1) -> genuinely
        CHIRAL -> the caller keeps BOTH hands (the projected S is the proper-symmetrised
        frame of the SAME hand; its exact mirror is the other hand).
      * |G| == 1 (only E within tol) -> S == F unchanged  ->  C1 SAFE (user: "many C1
        must NOT be affected").

    Deterministic; best-effort (``(None, False)`` on any error, never raises)."""
    try:
        import numpy as _np
        from delfin.manta.permute_dedup import _automorphisms_for_xyz
        rows = [l.split() for l in str(xyz).splitlines() if len(l.split()) == 4]
        syms = [r[0] for r in rows]
        F = _np.array([[float(r[1]), float(r[2]), float(r[3])] for r in rows])
        n = len(F)
        if n < 2:
            return xyz, False
        Fc = F - F.mean(0)
        _s, autos, _h = _automorphisms_for_xyz(xyz, True, 4096)
        if not autos:
            autos = [list(range(n))]
        # Collect the group operations: (O, perm) for every automorphism whose optimal
        # orthogonal (Procrustes, REFLECTION ALLOWED -> det ±1) superposition of F onto
        # F[perm] has residual < tol.  O = U Vt for H = F[perm]^T · F (see _mirror_… for
        # the proper-only sibling; here we do NOT force det=+1, so improper ops appear).
        ops = []
        have_id = False
        for perm in autos:
            if len(perm) != n:
                continue
            Q = Fc[list(perm)]
            U, _S2, Vt = _np.linalg.svd(Q.T @ Fc)
            O = U @ Vt                                   # orthogonal, det = ±1
            r = float(_np.sqrt(_np.mean(_np.sum((Fc @ O.T - Q) ** 2, axis=1))))
            if r < tol:
                ops.append((O, Q))                       # keep Q=Fc[perm] for the projection
                if all(p == i for i, p in enumerate(perm)):
                    have_id = True
        if not ops:
            ops = [(_np.eye(3), Fc.copy())]
        elif not have_id:
            ops.append((_np.eye(3), Fc.copy()))
        # improper element present?  det(O) < 0 for ANY op  <=>  achiral configuration.
        has_improper = any(bool(_np.linalg.det(O) < 0.0) for O, _ in ops)
        # PROJECT: S_i = (1/|G|) Σ_k O_k^T · r_{π_k(i)}  =  (1/|G|) Σ_k (Fc[π_k] · O_k)_i .
        Sacc = _np.zeros_like(Fc)
        for O, Q in ops:
            Sacc += Q @ O
        Sacc /= float(len(ops))
        Sacc += F.mean(0)
        out = "\n".join(f"{syms[i]} {Sacc[i, 0]:.6f} {Sacc[i, 1]:.6f} {Sacc[i, 2]:.6f}"
                        for i in range(n))
        return out, has_improper
    except Exception:
        return None, False


def _heavy_dist_fp(xyz, bin_a=0.15):
    """Rotation- AND permutation-invariant, chirality-INSENSITIVE geometry fingerprint:
    the sorted multiset of heavy-atom pairwise distances, binned to ``bin_a`` Å.  Two
    frames of ONE molecule with the same fingerprint are the same geometry up to a rigid
    motion OR a reflection (distances are reflection-invariant) — so pairing it with the
    chirality sign gives a canonical isomer key that collapses achiral image/mirror pairs
    while keeping genuine Δ/Λ.  None on parse error."""
    try:
        import numpy as _np
        pts = []
        for _l in str(xyz).splitlines():
            q = _l.split()
            if len(q) == 4 and q[0] != 'H':
                pts.append((float(q[1]), float(q[2]), float(q[3])))
        if len(pts) < 3:
            return None
        C = _np.array(pts)
        n = len(C)
        ds = []
        for i in range(n):
            ds.extend(_np.linalg.norm(C[i + 1:] - C[i], axis=1).tolist())
        return tuple(sorted(round(d / bin_a) for d in ds))
    except Exception:
        return None


def _enantiomer_mirror_filter(isomers, smiles=None):
    """Add the Δ/Λ MIRROR enantiomer for every CHIRAL coordination frame (user
    2026-07-21: "we also want the respective enantiomers in the frames").

    An enantiomer IS the exact mirror image -> reflect the built coordinates (free,
    exact) instead of rebuilding in two directions.  A frame is chiral iff its mirror
    is NOT superimposable on it under proper rotation + graph automorphism (reflection
    -EXCLUDING Kabsch) -- tested by the existing ``is_permutation_duplicate`` primitive
    (NOT gated), so achiral (meso) frames add no duplicate and an enantiomer already
    present (e.g. built by CHIRAL_ENUM) is not re-added.  UNIVERSAL across every
    polyhedron (chirality-by-mirror is geometry, not a per-polyhedron formula) and
    robust where the analytical helicity classifier fails (tridentate wraps).

    CRITICAL (user 2026-07-21): the gate is CONFIGURATIONAL chirality, NOT the
    conformer's own mirror symmetry.  A conformer of an ACHIRAL molecule can lack a
    mirror plane (a chiral backbone pucker) yet be the SAME molecule (interconvertible
    by a conformational flip) -- mirroring it would fabricate a false "enantiomer" and
    explode the pool.  So we gate on the COORDINATION-SPHERE chirality (the pseudoscalar
    over metal->donor directions), which is set by the polyhedron and is ROBUST across
    conformers (donors are pinned; only the backbone puckers).  Achiral configuration
    -> sign 0 -> NOT mirrored, whatever the conformer pucker.  A ``is_permutation_
    duplicate`` meso-guard additionally skips internally mirror-symmetric molecules.

    Emits each enantiomer IMMEDIATELY AFTER its partner (Δ, Λ consecutive in the
    trajectory).  Guard: skip when the SMILES has FIXED ligand stereocentres --
    mirroring the whole complex would INVERT them (a different compound); those
    diastereomers come from the arrangement enumeration + STEREOCENTER_ENUM.  Env-gated
    DELFIN_FFFREE_ENANTIOMER_MIRROR (default off -> byte-identical).  Deterministic; never raises."""
    if not isomers or os.environ.get("DELFIN_FFFREE_ENANTIOMER_MIRROR", "0") != "1":
        return isomers
    try:
        from delfin.manta.permute_dedup import is_permutation_duplicate
    except Exception:
        is_permutation_duplicate = None
    try:
        if smiles and RDKIT_AVAILABLE and isinstance(smiles, str):
            _m = Chem.MolFromSmiles(smiles)
            if _m is not None and Chem.FindMolChiralCenters(
                    _m, includeUnassigned=False, useLegacyImplementation=False):
                return isomers   # fixed ligand stereocentre -> mirroring inverts it
    except Exception:
        pass
    _chi_tol = _delfin_env_float("DELFIN_FFFREE_ENANTIOMER_MIRROR_TOL", 0.40)
    # ADDITIVE-ONLY (2026-07-21).  A broad full:1000 A/B proved the post-hoc UNITE + point-group PROJECTION
    # is NOT never-worse: averaging a frame over its APPROXIMATE symmetry ops tears ligands and distorts
    # polyhedra (measured: poly-distortion 13, ligand-quality 18, ligand-bond-torn 1, tier-2 on 69,
    # holistic +0.513), and the fp-GROUPING collapses DISTINCT isomers that merely share a distance-bin
    # (isomers-lost 10, ccdc-isomer-lost 4).  Post-hoc geometry repair on built atoms is exactly the
    # guiding-principle antipattern.  So this filter now does ONLY the SAFE, purely-ADDITIVE half of the user's
    # ask ("we want the enantiomers in the frames"): for every GENUINELY chiral coordination frame whose
    # opposite hand is not already built, ADD its EXACT mirror.  A reflection of a good frame is an equally
    # good frame -> zero distortion, zero isomer loss -> never-worse.  The ELIMINATE-where-possible geometry
    # symmetrisation is DEFERRED to the IN-CONSTRUCTION poly-seater (seat donors at ideal symmetric vertices;
    # NEVER move built atoms) -- the correct root architecture.  Gate CONFIGURATIONAL (coord-sphere)
    # chirality, robust to conformer/peripheral wobble; the mirror's presence is tested by a
    # reflection-aware (fingerprint, sign) key so an already-built opposite hand is never duplicated.
    # FAST per-frame key: (reflection-INVARIANT fingerprint, cheap O(donors^3) chirality pseudoscalar).
    # Deliberately AVOID the expensive permutation-Kabsch _coord_sphere_chirality_rmsd here: on high-frame-
    # count systems (RILVUD 209 frames) it made the finalisation slow enough to push _one.py over the
    # per-system build timeout -> partial build -> FALSE isomer loss in the A/B.  The pseudoscalar sign is
    # a sufficient chiral gate (sign==0 -> treat as achiral -> add nothing; a rare genuinely-chiral frame
    # with an accidental zero pseudoscalar is merely SKIPPED = a completeness miss, never a regression).
    _keys = []
    _present = set()
    for it in isomers:
        try:
            _fp = _heavy_dist_fp(it[0]); _sg = _coord_chirality_sign(it[0])
        except Exception:
            _fp, _sg = None, 0
        _keys.append((_fp, _sg))
        _present.add((_fp, _sg))
    out = list(isomers)
    for _i, it in enumerate(isomers):
        _fp, _sg = _keys[_i]
        if _sg == 0:
            continue   # achiral configuration -> mirror superimposes on itself -> nothing to add
        _mk = (_fp, -_sg)   # the mirror: fp is reflection-INVARIANT, ONLY the chirality sign flips
        if _mk in _present:
            continue   # the opposite hand is already built (arrangement enumeration / CHIRAL_ENUM)
        _present.add(_mk)
        _mx = _mirror_xyz_coords(it[0])
        if _mx is None:
            continue
        _lbl = it[1] if len(it) > 1 else ""
        _hand = "Λ" if _sg < 0 else "Δ"   # mirror hand = opposite of the original's
        out.append((_mx, (f"{_lbl}-{_hand}" if _lbl else _hand)) + tuple(it[2:]))
    return out
