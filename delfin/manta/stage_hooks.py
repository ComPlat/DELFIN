"""The additive post-construction stages of the MANTA constructor, each behind its own DELFIN_FFFREE_* switch and delegating to a delfin.manta module.

Moved verbatim from delfin/smiles_converter.py (MANTA split, 2026-10); every
statement is the original text, only the imports below are new.
"""
from __future__ import annotations

import os
from typing import List, Tuple

from delfin.common.logging import get_logger
from delfin.manta.converter_flags import (
    _class_conditional_flag,
    _delfin_env_float,
    _delfin_env_int,
    _trace_seating,
)
from delfin.manta.hapto_detect import (
    _classify_complex_class,
    _find_hapto_groups,
)
from delfin.manta.ml_tables import (
    _METAL_SET,
    _get_ml_bond_length,
)
from delfin.manta.topology_checks import (
    _count_xyz_clashes,
    _fragment_topology_ok,
)

logger = get_logger("delfin.smiles_converter")


def _apply_h_placement_if_enabled(results):
    """H placement for the FF-free path (19.08.2026, DELFIN_FFFREE_H_PLACEMENT).

    Needs NO ``mol``: all three stages work on the XYZ text, so the atom-order trap that
    the docstring of ``_ffree_shared_tail`` warns about does not apply.  That is exactly
    what ``FIX_F19`` fails on, although its motion class is actually the best in the tree --
    it reads RDKit hybridisation and needs ``frame_to_mol_map``.

    ⚠ THE MODULE CHECKS THE SWITCH ITSELF.  Deliberately NOT again here: a switch that is
    read in two places drifts.

    ⚠ WHY THIS IS SAFE: heavy atoms are pinned byte-exactly and checked at the end
    (``aborted_heavy_moved``); only H are moved.  Stage B keeps the DIRECTION and
    changes only the length, so it can flip neither an angle nor a sign.  For stages A
    and C there is a stereo gate: if in a frame a real candidate flips or a centre gets
    flattened, the WHOLE frame rolls back.  Measured: of 1076 real candidates exactly one
    flipped, and the gate costs nothing.

    ⚠ STAGE A REBUILDS NOTHING, it DELEGATES to ``_vsepr_repair`` -- the only one of the
    seven existing H mechanisms that is unrefuted and that, solely because of its position
    (3038 lines behind the FF-free return), never fired FF-free.
    """
    try:
        from delfin.manta import _h_placement as _hp
    except Exception:
        return results
    try:
        return _hp.apply_to_results(results)
    except Exception:
        return results


def _apply_mirror_enum_if_enabled(results):
    """Mirror completion of the manifold (17.08.2026).

    Appends to each frame its mirror image.  Needs NO ``mol``: a mirroring is a pure
    coordinate operation, so the atom-order trap that the docstring of
    ``_ffree_shared_tail`` warns about does not apply.

    ⚠ THE MODULE CHECKS THE SWITCH ITSELF (``_mirror_enum._is_enabled``).  Deliberately NOT
    again here: a switch that is read in two places drifts -- and the mix-up "setting value
    mistaken for module switch" (``TORSION_GRID`` versus ``TORSION_RELAX``, 17.08.) came
    from exactly such a duplication.

    ⚠ WHY THIS IS SAFE: a mirroring is an ISOMETRY.  Every bond length, every angle,
    every M-D distance stays EXACT (self-test: distance-matrix delta 0.0); only the signs
    flip.  No relaxation, no clash possible.

    ⚠ AND WHY THE FRAME SURVIVES: both dedups are determinant-corrected and explicitly
    forbid mirrorings (``permute_dedup._kabsch_rmsd_perm:197``,
    ``assemble_complex._kabsch_rot:63``) -- an enantiomer never aligns and is therefore
    NOT rejected as a duplicate.  Without this property the pass would have been a null
    test; it was checked in the source BEFORE building, not assumed.

    Default of the module switch is 0 -> byte-identical.
    """
    if not results:
        return results
    try:
        from delfin.manta._mirror_enum import expand_results as _mirror_expand
        return _mirror_expand(results)
    except Exception as _mx:
        try:
            logger.debug("MIRROR-ENUM skipped: %s", _mx)
        except Exception:
            pass
        return results


def _apply_me_bond_snap_if_enabled(results):
    """Set terminal M=E multiple bonds on the FF-FREE path (17.08.2026).

    Calls ``delfin.manta._me_bond_snap.snap_me_bonds`` per frame.  The module works
    exclusively on the frame (symbols + coordinates + geometric adjacency) and needs
    NO ``mol`` -- so the atom-order trap, on which the ring-pucker emitter once already
    became a null test at exactly this position, does not apply.

    ⚠ THE SWITCH IS THE SAME as for the legacy setting sites
    (``DELFIN_FFFREE_ME_BOND_LEN``), not a second one.  One axis, one switch --
    otherwise an A/B measures not the mechanism but the wiring.

    ⚠ WITHOUT A CALIBRATED TABLE NOTHING HAPPENS, and that is a license condition, not a
    convenience: ``_ml_me_band`` returns ``None`` when ``DELFIN_ML_ME_BANDS`` is missing
    or the pair is below the minimum sample size.  Then ``_get_ml_bond_length`` falls
    back to the sigma value -- the difference is zero, and ``_target`` returns ``None``.
    The comparison against the sigma value is thus the ONLY place where it is decided
    whether a pair has a band at all; no CCDC number stands here.

    Default 0 -> not called -> byte-identical.
    """
    if not results:
        return results
    if not _delfin_env_int("DELFIN_FFFREE_ME_BOND_LEN", 0):
        return results
    try:
        from delfin.manta._me_bond_snap import snap_me_bonds

        def _target(m_sym: str, d_sym: str):
            try:
                sigma = float(_get_ml_bond_length(m_sym, d_sym, "sigma"))
                me = float(_get_ml_bond_length(m_sym, d_sym, "me"))
            except Exception:
                return None
            if me <= 0.0 or abs(me - sigma) < 1e-9:
                return None          # no band for this pair -> do nothing
            return me

        new_results: List[Tuple[str, str]] = []
        for (xyz, label) in results:
            try:
                new_results.append((snap_me_bonds(xyz, _target), label))
            except Exception:
                new_results.append((xyz, label))
        return new_results
    except Exception as _me_exc:
        try:
            logger.debug("ME-BOND-SNAP skipped: %s", _me_exc)
        except Exception:
            pass
        return results


def _apply_coord_angle_fix_if_enabled(mol, results, dual_parse_done: bool):
    """Iter-13 Baustein 3 dispatch helper.

    Apply post-ETKDG/UFF coordination-angle correction to ``results`` if a
    metal is present, this is the outer (non dual-parse) call, and the
    DELFIN_FFFREE_COORD_ANGLE_FIX env-flag is set.  Fail-safe: any exception in the
    corrector returns ``results`` unchanged.

    Centralised here so every scaffold-path return point in
    ``smiles_to_xyz_isomers`` (mono σ, mono hapto, multi σ-σ, multi
    σ-hapto, multi hapto-hapto, fallback single-conformer) can share one
    insertion site without code drift.

    Bit-exact when ``DELFIN_FFFREE_COORD_ANGLE_FIX=0`` (default).  No effect on results
    that contain no metal.  ``mol`` may be ``None`` — the underlying
    corrector operates on XYZ text only.
    """
    if not results:
        return results
    if dual_parse_done:
        return results
    if not _delfin_env_int("DELFIN_FFFREE_COORD_ANGLE_FIX", 0):
        return results
    # Quick element scan — if no metal symbol appears in any result XYZ,
    # the corrector would no-op anyway.  Cheap pre-check avoids import
    # cost for organic/no-metal complexes that share this path.
    try:
        any_metal = False
        for entry in results:
            xyz_text = entry[0] if entry else ""
            if not isinstance(xyz_text, str):
                continue
            for ln in xyz_text.splitlines():
                tok = ln.strip().split(" ", 1)[0] if ln.strip() else ""
                if tok and tok[0].isalpha() and len(tok) <= 2:
                    # Defer real metal classification to corrector;
                    # break early once we see a plausibly metallic symbol.
                    if tok not in ("H", "C", "N", "O", "F", "P", "S",
                                   "Cl", "Br", "I", "B", "Si"):
                        any_metal = True
                        break
            if any_metal:
                break
        if not any_metal:
            return results
    except Exception:
        # Pre-check error → fall through to full corrector (it is fail-safe).
        pass
    try:
        from delfin.manta._coord_angle_corrector import correct_results as _b3_correct
        return _b3_correct(mol, results)
    except Exception as _b3_exc:
        try:
            logger.debug("Baustein 3 coord-angle correction skipped: %s", _b3_exc)
        except Exception:
            pass
        return results


def _apply_5j_a_cp_piano_stool_if_enabled(mol, results, dual_parse_done: bool):
    """Welle-5j Agent A dispatch helper — Cp piano-stool hapticity refinement.

    Welle-5i Agent C catalogued 34 hapto BROKEN-TO-BROKEN files; 83 % (28 / 34)
    were η⁵-cyclopentadienyl coordination misclassified as η⁶-arene by the
    downstream hapticity / coord-geometry detector.  This dispatch helper
    runs ``delfin.manta._cp_piano_stool.correct_results`` on every metal-bearing
    structure to snap M-ring-centroid distance + axial orientation to the
    ideal η⁵ piano-stool geometry, which makes the detector classify the
    ring as Cp (CN=5 polyhedron) instead of arene.

    Insertion order: AFTER ``_apply_coord_angle_fix_if_enabled`` (B3 rotates
    donor X-side) and BEFORE ``_apply_baustein4_if_enabled`` (B4 then
    re-projects ring-attached H onto the post-snap ring plane).

    Per-conformer / per-violation rollback inside the helper.  Bit-exact
    when ``DELFIN_5J_A_CP_PIANO_STOOL=0`` (default).  Skipped on inner
    dual-parse calls (matches B3 / B4 dispatch contract — heavy-atom
    signature dedup must see consistent coordinates).
    """
    if not results:
        return results
    if dual_parse_done:
        return results
    if not _delfin_env_int("DELFIN_5J_A_CP_PIANO_STOOL", 0):
        return results
    try:
        from delfin.manta._cp_piano_stool import correct_results as _cp_correct
        return _cp_correct(mol, results)
    except Exception as _cp_exc:
        try:
            logger.debug("5j-A Cp piano-stool refinement skipped: %s", _cp_exc)
        except Exception:
            pass
        return results


def _apply_baustein4_if_enabled(mol, results, dual_parse_done: bool):
    """Iter-14 Baustein 4 dispatch helper — RigidPiFragment H projection.

    Apply post-ETKDG/UFF ring-attached H projection to ``results`` if this
    is the outer (non dual-parse) call and the DELFIN_BAUSTEIN4 env-flag
    is set.  π-rigid-body invariant: every transformation that moves a
    π-frame must drag attached H atoms rigidly with the ring.  Iter-9 H1 and Iter-12/13 B3 already
    handle some cases; B4 is the universal final pass.

    Bit-exact when ``DELFIN_BAUSTEIN4=0`` (default).  Operates on XYZ text
    only (no RDKit dependency for the projection step).  ``mol`` may be
    ``None``.

    Insertion order: B3 first (rotates X-side around donor — moves H atoms
    on rotated bonds), B4 second (projects H onto current ring planes).
    """
    if not results:
        return results
    if dual_parse_done:
        return results
    if not _delfin_env_int("DELFIN_BAUSTEIN4", 0):
        return results
    try:
        from delfin.manta._pi_h_projector import correct_results as _b4_correct
        new_results = _b4_correct(mol, results)
        # Iter-6 T3 — B4-clash-rollback (preventive; B4 default OFF anyway).
        # If DELFIN_B4_CLASH_ROLLBACK=1: count heavy/H non-bonded clashes
        # before vs after B4; rollback per-frame on regression.  Bit-exact
        # to existing behaviour when the env-flag is 0 (default).
        if _delfin_env_int("DELFIN_B4_CLASH_ROLLBACK", 0):
            try:
                if isinstance(new_results, list) and len(new_results) == len(results):
                    rolled: List[Tuple[str, str]] = []
                    for (old_xyz, old_lbl), new_pair in zip(results, new_results):
                        try:
                            new_xyz, new_lbl = new_pair
                        except Exception:
                            rolled.append((old_xyz, old_lbl))
                            continue
                        try:
                            n_old = _count_xyz_clashes(old_xyz)
                            n_new = _count_xyz_clashes(new_xyz)
                        except Exception:
                            rolled.append((new_xyz, new_lbl))
                            continue
                        if n_new > n_old:
                            try:
                                logger.debug(
                                    "B4 rollback: clashes %d->%d, frame %s",
                                    n_old, n_new, new_lbl,
                                )
                            except Exception:
                                pass
                            rolled.append((old_xyz, old_lbl))
                        else:
                            rolled.append((new_xyz, new_lbl))
                    return rolled
            except Exception as _b4_roll_exc:
                try:
                    logger.debug("B4 clash-rollback skipped: %s", _b4_roll_exc)
                except Exception:
                    pass
        return new_results
    except Exception as _b4_exc:
        try:
            logger.debug("Baustein 4 H-projection skipped: %s", _b4_exc)
        except Exception:
            pass
        return results


def _apply_stereocenter_enum_if_enabled(mol, results, dual_parse_done: bool):
    """Stereocentre-fold completeness dispatch (env-gated DELFIN_STEREOCENTER_ENUM, default OFF).

    ADDITIVELY appends every buildable coordination-created X-H stereocentre fold (both R and S at
    each such donor — the completeness law) that the base manifold does not already contain, so the
    eye's ``ccdc_isomer_realized`` hard floor (is the crystal's N-H fold present?) can go
    FALSE->TRUE without ever dropping an existing frame.  Deterministic; bit-exact no-op when the
    env-flag is unset or the structure has no metal-created X-H stereocentre.  Skipped on the inner
    dual-parse call (matches the B3/B4 dispatch contract — the union dedup must see one consistent
    frame set).  ``mol`` is unused (the corrector operates on XYZ text only, like Baustein 3).
    """
    if not results:
        return results
    if dual_parse_done:
        return results
    try:
        from delfin.manta import _stereocenter_enum as _sc
        if not _sc._is_enabled():
            return results
        return _sc.expand_results(results)
    except Exception as _sc_exc:
        try:
            logger.debug("stereocentre-enum expansion skipped: %s", _sc_exc)
        except Exception:
            pass
        return results


def _apply_atropisomer_enum_if_enabled(mol, results, dual_parse_done: bool):
    """AXIAL completeness -- env-gated, default OFF (the switch name stands in exactly ONE
    place: `_atropisomer_enum._atrop_enabled`; here it is only asked, never guessed).

    Additively appends, for every stereogenic axis (biaryl AND mesomeric aryl amide), the
    MISSING sign, so that `ccdc_atropisomer_realized` can go from FALSE to TRUE without a
    frame ever being lost.

    WHY.  Measured on 16.08.2026 on 965 systems: the axis stands at **0 of 44** -- and
    the cause is NOT the geometry (mean twist 42.8 degrees over 145 axes,
    BIRVUW reaches 85.2 degrees), but completeness: **41 % of the axes carry only
    ONE handedness**, only 33 % both.  Since the axis demands that EVERY axis of a
    molecule finds its crystal sign, at ~3.3 axes per system a full hit is
    almost ruled out.

    Mirrors the contract of `_apply_stereocenter_enum_if_enabled`: additive, deterministic,
    bit-exact no-op when the switch is off or no stereogenic axis exists; skipped on the
    inner dual-parse call (the union dedup must see ONE consistent
    frame set).  ``mol`` is unused -- the corrector works on XYZ text only.
    """
    if not results:
        return results
    if dual_parse_done:
        return results
    try:
        from delfin.manta import _atropisomer_enum as _at
        if not _at._atrop_enabled():
            return results
        return _at.expand_atropisomers(results)
    except Exception as _at_exc:
        try:
            logger.debug("atropisomer-enum expansion skipped: %s", _at_exc)
        except Exception:
            pass
        return results


def _apply_bond_decollapse_if_enabled(mol, results, dual_parse_done: bool):
    """Iter-25 (2026-05-20) dispatch — final bond-decollapse corrector.

    Calibration breakthrough: ~79% of hapto structures emit collapsed ligands
    (heavy-heavy pairs 0.24-1.2 A — fused substituents / overlapping rings).
    This FINAL pass spring-relaxes the non-metal heavy graph to physical bond
    lengths + repels superimposed atoms, freezing metals + the coordination
    sphere (M-D invariant preserved).  Per-frame rollback: kept only if the
    collapsed-bond count strictly drops and no M-D bond breaks.  Geometry-only
    (atom-order independent).  Aggregate validated -20.5pp collapse on a
    150-hapto sample.  Class-cond default-ON {hapto, multi_hapto}; bit-exact
    when flag off / class excluded / no collapsed bond present.
    """
    if not results:
        return results
    # DIAGNOSTIC/ERDBEBEN: DELFIN_BOND_DECOLLAPSE_FORCE=1 bypasses BOTH the dual-parse skip and the
    # hapto-only class scope, so the spring-relax de-collapse runs on ANY class (e.g. sigma AQIBAE) -- to
    # test whether the existing pass fixes the metal-context cage collapse before building a re-embed pass.
    _bd_force = os.environ.get("DELFIN_BOND_DECOLLAPSE_FORCE", "0") == "1"
    if dual_parse_done and not _bd_force:
        return results
    if not _bd_force and not _class_conditional_flag(
        "DELFIN_BOND_DECOLLAPSE", mol, default=0,
        default_classes=["hapto", "multi_hapto"],
    ):
        return results
    try:
        from delfin.manta._bond_decollapse import correct_results as _bd_correct
        _bd_out = _bd_correct(mol, results)
        if os.environ.get("DELFIN_TRACE_SEATING", "0") == "1":
            _trace_seating("BOND_DECOLLAPSE ran: %d frames (force=%s dual_parse=%s)"
                           % (len(results), _bd_force, dual_parse_done))
        return _bd_out
    except Exception as _bd_exc:
        try:
            logger.debug("Iter-25 bond-decollapse skipped: %s", _bd_exc)
        except Exception:
            pass
        return results


def _apply_isolated_reseat_if_enabled(mol, results, dual_parse_done: bool):
    """ERDBEBEN dispatch — re-seat planar-COLLAPSED ligand fragments from a clean ISOLATED embed
    (delfin.manta._isolated_reseat).  The metal-context whole-complex ETKDG collapses rigid cages
    (AQIBAE) though the isolated fragment embeds 3D 20/20; the spring-relax de-collapse cannot pop the
    planar local minimum, only a fresh embed can.  Runs FIRST among the final passes so the subsequent
    aromatic-planarity / bond-length passes refine the re-seated fragment.  Gated
    DELFIN_FFFREE_ISOLATED_SEAT (default off -> byte-identical); per-frame rollback keeps it never-worse
    (collapse must drop, no M-D break, no worse clash).  Runs on ALL classes -- the collapse is not
    class-specific -- and regardless of dual_parse (a collapsed frame must be fixed either way)."""
    if not results:
        return results
    if os.environ.get("DELFIN_FFFREE_ISOLATED_SEAT", "0") != "1":
        return results
    try:
        from delfin.manta._isolated_reseat import correct_results as _ir_correct
        _ir_out = _ir_correct(mol, results)
        return _ir_out
    except Exception as _ir_exc:
        try:
            logger.debug("ERDBEBEN isolated-reseat skipped: %s", _ir_exc)
        except Exception:
            pass
        return results


def _apply_aromatic_planarity_if_enabled(mol, results, dual_parse_done: bool):
    """Iter-24 (2026-05-20) dispatch — post-UFF aromatic-ring flattening.

    Re-grounding forensic: hapto/multi_hapto true-aromatic rings (M_coord
    chelate rings excluded) pucker 72-75 % @ mean OOP 0.34-0.38 Å.  Flatten
    them onto their SVD best-fit plane, centroid-preserving (M-ring distance
    invariant) + ring-H dragged + per-frame M-D-invariant rollback.

    Class-conditional default-ON for {hapto, multi_hapto} via
    ``DELFIN_AROMATIC_PLANARITY`` (sigma 25 % is borderline/threshold-noise,
    multi_sigma 10 % already good — excluded).  Bit-exact when the flag is 0
    or the class is excluded or no qualifying ring is present.
    """
    if not results:
        return results
    if dual_parse_done:
        return results
    # Iter-24 WIP: default-OFF.  Pre-commit forensic showed the naive
    # atom-wise flatten is inconsistent (some pendant rings worsen; coordinated
    # η-rings trip the M-D invariant).  Needs redesign (centroid-based M-D for
    # coordinated rings + fix flatten/H-projection interaction) before default-ON.
    if not _class_conditional_flag(
        "DELFIN_AROMATIC_PLANARITY", mol, default=0,
    ):
        return results
    try:
        from delfin.manta._aromatic_ring_flattener import correct_results as _arom_correct
        return _arom_correct(mol, results)
    except Exception as _arom_exc:
        try:
            logger.debug("Iter-24 aromatic-planarity skipped: %s", _arom_exc)
        except Exception:
            pass
        return results


def _apply_arom_planarize_if_enabled(mol, results, dual_parse_done: bool):
    """Iter-33 (2026-06-19) dispatch — universal aromatic ring-SYSTEM
    planarisation (heteroaromatic + fused/polycyclic + coordinated).

    The Iter-24 ``_aromatic_ring_flattener`` deliberately skips coordinated
    rings and only fuses pendant rings, so M-bound heteroaromatic donor rings
    and partly-coordinated polycyclic chelates (the worst pucker offenders,
    e.g. AXAGOY's C5N donor rings @ 0.13-0.16 Å OOP) stay bent.  This pass
    flattens EVERY aromatic ring-system as one rigid unit onto its best-fit
    plane while preserving the M-D invariant: coordinated ring atoms are
    anchored (never moved), the plane is constrained through them, and only
    non-anchor atoms + ring-H + first substituents are projected.

    Default-OFF, byte-identical to the build commit when
    ``DELFIN_FFFREE_AROM_PLANARIZE`` is unset (or 0).  Universal (geometric
    aromatic-ring perception, no SMILES specialisation).  Per-system +
    frame-level never-worse + hard M-D-invariant rollback inside the module.
    """
    if not results:
        return results
    if dual_parse_done:
        return results
    if not _class_conditional_flag(
        "DELFIN_FFFREE_AROM_PLANARIZE", mol, default=0,
    ):
        return results
    try:
        from delfin.manta._arom_planarize import correct_results as _ap_correct
        return _ap_correct(mol, results)
    except Exception as _ap_exc:
        try:
            logger.debug("Iter-33 arom-planarize skipped: %s", _ap_exc)
        except Exception:
            pass
        return results


def _apply_arom_bond_length_if_enabled(mol, results, dual_parse_done: bool):
    """Resonance-aware aromatic bond-length equalisation (eye organic-bond-length
    axis).

    MANTA seats aromatic rings from covalent/distance-geometry priors, leaving the
    ring bonds drifting toward the SINGLE-bond covalent length and/or alternating
    (Kekulé-localised) rather than at the DELOCALISED, mesomerism-equalised length
    — the largest systematic organic-geometry defect on the champion pool.  This
    final pass reshapes each perceived aromatic ring system so every ring bond
    sits at its first-principles delocalised target (Pyykkö single↔double radius
    interpolation at the Hückel benzene fraction f = 2/3: C–C 1.393, C–N 1.333,
    …), preserving angles/planarity (in-plane PBD, metal-coordinated ring atoms
    anchored) and rigidly dragging substituents — LENGTHS ONLY.  Per-frame
    never-worse rollback (ring bond-length deviation must strictly drop, M–D
    invariant + clash guarded).

    Default-OFF, byte-identical to the build commit when
    ``DELFIN_FFFREE_AROM_BOND_LENGTH`` is unset (or 0).  Universal (geometric
    aromatic-ring perception + first-principles covalent-radius target, no SMILES
    specialisation).
    """
    if not results:
        return results
    if dual_parse_done:
        return results
    if not _class_conditional_flag(
        "DELFIN_FFFREE_AROM_BOND_LENGTH", mol, default=0,
    ):
        return results
    try:
        from delfin.manta._arom_bond_length import correct_results as _abl_correct
        return _abl_correct(mol, results)
    except Exception as _abl_exc:
        try:
            logger.debug("arom-bond-length skipped: %s", _abl_exc)
        except Exception:
            pass
        return results


def _apply_pi_coplanar_m_if_enabled(mol, results, dual_parse_done: bool):
    """Iter-34 (2026-06-19) dispatch — coordinated planar π-donor co-planar-M
    orienter (eye-flagged ABIZIW).

    A coordinated, internally-planar aromatic π-donor that binds through an
    in-plane sp2 σ lone pair (pyridine-type ring N, amidinate, …) must lie so
    its ring mean-plane CONTAINS the metal.  Iter-33 ``AROM_PLANARIZE`` flattens
    the ring internally but does NOT rotate the coordinated ring plane through
    the metal, so the ring sits tilted with M out of its plane (ABIZIW C5N donor
    rings: M out-of-plane 0.9-1.4 Å).  This pass rotates the RIGID ring about
    the FIXED donor so the donor's in-plane lone-pair points at M (M brought
    into the ring plane).  Hapto / η π-face donors (which bind perpendicular)
    are excluded by a geometric face-on test.  Donor + metal frozen → M-D bond
    preserved exactly; per-ring never-worse on M-out-of-plane and inter-ligand
    clash; hard M-D invariant rollback.

    Default-OFF, byte-identical to the build commit when
    ``DELFIN_FFFREE_PI_COPLANAR_M`` is unset (or 0).  Universal (geometric
    perception, no SMILES specialisation).
    """
    if not results:
        return results
    if dual_parse_done:
        return results
    if not _class_conditional_flag(
        "DELFIN_FFFREE_PI_COPLANAR_M", mol, default=0,
    ):
        return results
    try:
        from delfin.manta._pi_coplanar_m import correct_results as _pcm_correct
        return _pcm_correct(mol, results)
    except Exception as _pcm_exc:
        try:
            logger.debug("Iter-34 pi-coplanar-M skipped: %s", _pcm_exc)
        except Exception:
            pass
        return results


def _apply_hapto_clearance_if_enabled(mol, results, dual_parse_done: bool):
    """Iter-21 (2026-05-19): Welle-5f-F 81f8a1f-style M-X clearance final post-pass.

    For hapto + multi-hapto class systems, apply the radial M-X clearance
    push from ``_hapto_final_clearance.enforce_m_x_clearance_xyz`` as the
    LAST post-emit step (after B3 angle-corrector + B4 π-H projection),
    bridging the gap between the candidate-select-time gate and the
    actual emitted XYZ.

    Default-ON for hapto+multi_hapto via _class_conditional_flag.  Per
    cross-archive analysis CROSS_ARCHIVE_RERUN_2026_05_18.md: 81f8a1f
    is Champion in cshm_max_max (47.25 vs HEAD 77.66, +30pp gap) via
    piano-stool cone + inline ring + rigid-body multi-metal mechanism.
    Welle-5c-CV (2026-05-16): per-bond intact rate edge +1.34pp + per-
    file +4.33pp on 901-file hapto intersection.

    Universal-fundamental: graph-only metal + hapto-group detection,
    element symbols only, no SMILES regex.

    Bit-exact when DELFIN_5F_F_HAPTO_FINAL_CLEARANCE=0.  Operator override
    via DELFIN_5F_F_HAPTO_FINAL_CLEARANCE_CLASSES=csv overrides class list.
    """
    if not results:
        return results
    if dual_parse_done:
        return results
    if not _class_conditional_flag(
        "DELFIN_5F_F_HAPTO_FINAL_CLEARANCE", mol, default=0,
        default_classes=["hapto", "multi_hapto"],
    ):
        return results
    try:
        from delfin.manta._hapto_final_clearance import enforce_m_x_clearance_xyz
        _find_hg = _find_hapto_groups
    except Exception as _imp_exc:
        try:
            logger.debug("hapto-clearance import failed: %s", _imp_exc)
        except Exception:
            pass
        return results
    # Per-frame application
    new_results: List[Tuple[str, str]] = []
    n_modified = 0
    for (xyz, label) in results:
        try:
            new_xyz = enforce_m_x_clearance_xyz(xyz, mol, _find_hg)
            if new_xyz != xyz:
                n_modified += 1
            new_results.append((new_xyz, label))
        except Exception as _hf_exc:
            try:
                logger.debug("hapto-clearance frame error (kept input): %s", _hf_exc)
            except Exception:
                pass
            new_results.append((xyz, label))
    try:
        if n_modified > 0:
            logger.debug("hapto-clearance modified %d/%d frames", n_modified, len(results))
    except Exception:
        pass
    return new_results


def _apply_baustein5_if_enabled(mol, results, dual_parse_done: bool):
    """Baustein 5 dispatch helper — PBD post-UFF geometry corrector (v2).

    Only-fix-what-is-broken philosophy:
    - Topology hard-gate [0.93, 1.07] × M-D ideal (tightened from v1)
    - M-D drift gate 0.05Å — if optimization drifts good M-D bonds → fallback
    - Conservative defaults (max_iter=10, step=0.1)
    - Phase A.5 Hungarian symmetry projection DISABLED (caused wave8-b5 damage)
    - Per-frame fallback to input on any failure

    Bit-exact when DELFIN_BAUSTEIN5=0 (default).
    """
    if not results:
        return results
    if dual_parse_done:
        return results
    if not _delfin_env_int("DELFIN_BAUSTEIN5", 0):
        return results
    try:
        from delfin.manta._post_optimizer import post_optimize_geometry
        cls = "sigma"
        try:
            cls = _classify_complex_class(mol)
        except Exception:
            pass
        # Patch P-H-TRACK: rigid-H tracking in B5 v2 stages 1/3 (default OFF).
        # When DELFIN_B5_RIGID_H=1, bonded H atoms are dragged with their heavy
        # parent during single-atom moves so C-H / N-H / O-H bonds are not
        # stretched by the heavy-atom-only topology gate.
        #
        # Phase 3B per-class override (analogue of DELFIN_SIGMA_*_CLASSES):
        #   export DELFIN_B5_RIGID_H_CLASSES="sigma,multi_sigma"
        #     → enable rigid-H only for those classes (e.g. when pool-verdict
        #     shows rigid-H schadet hapto/multi_hapto coordination).
        # Empty _CLASSES env (default) → fall back to scalar DELFIN_B5_RIGID_H.
        b5_rigid_h = _class_conditional_flag("DELFIN_B5_RIGID_H", mol)
        new_results: List[Tuple[str, str]] = []
        for (xyz, label) in results:
            try:
                new_xyz, report = post_optimize_geometry(
                    xyz, mol, class_label=cls,
                    rigid_h=b5_rigid_h,
                )
                if report.get("topology_preserved", False):
                    new_results.append((new_xyz, label))
                else:
                    new_results.append((xyz, label))
            except Exception:
                new_results.append((xyz, label))
        return new_results
    except Exception as _b5_exc:
        try:
            logger.debug("Baustein 5 PBD post-optimizer skipped: %s", _b5_exc)
        except Exception:
            pass
        return results


def _apply_baustein6_if_enabled(mol, results, dual_parse_done: bool):
    """Baustein 6 dispatch helper — variational L-BFGS-B refiner + 4-tier symmetry.

    Per-frame: minimises 8-term U_total (bond + angle + clash + topology log-
    barrier + Hungarian coord-sphere + Morgan equivalence + fragment archetype
    + global PG) via SciPy L-BFGS-B with analytic gradient. Topology hard-gate
    inside the refiner enforces no M-D break; on failure the input frame is
    returned unchanged. Runs AFTER B5 in the converter (B5 fixes catastrophic
    moves first, B6 smoothly balances the residual forces).

    Env-flags (Phase 3B per-class override pattern):
      DELFIN_B6_WIRED=0             (default OFF — bit-exact when disabled)
      DELFIN_BAUSTEIN6=0            (alias for DELFIN_B6_WIRED, checked second)
      DELFIN_BAUSTEIN6_CLASSES=     (optional class allow-list, e.g. "sigma,no_metal")

    Determinism: L-BFGS-B + Hungarian + analytic gradients are all
    deterministic; no random seeds are consumed.
    """
    if not results:
        return results
    if dual_parse_done:
        return results
    if (not _delfin_env_int("DELFIN_B6_WIRED", 0)
            and not _class_conditional_flag("DELFIN_BAUSTEIN6", mol, default=0)):
        return results
    try:
        from delfin.manta._variational_refiner import variational_refine
        cls = "sigma"
        try:
            cls = _classify_complex_class(mol)
        except Exception:
            pass
        new_results: List[Tuple[str, str]] = []
        for (xyz, label) in results:
            try:
                new_xyz, report = variational_refine(
                    xyz, mol, class_label=cls,
                    max_iter=200,
                )
                if (report.get("topology_preserved", False)
                        and not report.get("fallback_used", True)):
                    new_results.append((new_xyz, label))
                else:
                    # VISIBILITY (2026-07-30).  B6 refuses silently: every bail-out path in the
                    # refiner returns the input XYZ with fallback_used=True, and the caller logged the
                    # reason at DEBUG only -- so `DELFIN_B6_WIRED=1` measured affected=0 on a 35-system
                    # probe and there was no way to tell WHY.  An 8-term functional that declines
                    # every frame and says nothing cannot be developed.  Only reachable when B6 is
                    # armed, so default behaviour is untouched.
                    import sys as _s6
                    print(f"[B6] declined {label}: err={report.get('error')!r} "
                          f"topo_in={report.get('topo_ok_input')} "
                          f"topo={report.get('topology_preserved')} "
                          f"fallback={report.get('fallback_used')} "
                          f"E={report.get('energy_initial')}->{report.get('energy_final')} "
                          f"it={report.get('iterations')} "
                          f"pg={report.get('global_pg')}", file=_s6.stderr)
                    new_results.append((xyz, label))
            except Exception as _b6_frame_exc:
                import sys as _s6
                print(f"[B6] raised on {label}: {_b6_frame_exc!r}", file=_s6.stderr)
                new_results.append((xyz, label))
        return new_results
    except Exception as _b6_exc:
        # SAME VISIBILITY HOLE one level up (2026-07-30): an ImportError here (scipy absent,
        # symbol renamed) silently returns the untouched frames, so the per-frame print below
        # never runs and the probe again reads as "affected=0" -- indistinguishable from "the
        # functional declined every frame".  Only reachable when B6 is armed.
        import sys as _s6
        print(f"[B6] dispatch failed entirely: {_b6_exc!r}", file=_s6.stderr)
        try:
            logger.debug("Baustein 6 variational refine skipped: %s", _b6_exc)
        except Exception:
            pass
        return results


def _apply_fixer_f19_if_enabled(mol, results, dual_parse_done: bool):
    """F19 sp3-H tetrahedrality fixer dispatch helper.

    Per-frame surgical post-pass: detect sp3 heavy atoms whose attached H
    atoms participate in (X-A-H) or (H-A-H) angles deviating from the ideal
    tetrahedral 109.5° by more than ``DELFIN_FIX_F19_TOL_DEG`` (default 10°),
    then repair by in-place rotation of ONLY the offending H atoms (never
    heavy atoms, never metals).  A-H bond lengths are preserved by
    construction; topology is unchanged.

    Env-flags:
        DELFIN_FIX_F19=0           (default OFF — bit-exact when disabled)
        DELFIN_FIX_F19_TOL_DEG=10  (tolerance in degrees)
        DELFIN_FIX_F19_CLASSES=    (optional class allow-list, Phase 3B pattern)

    Insertion order: AFTER B5 PBD post-optimizer.  B5 may slightly shift
    heavy atoms via Stage 1/2; F19 then re-aligns H atoms to the post-B5
    heavy frame.  Per-frame fallback to input on any failure.
    """
    if not results:
        return results
    if dual_parse_done:
        return results
    if not _class_conditional_flag("DELFIN_FIX_F19", mol, default=0):
        return results
    if mol is None:
        return results
    try:
        from delfin.manta._fix_sp3_h_tetrahedrality import fix_sp3_h_tetrahedrality
        tol = _delfin_env_float("DELFIN_FIX_F19_TOL_DEG", 10.0)
        new_results: List[Tuple[str, str]] = []
        for (xyz, label) in results:
            try:
                new_xyz, report = fix_sp3_h_tetrahedrality(
                    xyz, mol, tolerance_deg=tol,
                )
                # Topology guaranteed True by construction; defensive check
                # falls back to input on any unexpected failure.
                if report.get("topology_preserved", True):
                    new_results.append((new_xyz, label))
                else:
                    new_results.append((xyz, label))
            except Exception:
                new_results.append((xyz, label))
        return new_results
    except Exception as _f19_exc:
        try:
            logger.debug("Fixer F19 (sp3-H tetrahedrality) skipped: %s",
                         _f19_exc)
        except Exception:
            pass
        return results


def _apply_hydroxyl_geom_if_enabled(mol, results, dual_parse_done: bool):
    """#329 pendant-hydroxyl C-O-H angle fixer dispatch helper.

    Per-frame surgical post-pass: detect pendant (non-coordinating) hydroxyl O
    (bonded to exactly one H and one C, no metal contact) whose C-O-H angle is
    too wide (>115°) or too narrow (<100°), and rotate ONLY that H about O to
    108.5° (O-H length preserved; 1-3 H···C contact correctly excluded from the
    clash test).  Never moves heavy atoms or metals; per-hydroxyl clash
    rollback; topology unchanged by construction.

    Data (2026-06-23, V2R): C-O-H >115° 46%→2%, median 113→108.5°, ideal band
    33%→80%, real (1-3-aware) clashes 170→167 (none new), deterministic.  The
    C-O length leg is default-OFF (geometric hybridisation classification
    unreliable; much of the "C-O too short" signal was correct phenols at 1.36).

    Env-flag:
        DELFIN_FFFREE_HYDROXYL_GEOM=0   (default OFF — bit-exact when disabled)
    """
    if not results:
        return results
    if dual_parse_done:
        return results
    if os.environ.get("DELFIN_FFFREE_HYDROXYL_GEOM", "0") != "1":
        return results
    try:
        from delfin.manta._fix_hydroxyl_geometry import fix_hydroxyl_geometry
        new_results: List[Tuple[str, str]] = []
        for (xyz, label) in results:
            try:
                new_xyz, report = fix_hydroxyl_geometry(xyz, mol)
                if report.get("topology_preserved", True):
                    new_results.append((new_xyz, label))
                else:
                    new_results.append((xyz, label))
            except Exception:
                new_results.append((xyz, label))
        return new_results
    except Exception as _hyd_exc:
        try:
            logger.debug("Hydroxyl-geom fixer (#329) skipped: %s", _hyd_exc)
        except Exception:
            pass
        return results


def _apply_f19_to_fallback_xyz(xyz_content, mol):
    """Stream-B Fix 2 (batch-2 = FULL): repair sp3-H tetrahedrality on a RAW
    XYZ produced by ``_smiles_to_xyz_unsanitized_fallback`` (and the sigma
    return paths that feed off it).

    The unsanitized-fallback path strips H before embed, then ``AddHs(addCoords=
    True)`` places sp3 H octahedrally (90/180° "square-planar CH3"), and the
    embed/AddHs geometry is often degenerate (overlapping geminal H, so a
    purely geometric H-count over-counts and the centre is skipped).  The main
    sanitized path runs ``_apply_fixer_f19_if_enabled``, but the fallback path
    bypasses it entirely, so XULPUP-class complexes keep their broken methyls.

    Why batch-2 is FULL (root-of-partial fixed): the legacy per-angle rotation
    (``fix_sp3_h_tetrahedrality``) corrected only ONE (X-A-H) angle per H and
    never reconstructed a joint umbrella, and it derived H-count from the
    degenerate embedded geometry.  Batch-2 instead

      (1) ``fix_sp3_h_tetrahedrality_full`` — reads connectivity from the
          RDKit ``mol`` bond graph (true per-centre H-count, so a methylene C
          is known to carry exactly two H even when AddHs overlapped them) and
          places EVERY sp3 H jointly onto the ideal tetrahedral / pyramidal /
          bent sites (C/N/O/P/S/Si/B incl. protonated/charged centres), then
      (2) ``_h_vsepr_realism.correct_xyz`` — a geometry-based finishing pass
          (now seeing a clean, non-degenerate frame) plus the 5f-D inter-
          substituent rotamer relief that staggers crowded multi-methyl
          centres (neopentane / NMe3 / PMe3 / tBu-).

    Both stages are gated by the SAME ``DELFIN_FIX_F19`` flag (default 0 →
    byte-identical: nothing runs, raw ``xyz_content`` returned untouched).
    Topology-safe (only H move, A-H lengths preserved, metals never touched),
    deterministic (fixed tables/probes, sorted iteration), never non-finite
    (NaN guards + per-H rollback).  The mol's atom ordering matches the XYZ
    (both come from the same post-AddHs ``mol``), which the index-based stage-1
    relies on; stage-2 is index-free (pure geometry).

    Returns the (possibly repaired) XYZ string; falls back to the input on any
    failure.
    """
    if xyz_content is None or mol is None:
        return xyz_content
    try:
        # Cheap gate first → byte-identical OFF (no hybridization mutation, no
        # fixer call) when DELFIN_FIX_F19 is unset/0.
        if not _class_conditional_flag("DELFIN_FIX_F19", mol, default=0):
            return xyz_content
        # The unsanitized fallback leaves hybridization UNSPECIFIED; perceive it
        # so the sp3 detector can see the methyls.  Done on a copy so the
        # caller's mol is never mutated.  Topology-safe (perception only).
        try:
            from rdkit import Chem as _Chem
            mol_h = _Chem.Mol(mol)
            try:
                _Chem.GetSymmSSSR(mol_h)
            except Exception:
                pass
            _Chem.SetHybridization(mol_h)
        except Exception:
            mol_h = mol
        # Stage 1 — graph-driven full umbrella reconstruction.
        result = xyz_content
        try:
            from delfin.manta._fix_sp3_h_tetrahedrality import (
                fix_sp3_h_tetrahedrality_full,
            )
            tol = _delfin_env_float("DELFIN_FIX_F19_TOL_DEG", 10.0)
            new_xyz, _rep = fix_sp3_h_tetrahedrality_full(
                result, mol_h, tolerance_deg=tol,
            )
            if isinstance(new_xyz, str) and new_xyz:
                result = new_xyz
        except Exception as _stage1_exc:
            try:
                logger.debug("F19 stage-1 (full umbrella) skipped: %s",
                             _stage1_exc)
            except Exception:
                pass
        # Stage 2 — VSEPR geometry finishing + inter-substituent rotamer relief
        # (5f-D).  Enable 5f-D only for THIS call; restore the prior env value
        # deterministically afterwards so we never leak global state.
        try:
            from delfin.manta import _h_vsepr_realism as _vsepr
            _prev_5fd = os.environ.get("DELFIN_5F_D_ALKYL_ROTAMER")
            os.environ["DELFIN_5F_D_ALKYL_ROTAMER"] = "1"
            try:
                new_xyz = _vsepr.correct_xyz(result)
            finally:
                if _prev_5fd is None:
                    os.environ.pop("DELFIN_5F_D_ALKYL_ROTAMER", None)
                else:
                    os.environ["DELFIN_5F_D_ALKYL_ROTAMER"] = _prev_5fd
            if isinstance(new_xyz, str) and new_xyz:
                result = new_xyz
        except Exception as _stage2_exc:
            try:
                logger.debug("F19 stage-2 (VSEPR finish) skipped: %s",
                             _stage2_exc)
            except Exception:
                pass
        return result if (isinstance(result, str) and result) else xyz_content
    except Exception as _f19fb_exc:
        try:
            logger.debug("Fixer F19 on fallback XYZ skipped: %s", _f19fb_exc)
        except Exception:
            pass
        return xyz_content


def _apply_fixer_f25_if_enabled(mol, results, dual_parse_done: bool):
    """F25 sp3-N pyramidality fixer dispatch helper.

    Per-frame surgical post-pass: detect sp3 N atoms whose 3-neighbour-angle
    sum exceeds ``DELFIN_FIX_F25_THRESHOLD_DEG`` (default 348°, i.e. less
    than 12° pyramidalization), then bend H substituents out of the local
    trigonal plane to recover a target pyramidal apex
    (``DELFIN_FIX_F25_TARGET_DEG``, default 328°).  Only H atoms move;
    heavy chain stays rigid.  Skips metal-bonded N and amide-N.

    Env-flags:
        DELFIN_FIX_F25=0                (default OFF — bit-exact when disabled)
        DELFIN_FIX_F25_THRESHOLD_DEG=348
        DELFIN_FIX_F25_TARGET_DEG=328
        DELFIN_FIX_F25_CLASSES=         (optional class allow-list)

    Insertion order: AFTER F19.  F19 aligns local H angles; F25 corrects
    larger-scale planar-sp3-N flatness.  Per-frame fallback to input on
    any failure.
    """
    if not results:
        return results
    if dual_parse_done:
        return results
    if not _class_conditional_flag("DELFIN_FIX_F25", mol, default=0):
        return results
    if mol is None:
        return results
    try:
        from delfin.manta._fix_sp3_n_pyramidality import fix_sp3_n_pyramidality
        thresh = _delfin_env_float("DELFIN_FIX_F25_THRESHOLD_DEG", 348.0)
        target = _delfin_env_float("DELFIN_FIX_F25_TARGET_DEG", 328.0)
        new_results: List[Tuple[str, str]] = []
        for (xyz, label) in results:
            try:
                new_xyz, _report = fix_sp3_n_pyramidality(
                    xyz, mol,
                    planarity_threshold_deg=thresh,
                    target_sum_deg=target,
                )
                new_results.append((new_xyz, label))
            except Exception:
                new_results.append((xyz, label))
        return new_results
    except Exception as _f25_exc:
        try:
            logger.debug("Fixer F25 (sp3-N pyramidality) skipped: %s",
                         _f25_exc)
        except Exception:
            pass
        return results


def _apply_fixer_sp2n_planarize_if_enabled(mol, results, dual_parse_done: bool):
    """SP2N-PLANARIZE sp2-N planarisation fixer dispatch helper.

    Per-frame surgical post-pass: detect acyclic sp2 nitrogen that the
    assembly/relax has pyramidalised + desymmetrised — primarily the broken
    C-nitro group (N bonded to exactly 2 O + 1 C, total degree 3) where the
    AVIDAM-class defect appears (N pushed ~1 Å out of the C-O-O plane, N-O
    desymmetrised to ~1.43/1.12 Å, O-N-O collapsed to ~107°).  Rebuilds a
    symmetric **planar** NO2 (N + 2 O + C coplanar, N-O = ``no_target`` Å,
    O-N-O = ``ono_target`` °) while keeping the C skeleton + C-N bond rigid.

    Optional ``DELFIN_FIX_SP2N_AMIDE_IMINE=1`` extends to acyclic amide/imine
    N (geometry-only plane projection; no bond retarget).

    Env-flags:
        DELFIN_FFFREE_SP2N_PLANARIZE=0    (default OFF — bit-exact when disabled)
        DELFIN_FIX_SP2N_OOP_DEG_A=0.20    (out-of-plane trigger, Å)
        DELFIN_FIX_SP2N_NO_A=1.22         (target N-O bond length, Å)
        DELFIN_FIX_SP2N_ONO_DEG=125.0     (target O-N-O angle, °)
        DELFIN_FIX_SP2N_AMIDE_IMINE=0     (extend to amide/imine N)
        DELFIN_FFFREE_SP2N_PLANARIZE_CLASSES=  (optional class allow-list)

    Insertion order: AFTER F25 (sp3-N pyramidality).  F25 pyramidalises
    over-flat sp3 N; this is the inverse (it flattens pyramidalised sp2 N) —
    they touch disjoint N (RDKit sp3 vs sp2 / no-H vs H), so order is benign.
    Skips metal-bonded N; per-group rollback on new clash; per-frame fallback
    to input on any failure.  Bit-exact when the flag is 0.
    """
    if not results:
        return results
    if dual_parse_done:
        return results
    if not _class_conditional_flag("DELFIN_FFFREE_SP2N_PLANARIZE", mol,
                                   default=0):
        return results
    if mol is None:
        return results
    try:
        from delfin.manta._fix_sp2n_planarize import planarize_sp2_nitrogen
        oop_thr = _delfin_env_float("DELFIN_FIX_SP2N_OOP_DEG_A", 0.20)
        no_t = _delfin_env_float("DELFIN_FIX_SP2N_NO_A", 1.22)
        ono_t = _delfin_env_float("DELFIN_FIX_SP2N_ONO_DEG", 125.0)
        amide_imine = (os.environ.get("DELFIN_FIX_SP2N_AMIDE_IMINE", "0")
                       == "1")
        new_results: List[Tuple[str, str]] = []
        for (xyz, label) in results:
            try:
                new_xyz, _report = planarize_sp2_nitrogen(
                    xyz, mol,
                    oop_threshold_A=oop_thr,
                    no_target_A=no_t,
                    ono_target_deg=ono_t,
                    include_amide_imine=amide_imine,
                )
                new_results.append((new_xyz, label))
            except Exception:
                new_results.append((xyz, label))
        return new_results
    except Exception as _sp2n_exc:
        try:
            logger.debug("Fixer SP2N-PLANARIZE skipped: %s", _sp2n_exc)
        except Exception:
            pass
        return results


def _apply_fixer_sp2c_planarize_if_enabled(mol, results, dual_parse_done: bool):
    """SP2C-PLANARIZE sp2-CARBON planarisation fixer dispatch helper.

    Per-frame surgical post-pass: detect an acyclic sp2 CARBON (azomethine /
    imine / vinyl / enone =CH- or =CR-: 3-coordinate, RDKit-sp2, a double bond to
    N/C/O) that the force-field-free ETKDG embed has left PYRAMIDALISED (out of
    the plane of its three neighbours), and project it back into that plane
    (geometry-only; bond lengths shift only by the small residual).  Closes the
    exact gap the N-only ``_fix_sp2n_planarize`` left open: a backbone N=CH-C
    carbon (NAYKOQ C36: Walsh ~19°, angle-sum 324° = the sp3 fallback, in 17/30
    ETKDG folds) was planarised by NOTHING.  Ring C is skipped (owned by the
    aromatic ring passes); metal-bonded C is skipped (obeys coordination).

    Env-flags:
        DELFIN_FFFREE_SP2C_PLANARIZE=0    (default OFF — bit-exact when disabled;
                                           champion-ON via _CHAMPION_FLAGS)
        DELFIN_FIX_SP2C_OOP_A=0.20        (out-of-plane trigger, Å)
        DELFIN_FFFREE_SP2C_PLANARIZE_CLASSES=  (optional class allow-list)

    Insertion order: AFTER SP2N-PLANARIZE (they touch disjoint atoms: N vs C).
    Per-atom rollback on new clash / no-improvement; per-frame fallback to input
    on any failure.  Bit-exact when the flag is 0.
    """
    if not results:
        return results
    if dual_parse_done:
        return results
    if not _class_conditional_flag("DELFIN_FFFREE_SP2C_PLANARIZE", mol,
                                   default=0):
        return results
    if mol is None:
        return results
    try:
        from delfin.manta._fix_sp2c_planarize import planarize_sp2_carbon
        oop_thr = _delfin_env_float("DELFIN_FIX_SP2C_OOP_A", 0.20)
        # TOPO-GATE (DELFIN_FFFREE_SP2C_PLANARIZE_TOPOGATE): re-run construction's
        # OWN broken-frame gate (_fragment_topology_ok) on the PLANARISED frame and
        # roll the WHOLE frame back if the projection changed the connectivity.  The
        # fixers run AFTER _clean_gate_filter, so without this a projection that tears
        # a bond ships a broken frame (AGEPAH broken_regressed).  Re-gating culls the
        # damage instead of shipping it -- never-worse w.r.t. topology.
        _topo_gate = bool(_delfin_env_int(
            "DELFIN_FFFREE_SP2C_PLANARIZE_TOPOGATE", 0))
        _smi = None
        if _topo_gate:
            try:
                from rdkit import Chem as _Chem
                _smi = _Chem.MolToSmiles(mol)
            except Exception:
                _topo_gate = False
        new_results: List[Tuple[str, str]] = []
        for (xyz, label) in results:
            try:
                new_xyz, _report = planarize_sp2_carbon(
                    xyz, mol, oop_threshold_A=oop_thr)
                if (_topo_gate and _smi and new_xyz != xyz
                        and not _fragment_topology_ok(new_xyz, _smi)):
                    new_xyz = xyz          # rollback: fix broke the topology
                new_results.append((new_xyz, label))
            except Exception:
                new_results.append((xyz, label))
        return new_results
    except Exception as _sp2c_exc:
        try:
            logger.debug("Fixer SP2C-PLANARIZE skipped: %s", _sp2c_exc)
        except Exception:
            pass
        return results


def _apply_fixer_wuxqak_if_enabled(mol, results, dual_parse_done: bool):
    """WUXQAK sp3-C linear-collapse fixer dispatch helper.

    Per-frame surgical post-pass for M-CH2-X / M-CHR-R' patterns where an
    sp3 C bonded to a metal collapses to a near-linear M-C-X angle
    (>``DELFIN_FIX_WUXQAK_ANGLE_DEG``, default 150°).  Rotates only the
    rigid BFS sub-fragment of X around C (M and C fixed) to recover
    ``DELFIN_FIX_WUXQAK_TARGET_DEG`` (default 109.5°).  Per-violation
    rollback if new clash or any intact M-D bond drifts > 0.05 Å.

    Env-flags:
        DELFIN_FIX_WUXQAK=0              (default OFF — bit-exact when disabled)
        DELFIN_FIX_WUXQAK_ANGLE_DEG=150
        DELFIN_FIX_WUXQAK_TARGET_DEG=109.5
        DELFIN_FIX_WUXQAK_CLASSES=        (optional class allow-list)

    Skips frames without a metal atom (fast pre-check).  Per-frame fallback
    to input on any failure or on topology-broken result.
    """
    if not results:
        return results
    if dual_parse_done:
        return results
    if not _class_conditional_flag("DELFIN_FIX_WUXQAK", mol, default=0):
        return results
    if mol is None:
        return results
    # Cheap pre-check: skip entirely if no result XYZ contains a metal.
    try:
        any_metal = False
        for entry in results:
            xyz_text = entry[0] if entry else ""
            if not isinstance(xyz_text, str):
                continue
            for ln in xyz_text.splitlines():
                tok = ln.strip().split(" ", 1)[0] if ln.strip() else ""
                if tok and tok[0].isalpha() and len(tok) <= 2:
                    if tok not in ("H", "C", "N", "O", "F", "P", "S",
                                   "Cl", "Br", "I", "B", "Si"):
                        any_metal = True
                        break
            if any_metal:
                break
        if not any_metal:
            return results
    except Exception:
        pass
    try:
        from delfin.manta._fix_wuxqak_sp3_c_linear import fix_wuxqak_sp3_c_linear
        angle_thr = _delfin_env_float("DELFIN_FIX_WUXQAK_ANGLE_DEG", 150.0)
        target = _delfin_env_float("DELFIN_FIX_WUXQAK_TARGET_DEG", 109.5)
        new_results: List[Tuple[str, str]] = []
        for (xyz, label) in results:
            try:
                new_xyz, report = fix_wuxqak_sp3_c_linear(
                    xyz, mol,
                    angle_threshold_deg=angle_thr,
                    target_angle_deg=target,
                )
                if report.get("topology_preserved", True):
                    new_results.append((new_xyz, label))
                else:
                    new_results.append((xyz, label))
            except Exception:
                new_results.append((xyz, label))
        return new_results
    except Exception as _wux_exc:
        try:
            logger.debug("Fixer WUXQAK (sp3-C linear) skipped: %s", _wux_exc)
        except Exception:
            pass
        return results


def _apply_fixer_bridging_anion_if_enabled(mol, results, dual_parse_done: bool):
    """μ-X bridging-anion M-X-M angle fixer dispatch helper.

    Per-frame surgical post-pass for bimetallic complexes with bridging
    anionic donors X (Cl, Br, OH, OR, NR2, SR, CR3 — universal, graph-
    driven detection).  Detects μ-X atoms via the mol bond graph (any
    non-metal atom bonded to ≥2 metals).  For each, infers the chemistry-
    realistic M-X-M angle window from the local topology:

      * 4-ring motif ([M2(μ-X)2] dimer)             → target 92°  (window 80-105°)
      * sp2/oxo motif (no non-metal substituents)    → target 145° (window 125-175°)
      * sp3 single-bridge (μ-OH/μ-OR/μ-NR2/μ-CR3)    → target 110° (window 95-130°)

    When the current M-X-M angle falls outside its window, the bridging
    X (together with its non-metal substituent BFS fragment) is
    translated onto the perpendicular bisector of the M-M segment by an
    amount that places the new M-X-M angle on the target.  Both metals
    stay fixed; |M-X| is preserved by construction.  Per-violation
    rollback if a new heavy-atom clash appears or any intact M-D bond
    drifts > 0.05 Å.

    Env-flags:
        DELFIN_FIX_BRIDGING_ANION=0       (default OFF — bit-exact when disabled)
        DELFIN_FIX_BRIDGING_ANION_CLASSES=    (optional class allow-list, see
                                             ``_class_conditional_flag``)

    Insertion order: AFTER F19/F25/WUXQAK so that all upstream single-
    metal donor / sp3-N / sp3-C fixers have settled first.  Per-frame
    fallback to input on any failure or on topology-broken result.
    """
    if not results:
        return results
    if dual_parse_done:
        return results
    if not _class_conditional_flag("DELFIN_FIX_BRIDGING_ANION", mol, default=0):
        return results
    if mol is None:
        return results
    # Cheap pre-check: skip frames without a metal symbol.  Bridging
    # anions are bimetallic-only, so an additional cheap test would be
    # "≥2 metal atoms in the frame", but the per-frame fixer already
    # no-ops if no μ-X bridge is detected, so the simple metal-presence
    # gate keeps the fast-path identical to F25/WUXQAK.
    try:
        any_metal = False
        for entry in results:
            xyz_text = entry[0] if entry else ""
            if not isinstance(xyz_text, str):
                continue
            for ln in xyz_text.splitlines():
                tok = ln.strip().split(" ", 1)[0] if ln.strip() else ""
                if tok and tok[0].isalpha() and len(tok) <= 2:
                    if tok not in ("H", "C", "N", "O", "F", "P", "S",
                                   "Cl", "Br", "I", "B", "Si"):
                        any_metal = True
                        break
            if any_metal:
                break
        if not any_metal:
            return results
    except Exception:
        pass
    try:
        from delfin.manta._fix_bridging_anion import fix_bridging_anion_angles
        new_results: List[Tuple[str, str]] = []
        for (xyz, label) in results:
            try:
                new_xyz, report = fix_bridging_anion_angles(xyz, mol)
                if report.get("topology_preserved", True):
                    new_results.append((new_xyz, label))
                else:
                    new_results.append((xyz, label))
            except Exception:
                new_results.append((xyz, label))
        return new_results
    except Exception as _bri_exc:
        try:
            logger.debug("Fixer bridging-anion (μ-X M-X-M) skipped: %s",
                         _bri_exc)
        except Exception:
            pass
        return results


def _apply_xtb_cascade_if_enabled(mol, results, dual_parse_done: bool):
    """Welle-5l Track-2 dispatch helper — xTB-cascade post-UFF refinement.

    Per-frame GFN2-xTB optimization with M-D-invariant rollback gate.  This
    is the OPT-IN cascade-stage that wraps ``delfin.manta._xtb_refiner.refine_with_xtb``.
    Runs AFTER all UFF-side post-processing (B3 / B4 / B5 / B6 / F19 / F25 /
    WUXQAK / bridging-anion / 5b-A / 5b-B / 5f-C / 5f-D / 5j-A) so the cascade
    refines the FINAL pipeline geometry, not an intermediate one.

    Architecture (per ``project_core_swap_decision`` 2026-05-14):
        scaffold → ETKDG → UFF (+ all bandages) → xTB cascade [optional]

    Env-flags (Phase 3B per-class override pattern, ``_class_conditional_flag``):
        DELFIN_CASCADE_REFINER=0                (default OFF — bit-exact when 0)
        DELFIN_CASCADE_REFINER_CLASSES=         (csv class allow-list, e.g.
                                                 "sigma,multi_sigma")

    Resolution (matches ``_class_conditional_flag`` precedence):
        1. If ``DELFIN_CASCADE_REFINER_CLASSES`` is set → cascade fires iff
           ``_classify_complex_class(mol)`` is in the csv list, regardless of
           the scalar value.  This is the recommended deployment pattern.
        2. Else fall back to scalar ``DELFIN_CASCADE_REFINER`` (0 / 1).

    Per-frame contract:
        - Charge: ``Chem.GetFormalCharge(mol)`` (sum of explicit formal charges).
        - Multiplicity (uhf = n_unpaired_e-): 0 if total electron count is even
          (closed-shell singlet), 1 if odd (doublet).  No spin-state search
          (per project_core_swap_decision xtb spin-state policy 2026-05-14).
        - Subprocess timeout: 30s per frame (matches iter16 config).
        - Rollback gate: bond > 1.50 × Σr_cov in refined geometry → revert
          to pre-cascade XYZ (gate lives inside ``refine_with_xtb`` —
          M-D bonds detected with 1.30× factor in input, break threshold
          1.50×, both relative to ``Σr_cov``).
        - Element-count / atom-order mismatch → revert (xtb output corrupt).
        - Any exception → revert to pre-cascade XYZ (fail-safe).

    Bit-exact when ``DELFIN_CASCADE_REFINER=0`` and ``DELFIN_CASCADE_REFINER_CLASSES``
    is unset (the default).  Skipped on inner dual-parse calls so the heavy-
    atom signature dedup sees consistent UFF coordinates.

    Args:
        mol: the parent RDKit ``Mol`` (used for class classification +
            total charge + electron count); may be ``None`` (skip).
        results: list of ``(xyz_str, label)`` tuples produced by the
            up-pipeline.
        dual_parse_done: ``True`` on the inner dual-parse call → skip.

    Returns:
        Same shape as ``results``.  Frames where xTB succeeded are
        replaced with the refined XYZ; everything else is unchanged.
    """
    if not results:
        return results
    if dual_parse_done:
        return results
    if mol is None:
        return results
    # Class-conditional gate (matches B5 rigid-H / F19 / etc. pattern).
    # Iter-17b 2026-05-18: revert Iter-17a default_classes=["multi_sigma"]
    # per user directive — "the structures must be as good as possible before
    # going to xtb, so cascade comes at the very end!"  Pre-xTB pipeline
    # (scaffold + ETKDG + UFF + bandages) must be optimized first.
    # Cascade is reserved for FINAL phase after all pre-xTB optimizations
    # are stable.  Env-flag opt-in preserved for testing.
    if not _class_conditional_flag("DELFIN_CASCADE_REFINER", mol, default=0):
        return results
    # Metal-only: no point spending xtb seconds on organic-only ligands
    # (UFF is parametrised for them).  Cheap pre-check.
    try:
        any_metal = any(a.GetSymbol() in _METAL_SET for a in mol.GetAtoms())
    except Exception:
        any_metal = True  # fall through on probe error
    if not any_metal:
        return results
    try:
        from delfin.manta._xtb_refiner import refine_with_xtb as _xtb_refine
    except Exception as _imp_exc:
        try:
            logger.debug("xTB cascade import failed: %s", _imp_exc)
        except Exception:
            pass
        return results

    # Total charge from RDKit formal charges (per spin-state policy).
    try:
        total_charge = sum(int(a.GetFormalCharge() or 0) for a in mol.GetAtoms())
    except Exception:
        total_charge = 0

    # Minimum multiplicity from electron-count parity (per project_core_swap_decision
    # xtb spin-state policy 2026-05-14): even → uhf=0 (closed-shell singlet),
    # odd → uhf=1 (doublet, one unpaired electron).  No spin-state scan.
    try:
        n_electrons = (
            sum(int(a.GetAtomicNum()) for a in mol.GetAtoms())
            - total_charge
        )
        uhf = int(n_electrons & 1)
    except Exception:
        uhf = 0

    timeout_s = _delfin_env_float("DELFIN_CASCADE_REFINER_TIMEOUT_S", 30.0)
    max_iter = _delfin_env_int("DELFIN_CASCADE_REFINER_MAX_ITER", 200)
    gfn = _delfin_env_int("DELFIN_CASCADE_REFINER_GFN", 2)

    n_refined = 0
    n_reverted = 0
    new_results: List[Tuple[str, str]] = []
    for (xyz, label) in results:
        try:
            refined = _xtb_refine(
                xyz,
                charge=total_charge,
                gfn=gfn,
                max_iter=max_iter,
                timeout_s=timeout_s,
                uhf=uhf,
            )
        except Exception as _rf_exc:
            try:
                logger.debug("xTB cascade frame error (kept input): %s", _rf_exc)
            except Exception:
                pass
            refined = xyz
        if refined is xyz or refined == xyz:
            n_reverted += 1
            new_results.append((xyz, label))
        else:
            n_refined += 1
            new_results.append((refined, label))
    try:
        logger.debug(
            "xTB cascade: %d refined / %d reverted (n_total=%d, charge=%d, uhf=%d)",
            n_refined, n_reverted, len(results), total_charge, uhf,
        )
    except Exception:
        pass
    return new_results
