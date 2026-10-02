"""Switches, tunables, timeouts, seed schedule and quality profiles of the MANTA constructor (every module-level env read stays at import time, every helper reads at call time).

Moved verbatim from delfin/smiles_converter.py (MANTA split, 2026-10); every
statement is the original text, only the imports below are new.
"""
from __future__ import annotations

import hashlib
import os
import random
import threading
from typing import Dict, List, Optional, Tuple

from delfin.common.logging import get_logger
from delfin.manta.hapto_detect import (
    _classify_complex_class,
)
from delfin.manta.ml_tables import (
    _ALL_POLYHEDRA_BY_CN,
    _PREFERRED_CN4_GEOMETRY,
)

logger = get_logger("delfin.smiles_converter")


def _preferred_cn4_for(metal_symbol: str, mol=None, metal_idx=None) -> str:
    """Preferred CN4 geometry -- with d-count, when it can be derived.

    THE MEASURED DEFECT (16.08.2026).  The table above lists `'Cu': 'TH'` -- tetrahedron.
    That is right for **Cu(I) d10** and wrong for **Cu(II) d9**, which is square-planar
    up to JT-elongated octahedral -- and Cu(II) is by far the more frequent case.  The same
    row for both oxidation states.  (Likewise `'Au': 'SQ'`: right for Au(III) d8, wrong for
    Au(I) d10, which is LINEAR -- that stays untouched here for now, because CN2 is a
    different site.)

    The builder already thinks in d-counts at this point -- the comments of the table say
    verbatim "d8 metals", "d6 low-spin", "d0-d5, d7, d10".  Only the result stands there as
    a hard-wired ELEMENT LIST, and an element has several oxidation states.

    d8  -> SQ (square-planar, ligand-field stabilisation)
    d9  -> SQ (JT-elongated; the four short bonds form the plane)
    d10 -> TH (no LFSE, sterically determined)

    ⚠ ONLY when the oxidation state is CERTAIN.  `oxidation_state` returns `None` as soon
    as a donor is not classifiable or more than one metal is present -- then it stays with
    the element table, i.e. with today's behaviour.  Default OFF -> byte-identical.
    """
    base = _PREFERRED_CN4_GEOMETRY.get(metal_symbol, 'SQ')
    if mol is None or metal_idx is None:
        return base
    if not _delfin_env_int("DELFIN_FFFREE_DN_GEOMETRY", 0):
        return base
    try:
        from delfin.manta._oxidation_state import oxidation_state as _ox
        r = _ox(mol, int(metal_idx))
    except Exception:
        return base
    if r is None:
        return base
    _os, d = r
    if d in (8, 9):
        return 'SQ'
    if d == 10:
        return 'TH'
    return base


def _all_polyhedra_codes(n_coord: int, metal_symbol: str,
                         legacy_codes: List[str]) -> List[str]:
    """Iter-2 Subagent A - polyhedron-completeness dispatch.

    When DELFIN_ALL_POLYHEDRA=1 (default ON): return the complete polyhedra
    set for ``n_coord`` with the metal's preferred polyhedron first
    (stable order, additive vs ``legacy_codes``).
    When DELFIN_ALL_POLYHEDRA=0: return ``legacy_codes`` unchanged.

    Coordination-isomer completeness contract: DELFIN must enumerate ALL
    polyhedra possible for a given coordination number, not pick one
    preferred candidate via a metal-identity heuristic.  Picking by
    metal-identity hides isomers the chemistry actually permits and that
    DFT-ranking downstream needs to discriminate.
    """
    # Iter-5: default flipped 1→0 (was Iter-2A default-ON, but pumps low-fidelity
    # polyhedra into pool — see iter5_polyhedron_md_forensik.md).  Opt-in via
    # DELFIN_ALL_POLYHEDRA=1.
    if not _delfin_env_int("DELFIN_ALL_POLYHEDRA", 0):
        return list(legacy_codes)
    full = _ALL_POLYHEDRA_BY_CN.get(int(n_coord))
    if not full:
        return list(legacy_codes)
    seen: set = set()
    ordered: List[str] = []
    for code in list(legacy_codes) + list(full):
        if code in seen:
            continue
        seen.add(code)
        ordered.append(code)
    return ordered


# ---------------------------------------------------------------------------
# Timeout-guarded embedding
# ---------------------------------------------------------------------------
# ============================================================================
# DELFIN SMILES→XYZ pipeline — user-tunable configuration
# ----------------------------------------------------------------------------
# Every knob below has a sensible default, a short docstring, and an
# environment-variable override of the form ``DELFIN_<NAME>`` so the end user
# can tune the pipeline without editing code.  Defaults target
# DFT-rankable output on a 32-core workstation.
#
# Quality profiles
#     ``max``    (default) — thickest pool of candidates, slowest, best
#                geometry quality.  40 ETKDG seeds, 3 chelate conformer
#                ranks, 3 ranked templates, pre-UFF cap 5·max_isomers.
#     ``normal`` — middle ground for the dashboard.  20 seeds, 2 ranks,
#                2 templates, cap 4·max_isomers.
#     ``fast``   — quickest feedback, smallest variety.  12 seeds, 1
#                rank, 1 template, cap 3·max_isomers.
# The profile applies to ``smiles_to_xyz_isomers(quality_mode=...)``; passing
# ``None`` falls back to the module defaults below.
# ============================================================================


def _delfin_env_int(name: str, default: int) -> int:
    try:
        return int(os.environ.get(name, str(default)))
    except Exception:
        return default


def _deterministic_mode() -> bool:
    """Master determinism switch (``DELFIN_DETERMINISTIC=1``, default OFF).

    When set, the structure generator is forced bit-identical CROSS-PROCESS
    (independent of ``PYTHONHASHSEED``, CPU load and wall-clock timing).  The
    switch IMPLIES: deterministic OpenBabel/embed paths, sorted enumeration
    (``DELFIN_FFFREE_DETERMINISTIC_ENUM=1``), zero multinuclear wall budgets,
    and disabled embed timeouts (the deterministic path is always taken).
    Unset/0 → byte-identical to legacy behaviour.
    """
    return os.environ.get("DELFIN_DETERMINISTIC", "0") == "1"


def _hapto_scaffold_primary_enabled() -> bool:
    """Hapto scaffold-primary routing (DELFIN_FFFREE_HAPTO_SCAFFOLD_PRIMARY=1,
    default OFF -> byte-id).  For η-coordinated complexes the FF-free RIGID_HAPTO
    path collapses the η-face onto the metal (measured η6 ~27% geometry-clean,
    M-centroid ~0.66-1.18 A vs the ~1.6 A target).  The legacy analytical hapto
    scaffold (smiles_to_xyz hapto branch) places the η-ring correctly (measured
    η6 ~84% clean, η4 100%).  With this flag ON the FF-free builder does NOT
    short-circuit a hapto complex: the scaffold path is tried FIRST (correct
    geometry) and the FF-free RIGID result is kept only as a FALLBACK for the
    cases the scaffold cannot build (never-worse on build-rate, strictly better
    geometry).  Confirmed by the 2026-05-20 per-group collapse map (η6-arene/
    sandwich build well; RIGID_HAPTO was a coverage-grab that regressed geometry)."""
    return os.environ.get("DELFIN_FFFREE_HAPTO_SCAFFOLD_PRIMARY", "0") == "1"


def _union_prepend_ffree(results_hapto, ffree_union):
    """DELFIN_FFFREE_UNION_HAPTO -- default 0, i.e. byte-identical OFF.

    THE CONTRACT BREACH THIS SWITCH CLOSES.  `DELFIN_FFFREE_UNION` promises
    *ADD, never replace*: the FF-free builder no longer returns early, its
    frames are parked in ``_ffree_union`` and PREPENDED at the very end of the
    legacy pipeline.  The hapto branch, however, returns ~2700 lines BEFORE
    this merge.  On every system with η-coordination ``_ffree_union`` is thus
    never read -- UNION does not merge there, it HANDS the system OVER to the
    hapto path, and that one replaces the enumerated manifold with a single
    repeated η label.

    MEASURED 28.08.2026 (`harness/union_additiv_zensus.py --hapto`), label
    multiset over 9720 systems of both arms of `unioniso10k`:
        ADDITIVE 9684 (99.6 %) · SWITCHERS 36 (0.4 %)
        WINNERS (topo F->T)  204 -> 204 additive  (100 %)
        LOSERS  (topo T->F)    3 ->   3 switchers (100 %)
        of the 36 switchers, 31 are ENTIRELY η-labelled in the ON arm,
        among them ALL THREE losers: CAZMEX, JACVES, XIDKON.
    Evidence case HACVOZ: OFF `T-4-hapto-iso1 / iso1-r1 / iso2 / iso2-r1 / iso3`
    (enumerated, 5 labels) -> ON `η6-arene` 50x (ONE label).  Same frame
    count, every OFF label gone.  That is not a deletion, that is a
    hand-over -- and exactly that is what the contract is supposed to prevent.

    ⚠️ WHAT THIS SWITCH DOES NOT DO.  It does NOT judge which builder delivers
    the better η geometry; the analytical scaffold is measured better there
    (η6 ~84 % clean against ~27 % for RIGID_HAPTO, see
    `_hapto_scaffold_primary_enabled`).  It only establishes what UNION
    promises: BOTH manifolds, neither replaced.  Frame 0 stays the FF-free
    construction, the complete hapto manifold follows behind it.

    ⚠️ KNOWN GAP, deliberately left open.  If
    `DELFIN_FFFREE_HAPTO_SCAFFOLD_PRIMARY` is on, ``_ffree_union`` stays None
    (the FF-free frames then sit in ``_hapto_ff_fallback``) and this switch
    does nothing.  The flag is NOT in the champion; the combination is
    unmeasured and is not silently repaired along the way here.

    Dedup is done as in the merge itself -- over the exact XYZ text --, so
    that no frame appears twice in the manifold."""
    if not ffree_union or _delfin_env_int("DELFIN_FFFREE_UNION_HAPTO", 0) != 1:
        return results_hapto
    try:
        _seen = {x for x, _l in ffree_union}
        _extra = [(x, l) for x, l in (results_hapto or []) if x not in _seen]
        return list(ffree_union) + _extra
    except Exception:
        return results_hapto


def _hapto_seat_rigid_enabled() -> bool:
    """DELFIN_FFFREE_HAPTO_SEAT_RIGID -- default 0, i.e. byte-identical OFF.

    THE MEASUREMENT (19.08.2026, candidate census on the 42 eta systems with a
    flat stereocentre).  The hapto seed comes from
    ``smiles_to_xyz(hapto_approx=True)``, and there
    ``_select_best_hapto_candidate`` chooses between two build modes:

      ``scaffold``  -- ``_build_hapto_scaffold``: the eta ring is seated as a regular
                       polygon (correct and RIGID), but EVERYTHING else ATOM BY ATOM
                       via BFS-VSEPR out of ONE parent atom, without any clash
                       check (lines 16119-16219).
      ``hybrid*``   -- ``_build_hybrid_hapto_complex``: the same scaffold, but every
                       ligand branch is embedded SEPARATELY (ETKDG + UFF) and afterwards
                       rotated/shifted onto its anchor points as a RIGID BODY
                       (``_embed_hybrid_fragment`` + ``_align_hybrid_fragment_onto_scaffold``).
                       The internal geometry stays untouched in the process.

    Measured (collapsed bonds of the SEED, ``_bond_decollapse`` floor):
        AHIGUW  rigid  1 : atom-wise 21      ALEKIM  rigid  3 : atom-wise 10
        ALEKOS  rigid  2 : atom-wise 12      BEXVAC  rigid  6 : atom-wise 64
        ABEWUZ  rigid 10 : atom-wise 47      BAKLAB  rigid 54 : atom-wise 51
    In 5 of 6 cases the RIGID build carries a multiple LESS collapse -- and still
    loses the selection, because ``_hapto_candidate_quality_score`` has no collapse
    term at all: its only overlap test ``_has_atom_clash`` measures against
    0.80 A, while the collapse floor for C-C sits at 0.82*1.52 = 1.25 A.  The score
    distinguishes the two builds by ~0.4 %, the geometry by three- to
    tenfold.

    WHY SELECTION AND NOT A REBUILD.  The cost law measures ordering/selection at ~0,
    re-embedding at +11.9 pp.  The rigid body is already built; building it a
    second time would be the most expensive class for the same result.  This
    line merely lets it arrive.

    RESULT ON (A/B against the state before the intervention, champion environment):
        SEED, 42 systems:  collapse sum 677 -> 498, median 11 -> 9,
                           11 systems better, 0 worse, 31 untouched.
        FRAMES, 12 systems / 330 frames: collapse sum 5560 -> 2048 (-63 %),
                           6 better, 0 worse.
        OFF: byte-identical 42/42 (seed) and 12/12 (frames).

    ⚠ THE REASON WHY THIS SWITCH MUST NOT LAND LIKE THIS -- THE FRAMES BECOME
    FEWER.  On the same 12 systems the frame count drops from 320 to 230, and
    the lost ones are one and all the ``hd-ta-*`` frames of the OB rotor
    branch (AHIGUW 24 -> 4, EQEYIJ 22 -> 1).  The rigid seating thus delivers
    CLEANER frames, but FEWER -- and completeness is the more sacred
    axis.  The collapse gain is therefore partly a counter artefact (fewer
    frames carry less collapse); PER FRAME the gain nevertheless
    holds (AHIGUW 18.4 -> 5.0; BEXVAC 65.5 -> 10.7; ALEKOS 7.7 -> 1.5).
    WHY the rotor frames drop out is NOT measured -- that is the next
    question, not a settled one.  As long as it is open, the right
    blueprint belongs in the ADDITIVE form: build both seeds and keep BOTH frame
    sets (the UNION pattern), instead of swapping one for the other.
    """
    return _delfin_env_int("DELFIN_FFFREE_HAPTO_SEAT_RIGID", 0) == 1


def _hapto_candidate_collapsed_bonds(cand_mol) -> Optional[int]:
    """Collapsed bonds of a hapto candidate -- ONE source for the floor.

    The floor (``COLLAPSE_FLOOR``, default 0.82) and the geometric bond
    perception come unchanged from ``delfin.manta._bond_decollapse``; NOTHING is
    rebuilt here and no second threshold is invented.  The import sits inside the
    function so that the switched-off path never pays for it.

    Return ``None`` means "not measurable" -- the caller must treat that like "no
    verdict", never like "zero collapse".
    """
    if cand_mol is None:
        return None
    try:
        import numpy as np
        from delfin.manta import _bond_decollapse as _BD
        conf = cand_mol.GetConformer(0)
        n = cand_mol.GetNumAtoms()
        syms = [cand_mol.GetAtomWithIdx(i).GetSymbol() for i in range(n)]
        pts = []
        for i in range(n):
            p = conf.GetAtomPosition(i)
            pts.append([p.x, p.y, p.z])
        P = np.asarray(pts, dtype=float)
        if P.shape[0] == 0:
            return None
        bonds = _BD._geometric_bonds(syms, P)
        return int(_BD._count_collapsed(syms, P, bonds))
    except Exception as _exc:
        logger.debug("Kollapszaehlung am Hapto-Kandidaten fehlgeschlagen: %s", _exc)
        return None


def _apply_uff_jitter(
    xyz_delfin: str,
    atom_indices,
    smiles: Optional[str],
    magnitude: float = 0.05,
) -> str:
    """Apply deterministic hash-seeded micro-jitter to a subset of atoms.

    Iter-8.7 (2026-05-12): the topology-enumerator path places metal +
    monodentate donors at the ideal symmetric ``_TOPO_GEOMETRY_VECTORS``
    positions (Td/Oh/TPR/ICOS/CUBO) and then freezes them during OB-UFF.
    The frozen polyhedron is a stationary point of UFF, so identical
    template-perfect bond lengths (e.g. Zn-Br = 2.3500 Å exact, Ru-N =
    2.060/2.400 Å exact) survive the optimization byte-for-byte and the
    surrounding ligand internals can collapse symmetrically.

    Adding a small position offset (±``magnitude`` Å, default 0.05) to
    each frozen atom BEFORE UFF breaks the symmetric stationary point
    while remaining within the M-D distance tolerance window of the
    downstream ``_verify_metal_connectivity`` check.  Per-coordinate
    offsets are drawn from a SMILES-hash-seeded RNG so the same SMILES
    always produces the same jitter (determinism for cache/dedup).

    Args:
        xyz_delfin: DELFIN-format coordinate block ("symbol x y z" per line,
            no header).
        atom_indices: iterable of 0-based atom indices to perturb.  Atoms
            not in this set keep their original coordinates exactly.
        smiles: SMILES string used as the hash seed.  ``None`` falls back
            to a fixed seed (still deterministic).
        magnitude: half-range of the uniform offset distribution in Å.

    Returns:
        Modified DELFIN-format XYZ string.  Returns ``xyz_delfin``
        unchanged on any parse error.
    """
    try:
        seed_src = (smiles or "").encode("utf-8")
        seed_int = int.from_bytes(
            hashlib.sha256(seed_src).digest()[:8], "big", signed=False
        )
        rng = random.Random(seed_int)

        idx_set = set(int(i) for i in atom_indices)
        if not idx_set:
            return xyz_delfin

        out_lines = []
        for i, raw in enumerate(xyz_delfin.splitlines()):
            if not raw.strip():
                out_lines.append(raw)
                continue
            parts = raw.split()
            if len(parts) < 4:
                out_lines.append(raw)
                continue
            sym = parts[0]
            try:
                x = float(parts[1])
                y = float(parts[2])
                z = float(parts[3])
            except ValueError:
                out_lines.append(raw)
                continue
            # Draw three offsets in a deterministic order keyed by atom
            # index so adding/removing atoms doesn't reshuffle the whole
            # sequence.
            dx = rng.uniform(-magnitude, magnitude)
            dy = rng.uniform(-magnitude, magnitude)
            dz = rng.uniform(-magnitude, magnitude)
            if i in idx_set:
                x += dx
                y += dy
                z += dz
            out_lines.append(f"{sym:4s} {x:12.6f} {y:12.6f} {z:12.6f}")
        return "\n".join(out_lines) + ("\n" if xyz_delfin.endswith("\n") else "")
    except Exception:
        return xyz_delfin


def _class_conditional_flag(name: str, mol, default: int = 0,
                            default_classes=None) -> bool:
    """Class-conditional env-flag (Iter-8.6c pattern, generalized Phase 3B).

    Precedence:
      1. If ``DELFIN_<name>_CLASSES`` is set (comma-separated class list)
         → flag is True iff ``_classify_complex_class(mol)`` is in that list.
      2. Else if ``default_classes`` is provided (Wave-4 extension)
         → flag is True iff ``_classify_complex_class(mol)`` is in that list,
         AND ``DELFIN_<name>`` is not explicitly set to ``0``.
      3. Else ``DELFIN_<name>`` integer (default ``default``).

    Wave-4 extension: ``default_classes`` lets us ship class-targeted
    default-ON behaviour without forcing the operator to set an env-var.
    Operator can still override via DELFIN_<name>_CLASSES (different list)
    or DELFIN_<name>=0 (disable entirely).
    """
    classes_env = f"{name}_CLASSES"
    env_classes_raw = os.environ.get(classes_env, "") or ""
    classes = {x.strip() for x in env_classes_raw.split(",") if x.strip()}
    if classes:
        try:
            return _classify_complex_class(mol) in classes
        except Exception:
            return False
    if default_classes:
        if _delfin_env_int(name, 1) == 0:
            return False
        try:
            return _classify_complex_class(mol) in set(default_classes)
        except Exception:
            return False
    return bool(_delfin_env_int(name, default))


def _every_append_gate_enabled(mol) -> bool:
    """Iter-20 (2026-05-19) wrapper for DELFIN_EVERY_APPEND_GATE.

    Default-flipped 0 → 1 for sigma class only (default_classes=["sigma"]).
    Per cross-archive analysis CROSS_ARCHIVE_RERUN_2026_05_18.md, 123a130
    is Champion in F3_bond (28.79%), bvs (72.75%), lig_realistic (26.22%),
    F23_funcgrp_geom (0.51%), cshm_mean (4.16) via per-bond emit-gate.

    Restricted to sigma class to avoid Co+2 M1B determinism issue
    (feedback_head_co_baseline_latent_nondet: EVERY_APPEND_GATE port drops
    Co+2 determinism from 80% to 60% via Chem.RWMol + AddConformer side
    effects).  Hapto / multi_hapto unchanged.

    Operator override via DELFIN_EVERY_APPEND_GATE=0 disables entirely,
    DELFIN_EVERY_APPEND_GATE_CLASSES=csv overrides the class allow-list.
    """
    return _class_conditional_flag(
        "DELFIN_EVERY_APPEND_GATE", mol, default=0,
        default_classes=["sigma"],
    )


def _multihapto_etkdg_fallback_enabled(mol) -> bool:
    """Iter-23 (2026-05-20) wrapper for DELFIN_MULTIHAPTO_ETKDG_FALLBACK.

    Re-enables the "fallback-as-feature" that commit 81f8a1f had by accident
    (WAVE7_Q archeology): when the analytical hapto scaffold yields only
    topology-broken candidates, fall back to ``_try_multiple_strategies``
    (raw stk + ETKDG seed=42).  fdeb9cb silently killed this when it fixed an
    orphan NameError that used to trigger the outer try/except fallback.

    Default-ON for multi_hapto only (where the scaffold most often produces
    broken Sn-bridge topology); zero-downside guardrail at the call site only
    swaps to the fallback when the fallback XYZ is *strictly* topology-OK.

    Operator override: DELFIN_MULTIHAPTO_ETKDG_FALLBACK=0 disables entirely,
    DELFIN_MULTIHAPTO_ETKDG_FALLBACK_CLASSES=csv overrides the class allow-list.
    """
    return _class_conditional_flag(
        "DELFIN_MULTIHAPTO_ETKDG_FALLBACK", mol, default=0,
        default_classes=["multi_hapto"],
    )


def _delfin_env_float(name: str, default: float) -> float:
    try:
        return float(os.environ.get(name, str(default)))
    except Exception:
        return default


# --- Iter-3 conformer-restore (post-123a130 multi-conformer regression) -----
# Default-ON env flag that lifts ``DELFIN_CHELATE_RANK_VARIANTS``,
# ``DELFIN_PRE_UFF_CAP_MULTIPLIER`` and ``DELFIN_TOPO_TEMPLATE_TOP_K`` from
# their post-Iter-3.5 minima (3,5,3) to (8,10,5) so that the per-(CF,perm)
# variant counter that emits ``-conf2`` / ``-conf3`` / ``-conf4`` labels
# downstream sees a deeper pool of distinct ranked-template / chelate-rank
# / scratch-builder candidates.  Smoke (115-SAKYIM 17→25, X10-JEMYOS 9→13,
# D-BURBOI 28→30) on Re/Mn/Fe CN-5/6 systems.  Time impact +0.5..+1.0 s on
# small systems, negligible on >50-atom complexes (UFF dominates).
# Set ``DELFIN_RESTORE_NUMCONFS=0`` for bit-exact pre-patch output.
_RESTORE_NUMCONFS = bool(_delfin_env_int("DELFIN_RESTORE_NUMCONFS", 1))


_RESTORE_RANKS    = 8 if _RESTORE_NUMCONFS else 3


_RESTORE_CAP_MULT = 10 if _RESTORE_NUMCONFS else 5


_RESTORE_TOPK     = 5 if _RESTORE_NUMCONFS else 3


# --- Conformer sampling & topology enumeration -----------------------------
DELFIN_TOP_LEVEL_SEED_COUNT: int = _delfin_env_int(
    "DELFIN_TOP_LEVEL_SEED_COUNT", 20
)
"""Number of ETKDG seeds for the top-level conformer sampling pool.

Default 20 is the best tradeoff between pool depth and stability of
the downstream topo-isomer builder:
beyond ~24 seeds, ``_rank_template_conformers`` routinely picks only
geometry-favoured SP-like templates for CN-5 systems (Fe(CO)3(NHC)2)
and suppresses the TBP canonical-form labels that the enumerator
otherwise emits.  Set ``DELFIN_TOP_LEVEL_SEED_COUNT=40`` or
``quality_mode='max'`` to widen the pool for extremely difficult
large ligand systems where more variety outweighs the template-bias
effect."""


# --- Class-aware ETKDG seed-count (opt-in, env-gated, default OFF) ---------
# Different coordination classes have different conformational dimensionality:
#
#   * sigma         : flexible monodentate donors -> ~20 seeds (default)
#   * hapto         : eta-ring constraints fix donor face -> fewer seeds
#   * multi_sigma   : bimetallic sigma -> bridge geometry adds DOF -> more
#   * multi_hapto   : bimetallic eta -> ring constraints + bridge DOF -> mid
#   * no_metal      : organic ligand alone -> same as sigma default
#
# Counts below are chemistry-motivated; pool-evaluator gates them.  Operator
# can override per class via ``DELFIN_CLASS_AWARE_SEEDS_<CLASS>=N`` (e.g.
# ``DELFIN_CLASS_AWARE_SEEDS_HAPTO=8`` to try an even tighter slice for the
# hapto class without touching the others).  ``DELFIN_CLASS_AWARE_SEEDS=0``
# (default) keeps the bit-exact pre-patch pipeline behaviour using the
# scalar ``DELFIN_TOP_LEVEL_SEED_COUNT``.
DELFIN_CLASS_AWARE_SEEDS: int = _delfin_env_int(
    "DELFIN_CLASS_AWARE_SEEDS", 0
)
"""Env-flag toggling class-aware ETKDG seed counts.  Default 0 (= disabled,
bit-exact pre-patch behaviour).  Set to 1 to enable per-class seed counts
from ``_class_aware_seed_count`` at every top-level ETKDG embedding site."""


# Chemistry-motivated per-class defaults (see module-level comment above).
# These values are deliberately documented in code so the operator can
# audit them without grepping; per-class env overrides go through
# ``DELFIN_CLASS_AWARE_SEEDS_<CLASS>``.
_CLASS_AWARE_SEED_DEFAULTS: Dict[str, int] = {
    "sigma":       20,
    "hapto":       12,
    "multi_sigma": 30,
    "multi_hapto": 15,
    "no_metal":    20,
}


def _class_aware_seed_count(class_label: Optional[str]) -> int:
    """Return the per-class ETKDG seed count.

    Resolution precedence:
      1. ``DELFIN_CLASS_AWARE_SEEDS_<CLASS>`` env-var (per-class override).
      2. ``_CLASS_AWARE_SEED_DEFAULTS[class_label]`` (chemistry defaults).
      3. ``DELFIN_TOP_LEVEL_SEED_COUNT`` (global scalar fallback) for any
         unknown / ``None`` class label.

    The minimum returned value is ``1``; callers may rely on
    ``_PIPELINE_SEEDS[:n]`` returning at least one seed.

    Side-effect-free: reads only ``os.environ`` and module constants.
    """
    fallback = max(1, int(DELFIN_TOP_LEVEL_SEED_COUNT))
    if not class_label or not isinstance(class_label, str):
        return fallback
    key = class_label.strip().lower()
    default = _CLASS_AWARE_SEED_DEFAULTS.get(key)
    if default is None:
        return fallback
    env_name = "DELFIN_CLASS_AWARE_SEEDS_" + key.upper()
    return max(1, int(_delfin_env_int(env_name, int(default))))


def _resolve_top_level_seed_count(mol) -> int:
    """Return the ETKDG seed count to use for ``mol`` at top-level sites.

    Honours ``DELFIN_CLASS_AWARE_SEEDS``: when 0 (default) returns the
    global ``DELFIN_TOP_LEVEL_SEED_COUNT``; when 1 classifies ``mol`` via
    ``_classify_complex_class`` and routes through
    ``_class_aware_seed_count``.

    Fail-safe: any classification exception falls back to the scalar so
    the patch can never regress conformer-pool depth for SMILES that
    happen to trip the classifier.
    """
    fallback = max(1, int(DELFIN_TOP_LEVEL_SEED_COUNT))
    if not _delfin_env_int("DELFIN_CLASS_AWARE_SEEDS", 0):
        n = fallback
    else:
        try:
            cls = _classify_complex_class(mol) if mol is not None else None
            n = _class_aware_seed_count(cls)
        except Exception:
            n = fallback
    # OPT-IN speed/quality knob (2026-07-06, default OFF -> full quality, byte-identical): fewer ETKDG
    # seeds LOWERS conformer quality, so we NEVER cut seeds automatically — the user decides the trade
    # (this cap is off unless DELFIN_ADAPTIVE_SEED_CAP=1).  The correct same-quality-faster path is
    # PARALLELISING all seeds across cores (see the embed pipeline), not reducing them.  When explicitly
    # enabled, scales seeds down for very large molecules (>120 heavy) as a last-resort finish guarantee.
    try:
        if _delfin_env_int("DELFIN_ADAPTIVE_SEED_CAP", 0) and mol is not None:
            n_heavy = sum(1 for a in mol.GetAtoms() if a.GetAtomicNum() > 1)
            if n_heavy > 120:
                cap = 3 if n_heavy > 200 else (4 if n_heavy > 160 else 6)
                n = min(n, cap)
    except Exception:
        pass
    return n


# --- Multi-sigma path V2 (Iter-multi_sigma audit, class-conditional default-ON) ---
# Forensics 2026-05-13: multi_sigma class shows 30.3%/36.4% (A/B) coverage,
# the lowest of any class (n=33).  Root cause is wall-clock budget exhaustion:
# 23/33 SMILES (70%) hit the external 600 s timeout while the embedding
# pipeline keeps churning through 20+ ETKDG seeds × 25 s _MULTIEMBED_TIMEOUT
# on 50-125-atom bimetallic systems.  Per-seed embedding alone consumes
# 10-12 s for a 100-atom Os-Sn complex (measured), so 20 top-level seeds +
# 20 multi-metal augmentation seeds × per-call timeout easily exceeds the
# pool-evaluator budget.  This V2 path tightens seed counts and adds a
# wall-clock budget for the augmentation block, without touching small
# multi-sigma molecules (which already work in ≤60 s).
#
# Wire-in 2026-05-15 — Phase 1.5: class-conditional default-ON.
# The patch is now active by default for SMILES classified as ``multi_sigma``;
# every other class is bit-exact pre-patch.  This recovers the 70% multi_sigma
# timeout-driven coverage loss documented above without affecting sigma/hapto
# pools that already converge well below the pool budget.
#
# Resolution semantics (mirrors ``_class_conditional_flag``):
#   * ``DELFIN_MULTI_SIGMA_PATH_V2`` unset    → active iff class is in
#     ``DELFIN_MULTI_SIGMA_PATH_V2_DEFAULT_CLASSES`` (multi_sigma).
#   * ``DELFIN_MULTI_SIGMA_PATH_V2=0``        → globally OFF (escape hatch).
#   * ``DELFIN_MULTI_SIGMA_PATH_V2=1``        → active iff class is in
#     the default class list (same as unset; kept for explicitness).
#   * ``DELFIN_MULTI_SIGMA_PATH_V2_CLASSES="…"`` → operator-supplied class
#     whitelist takes precedence over both of the above.
DELFIN_MULTI_SIGMA_PATH_V2_DEFAULT_CLASSES: Tuple[str, ...] = ("multi_sigma",)


def _multi_sigma_v2_active(mol) -> bool:
    """Return True iff the multi_sigma V2 path is enabled for *mol*.

    Resolution precedence (mirrors ``_class_conditional_flag`` semantics):
      1. ``DELFIN_MULTI_SIGMA_PATH_V2_CLASSES`` env-var (comma-separated
         class whitelist) — operator override, wins over the scalar gate.
      2. Else if ``DELFIN_MULTI_SIGMA_PATH_V2=0`` is explicitly set:
         return False (operator escape hatch — global opt-out).
      3. Else (scalar unset or ``=1``): True iff
         ``_classify_complex_class(mol)`` is in
         ``DELFIN_MULTI_SIGMA_PATH_V2_DEFAULT_CLASSES`` (i.e. multi_sigma).

    Fail-safe: classification exceptions return False so the V2 path can
    never regress for SMILES that happen to trip the classifier.
    """
    if mol is None:
        return False
    classes_env = os.environ.get("DELFIN_MULTI_SIGMA_PATH_V2_CLASSES", "") or ""
    classes = {x.strip() for x in classes_env.split(",") if x.strip()}
    if classes:
        try:
            return _classify_complex_class(mol) in classes
        except Exception:
            return False
    # Class-conditional default-ON: env unset or =1 → activate for the
    # default class set.  Explicit =0 is the global opt-out.
    if _delfin_env_int("DELFIN_MULTI_SIGMA_PATH_V2", 1) == 0:
        return False
    try:
        return _classify_complex_class(mol) in set(
            DELFIN_MULTI_SIGMA_PATH_V2_DEFAULT_CLASSES
        )
    except Exception:
        return False


def _multi_sigma_v2_budget(n_heavy: int) -> Dict[str, float]:
    """Return the seed/timeout budget for the multi_sigma V2 path.

    Scales with heavy-atom count.  Targets a per-SMILES wall-clock
    ≲ 200 s for the embedding pipeline (leaves room for topology
    enumeration, OB UFF and the final geometry gate).

    Returns a dict with:
      * ``seeds_top``: top-level ETKDG seed count cap.
      * ``seeds_mm_aug``: multi-metal augmentation seed count cap.
      * ``embed_timeout``: per-seed ``EmbedMultipleConfs`` timeout (s).
      * ``mm_walltime``: wall-clock budget for the augmentation block (s).

    Small multi-sigma SMILES (≤ 40 atoms) keep the pre-patch seed
    schedule — they already converge in ≤ 60 s.
    """
    n = max(1, int(n_heavy))
    if n <= 40:
        # Small bimetallic systems: no regression risk, keep defaults.
        return {
            "seeds_top": 20,
            "seeds_mm_aug": 10,
            "embed_timeout": 25.0,
            "mm_walltime": 120.0,
        }
    if n <= 60:
        return {
            "seeds_top": 8,
            "seeds_mm_aug": 6,
            "embed_timeout": 15.0,
            "mm_walltime": 90.0,
        }
    if n <= 90:
        return {
            "seeds_top": 5,
            "seeds_mm_aug": 4,
            "embed_timeout": 12.0,
            "mm_walltime": 60.0,
        }
    # Very large multi-metal (Fe4-tetranuclear, Pt-cyclam-OTf3, …).
    return {
        "seeds_top": 3,
        "seeds_mm_aug": 3,
        "embed_timeout": 10.0,
        "mm_walltime": 45.0,
    }


DELFIN_CHELATE_RANK_VARIANTS: int = _delfin_env_int(
    "DELFIN_CHELATE_RANK_VARIANTS", 3
)
"""Distinct chelate-conformer puckers tried per (CF, permutation).
Downstream dedup prunes duplicates; larger values widen the survivor
pool for the stricter gate."""


DELFIN_TOPO_TEMPLATE_TOP_K: int = _delfin_env_int(
    "DELFIN_TOPO_TEMPLATE_TOP_K", 3
)
"""Number of ranked ETKDG templates used per permutation by
``_build_topology_xyz_from_template``."""


DELFIN_PRE_UFF_CAP_MULTIPLIER: int = _delfin_env_int(
    "DELFIN_PRE_UFF_CAP_MULTIPLIER", 5
)
"""Pre-UFF candidate cap as multiple of ``max_isomers`` in
``_generate_topological_isomers``.  Larger ⇒ more UFF work but
deeper pool for dedup and gate."""


# --- Parallelism caps (keep modest so the pipeline runs on shared nodes) ---
DELFIN_MAX_THREAD_WORKERS: int = _delfin_env_int(
    "DELFIN_MAX_THREAD_WORKERS", 64
)
"""Upper bound for ThreadPoolExecutor workers (ETKDG sampling,
classification, topology build grid).  Scales to ``os.cpu_count()``
or this cap, whichever is smaller.  Default raised from 32 → 64 to
use more of the available CPU on batch / research runs; the per-call
``num_confs`` and ``cap_mult`` quality-profile knobs keep total work
bounded."""


DELFIN_MAX_PROCESS_WORKERS: int = _delfin_env_int(
    "DELFIN_MAX_PROCESS_WORKERS", 64
)
"""Upper bound for ProcessPoolExecutor workers (batch UFF).  OB
holds the GIL so true parallelism requires processes.  Default 64
matches the thread cap; lower via env var on shared login nodes
where 64 concurrent OB processes would saturate RAM."""


# --- Topology-gate thresholds (``_verify_topology_from_graph``) -----------
DELFIN_RULE4_PI_PLANAR_TOL_FRAC: float = _delfin_env_float(
    "DELFIN_RULE4_PI_PLANAR_TOL_FRAC", 0.25
)
"""Rule 4: π-ring planarity tolerance, expressed as a fraction of the
mean in-ring bond length."""


DELFIN_RULE5_INTERFRAG_COV_FACT: float = _delfin_env_float(
    "DELFIN_RULE5_INTERFRAG_COV_FACT", 1.15
)
"""Rule 5: minimum inter-fragment heavy-atom separation, expressed as a
multiple of the pair covalent-radius sum."""


DELFIN_RULE5_INNER_SPHERE_FACT: float = _delfin_env_float(
    "DELFIN_RULE5_INNER_SPHERE_FACT", 1.00
)
"""Rule 5 (softened): minimum separation for pairs where both atoms
sit in the inner coordination sphere of a multi-metal cluster."""


DELFIN_RULE6_METALLACYCLE_MAX_DEV: float = _delfin_env_float(
    "DELFIN_RULE6_METALLACYCLE_MAX_DEV", 0.40
)
"""Rule 6: maximum out-of-plane deviation (Å) of any atom of an
sp²-chelate metallacycle from the cycle's best-fit plane.  Enforces
the "metal in π-plane" rule at the gate level; tightened from 0.60 Å
to catch ~10° pyramidalisation of cyclometallated donor carbons."""


DELFIN_RULE7_SP2_OOP_MAX: float = _delfin_env_float(
    "DELFIN_RULE7_SP2_OOP_MAX", 0.35
)
"""Rule 7 sp² (purely organic): maximum out-of-plane distance (Å) of
a trigonal-planar atom from the plane of its three organic heavy
neighbours."""


DELFIN_RULE7_SP2_OOP_MAX_METAL: float = _delfin_env_float(
    "DELFIN_RULE7_SP2_OOP_MAX_METAL", 0.55
)
"""Rule 7 sp² (metal-bonded): looser out-of-plane budget for sp²
atoms whose neighbour set includes the metal (NHC carbenes,
cyclometallated donors, carbonyl carbons).  UFF without metal-specific
parameters routinely pyramidalises these by 0.3–0.5 Å without
distorting the rest of the topology, so the gate lets them through
while Rule 6 and the UFF metallacycle-torsion constraint keep the
average metal-in-π-plane behaviour enforced."""


DELFIN_RULE7_SP2_ANGLE_MIN: float = _delfin_env_float(
    "DELFIN_RULE7_SP2_ANGLE_MIN", 90.0
)
"""Rule 7 sp² lower angle bound (deg); 120° ± 30°."""


DELFIN_RULE7_SP2_ANGLE_MAX: float = _delfin_env_float(
    "DELFIN_RULE7_SP2_ANGLE_MAX", 150.0
)
"""Rule 7 sp² upper angle bound (deg)."""


DELFIN_RULE7_SP_MIN_ANGLE_DEG: float = _delfin_env_float(
    "DELFIN_RULE7_SP_MIN_ANGLE_DEG", 150.0
)
"""Rule 7 sp: minimum X-A-Y angle (deg) for a 2-coordinate atom with a
triple / cumulated-double bond.  Linear geometry target."""


DELFIN_RULE7_SP3_MIN_ANGLE_DEG: float = _delfin_env_float(
    "DELFIN_RULE7_SP3_MIN_ANGLE_DEG", 80.0
)
"""Rule 7 sp³: minimum X-A-Y angle (deg) for a 4-coordinate saturated
atom.  Tolerates cyclopropane-style strain down to ~84°."""


# --- Chelate conformer search ---------------------------------------------
DELFIN_CHELATE_N_TRIALS: int = _delfin_env_int(
    "DELFIN_CHELATE_N_TRIALS", 40
)
"""Number of ETKDG seeds explored in ``_chelate_conformer_candidates``."""


DELFIN_CHELATE_ACCEPT_DELTA: float = _delfin_env_float(
    "DELFIN_CHELATE_ACCEPT_DELTA", 0.10
)
"""Native-bite-match tolerance (Å) for a chelate conformer to be
accepted as a "good fit"."""


DELFIN_CHELATE_REJECT_DELTA: float = _delfin_env_float(
    "DELFIN_CHELATE_REJECT_DELTA", 0.50
)
"""Native-bite-match tolerance (Å) above which a chelate conformer is
outright rejected."""


DELFIN_CHELATE_CAP_30: int = _delfin_env_int("DELFIN_CHELATE_CAP_30", 12)
"""n_trials cap in ``_chelate_conformer_candidates`` for fragments with
>30 atoms.  Default 12 catches medium-size terdentate systems (terpy,
tpy-phenyl, tpy-NMe2, cyclam derivatives) that time out at the full
40 trials even though they are not >60 atoms.  Each seed costs ~6 s
on this fragment size, so 12 × 6 s = 72 s worst case per chelate
attempt, leaving headroom for multi-chelate multi-arrangement runs
inside the 900 s outer budget.  Verified direct tests:
  25-Os(terpy)(py-Ph)       n=3-5 in 80-120 s
  27-Fe(terpy)(py-p-tolyl)  n=4   in 50 s
  28-Fe(terpy-NMe2)2        n=3   in 35 s (DELFIN_CHELATE_N_TRIALS=10)"""


DELFIN_CHELATE_CAP_60: int = _delfin_env_int("DELFIN_CHELATE_CAP_60", 15)
"""n_trials cap in ``_chelate_conformer_candidates`` for fragments with
>60 atoms.  Default 15 is the conservative post-b32856f value; raise
to widen the chelate-pose search for heavy polydentate ligands at the
cost of wall-time (each seed costs ~6 s on >60-atom fragments)."""


DELFIN_CHELATE_CAP_90: int = _delfin_env_int("DELFIN_CHELATE_CAP_90", 8)
"""n_trials cap in ``_chelate_conformer_candidates`` for fragments with
>90 atoms.  Default 8 is the conservative post-b32856f value; raise
for regressions where the correct chelate pose never makes it into
the cap-8 subset."""


# --- Chelate-rank class-aware donor priority -------------------------------
# Default-OFF env-flag that augments the chelate-conformer ranking inside
# ``_chelate_conformer_candidates`` with a class-specific secondary score.
# Currently the candidates are ordered purely by donor-donor distance fit
# (``delta``), which on flexible polydentate ligands produces N nearly
# identical top conformers (same backbone pucker, same donor orientation).
# When the operator opts in, the secondary score breaks those near-ties in
# favour of conformers whose donor-element ordering matches the complex
# class — broadening the chelate-rank survivor pool and the downstream
# named-isomer coverage (see ``feedback_named_isomer_coverage`` memory).
#
# Score model (all values in Å-equivalent so the composite key is
# ``delta + alpha * class_penalty``, alpha chosen so class re-ranking only
# ever flips conformers whose ``delta`` differs by < ~0.05 Å; the strict
# fit gate is preserved):
#   class_penalty = sum(W[class][donor_element] for donor in fragment) /
#                   max(1, n_donors)
#   composite     = delta + DELFIN_CHELATE_RANK_CLASS_ALPHA * class_penalty
#
# Per-element weights encode chemistry heuristics validated against
# domain references (Lever AOM σ/π donor scales, Persson HSAB
# affinity table) without any SMILES-specific shortcuts:
#   * sigma         — favour N/O/S strong-σ chelating donors, mildly
#                     deprioritise halides (Cl/Br/I) and bulky P which
#                     historically over-rank by ETKDG distance fit alone.
#   * hapto         — neutral on σ donors, lift carbon donors (cyclometal,
#                     η-anchor) since the hapto class encodes
#                     mixed-η/σ topology.
#   * multi_sigma   — bridging-friendly donors (μ-OR, μ-Cl) up-weighted to
#                     keep bridge poses in the pool for the second metal.
#   * multi_hapto   — same as hapto, plus carbon-donor bonus to feed
#                     mixed-η bridging arrangements.
#   * no_metal      — never re-ranks (no_metal SMILES have no chelates).
#
# Negative weight = lower penalty = better (composite score smaller).
# Set ``DELFIN_CHELATE_RANK_CLASS_AWARE=1`` to enable.
# Operator override via ``DELFIN_CHELATE_RANK_CLASS_AWARE_CLASSES=cls1,cls2``
# limits the re-rank to the listed classes (e.g. "sigma" only).
DELFIN_CHELATE_RANK_CLASS_AWARE: int = _delfin_env_int(
    "DELFIN_CHELATE_RANK_CLASS_AWARE", 0
)
"""1 = enable class-aware chelate-conformer re-ranking; 0 = HEAD baseline
(pure donor-distance fit).  Default 0 keeps bit-exact behaviour."""


DELFIN_CHELATE_RANK_CLASS_ALPHA: float = _delfin_env_float(
    "DELFIN_CHELATE_RANK_CLASS_ALPHA", 0.04
)
"""Mixing coefficient for the class-aware penalty term.  Default 0.04 Å so
two conformers tying on ``delta`` to within ~0.04 Å can be re-ordered by
the class score, but a conformer with markedly better distance fit
(>0.05 Å advantage) always wins regardless of class preference."""


# Per-element class-aware weights: lower is better.  Conservative range
# [-0.5, +0.5] keeps the per-donor contribution well below ``alpha`` after
# averaging, so the composite never overwhelms the primary distance fit.
_CHELATE_CLASS_DONOR_WEIGHTS: Dict[str, Dict[str, float]] = {
    "sigma": {
        "N":  -0.40,  # strong-σ amine/imine/pyridine — preferred
        "O":  -0.30,  # carboxylate / phenolate / alkoxide
        "S":  -0.20,  # thiolate / thioether
        "P":   0.10,  # phosphine — over-ranked by distance alone
        "C":   0.00,  # NHC carbene / cyclometal — neutral
        "Cl":  0.30,  # halide — usually monodentate, deprioritise
        "Br":  0.30,
        "I":   0.30,
        "F":   0.20,
    },
    "hapto": {
        "N":  -0.20,
        "O":  -0.20,
        "S":  -0.10,
        "P":   0.00,
        "C":  -0.30,  # η-anchor / cyclometal carbon — preferred
        "Cl":  0.20,
        "Br":  0.20,
        "I":   0.20,
        "F":   0.10,
    },
    "multi_sigma": {
        "N":  -0.30,
        "O":  -0.40,  # bridging μ-OR / μ-O common in dimers — strong bonus
        "S":  -0.30,  # μ-S bridges
        "P":   0.05,
        "C":   0.00,
        "Cl": -0.20,  # μ-Cl bridges — common, lift them
        "Br": -0.10,
        "I":   0.00,
        "F":   0.10,
    },
    "multi_hapto": {
        "N":  -0.20,
        "O":  -0.30,
        "S":  -0.20,
        "P":   0.00,
        "C":  -0.30,
        "Cl": -0.10,
        "Br":  0.00,
        "I":   0.10,
        "F":   0.10,
    },
    "no_metal": {},   # no re-rank applied
}


def _chelate_class_donor_penalty(
    mol_or_frag,
    donor_atom_indices,
    class_label: str,
) -> float:
    """Compute the average class-aware donor-priority penalty for one
    chelate conformer.

    Returns the **mean** of per-donor weights from
    ``_CHELATE_CLASS_DONOR_WEIGHTS[class_label]``.  Unknown elements
    contribute 0.0 (neutral), so the helper is robust to exotic donors
    (Se, Te, Si, B) that the table does not enumerate.

    Returns 0.0 unconditionally for ``no_metal`` and unknown classes so
    the function is bit-exact when the env-flag is off, when the class
    table is empty, or when the fragment has no donors.
    """
    weights = _CHELATE_CLASS_DONOR_WEIGHTS.get(class_label or "", {})
    if not weights or not donor_atom_indices:
        return 0.0
    if mol_or_frag is None:
        return 0.0
    total = 0.0
    counted = 0
    for idx in donor_atom_indices:
        try:
            atom = mol_or_frag.GetAtomWithIdx(int(idx))
        except Exception:
            continue
        try:
            sym = atom.GetSymbol()
        except Exception:
            continue
        if sym in weights:
            total += float(weights[sym])
            counted += 1
        else:
            # Unknown donor element → neutral contribution (no penalty,
            # no bonus).  Counted so the average stays meaningful.
            counted += 1
    if counted == 0:
        return 0.0
    return total / counted


# --- Final-result geometry gate (smiles_to_xyz_isomers) -------------------
DELFIN_SEVERE_DIST_MAX_ABS: float = _delfin_env_float(
    "DELFIN_SEVERE_DIST_MAX_ABS", 2.4
)
"""Absolute upper bound (Å) on any non-metal covalent bond length used
by ``_has_severe_covalent_distortion``.  Raise for SMILES with legitimate
long bonds (I-I ~2.67, Te-Te ~2.72) that the default falsely rejects."""


DELFIN_SEVERE_DIST_MAX_SCALE: float = _delfin_env_float(
    "DELFIN_SEVERE_DIST_MAX_SCALE", 1.8
)
"""Maximum bond length relative to the sum of covalent radii, used by
``_has_severe_covalent_distortion``.  Raise to tolerate UFF-overstretched
bonds that still describe the correct topology."""


DELFIN_FFFREE_COLLAPSE_REJECT: int = _delfin_env_int("DELFIN_FFFREE_COLLAPSE_REJECT", 0)
"""ERDBEBEN collapse-reject (default 0 -> byte-identical).  When 1,
``_has_severe_covalent_distortion`` ALSO rejects a conformer that has
PLANAR-COLLAPSED -- a topology that must be 3D squashed into a plane by a
distance-geometry degeneracy (the eye's AQIBAE blind-spot).  It gives the
existing gate the symmetric LOWER bound it lacks, via two universal,
silent-on-clean signals that mirror the WEDDELL ``find_planar_collapse``
axis: (A) an sp3 centre whose four nearest bonded neighbours have a
tetra-volume below ``DELFIN_COLLAPSE_TETRA_VOL_MIN`` A^3 (a real sp3 centre
is ~1.5-2.5, clean-CCDC floor 0.81, even in a flat ring); (B) any covalent
bond compressed below ``DELFIN_COLLAPSE_BOND_MIN_SCALE`` x the covalent-radius
sum (below every real bond order -- clean shortest single 0.73x, triple N#N
~0.77x).  The existing embed-retry then re-embeds the conformer; one that
stays collapsed is dropped (it is physically impossible)."""


DELFIN_COLLAPSE_TETRA_VOL_MIN: float = _delfin_env_float(
    "DELFIN_COLLAPSE_TETRA_VOL_MIN", 0.40
)
"""sp3 tetra-volume (A^3) floor for the collapse-reject; below this the centre
is planar-collapsed (clean-CCDC min 0.81, so 0.40 is a 2x silent margin)."""


DELFIN_COLLAPSE_BOND_MIN_SCALE: float = _delfin_env_float(
    "DELFIN_COLLAPSE_BOND_MIN_SCALE", 0.65
)
"""Lower bound on a covalent bond length relative to the covalent-radius sum
for the collapse-reject (clean-CCDC shortest bond 0.73x -> 0.65 is silent)."""


DELFIN_FINAL_GATE_ENABLED: int = _delfin_env_int("DELFIN_FINAL_GATE_ENABLED", 1)
"""When 1 (default), the three ``_xyz_passes_final_geometry_checks``
guards in ``smiles_to_xyz_isomers`` (main conformer loop, linkage
isomer append, alt-binding append) reject structures that pass fragment-
topology but fail the severe-distortion check.  Set to 0 to bypass the
gate when diagnosing regressions where it is over-rejecting."""


DELFIN_TOPOLOGY_STRICT_MODE: int = _delfin_env_int("DELFIN_TOPOLOGY_STRICT_MODE", 0)
"""Iter-5 forensik (123a130 vs HEAD): default 0 = HEAD bit-identical.
When 1, short-circuits the 5 NEW frame-emitter paths added Apr 28-May 6
(HD-TA, GIE, _emit_chelate_pucker_variants, _emit_all_trans_by_type_arrangements,
H1 universal aromatic-H projection over all results) that pass HEAD's strict
gate at 1.15x cov-sum but fail the downstream covalent-radius bond-multiset
test at 1.26x cov-sum effective.  Restores 123a130-style frame mix while
keeping all other Iter-1...4 chemistry improvements (Lambda/Delta-universal,
Burnside, conformer-mult, B3/B4 corrector, in-pipeline pi-H projection at
ring-snap call sites).

Forensik report: agent_workspace/quality_framework/results/iter5_topology_forensik.md
Acceptance target: topo_pct_match >= 65 percent, topo_pct_extra_fragment <= 1.0 percent."""


# --- Symmetry / ideal-polyhedron scoring ---------------------------------
DELFIN_SYMMETRY_WEIGHT: float = _delfin_env_float(
    "DELFIN_SYMMETRY_WEIGHT", 1.0
)
"""Multiplier on the symmetry-bonus term inside
``_geometry_quality_score``.  Default 1.0 preserves legacy scoring;
raise (e.g. 2.0) to promote highly symmetric polyhedra above distorted
lookalikes.  Kept conservative by default so sweep data stays
comparable across commits; raise per-run for targeted ranking."""


DELFIN_IDEAL_POLYHEDRON_MAX_DEV: float = _delfin_env_float(
    "DELFIN_IDEAL_POLYHEDRON_MAX_DEV", 0.0
)
"""Maximum single-angle deviation (degrees) from the nearest ideal
polyhedron (Td / Oh / TBP / SP / PBP / SAP / TPR) that a conformer
may exhibit before ``_has_bad_geometry`` rejects it.  Default 0.0
disables the additional polyhedron-fidelity gate; set to e.g. 20.0
so a Jahn-Teller distortion (~5-10°) passes but a 25° broken
octahedron is rejected.  Measured per metal — any metal that exceeds
the threshold marks the whole structure bad."""


# --- Timeouts (seconds, per-call) -----------------------------------------
_EMBED_TIMEOUT: float = _delfin_env_float("DELFIN_EMBED_TIMEOUT", 10.0)
"""Per-call timeout for RDKit ``EmbedMolecule`` / stk embedding.

Large metal complexes with highly connected ring systems can cause ETKDG
distance-bounds calculation to hang indefinitely.  Override via env var
(e.g. ``DELFIN_EMBED_TIMEOUT=20``) on heavy macrocycles if needed.
"""


_OB_ROTOR_TIMEOUT: float = _delfin_env_float("DELFIN_OB_ROTOR_TIMEOUT", 15.0)
"""Per-call timeout for Open Babel rotor search."""


_MULTIEMBED_TIMEOUT: float = _delfin_env_float(
    "DELFIN_MULTIEMBED_TIMEOUT", 25.0
)
"""Per-call timeout for ``EmbedMultipleConfs``.

Running multiple seeds in parallel amortises the cost; a single slow
seed should not strangle the whole pool.  Increase via env var on
very large molecules.
"""


_CHELATE_EMBED_TIMEOUT: float = _delfin_env_float(
    "DELFIN_CHELATE_EMBED_TIMEOUT", 6.0
)
"""Per-seed timeout for the chelate conformer search."""


# Thread-local override for ``_MULTIEMBED_TIMEOUT``.  Used by the multi-sigma
# V2 path (and any future class-conditional patch) to push a tighter per-seed
# timeout into ``_embed_multiple_confs_with_timeout`` without changing the
# global default or thread-stamping every callsite.  Inactive when None.
_MULTIEMBED_TIMEOUT_OVERRIDE = threading.local()


# --- Deterministic seed schedule (shared across all stages) ---------------
def _generate_pipeline_seeds(n: int) -> Tuple[int, ...]:
    """Generate ``n`` distinct, deterministic seeds for ETKDG.

    The first 40 values match the hand-picked low-prime schedule that
    earlier regressions relied on (mono-metal determinism tests pin
    these exact integers).  Remaining values are drawn from a simple
    primes-via-trial-division generator starting at the largest
    hand-picked prime (26821), so every seed is a distinct positive
    integer < 2³¹ that changes across slots.
    """
    base = [
        31, 42, 7, 97, 13, 61, 83, 127, 211, 307,
        401, 503, 1009, 1619, 2027, 2531, 3181, 3847,
        4547, 5323, 6199, 7069, 8101, 9203, 10253,
        11329, 12409, 13499, 14591, 15683, 16787,
        17891, 19001, 20113, 21227, 22343, 23459,
        24571, 25703, 26821,
    ]
    if n <= len(base):
        return tuple(base[:n])
    seen = set(base)
    out = list(base)
    candidate = base[-1] + 2  # odd numbers upward
    while len(out) < n:
        is_prime = True
        r = int(candidate ** 0.5) + 1
        for d in range(3, r, 2):
            if candidate % d == 0:
                is_prime = False
                break
        if is_prime and candidate not in seen:
            seen.add(candidate)
            out.append(candidate)
        candidate += 2
    return tuple(out)


_PIPELINE_SEEDS: Tuple[int, ...] = _generate_pipeline_seeds(1024)
"""Single deterministic seed schedule shared by every stage of the
metal isomer pipeline (top-level conformer sampling, chelate
conformer search, fragment embedding, legacy fallbacks).  Keeping
one canonical source makes determinism guarantees hold end-to-end
and lets cross-stage caches stay coherent."""


_TOP_LEVEL_SEEDS: Tuple[int, ...] = _PIPELINE_SEEDS[
    :max(1, DELFIN_TOP_LEVEL_SEED_COUNT)
]
"""Resolved seed tuple actually used for top-level conformer sampling;
a slice of ``_PIPELINE_SEEDS`` controlled by
``DELFIN_TOP_LEVEL_SEED_COUNT``."""


# --- Quality-profile presets ----------------------------------------------
# ``alt_tries`` caps the template-iteration count inside
# ``_generate_alternative_binding_modes`` / linkage isomer generation so
# the dashboard doesn't spend 30+ s churning through 8 templates for
# alt-binding modes that never validate on cyclometallated / conjugated
# ligand systems like Ir(ppy)2.
_DELFIN_PROFILES: Dict[str, Dict[str, int]] = {
    "fast":    {"seeds": 12, "ranks": 3, "topk": 1, "cap_mult": 3,  "alt_tries": 0},
    "normal":  {"seeds": 20, "ranks": 3, "topk": 2, "cap_mult": 4,  "alt_tries": 4},
    "max":     {"seeds": 40, "ranks": 4, "topk": 4, "cap_mult": 8,  "alt_tries": 8},
    "extreme": {"seeds": 60, "ranks": 5, "topk": 5, "cap_mult": 12, "alt_tries": 12},
}


# ``extreme`` is for research / benchmark runs on a compute node with
# spare CPU and RAM — it widens every knob so the dashboard's quality
# cap does not stop the pipeline from covering every constitutional
# isomer + every backbone pucker an LHC-size candidate pool might
# turn up.  Stays fully deterministic (fixed seed schedule).
# ``ranks=3`` is shared across all profiles because ligand-conformer
# variety (chair / boat / twist puckers on flexible chelates) is part
# of the constitutional-isomer output the user needs in every mode.
# ``alt_tries=0`` disables alt-binding-mode exploration entirely for
# ``fast`` because iterating rewired donor candidates (e.g. ppy phenyl
# meta/para carbons) is almost always rejected by the fragment-topology
# gate on cyclometallated systems and costs 40-60 s per Ir(ppy)2-type
# SMILES without adding output variety.


def _resolve_quality_profile(name: Optional[str]) -> Dict[str, int]:
    """Return the profile dict for ``name`` or the module defaults."""
    if name is None:
        return {
            "seeds": DELFIN_TOP_LEVEL_SEED_COUNT,
            "ranks": DELFIN_CHELATE_RANK_VARIANTS,
            "topk":  DELFIN_TOPO_TEMPLATE_TOP_K,
            "cap_mult": DELFIN_PRE_UFF_CAP_MULTIPLIER,
            "alt_tries": 8,
        }
    key = name.lower().strip()
    if key not in _DELFIN_PROFILES:
        raise ValueError(
            f"quality_mode must be one of {sorted(_DELFIN_PROFILES)}, "
            f"got {name!r}"
        )
    return dict(_DELFIN_PROFILES[key])


# Default seed used by every embed call that does not receive a SMILES-derived
# seed.  RDKit's own default (``randomSeed = -1``) is wall-clock based and makes
# the embed non-deterministic; this fixed value restores reproducibility.
_DEFAULT_EMBED_SEED: int = 42


def _deterministic_embed_seed(smiles: Optional[str] = None) -> int:
    """Return a fixed, reproducible RDKit ``randomSeed`` for an embed call.

    DELFIN must be deterministic: the same SMILES + same code + same env has to
    produce bit-identical XYZ output.  RDKit's default ``randomSeed = -1`` seeds
    the embedder from the wall clock, so any embed call that forgets an explicit
    seed (notably the exception-fallback paths for KekulizeException /
    AtomValenceException) becomes non-reproducible.

    When *smiles* is given, the seed is derived from a BLAKE2b hash of the
    string so different molecules still get diverse-but-fixed seeds.  The
    derivation is purely a function of the input string, so it is universal
    (no per-SMILES special-casing) and stable across runs, machines and Python
    process restarts (unlike the salted built-in ``hash()``).  When *smiles* is
    ``None`` or empty the module-level :data:`_DEFAULT_EMBED_SEED` is returned.

    Args:
        smiles: Canonical (or raw) SMILES string, or ``None``.

    Returns:
        A non-negative ``int`` suitable for ``EmbedParameters.randomSeed`` /
        the ``randomSeed=`` keyword of ``AllChem.EmbedMolecule``.
    """
    if not smiles:
        return _DEFAULT_EMBED_SEED
    digest = hashlib.blake2b(smiles.encode("utf-8"), digest_size=4).digest()
    # Keep it inside the positive 31-bit range RDKit expects for a seed.
    return int.from_bytes(digest, "big") & 0x7FFFFFFF


def _trace_seating(*parts):
    """DIAGNOSTIC (env DELFIN_TRACE_SEATING=1, default OFF -> byte-identical): log which seating
    branch a donor takes, to stderr.  Pure OBSERVABILITY -- changes NO geometry.  Surfaced per system
    by loop.py --debug <rid> (which captures the build subprocess stderr)."""
    if os.environ.get("DELFIN_TRACE_SEATING", "0") == "1":
        try:
            import sys as _sys
            print("[SEATING]", *parts, file=_sys.stderr, flush=True)
        except Exception:
            pass
