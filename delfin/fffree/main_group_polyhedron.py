"""delfin.fffree.main_group_polyhedron — VSEPR-aware polyhedron picker for
main-group / p-block metals with a stereochemically-active lone pair
(Mogul-PRIMARY extension, 2026-06-07, hmaximilian).

Background
==========

The legacy polyhedron picker (``delfin.fffree.polyhedra.geometries_for_cn``)
is calibrated against d-block transition metals.  Main-group / p-block
metals routinely classified by ``delfin._bond_decollapse._is_metal`` as
"metals" (In, Sn, Sb, Pb, Bi, Tl, Ge, Po; and lighter p-block: Al, Ga)
nominally fit the same CN → polyhedron lookup, but the **lower oxidation
states** carry a **stereochemically-active lone pair** that occupies one
vertex of an *effective* (CN+1) polyhedron:

   *  Sn(II), Pb(II), Sb(III), Bi(III), Ge(II), As(III), Tl(I), In(I)
      → lone pair active → effective CN = observed CN + 1
   *  Sn(IV), Pb(IV), Sb(V), Bi(V), Ge(IV), As(V), Tl(III), In(III), Al(III), Ga(III)
      → no LP → effective CN = observed CN (legacy / TM picker is fine)

Consequence on the build:

   *  Sn(II) with **CN = 3** σ-donors does NOT want a trigonal-planar
      polyhedron — it wants a *pyramidal-of-Td* (NH₃-like), i.e. the
      donors live on three of four T-4 vertices and the LP sits at the
      fourth.  The legacy picker selects ``SP-3 trigonal planar`` (flat
      120°), forcing the three donors into a plane that the real
      geometry does NOT span.
   *  Sn(II) with **CN = 4** σ-donors wants a **see-saw / disphenoidal**
      (TBP minus an equatorial = SF₄-like), not a regular tetrahedron.
   *  In(I) with **CN = 2** wants a **bent / V-shape** (SP-3 with one
      vertex missing = SO₂-like, ~96°), not a linear D∞h pair.

This module provides:

   *  ``has_stereo_active_lp(metal_sym, oxidation_state)`` — universal
      VSEPR rule (period + group + ox state → LP presence).
   *  ``main_group_polyhedron_for(metal_sym, observed_cn, ox_state)`` —
      returns the EFFECTIVE polyhedron name + the donor-only vertex
      indices (so the assembler still places exactly ``observed_cn``
      donors; the LP vertex is left empty).
   *  ``effective_ref_vectors_main_group(geometry)`` — vertex unit
      vectors for the donor-only sub-set of the LP-aware polyhedron
      (drop-in for ``polyhedra.ref_vectors``).

Design contract
===============

   *  **NO SMILES patterns.**  Universal lookup from periodic-table
      Z + group + ox-state → LP presence.
   *  **NO per-class branches.**  The donor-only vertex array is
      derived from the same registered polyhedra (T-4, TBP-5, SPY-5,
      OC-6) by *dropping* the LP vertex; no new geometry definitions
      are needed below CN 4 (we re-use T-4 / TBP-5 with one vertex
      reserved).
   *  **Default-OFF byte-identical.**  When
      ``DELFIN_FFFREE_MAIN_GROUP_LP=0`` (or MOGUL_PRIMARY off and the
      gate is unset) the picker is a no-op and callers fall through to
      the legacy ``polyhedra.geometries_for_cn`` rule.
   *  **Auto-ON under MOGUL_PRIMARY.**  When MOGUL_PRIMARY=1 and the
      LP flag is unset, the LP-aware picker is active so the affected
      main-group cases (In CN=1/2, Sn CN=1/3, Sb/Bi CN=3) build the
      correct distorted-TM geometry by default.

The picker only fires when:

   1. The metal element is one of the LP-bearing p-block metals
      tabulated in ``_LP_PRESENT_BY_METAL`` AND
   2. The oxidation state matches the LP-active row in the table
      (e.g. Sn²⁺ but NOT Sn⁴⁺), AND
   3. The observed CN < formal capacity (so the LP vertex is "needed").

When any condition fails the picker returns ``None`` and callers
fall through to the legacy polyhedron list.

Author: hmaximilian <hmaximilian496@gmail.com>
Branch: feat-mogul-primary-2026-06-07
"""
from __future__ import annotations

import math
import os
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np

__all__ = [
    "main_group_lp_enabled",
    "is_main_group_lp_metal",
    "has_stereo_active_lp",
    "effective_cn_with_lp",
    "main_group_polyhedron_for",
    "effective_ref_vectors_main_group",
    "MAIN_GROUP_GEOM_BY_CN",
]


# ---------------------------------------------------------------------------
# Env gates
# ---------------------------------------------------------------------------


_FLAG_MGLP = "DELFIN_FFFREE_MAIN_GROUP_LP"
_FLAG_MOGUL_PRIMARY = "DELFIN_FFFREE_MOGUL_PRIMARY"
_FLAG_PT3 = "DELFIN_FFFREE_PURE_TRACK3"


def main_group_lp_enabled() -> bool:
    """True iff the LP-aware main-group picker is wired on.

    Default ON when MOGUL_PRIMARY=1 (the affected main-group cases need
    the LP correction by default in the production stack); default OFF
    otherwise so the legacy non-MOGUL path stays byte-identical.

    Users can switch back to the legacy picker by setting
    ``DELFIN_FFFREE_MAIN_GROUP_LP=0`` explicitly.
    """
    raw = os.environ.get(_FLAG_MGLP, "").strip()
    if raw != "":
        return raw == "1"
    # Default ON under MOGUL_PRIMARY (or PURE_TRACK3); OFF otherwise.
    return (
        os.environ.get(_FLAG_MOGUL_PRIMARY, "0") == "1"
        or os.environ.get(_FLAG_PT3, "0") == "1"
    )


# ---------------------------------------------------------------------------
# Element classification: which main-group metals can carry a stereo-active LP
# ---------------------------------------------------------------------------


# Per-element rule:  ``oxidation_state in lp_active_ox_states`` → LP active.
#
# Universal periodic-table look-up.  Each row covers ONE element with the
# oxidation states whose valence-shell holds an unshared pair (group ox −
# 2 in the heavier p-block).  No SMILES, no per-class chemistry.
#
#  Group 13 (s² p¹): low ox = +1 carries one LP (only meaningful for
#       the heavier elements In, Tl; B/Al/Ga +1 are rare/exotic).
#  Group 14 (s² p²): low ox = +2 carries one LP (Ge, Sn, Pb).
#  Group 15 (s² p³): mid ox = +3 carries one LP (As, Sb, Bi).
#  Group 16 (s² p⁴): mid ox = +4 carries one LP (Se, Te, Po) — these
#       are classified as "metals" by ``_is_metal`` only for Po.
#
# Empty oxidation list  →  element NEVER carries a stereo-active LP at
# normal coordination (e.g. high-ox-state Sn⁴⁺ behaves as a d⁰ TM).
_LP_PRESENT_BY_METAL: Dict[str, Tuple[int, ...]] = {
    # Group 13
    "In":  (1,),                # In(I) carries LP, In(III) does not
    "Tl":  (1,),                # Tl(I) carries LP, Tl(III) does not
    # Group 14
    "Ge":  (2,),                # Ge(II) carries LP, Ge(IV) does not
    "Sn":  (2,),                # Sn(II) carries LP, Sn(IV) does not
    "Pb":  (2,),                # Pb(II) carries LP, Pb(IV) does not
    # Group 15  (also classified as "metals" by _is_metal for Sb, Bi)
    "As":  (3,),                # As(III) carries LP, As(V) does not
    "Sb":  (3,),                # Sb(III) carries LP, Sb(V) does not
    "Bi":  (3,),                # Bi(III) carries LP, Bi(V) does not
    # Group 16
    "Po":  (4,),                # Po(IV) carries LP, Po(VI) does not
    # NOTE: Al(III), Ga(III) — NEVER an LP under chemically relevant ox states.
}


def is_main_group_lp_metal(metal_sym: str) -> bool:
    """True iff ``metal_sym`` is a p-block metal that CAN carry a
    stereo-active lone pair at some oxidation state.

    Pure look-up — does not check the actual ox state of the input.
    See :func:`has_stereo_active_lp` for the full predicate.
    """
    return str(metal_sym) in _LP_PRESENT_BY_METAL


def has_stereo_active_lp(metal_sym: str, oxidation_state: Optional[int]) -> bool:
    """True iff the metal at the given ox state carries a stereo-active LP.

    Returns False when:
      * The metal is not in the LP-bearing table (every d-block, f-block,
        s-block metal, plus Al/Ga).
      * ``oxidation_state`` is None — we DON'T guess; the caller must
        supply an inferred ox state or accept the no-LP fall-through.
      * The metal IS in the table but the ox state is not in its
        LP-active list (e.g. Sn⁴⁺, Bi⁵⁺).

    Universal: no SMILES, no per-class branches.
    """
    rule = _LP_PRESENT_BY_METAL.get(str(metal_sym))
    if rule is None:
        return False
    if oxidation_state is None:
        return False
    return int(oxidation_state) in rule


def effective_cn_with_lp(observed_cn: int) -> int:
    """Effective polyhedron CN when one vertex is occupied by a lone pair.

    Pure arithmetic — no SMILES, no chemistry-specific branches.
    ``observed_cn`` is the number of σ-donors actually attached to the
    metal (as supplied by ``decompose``); the LP adds one pseudo-vertex.
    """
    return int(observed_cn) + 1


# ---------------------------------------------------------------------------
# Donor-only polyhedra (LP-aware)
# ---------------------------------------------------------------------------
#
# Convention: the LP-aware polyhedron has CN_eff = N_donors + 1 vertices in
# its "parent" T-4 / TBP-5 / SPY-5 / OC-6 frame.  We RESERVE one vertex for
# the LP and return the remaining N_donors vertex unit vectors as the donor
# positions.  The chosen "parent" polyhedron and the reserved vertex are
# tabulated below per (group, observed CN).
#
# The reserved-vertex choice is the classic VSEPR result:
#
#  observed CN = 1, LP active        → linear (parent SP-3 trigonal-planar
#                                       with 2 LPs); donors at vertex 0.
#                                       This handles diatomic R-M: pairs.
#  observed CN = 2, LP active        → bent / V-shape (parent SP-3 with one
#                                       LP at vertex 2); ~120° at the metal,
#                                       not 180° linear.
#  observed CN = 3, LP active        → pyramidal NH₃-like (parent T-4 with
#                                       LP at vertex 3); ~107° at metal.
#  observed CN = 4, LP active        → see-saw SF₄-like (parent TBP-5 with
#                                       LP at one equatorial vertex; donors
#                                       on 2 axial + 2 equatorial).
#  observed CN = 5, LP active        → square pyramid BrF₅-like (parent
#                                       OC-6 with LP at axial vertex).
#  observed CN >= 6, LP active       → octahedral with LP delocalised (no
#                                       geometric distortion at this CN);
#                                       picker returns None → legacy path.
#
# The vertex INDEX ORDER used here matches the legacy polyhedron tables in
# ``polyhedra.py`` so the donor-donor distance bounds populated by
# ``mogul_bounds.build_bounds_matrix`` against the donor-only sub-set are
# consistent with the rest of the pipeline.


# Canonical name → (parent_polyhedron_name, lp_vertex_idx, donor_vertex_idxs).
#
# When the canonical name is used as a ``geometry`` argument anywhere in the
# code path, the assembler picks the parent's vertex array and indexes it by
# ``donor_vertex_idxs`` to place exactly ``observed_cn`` donors.  The
# missing vertex (LP slot) is left as an unfilled coordinate sink that the
# DG embed implicitly handles (no bond, no donor → no constraint).
_MG_POLYHEDRA: Dict[str, Tuple[str, int, Tuple[int, ...]]] = {
    # observed CN = 1  (linear R-M:)
    "LP1-1 linear-mono lone-pair":
        ("SP-3 trigonal planar", 1, (0,)),
    # observed CN = 2  (bent / V-shape, ~120° SO₂-like)
    "LP2-2 bent lone-pair":
        ("SP-3 trigonal planar", 2, (0, 1)),
    # observed CN = 3  (pyramidal NH₃-like, ~107° at metal)
    "LP3-3 pyramidal lone-pair":
        ("T-4 tetrahedron", 3, (0, 1, 2)),
    # observed CN = 4  (see-saw / SF₄ disphenoidal)
    # Parent TBP-5: vertex order [+z axial, -z axial, eq_0, eq_1, eq_2].
    # The LP sits at ONE equatorial vertex (vertex 4) so donors occupy
    # the 2 axials + 2 equatorials — the empirical see-saw.
    "LP4-4 see-saw lone-pair":
        ("TBP-5 trigonal bipyramid", 4, (0, 1, 2, 3)),
    # observed CN = 5  (square pyramid BrF₅-like)
    # Parent OC-6: vertex order [+x, -x, +y, -y, +z, -z] in REFS.
    # LP at -z (vertex 5) so donors occupy the +z apex + the equatorial
    # square — the empirical square pyramid (C4v).
    "LP5-5 square-pyramid lone-pair":
        ("OC-6 octahedron", 5, (0, 1, 2, 3, 4)),
}

# Aliases (shorthand → canonical full name).
_MG_ALIASES: Dict[str, str] = {
    "LP1-1": "LP1-1 linear-mono lone-pair",
    "LP2-2": "LP2-2 bent lone-pair",
    "LP3-3": "LP3-3 pyramidal lone-pair",
    "LP4-4": "LP4-4 see-saw lone-pair",
    "LP5-5": "LP5-5 square-pyramid lone-pair",
}


# CN_observed → list of LP-aware polyhedra (most-common first).  The
# legacy first-candidate rule applies inside this list.
MAIN_GROUP_GEOM_BY_CN: Dict[int, List[str]] = {
    1: ["LP1-1 linear-mono lone-pair"],
    2: ["LP2-2 bent lone-pair"],
    3: ["LP3-3 pyramidal lone-pair"],
    4: ["LP4-4 see-saw lone-pair"],
    5: ["LP5-5 square-pyramid lone-pair"],
    # CN >= 6: LP is delocalised → no geometric distortion → fall through
    # to the legacy OC-6 picker (caller sees an empty list and uses
    # ``polyhedra.geometries_for_cn`` as before).
}


def is_main_group_lp_geometry(geometry: str) -> bool:
    """True iff ``geometry`` is a registered LP-aware main-group polyhedron name."""
    g = str(geometry)
    return g in _MG_POLYHEDRA or g in _MG_ALIASES


def _resolve_canonical(geometry: str) -> Optional[str]:
    g = str(geometry)
    if g in _MG_POLYHEDRA:
        return g
    if g in _MG_ALIASES:
        return _MG_ALIASES[g]
    return None


def effective_ref_vectors_main_group(geometry: str) -> np.ndarray:
    """Donor-only unit vectors for an LP-aware polyhedron.

    Drop-in for :func:`polyhedra.ref_vectors`: returns an ``(N_donors, 3)``
    array of unit-normalised donor positions for the parent polyhedron
    with the LP vertex removed.

    Raises ``KeyError`` on unknown name so the caller can chain to
    :func:`polyhedra.ref_vectors` via try/except.
    """
    canonical = _resolve_canonical(geometry)
    if canonical is None:
        raise KeyError(geometry)
    parent, _lp_idx, donor_idxs = _MG_POLYHEDRA[canonical]
    try:
        from delfin.fffree import polyhedra as _polyhedra
    except ImportError as e:
        raise KeyError(geometry) from e
    parent_v = _polyhedra.ref_vectors(parent)
    if parent_v is None:
        raise KeyError(geometry)
    parent_v = np.asarray(parent_v, dtype=float)
    if parent_v.ndim != 2 or parent_v.shape[1] != 3:
        raise KeyError(geometry)
    if max(donor_idxs) >= parent_v.shape[0]:
        raise KeyError(geometry)
    sub = parent_v[list(donor_idxs)]
    norms = np.linalg.norm(sub, axis=1, keepdims=True)
    norms = np.where(norms < 1e-9, 1.0, norms)
    return sub / norms


def main_group_polyhedron_for(
    metal_sym: str,
    observed_cn: int,
    oxidation_state: Optional[int],
) -> Optional[str]:
    """Return the LP-aware polyhedron name for ``(metal, CN, ox state)``.

    Returns ``None`` when:
      * The LP-aware picker is disabled (env-gate OFF).
      * The metal does not carry a stereo-active LP at this ox state.
      * ``observed_cn`` is outside the LP-aware table (CN < 1 or CN > 5).
      * No ox state was supplied AND the LP-active and LP-inactive cases
        cannot be disambiguated (returns None → caller falls through to
        legacy).

    When ``oxidation_state`` is ``None`` and the metal has a UNIQUE
    LP-active ox state (Sn → only 2 carries LP; In → only 1; Sb → only
    3), and ``observed_cn`` is *consistent* with that LP-active state
    (i.e. CN ≤ formal-coord-cap − 1), the picker assumes LP-active.
    For metals with multiple ox-state options (none in the current
    table, all p-block metals here have exactly one LP-active state)
    the rule is unambiguous.

    Universal: pure look-up + arithmetic, no SMILES patterns.
    """
    if not main_group_lp_enabled():
        return None
    metal = str(metal_sym)
    cn = int(observed_cn)
    if cn < 1:
        return None
    if metal not in _LP_PRESENT_BY_METAL:
        return None

    # Decide LP-active.
    lp_active = False
    if oxidation_state is not None:
        lp_active = has_stereo_active_lp(metal, int(oxidation_state))
    else:
        # No explicit ox state supplied.  The conservative inference: if
        # the metal has EXACTLY ONE LP-active ox state AND ``observed_cn``
        # is low enough to be consistent with that ox state, assume LP.
        # All entries in ``_LP_PRESENT_BY_METAL`` currently have exactly
        # one LP-active ox state, so the rule reduces to: "low CN → LP".
        # The CN cap (CN ≤ 4 here) is the universal VSEPR cap for
        # the relevant lone-pair-bearing distorted geometries; CN ≥ 5
        # in the absence of an ox-state hint is interpreted as the
        # high-ox-state (LP-inactive) regime.
        active_states = _LP_PRESENT_BY_METAL[metal]
        if len(active_states) == 1 and cn <= 4:
            lp_active = True

    if not lp_active:
        return None
    # CN >= 6: LP is delocalised under the chemistry of these metals
    # at the observed coordination; legacy OC-6 / TPR-6 are fine.
    geoms = MAIN_GROUP_GEOM_BY_CN.get(cn)
    if not geoms:
        return None
    # First-candidate (most-common) rule inside the LP list — matches the
    # legacy ``geometries_for_cn`` convention.
    return geoms[0]


# ---------------------------------------------------------------------------
# Convenience: parent-polyhedron lookup (for callers that need to know
# which legacy polyhedron the LP-aware geometry is a sub-set of, e.g. for
# Pólya-key resolution downstream)
# ---------------------------------------------------------------------------


def parent_polyhedron(geometry: str) -> Optional[str]:
    """Return the legacy polyhedron the LP-aware geometry is derived from."""
    canonical = _resolve_canonical(geometry)
    if canonical is None:
        return None
    return _MG_POLYHEDRA[canonical][0]


def donor_vertex_indices(geometry: str) -> Optional[Tuple[int, ...]]:
    """Return the donor-vertex indices (into the parent polyhedron) for
    an LP-aware geometry, or ``None`` on unknown name.
    """
    canonical = _resolve_canonical(geometry)
    if canonical is None:
        return None
    return _MG_POLYHEDRA[canonical][2]


def lp_vertex_index(geometry: str) -> Optional[int]:
    """Return the reserved LP vertex index in the parent polyhedron."""
    canonical = _resolve_canonical(geometry)
    if canonical is None:
        return None
    return _MG_POLYHEDRA[canonical][1]
