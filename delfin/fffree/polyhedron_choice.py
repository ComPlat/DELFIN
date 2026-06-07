"""Polyhedron-CHOICE table for ambivalent coordination numbers (Task #93,
2026-06-07).

Background
==========
The legacy FF-free dispatcher (``decompose._default_geometry``) picks ONE
polyhedron per ``(metal, CN)`` pair:

  * CN4: ``SP-4`` if metal in ``{Pt, Pd, Ni, Au, Rh, Ir}``  else ``T-4``
  * CN5: always ``TBP-5``
  * CN6: always ``OC-6``

For genuinely ambivalent ``(metal, CN)`` combinations this loses isomers:

  * Co²⁺ d⁷ CN4 (e.g. CETUCT ``[Co]Cl₂(thiourea)₂``) — high-spin prefers
    tetrahedral T-4, low-spin prefers square-planar SP-4.  The decompose
    default picks T-4 only, so SP-4 is never enumerated.
  * Cu²⁺ d⁹ CN4 — Jahn-Teller flattens T-4 into SP-4-like; both should be
    emitted.
  * Ni²⁺ d⁸ CN4 — usually SP-4 but tetrahedral Ni(II) does exist with weak-
    field bulky ligands (e.g. NiBr₂(PPh₃)₂); both should be emitted.
  * Fe²⁺ d⁶ CN4 — both T-4 (HS, Fe(II)) and SP-4 (LS or radical / model
    geometries) are real.
  * CN5: TBP-5 ↔ SPY-5 are interconverted by the Berry pseudorotation; real
    CN5 complexes distribute across both for almost every TM.
  * CN6: OC-6 ↔ TPR-6 (early-TM Mo / W d⁰ / d¹).

Strict-choice metals (single polyhedron, no enumeration):

  * Pd²⁺, Pt²⁺ CN4: strict SP-4 (low-spin d⁸, kinetically inert).
  * Zn²⁺ CN4: strict T-4 (d¹⁰, ligand-field-symmetric).

Module contract
===============
The table ``_POLYHEDRON_CHOICES_BY_METAL_CN`` maps ``(metal, cn)`` to an
ORDERED list of geometry names (the first entry is the "preferred" default,
mirroring the order in the existing ``_default_geometry`` dispatch).

``get_polyhedron_choices(metal, cn)`` returns this list, or ``[]`` for any
``(metal, cn)`` not in the table.  An empty list signals "no ambivalence —
fall back to the single ``_default_geometry`` choice".

``polyhedron_choice_active()`` returns ``True`` when
``DELFIN_FFFREE_POLYHEDRON_CHOICE=1`` (or ``DELFIN_FFFREE_PURE_TRACK3=1``).
Default OFF byte-identical to HEAD.

``polyhedron_choice_max_geometries()`` caps the per-SMILES geometry count
(default ``2``) so even a future expansion to 3+ polyhedra per ``(metal, cn)``
stays bounded.

Determinism: all returns are lex-ordered + immutable; the function is pure
(no I/O, no env lookups except the active gate); ``PYTHONHASHSEED=0`` is
honoured by definition.

This module is GEOMETRY-SELECTION ONLY.  The actual polyhedron vertex
arrays + Pólya isomer enumeration + geometric assembly live in
``polyhedra.py``, ``polya_isomer_count.py`` and ``assemble_complex.py`` /
``converter_backend.py``.  The choice table is the SINGLE source of truth
for "which polyhedra should we enumerate for this metal+CN".
"""
from __future__ import annotations

import os
from typing import Dict, List, Tuple


# Canonical geometry NAMES (matching ``polyhedra.REFS`` / ``GEOM_BY_CN``).
# Keep the strings verbatim — they are used as keys throughout the assembler
# and converter backend.
_T4 = "T-4 tetrahedron"
_SP4 = "SP-4 square planar"
_TBP5 = "TBP-5 trigonal bipyramid"
_SPY5 = "SPY-5 square pyramid"
_OC6 = "OC-6 octahedron"
_TPR6 = "TPR-6 trigonal prism"
_SP3 = "SP-3 trigonal planar"
_T3 = "T-3 T-shape"


# Ordered polyhedron list per (metal, CN).  Order matters: the FIRST entry is
# treated as the "preferred default" (mirrors the existing one-shot dispatch
# in ``decompose._default_geometry``); subsequent entries are emitted as
# additional isomers when the env-gate is on.
#
# Empty list (or absence from table) ⇒ "no ambivalence" ⇒ the converter
# backend falls back to the single ``_default_geometry`` choice for that
# (metal, cn).  Default OFF byte-identical to HEAD.
_POLYHEDRON_CHOICES_BY_METAL_CN: Dict[Tuple[str, int], List[str]] = {
    # ---------------------------------------------------------------- CN4 --
    # d⁸ metals (Pt/Pd/Ni/Au/Rh/Ir): SP-4 is the legacy default.
    # Pt²⁺, Pd²⁺ are kinetically inert + strictly square-planar.
    ("Pd", 4): [_SP4],
    ("Pt", 4): [_SP4],
    # Ni²⁺: predominantly SP-4 (low-spin d⁸) but T-4 is real for weak-field
    # bulky ligands (e.g. NiBr₂(PPh₃)₂).
    ("Ni", 4): [_SP4, _T4],
    # Au¹⁺/Au³⁺: SP-4 dominant; T-4 is uncommon but exists in some clusters.
    ("Au", 4): [_SP4, _T4],
    # Rh¹⁺/Ir¹⁺: SP-4 dominant; tetrahedral analogues for d¹⁰-like configs.
    ("Rh", 4): [_SP4, _T4],
    ("Ir", 4): [_SP4, _T4],
    # Cu²⁺ d⁹: Jahn-Teller flattens T-4 toward SP-4; both exist.
    ("Cu", 4): [_SP4, _T4],
    # Co²⁺ d⁷: HS T-4 (common with halides / thio donors) vs LS SP-4 (rarer).
    # CETUCT case — Co(II) + 2 Cl + 2 S-thiourea — currently emits only T-4.
    ("Co", 4): [_T4, _SP4],
    # Fe²⁺ d⁶ / Fe³⁺ d⁵: predominantly T-4 with weak donors; SP-4 known too.
    ("Fe", 4): [_T4, _SP4],
    # Zn²⁺ d¹⁰: strict T-4 (ligand-field-symmetric; SP-4 essentially absent).
    ("Zn", 4): [_T4],
    # Cd²⁺ d¹⁰: strict T-4 (same rationale as Zn).
    ("Cd", 4): [_T4],
    # Hg²⁺ d¹⁰: strict T-4.
    ("Hg", 4): [_T4],
    # Mn²⁺ d⁵: T-4 dominant (halide / thio); SP-4 known in tetrathiolates.
    ("Mn", 4): [_T4, _SP4],
    # ---------------------------------------------------------------- CN5 --
    # TBP-5 ↔ SPY-5 are interconverted by the Berry pseudorotation; for every
    # ambivalent CN5 metal we emit BOTH.  Order = TBP-5 first (legacy default).
    ("V",  5): [_TBP5, _SPY5],
    ("Cr", 5): [_TBP5, _SPY5],
    ("Mn", 5): [_TBP5, _SPY5],
    ("Fe", 5): [_TBP5, _SPY5],
    ("Co", 5): [_TBP5, _SPY5],
    ("Ni", 5): [_TBP5, _SPY5],
    ("Cu", 5): [_TBP5, _SPY5],
    ("Zn", 5): [_TBP5, _SPY5],
    ("Mo", 5): [_TBP5, _SPY5],
    ("W",  5): [_TBP5, _SPY5],
    ("Ru", 5): [_TBP5, _SPY5],
    ("Rh", 5): [_TBP5, _SPY5],
    ("Pd", 5): [_TBP5, _SPY5],
    ("Os", 5): [_TBP5, _SPY5],
    ("Ir", 5): [_TBP5, _SPY5],
    ("Pt", 5): [_TBP5, _SPY5],
    # ---------------------------------------------------------------- CN6 --
    # Early-TM Mo / W (d⁰ / d¹) sometimes prefer TPR-6 over OC-6.  Other CN6
    # metals are strict OC-6 — the table omits them so the legacy single
    # OC-6 dispatch is byte-identical.
    ("Mo", 6): [_OC6, _TPR6],
    ("W",  6): [_OC6, _TPR6],
}


def get_polyhedron_choices(metal: str, cn: int) -> List[str]:
    """Return the ordered list of polyhedron names to enumerate for the
    given ``(metal, cn)`` pair.

    An empty list signals "no ambivalence" — the caller falls back to the
    single :func:`delfin.fffree.decompose._default_geometry` choice.

    Pure function: no env lookups, no I/O, deterministic.  Defensive copy of
    the table value so callers cannot mutate the canonical list.

    Parameters
    ----------
    metal : str
        Element symbol (e.g. ``"Co"``).  Case-sensitive — matches RDKit /
        decompose conventions.  Unknown metals return ``[]``.
    cn : int
        Coordination number.  Out-of-range values return ``[]``.
    """
    if not isinstance(metal, str) or not metal:
        return []
    try:
        cn_int = int(cn)
    except (TypeError, ValueError):
        return []
    choices = _POLYHEDRON_CHOICES_BY_METAL_CN.get((metal, cn_int), [])
    # Defensive copy: caller cannot mutate the canonical table.
    return list(choices)


def polyhedron_choice_active() -> bool:
    """True when ``DELFIN_FFFREE_POLYHEDRON_CHOICE=1`` (or
    ``DELFIN_FFFREE_PURE_TRACK3=1``).  Default OFF — byte-identical to
    HEAD when both are unset."""
    return (
        os.environ.get("DELFIN_FFFREE_POLYHEDRON_CHOICE", "0") == "1"
        or os.environ.get("DELFIN_FFFREE_PURE_TRACK3", "0") == "1"
    )


def polyhedron_choice_max_geometries() -> int:
    """Cap on the number of polyhedra to enumerate per (metal, CN).

    Read from ``DELFIN_FFFREE_POLYHEDRON_CHOICE_MAX_GEOMETRIES`` (default 2)
    so even a future expansion to 3+ polyhedra per pair stays bounded.
    """
    raw = os.environ.get("DELFIN_FFFREE_POLYHEDRON_CHOICE_MAX_GEOMETRIES", "2")
    try:
        n = int(raw)
    except (TypeError, ValueError):
        return 2
    return max(1, n)


def additional_polyhedra(metal: str, cn: int, primary: str) -> List[str]:
    """Return the polyhedra to ENUMERATE IN ADDITION to ``primary`` for the
    given ``(metal, cn)`` pair, respecting the max-geometries cap.

    Used by ``converter_backend`` to know which extra polyhedra to feed into
    ``_enumerate_geometry`` when the env-gate is on.  Returns ``[]`` when:

    * The gate is OFF.
    * ``(metal, cn)`` is absent from the table (no ambivalence).
    * The table lists only ``[primary]`` (strict-choice metal).
    * Every alternative already equals ``primary``.

    Parameters
    ----------
    metal, cn : metal symbol + coordination number.
    primary : geometry name already emitted by the main loop (will be
        EXCLUDED from the returned list to avoid duplicates).
    """
    if not polyhedron_choice_active():
        return []
    choices = get_polyhedron_choices(metal, cn)
    if not choices:
        return []
    cap = polyhedron_choice_max_geometries()
    # Lex-deterministic: preserve the table order (first entry = preferred).
    # Filter out ``primary`` (already emitted) + dedup while preserving order.
    seen = {primary}
    extra: List[str] = []
    for g in choices:
        if g in seen:
            continue
        seen.add(g)
        extra.append(g)
        # ``cap`` includes ``primary`` itself, so we can emit at most cap-1
        # additional polyhedra.
        if len(extra) >= max(0, cap - 1):
            break
    return extra
