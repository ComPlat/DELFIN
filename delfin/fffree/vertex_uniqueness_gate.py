"""Vertex-uniqueness gate (build-time atom-collapse fail-safe).

Universal post-build defence against the converter-pipeline pathology where
two or more heavy atoms end up sharing (or very nearly sharing) the same 3D
position — most visibly when a polydentate ligand's donor atoms or two
monodentate donor atoms get routed to the same polyhedron vertex, or when
the legacy ETKDG path is reached for a charged metal SMILES that
``contains_metal()`` failed to recognise (``[FeH2-N]``, ``[IrH-N]``,
``[RhH-N]`` etc.).  The forensik scan (1000 random voll-pool XYZ files)
counted 4.8 % of structures with at least one heavy-heavy pair below
1.0 A, with the worst case ``X10-CIDCOJ_4d_Y_CN7_hetero.xyz`` having three
O-pair near-duplicates at 0.002 - 0.009 A apart (catastrophic atomic
overlap).

Gate behaviour
--------------
``is_active()``
    Returns ``True`` iff one of the opt-in environment flags is set:

    - ``DELFIN_FFFREE_VERTEX_UNIQUENESS=1`` (the dedicated knob), or
    - ``DELFIN_FFFREE_PURE_TRACK3=1``      (the universal Track-3 super-gate).

    Default OFF -> **byte-identical** to HEAD when the flags are unset
    (`drop_collapsed_isomers` returns its input list unchanged before any
    parsing happens).

``drop_collapsed_isomers(results, *, hard_threshold=0.5, soft_threshold=1.0)``
    Read-only filter over a ``[(xyz_str, label), ...]`` list (the
    ``smiles_to_xyz_isomers`` public-API output shape).  Drops isomers
    whose XYZ contains ANY heavy-heavy pair below ``hard_threshold``
    (default 0.5 A — clearly inside atomic overlap, no realistic
    geometry produces this).  ``soft_threshold`` is used only for the
    counter that the caller may log (no isomer is dropped on the soft
    threshold; the gate is conservative-by-default so it never trims
    a legitimate but cramped construction).

The gate runs at the public entry point right next to
``_filter_nonfinite_isomers``, so EVERY build path (legacy, hapto,
fffree-construction, GRIP-ensemble) is covered — no per-path patching.

Universality + determinism
--------------------------
Pure geometric criterion (any element pair, any ligand class).  No SMILES
inspection, no RDKit calls, no RNG.  Two runs with ``PYTHONHASHSEED=0``
produce bit-identical outputs because dictionary iteration is not used:
isomers are processed in their input order and either kept or dropped
based on a fixed-threshold comparison.

The function silently tolerates malformed XYZ lines (treated as
non-heavy-atom rows) so the gate cannot itself raise.
"""
from __future__ import annotations

import os
from typing import Iterable, List, Sequence, Tuple


# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------
_FLAG_PRIMARY = "DELFIN_FFFREE_VERTEX_UNIQUENESS"
_FLAG_SUPERGATE = "DELFIN_FFFREE_PURE_TRACK3"

# Default catastrophic-collapse threshold (Angstrom).
# Two heavy atoms below this distance are considered unrecoverable:
# even the shortest realistic heavy-heavy contact (a triple bond, e.g.
# acetylenic C#C at ~1.20 A) is comfortably above 0.5 A.  All near-zero
# duplicates from the forensik (CIDCOJ O-O 0.002 A, AQAVOF Fe-P 0.012 A,
# EGUKUN Cl-Cl 0.016 A, RUQVAC Cl-Cl 0.143 A) fall under this threshold
# by an order of magnitude.
DEFAULT_HARD_THRESHOLD = 0.5

# Soft threshold used purely for logging / counting (no isomer dropped).
DEFAULT_SOFT_THRESHOLD = 1.0


def is_active() -> bool:
    """Return ``True`` iff the gate should run.

    The gate is opt-in: default OFF -> byte-identical to HEAD.  Two
    environment flags activate it: the dedicated knob
    ``DELFIN_FFFREE_VERTEX_UNIQUENESS=1`` and the universal Track-3
    super-gate ``DELFIN_FFFREE_PURE_TRACK3=1``.  Any truthy value
    (``"1"`` is the canonical) suffices; ``"0"`` / unset / empty leaves
    the gate inactive.
    """
    for flag in (_FLAG_PRIMARY, _FLAG_SUPERGATE):
        v = os.environ.get(flag, "").strip()
        if v and v != "0":
            return True
    return False


# ---------------------------------------------------------------------------
# XYZ parsing helpers (tolerant; never raise)
# ---------------------------------------------------------------------------
def _parse_heavy_coords(xyz_str: str) -> Tuple[List[str], List[Tuple[float, float, float]]]:
    """Extract ``(syms, coords)`` for non-hydrogen atoms from an XYZ string.

    Tolerant: ignores the optional XYZ header (atom count + comment), skips
    lines that do not parse as ``<symbol> <x> <y> <z>``.  Used because the
    converter's public-API output shape is not strictly canonical (some
    paths emit headerless body-only XYZ; see ``_organic_conformer_pool``).
    """
    syms: List[str] = []
    coords: List[Tuple[float, float, float]] = []
    if not xyz_str:
        return syms, coords
    for ln in xyz_str.splitlines():
        parts = ln.split()
        if len(parts) < 4:
            continue
        sym = parts[0]
        # Skip pure H (case-sensitive: XYZ uses element symbols)
        if sym == "H":
            continue
        try:
            x = float(parts[1]); y = float(parts[2]); z = float(parts[3])
        except (ValueError, IndexError):
            continue
        syms.append(sym)
        coords.append((x, y, z))
    return syms, coords


def min_heavy_heavy_distance(xyz_str: str) -> float:
    """Smallest pairwise heavy-heavy distance in ``xyz_str``.

    Returns ``float("inf")`` if fewer than two heavy atoms could be
    parsed (empty XYZ, headerless single-atom fragment, parsing
    failure).  Pure Python (no numpy) so the helper is import-cheap
    and side-effect-free for downstream callers.
    """
    _syms, coords = _parse_heavy_coords(xyz_str)
    n = len(coords)
    if n < 2:
        return float("inf")
    best = float("inf")
    for i in range(n):
        xi, yi, zi = coords[i]
        for j in range(i + 1, n):
            xj, yj, zj = coords[j]
            dx = xi - xj; dy = yi - yj; dz = zi - zj
            d2 = dx * dx + dy * dy + dz * dz
            if d2 < best:
                best = d2
    return best ** 0.5


def has_collapsed_pair(xyz_str: str, threshold: float = DEFAULT_HARD_THRESHOLD) -> bool:
    """``True`` iff ``xyz_str`` contains any heavy-heavy pair below ``threshold``.

    Early-out: stops at the first violating pair (cheaper than computing
    the global minimum when a structure is heavily collapsed).
    """
    _syms, coords = _parse_heavy_coords(xyz_str)
    n = len(coords)
    if n < 2:
        return False
    thr2 = threshold * threshold
    for i in range(n):
        xi, yi, zi = coords[i]
        for j in range(i + 1, n):
            xj, yj, zj = coords[j]
            dx = xi - xj; dy = yi - yj; dz = zi - zj
            if dx * dx + dy * dy + dz * dz < thr2:
                return True
    return False


# ---------------------------------------------------------------------------
# Public gate
# ---------------------------------------------------------------------------
def drop_collapsed_isomers(
    results: Sequence,
    *,
    hard_threshold: float = DEFAULT_HARD_THRESHOLD,
    soft_threshold: float = DEFAULT_SOFT_THRESHOLD,
):
    """Drop isomers with catastrophic heavy-heavy collapse (gate-active only).

    Parameters
    ----------
    results
        ``[(xyz_str, label), ...]`` (the ``smiles_to_xyz_isomers``
        public-API output shape).  ``None`` and empty inputs are
        returned unchanged.
    hard_threshold
        Heavy-heavy distance below which an isomer is dropped.
    soft_threshold
        Logging-only threshold (no isomer dropped); kept in the API for
        downstream metric collection.

    Returns
    -------
    list
        The filtered list.  When ``is_active()`` is ``False``, returns
        ``results`` unchanged (literally the same object) so the gate
        is byte-identical to HEAD on the disable path.
    """
    if not is_active() or not results:
        return results
    kept = []
    for item in results:
        if isinstance(item, (tuple, list)) and item:
            xyz = item[0]
        else:
            xyz = item
        if not isinstance(xyz, str):
            kept.append(item)
            continue
        if has_collapsed_pair(xyz, hard_threshold):
            continue
        kept.append(item)
    # Soft-threshold counter is informational only.  No isomers are
    # dropped from the soft tier; we leave the variable assignment in
    # for callers that want to wire telemetry later.
    _ = soft_threshold
    return kept


def contains_metal_extended(smiles: str, metals: Iterable[str]) -> bool:
    """Charged-H-tolerant metal detection.

    The legacy ``contains_metal()`` regex (``\\[<M>[+\\-\\d\\]@]``) requires
    the character immediately after the metal symbol in a bracketed atom
    to be one of ``+``, ``-``, digit, ``]`` or ``@``.  That rejects the
    common charged-hydride / charged-coordinated organometallic SMILES
    where an explicit ``H`` follows the metal:

    - ``[FeH2-4]`` (Fe with 2H and -4 charge, e.g. AQAVOF)
    - ``[IrH-3]`` (Ir with 1H and -3 charge, e.g. EGUKUN, HAZTIN)
    - ``[RhH-3]`` (Rh with 1H and -3 charge, e.g. RUQVAC)
    - ``[CrH-3]``, ``[CoH-4]`` ... (any d-block charged hydride)

    These SMILES were silently classified as **non-metal** and routed
    through ETKDG's organic conformer pool, which produces catastrophic
    metal-on-donor overlap (the very 4 % build-collapse population that
    motivates this module).  The extended detector adds the
    ``\\[<M>H`` pattern so any bracketed metal-symbol followed by
    an explicit hydrogen count is recognised as a metal SMILES.

    Pure regex; no RDKit; deterministic.
    """
    import re
    s = smiles or ""
    for m in metals:
        esc = re.escape(m)
        # original window: +/-/digit/]/@
        if re.search(rf"\[{esc}[+\-\d\]@]", s, re.IGNORECASE):
            return True
        # explicit-H window (charged hydride / explicit-H heavy-atom block)
        if re.search(rf"\[{esc}H", s, re.IGNORECASE):
            return True
        # bare-bracket (`[Fe]`)
        if re.search(rf"\[{esc}\]", s, re.IGNORECASE):
            return True
    return False
