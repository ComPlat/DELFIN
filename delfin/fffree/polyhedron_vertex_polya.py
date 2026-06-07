"""delfin.fffree.polyhedron_vertex_polya — Polyhedron-Vertex Pólya Enumeration.

Fundamental architecture extension (2026-06-07): enumerate the distinct
cis / trans / fac / mer / Δ-Λ donor arrangements at a metal coordination
polyhedron BEFORE 3D embedding, so the conformer ensemble emitted by
``embed_fallback`` covers all geometric isomers (orbits under the polyhedron's
proper rotation group), not just a single donor-arrangement pattern.

Motivation (BEYRAY-class bug, voll-pool 2026-06-06)
---------------------------------------------------
SMILES  ``O=C1[O][Cu-4]2([OH2+])([OH2+])([O]C(=O)C3=CN=CC=[N+]32)[N+]2=CC=NC=C12``
(Cu²⁺, 2 bidentate pyridyl-carboxylate chelates + 2 H₂O, OC-6 with donor
multi-set {N:2, O:2, OH₂:2}) currently emits 16 frames labelled
``embed-conf0-grip`` … ``embed-conf15-grip``: all 16 are ETKDG conformers of a
SINGLE vertex-coloring.  The textbook count of distinct cis / trans isomers
for [M(N-O)₂(H₂O)₂] with N-O chelating is ≥ 4 (cis-N₂ / trans-N₂ / cis-O₂ /
trans-O₂ × Δ/Λ).  The binding-mode enumerator (κ/η/linkage) handles binding-
mode diversity but NOT the vertex-assignment diversity within a fixed
polyhedron, so the BEYRAY-class population is systematically under-sampled.

Mathematical foundation (Burnside)
----------------------------------
For a metal with coordination number ``n`` and donor multi-set
``M = {d₁: c₁, d₂: c₂, …}``::

    |orbits|  =  (1/|G|) Σ_{g ∈ G} |Fix(g)|

with ``G`` = proper-rotation group of the coordination polyhedron (Oh for
OC-6, Td for T-4, D₃h for TBP-5, …) acting on the vertex set, and ``Fix(g)``
the colorings invariant under ``g`` (every cycle of ``g`` is monochromatic).

The orbit enumeration itself is a thin wrapper around
:func:`delfin.fffree.polya_isomer_count.enumerate_isomers` — the existing
module already encodes every polyhedron's rotation group and is unit-tested
against the textbook isomer counts (OC-6 A₂B₂C₂ = 6, OC-6 A₃B₃ = 2 fac/mer,
SP-4 A₂B₂ = 2 cis/trans, T-4 ABCD = 2, …).  This module adds:

1. A polyhedron name → group key bridge (e.g. ``"OC-6 octahedron"`` →
   ``"octahedron"``), so callers can pass the canonical geometry string used
   in :mod:`delfin.fffree.polyhedra`.

2. **Chelate-constraint filtering** — chelating ligands force their two
   donors onto an adjacent (cis) vertex pair.  Orbit representatives that
   violate any chelate-cis constraint are filtered out.  The cis-edge map
   is derived geometrically from the polyhedron's ideal vertex set (same
   100 deg cut-off as :mod:`polya_isomer_count`).

3. **Public top-level API** — ``enumerate_orbits_with_chelates`` operates
   on the canonical polyhedron string + a flat donor-type list (one per
   vertex) + chelate-pair indices into the donor list.

4. **Cost control** — a ``max_orbits`` cap (env-overridable via
   ``DELFIN_FFFREE_POLYHEDRON_VERTEX_POLYA_MAX_ORBITS``) prevents
   exponential blow-up on rare CN9+ mixed-donor cases.

Determinism
-----------
Pure-Python, no random sampling, no hash-iteration on the critical path —
all internal containers are explicitly sorted before being returned.  Two
runs with the same input produce bit-identical orbit lists with
``PYTHONHASHSEED=0`` and idempotent across machines.

Default OFF, byte-identical
---------------------------
Nothing in this module mutates state.  It is callable from any layer, but
the integration point (:mod:`delfin.fffree.embed_fallback`) only consults
the orbit list when ``DELFIN_FFFREE_POLYHEDRON_VERTEX_POLYA=1`` is set.
"""
from __future__ import annotations

import math
import os
from collections import Counter
from typing import Dict, List, Optional, Sequence, Tuple

# ---------------------------------------------------------------------------
# Polyhedron string → Pólya group key (matches polya_isomer_count._GROUPS).
#
# Canonical geometry strings are produced by ``polyhedra.GEOM_BY_CN`` /
# ``polyhedra.geometries_for_cn``; we accept BOTH the full canonical form
# (e.g. ``"OC-6 octahedron"``) and the leading-token alias (e.g. ``"OC-6"``)
# so callers can pass whichever they have available.  Unknown strings raise
# ``KeyError`` (deterministic; caller catches and silently disables the
# orbit enumeration for that SMILES, falling back to the legacy 16-conformer
# ETKDG output -- env-default-OFF byte-identical contract).
# ---------------------------------------------------------------------------
_SHAPE_TO_GROUP: Dict[str, str] = {
    "L-2 linear": "linear",
    "SP-3 trigonal planar": "trigonal_planar",
    "T-3 T-shape": "tshape",
    # Y-3 (Y-shape, 2026-06-07 CN3/CN5 completeness extension): C2v point group,
    # same proper-rotation group (C2, order 2) as T-shape, so the *orbit*
    # structure is identical — Burnside counts are equal for any donor
    # multi-set.  Y-shape differs from T-shape only in the angular geometry
    # (~150°/150°/60° vs 90°/90°/180° at the metal); placement uses the same
    # T-3 vertex array (the ETKDG/relaxation phase reshapes geometry from
    # the proper-rotation skeleton).  This keeps :mod:`polyhedra` as the
    # single source of truth for vertex coordinates while still allowing
    # Y-3 to be enumerated as a distinct polyhedron name when callers
    # request it.
    "Y-3 Y-shape": "tshape",
    "T-4 tetrahedron": "tetrahedron",
    "SP-4 square planar": "square_planar",
    "TBP-5 trigonal bipyramid": "trigonal_bipyramid",
    "SPY-5 square pyramid": "square_pyramid",
    "OC-6 octahedron": "octahedron",
    "TPR-6 trigonal prism": "trigonal_prism",
    "PB-7 pentagonal bipyramid": "pentagonal_bipyramid",
    "SQAP-8 square antiprism": "square_antiprism",
    "TTP-9 tricapped trigonal prism": "tricapped_trigonal_prism",
}

# Leading-token aliases (e.g. "OC-6" → "OC-6 octahedron").
_ALIAS_TO_SHAPE: Dict[str, str] = {
    "L-2": "L-2 linear",
    "SP-3": "SP-3 trigonal planar",
    "TP-3": "SP-3 trigonal planar",       # textbook "TP-3" alias for trigonal-planar
    "T-3": "T-3 T-shape",
    "Y-3": "Y-3 Y-shape",                 # Y-shape CN3 (2026-06-07 extension)
    "T-4": "T-4 tetrahedron",
    "SP-4": "SP-4 square planar",
    "TBP-5": "TBP-5 trigonal bipyramid",
    "SPY-5": "SPY-5 square pyramid",
    "OC-6": "OC-6 octahedron",
    "TPR-6": "TPR-6 trigonal prism",
    "PB-7": "PB-7 pentagonal bipyramid",
    "PBP-7": "PB-7 pentagonal bipyramid",
    "SQAP-8": "SQAP-8 square antiprism",
    "SAP-8": "SQAP-8 square antiprism",
    "TTP-9": "TTP-9 tricapped trigonal prism",
}


def _canonical_shape(polyhedron: str) -> str:
    """Normalise ``polyhedron`` to the canonical key in :data:`_SHAPE_TO_GROUP`.

    Accepts canonical (``"OC-6 octahedron"``), leading-token (``"OC-6"``), and
    stripped lowercase variants.  Raises ``KeyError`` for unknown values.
    """
    if not isinstance(polyhedron, str):
        raise KeyError(polyhedron)
    p = polyhedron.strip()
    if p in _SHAPE_TO_GROUP:
        return p
    # Leading token (first space-separated component).
    head = p.split()[0] if p else ""
    if head in _ALIAS_TO_SHAPE:
        return _ALIAS_TO_SHAPE[head]
    # Case-insensitive direct match.
    for k in _SHAPE_TO_GROUP:
        if k.lower() == p.lower():
            return k
    raise KeyError(polyhedron)


def _group_key(polyhedron: str) -> str:
    """Return the :mod:`polya_isomer_count` group key for ``polyhedron``."""
    return _SHAPE_TO_GROUP[_canonical_shape(polyhedron)]


# ---------------------------------------------------------------------------
# Cis-edge geometry (derived from the SAME ideal vertex array the assembler
# places into, see :mod:`delfin.fffree.polyhedra`).  Bidentate chelating
# donors must land on a cis pair; orbit representatives with any chelate pair
# on a non-cis edge are filtered out.
# ---------------------------------------------------------------------------
_CHELATE_CIS_MAX_DEG = 115.0   # matches polya_isomer_count.CHELATE_CIS_MAX_DEG


def cis_edges_for(polyhedron: str) -> List[Tuple[int, int]]:
    """Return all (i, j) vertex pairs with angular separation ≤ 115 deg.

    Single source of truth: the same ideal vertex array used by
    :mod:`delfin.fffree.polyhedra.ref_vectors`, with the same 115 deg cut-off
    as :mod:`polya_isomer_count.CHELATE_CIS_MAX_DEG`.  Returns sorted tuples
    (lex by i, then j) for determinism.
    """
    try:
        import numpy as np
        from delfin.fffree import polyhedra as _PLY
    except Exception:
        return []
    shape = _canonical_shape(polyhedron)
    try:
        V = _PLY.ref_vectors(shape)
    except KeyError:
        # Y-3 fallback: same geometric skeleton as T-3 for cis-edge
        # derivation (angular relaxation is downstream of orbit selection).
        if shape == "Y-3 Y-shape":
            try:
                V = _PLY.ref_vectors("T-3 T-shape")
            except KeyError:
                return []
        else:
            return []
    cos_max = math.cos(math.radians(_CHELATE_CIS_MAX_DEG))
    n = len(V)
    out: List[Tuple[int, int]] = []
    for i in range(n):
        for j in range(i + 1, n):
            c = float(V[i] @ V[j])
            if c > cos_max:
                out.append((i, j))
    return sorted(out)


# ---------------------------------------------------------------------------
# Core orbit enumeration.
# ---------------------------------------------------------------------------
def cn_for(polyhedron: str) -> int:
    """Coordination number (= number of polyhedron vertices) for ``polyhedron``."""
    from delfin.fffree.polya_isomer_count import _get_group as _gg
    _, n = _gg(_group_key(polyhedron))
    return int(n)


def enumerate_orbits(
    polyhedron: str,
    donor_types: Sequence[str],
) -> List[Tuple[str, ...]]:
    """Enumerate distinct vertex-coloring orbits for ``polyhedron`` + ``donor_types``.

    Parameters
    ----------
    polyhedron : str
        Canonical geometry string (``"OC-6 octahedron"``, ``"SP-4"`` etc.).
    donor_types : sequence of str
        Donor-type label per vertex slot (length must equal CN of polyhedron).
        Labels are arbitrary strings; the orbit count depends only on their
        multiplicities (Burnside).

    Returns
    -------
    list of tuple of str
        One canonical representative per orbit, each a length-CN tuple of
        donor labels (vertex index → assigned donor type).  Lex-sorted for
        determinism.
    """
    from delfin.fffree.polya_isomer_count import enumerate_isomers
    gkey = _group_key(polyhedron)
    n = cn_for(polyhedron)
    types = list(donor_types)
    if len(types) != n:
        raise ValueError(
            f"donor_types length {len(types)} != CN {n} for polyhedron {polyhedron!r}"
        )
    spec = dict(Counter(types))
    orbits = enumerate_isomers(gkey, spec)
    # Deterministic order: lex sort.
    return sorted(orbits)


def _orbit_satisfies_chelates(
    orbit: Sequence[str],
    donor_types: Sequence[str],
    chelate_pairs: Sequence[Tuple[int, int]],
    cis_edges: Sequence[Tuple[int, int]],
) -> bool:
    """Return True iff ``orbit`` admits an assignment of the chelate pairs to
    cis-edges that is consistent with the per-vertex donor labels.

    The orbit is a tuple of donor labels per VERTEX; the chelate_pairs are
    indices into the DONOR list.  We need to check that some bijection from
    donor-slots to vertices (respecting the orbit's vertex labels) maps every
    chelate pair onto a cis-edge.

    Algorithm: per-orbit, treat the multiplicities of each donor label as
    interchangeable; for each chelate pair (i, j) with donor labels
    (a, b), check that there EXIST two cis vertices (u, v) carrying labels
    (a, b) (in some order).  Since each chelate consumes one (i, j) slot,
    we run a small backtrack that tries each cis-edge assignment and marks
    the consumed vertices.
    """
    cis_set = {(min(a, b), max(a, b)) for a, b in cis_edges}
    types_list = list(donor_types)
    chelate_list = [(int(i), int(j)) for i, j in chelate_pairs]
    n = len(orbit)

    # Per-label list of vertex indices carrying that label in the orbit.
    label_to_vertices: Dict[str, List[int]] = {}
    for v, lab in enumerate(orbit):
        label_to_vertices.setdefault(lab, []).append(v)

    def backtrack(idx: int, used: frozenset) -> bool:
        if idx == len(chelate_list):
            return True
        i, j = chelate_list[idx]
        a = types_list[i]
        b = types_list[j]
        # Candidate vertex pairs: u carries a, v carries b, neither used,
        # (u, v) (canonical sort) is a cis-edge.
        cand_a = [v for v in label_to_vertices.get(a, []) if v not in used]
        cand_b = [v for v in label_to_vertices.get(b, []) if v not in used]
        # Enumerate unordered pairs (u, v) with u carrying a and v carrying b
        # OR u carrying b and v carrying a (handles a == b too).
        tried: set = set()
        for u in cand_a:
            for w in cand_b:
                if u == w:
                    continue
                pair = (min(u, w), max(u, w))
                if pair in tried:
                    continue
                tried.add(pair)
                if pair not in cis_set:
                    continue
                if backtrack(idx + 1, used | {u, w}):
                    return True
        # Symmetric branch is already covered when a != b because we
        # iterate cand_a x cand_b; when a == b the same pair will be tried
        # in both orders but ``tried`` dedups.
        return False

    return backtrack(0, frozenset())


def _max_orbits_cap() -> int:
    """Env-overridable cap on the number of orbits emitted per SMILES."""
    raw = os.environ.get(
        "DELFIN_FFFREE_POLYHEDRON_VERTEX_POLYA_MAX_ORBITS", ""
    ).strip()
    if raw:
        try:
            v = int(raw)
            if v > 0:
                return v
        except (TypeError, ValueError):
            pass
    return 10


def enumerate_orbits_with_chelates(
    polyhedron: str,
    donor_types: Sequence[str],
    chelate_pairs: Sequence[Tuple[int, int]],
    *,
    max_orbits: Optional[int] = None,
) -> List[Tuple[str, ...]]:
    """Enumerate distinct vertex-coloring orbits, filtering for chelate-cis.

    Parameters
    ----------
    polyhedron : str
        Canonical geometry string.
    donor_types : sequence of str
        Length-CN flat list of donor-type labels (one per vertex slot).
    chelate_pairs : sequence of (i, j) int pairs
        Donor-list indices that MUST be cis (each pair = one bidentate
        chelating ligand whose two donors are connected in the ligand
        backbone).
    max_orbits : int, optional
        Cap on the number of orbit representatives returned.  Defaults to
        the value of ``DELFIN_FFFREE_POLYHEDRON_VERTEX_POLYA_MAX_ORBITS``
        (or 10 if unset).

    Returns
    -------
    list of tuple of str
        Canonical orbit representatives (lex-sorted) that admit a chelate-
        cis-compatible assignment.  Capped at ``max_orbits``.
    """
    all_orbits = enumerate_orbits(polyhedron, donor_types)
    if not chelate_pairs:
        out = all_orbits
    else:
        cis = cis_edges_for(polyhedron)
        out = [
            o for o in all_orbits
            if _orbit_satisfies_chelates(o, donor_types, chelate_pairs, cis)
        ]
    cap = int(max_orbits) if max_orbits is not None else _max_orbits_cap()
    return out[:cap]


# ---------------------------------------------------------------------------
# Multi-polyhedron orbit enumeration (CN3 + CN5 completeness, 2026-06-07).
#
# Motivation (ALESUH-class bug):
# Cu(II) CN5 with donor multiset {4N, 1Br} is ambivalent between TBP-5 and
# SPY-5 — both are real Berry-pseudorotation-related minima.  Picking only
# the FIRST polyhedron from :data:`delfin.fffree.polyhedra.GEOM_BY_CN`
# (as :func:`_default_polyhedron` does) silently drops half the geometric-
# isomer ensemble.  This block enumerates orbits across *every* polyhedron
# registered for the CN and returns a list of (polyhedron, orbit) tuples so
# downstream embedding can dispatch one ETKDG seed per (polyhedron, orbit).
#
# Same Burnside engine, same chelate-cis filter, same determinism contract.
# Scope: CN3 ({SP-3, T-3} ± Y-3) and CN5 ({TBP-5, SPY-5}) -- the two CNs
# where ambivalent geometric minima are textbook chemistry.  For other CNs
# the single-polyhedron path is preserved bit-exact (the multi-polyhedron
# helper is a strict addition).
# ---------------------------------------------------------------------------

# CNs for which multi-polyhedron orbit enumeration is enabled.  Conservative
# scope: CN3 and CN5 only.  Extending to CN4 (T-4 vs SP-4) and CN6 (OC-6 vs
# TPR-6) is a future iter and would expand the ensemble for d8 Pd/Pt CN4 and
# early-TM Mo/W CN6.
_MULTI_POLY_CNS: Tuple[int, ...] = (3, 5)


def polyhedra_for_cn(cn: int) -> List[str]:
    """Return canonical polyhedron strings for ``cn``.

    Pulls from :data:`delfin.fffree.polyhedra.GEOM_BY_CN` and intersects with
    the polyhedra registered in :data:`_SHAPE_TO_GROUP` (so callers don't get
    polyhedra without a Pólya group attached).  Y-3 (Y-shape) is appended to
    the CN3 list when the env-flag is on -- it's not in ``GEOM_BY_CN`` but
    shares the C2 proper-rotation group with T-3 and is a real coordination
    geometry for d8 / d10 CN3 (e.g. Au(I) CN3 with one wide donor angle).

    Deterministic: returns the polyhedra in the order they're registered in
    ``GEOM_BY_CN``, with Y-3 appended last for CN3 (so the existing CN3
    behavior is preserved when the caller iterates in order).
    """
    try:
        from delfin.fffree.polyhedra import GEOM_BY_CN
    except Exception:
        return []
    cands = list(GEOM_BY_CN.get(int(cn), []))
    # Keep only polyhedra with a registered Pólya group.
    out: List[str] = [p for p in cands if p in _SHAPE_TO_GROUP]
    # CN3: also expose Y-3 (Y-shape), which shares the T-3 proper-rotation
    # group but is geometrically distinct.  Append last so existing single-
    # polyhedron callers get the textbook behavior unchanged.
    if int(cn) == 3 and "Y-3 Y-shape" in _SHAPE_TO_GROUP:
        out.append("Y-3 Y-shape")
    return out


def enumerate_multi_polyhedron_orbits(
    cn: int,
    donor_types: Sequence[str],
    chelate_pairs: Sequence[Tuple[int, int]] = (),
    *,
    polyhedra: Optional[Sequence[str]] = None,
    max_orbits_per_polyhedron: Optional[int] = None,
) -> List[Tuple[str, Tuple[str, ...]]]:
    """Enumerate (polyhedron, orbit) pairs across every polyhedron for ``cn``.

    Parameters
    ----------
    cn : int
        Coordination number (= length of ``donor_types``).
    donor_types : sequence of str
        Length-CN donor-type label list (one per donor slot).
    chelate_pairs : sequence of (i, j) int pairs
        Donor-list-index pairs that must sit on a cis-edge of every
        polyhedron under consideration.  Empty by default.
    polyhedra : sequence of str, optional
        Explicit polyhedron list to enumerate over.  When ``None`` (default),
        :func:`polyhedra_for_cn` provides the list.  Each entry must be a
        canonical string in :data:`_SHAPE_TO_GROUP`.
    max_orbits_per_polyhedron : int, optional
        Cap on the number of orbit representatives per polyhedron.  Defaults
        to ``DELFIN_FFFREE_POLYHEDRON_VERTEX_POLYA_MAX_ORBITS`` (10 if unset).

    Returns
    -------
    list of (polyhedron, orbit)
        One tuple per (polyhedron, orbit) pair.  Sorted lex by
        (polyhedron, orbit) for determinism.  Each ``orbit`` is a length-CN
        tuple of donor-type labels per vertex.

    Notes
    -----
    * Orbit representatives are LOCAL to each polyhedron's vertex set
      (i.e. an orbit for TBP-5 and an orbit for SPY-5 are NOT directly
      comparable as vertex-label tuples even when they look identical --
      the vertex INDEX has a polyhedron-specific meaning, encoded in
      :data:`polyhedra.ref_vectors`).  The polyhedron string in the
      returned tuple disambiguates this.
    * Polyhedra without a Pólya group registration, or with an invalid
      donor-count match, are silently skipped (the function never raises).
    """
    if int(cn) != len(donor_types):
        raise ValueError(
            f"donor_types length {len(donor_types)} != CN {int(cn)}"
        )
    cand_polys = list(polyhedra) if polyhedra is not None else polyhedra_for_cn(int(cn))
    out: List[Tuple[str, Tuple[str, ...]]] = []
    for poly in cand_polys:
        try:
            orbits = enumerate_orbits_with_chelates(
                poly, list(donor_types), list(chelate_pairs),
                max_orbits=max_orbits_per_polyhedron,
            )
        except (KeyError, ValueError):
            continue
        for orbit in orbits:
            out.append((poly, tuple(orbit)))
    # Deterministic ordering: (polyhedron, orbit) lex sort.
    return sorted(out)


# ---------------------------------------------------------------------------
# SMILES → (polyhedron, donor_types, chelate_pairs) bridge.
#
# Graph-only, deterministic, universal.  Uses RDKit ONLY to parse the SMILES
# into a graph; everything else is plain Python.  No SMILES-pattern matching
# (the doctrine: "graph-only, never SMILES-specific shortcuts").
# ---------------------------------------------------------------------------
def _default_polyhedron(cn: int, metal: str) -> Optional[str]:
    """Pick a sensible default coordination polyhedron for a given CN.

    Returns the FIRST entry from :data:`delfin.fffree.polyhedra.GEOM_BY_CN`
    (the same choice the assembler defaults to before constructive
    enumeration).  Falls back to ``None`` for CNs that have no registered
    polyhedron.
    """
    try:
        from delfin.fffree.polyhedra import GEOM_BY_CN
    except Exception:
        return None
    cands = GEOM_BY_CN.get(int(cn), [])
    if not cands:
        return None
    return cands[0]


def _all_polyhedra(cn: int, metal: str) -> List[str]:
    """Return ALL registered polyhedra for ``cn`` (multi-poly dispatch).

    For ``cn in polyhedra._MULTI_POLY_CNS`` (currently 3, 4, 5, 8, 9, 10,
    11, 12): returns every geometry in ``GEOM_BY_CN[cn]`` whose Pólya
    group is registered in :data:`_SHAPE_TO_GROUP` (so orbit enumeration
    can actually run on it).  For other CNs returns the single default.

    Used by the multi-poly path in :func:`enumerate_orbits_for_smiles`
    to emit one ETKDG embed per (polyhedron, orbit) pair, rather than
    a single polyhedron's orbits only.

    Universal, no per-metal lookups.  Deterministic order
    = ``GEOM_BY_CN[cn]`` order.
    """
    try:
        from delfin.fffree.polyhedra import (
            GEOM_BY_CN, is_multi_poly_cn,
        )
    except Exception:
        return []
    cands = list(GEOM_BY_CN.get(int(cn), []))
    if not cands:
        return []
    if not is_multi_poly_cn(int(cn)):
        return cands[:1]
    # Keep only polyhedra whose Pólya group is implemented (so orbit
    # enumeration can run); drop any unimplemented variants silently.
    return [g for g in cands if g in _SHAPE_TO_GROUP] or cands[:1]


def _donor_label(symbol: str) -> str:
    """Canonical donor-type label (element-only).  Kept for callers that
    don't have a graph context handy.  When a graph is available, prefer
    :func:`_donor_label_from_atom` which can distinguish e.g. water-O
    from carboxylate-O via heavy-atom-neighbour count."""
    return str(symbol)


def _donor_label_from_atom(atom, metal_idx: int) -> str:
    """Graph-context-aware donor label.

    Two donors of the same element are still distinguishable as DIFFERENT
    coordination types when their non-metal-non-H neighbour count differs.
    This captures the textbook chemistry:

    * water O (no heavy non-metal neighbour)                  -> ``"O0"``
    * hydroxide / alkoxide / carbonyl O (1 heavy neighbour)   -> ``"O1"``
    * bridging / ether / carboxylate sp2 O (2 heavy neighbours)-> ``"O2"``
    * pyridyl / amine N is treated the same way for N         -> ``"N1"`` ...

    The same rule for any element keeps the function fully universal — it
    encodes "how the donor is anchored in its ligand" without any SMILES-
    pattern lookup.  ``metal_idx`` is excluded from the neighbour count
    (it's the bond we're characterising, not part of the donor's
    "environment").

    Universal, graph-only, deterministic.  Independent of formal charges
    and aromaticity flags (which can be representation-dependent).
    """
    try:
        from delfin._bond_decollapse import _is_metal
    except Exception:
        _is_metal = lambda s: False  # noqa: E731
    n_heavy = 0
    for nbr in atom.GetNeighbors():
        try:
            sym = nbr.GetSymbol()
        except Exception:
            continue
        if int(nbr.GetIdx()) == int(metal_idx):
            continue
        if _is_metal(sym):
            # Bridging via another metal: still counts as a heavy anchor for
            # the donor's coordination environment.
            n_heavy += 1
            continue
        if sym == "H":
            continue
        n_heavy += 1
    return f"{atom.GetSymbol()}{int(n_heavy)}"


def detect_from_smiles(
    smiles: str,
) -> Optional[Tuple[str, List[str], List[Tuple[int, int]], int]]:
    """Parse ``smiles`` and extract (polyhedron, donor_types, chelate_pairs, metal_idx).

    Graph-only: identifies the metal atom (first metal in the molecule),
    enumerates its direct neighbours as donors (donor_types is a per-donor
    list, ordered by RDKit atom index), and identifies chelate pairs as
    donor pairs whose two atoms are connected by a short non-metal path
    (length ≤ 4 bonds) in the metal-removed graph.

    Returns ``None`` when:
    * RDKit can't parse the SMILES,
    * no recognised metal is present,
    * the metal has < 2 donors (no polyhedron sensible),
    * no default polyhedron is registered for the metal's CN.

    The function never raises; failure paths return ``None`` so the caller
    can silently fall through to the legacy ETKDG path.
    """
    try:
        from rdkit import Chem
    except Exception:
        return None
    try:
        from delfin._bond_decollapse import _is_metal
    except Exception:
        return None
    if not smiles:
        return None
    try:
        mol = Chem.MolFromSmiles(smiles)
    except Exception:
        return None
    if mol is None:
        return None
    # First metal atom (deterministic: lowest atom index).
    metal_idx: Optional[int] = None
    metal_sym: str = ""
    for atom in mol.GetAtoms():
        sym = atom.GetSymbol()
        if _is_metal(sym):
            metal_idx = int(atom.GetIdx())
            metal_sym = sym
            break
    if metal_idx is None:
        return None
    metal_atom = mol.GetAtomWithIdx(metal_idx)
    donor_idxs: List[int] = []
    donor_types: List[str] = []
    for nbr in metal_atom.GetNeighbors():
        nsym = nbr.GetSymbol()
        if _is_metal(nsym):
            continue
        donor_idxs.append(int(nbr.GetIdx()))
        donor_types.append(_donor_label_from_atom(nbr, metal_idx))
    if len(donor_idxs) < 2:
        return None
    polyhedron = _default_polyhedron(len(donor_idxs), metal_sym)
    if polyhedron is None:
        return None
    chelate_pairs = _detect_chelate_pairs(mol, metal_idx, donor_idxs)
    return polyhedron, donor_types, chelate_pairs, metal_idx


def _detect_chelate_pairs(
    mol,
    metal_idx: int,
    donor_idxs: Sequence[int],
    *,
    max_path_len: int = 4,
) -> List[Tuple[int, int]]:
    """Return donor-list-index pairs that are chelate partners.

    A pair (i, j) of donors is a chelate iff there exists a path between
    ``donor_idxs[i]`` and ``donor_idxs[j]`` of length ≤ ``max_path_len`` in
    the metal-removed graph.  The path length cap (4) keeps the enumeration
    tight: typical bidentate chelating ligands (en, bpy, acac, glycinate,
    oxalate, pyridyl-carboxylate, …) form 5- or 6-membered metallacycles =
    path length 3 or 4 in the donor-donor backbone.  Longer paths
    (medium-ring chelates) are rare and currently treated as non-chelating
    — they will simply get fewer cis constraints, and the assembler's
    metallacycle embed handles them as a separate code path.

    Pure graph-BFS, deterministic, no RDKit substructure search.

    Returns list of (i, j) integer pairs (indices INTO ``donor_idxs``,
    NOT atom indices), sorted lex.
    """
    try:
        from delfin._bond_decollapse import _is_metal
    except Exception:
        return []
    n_donors = len(donor_idxs)
    if n_donors < 2:
        return []
    # Build adjacency excluding any metal atom.  We exclude ALL metals so
    # multi-metal cages still produce sensible chelate detection per metal.
    metal_atoms = {a.GetIdx() for a in mol.GetAtoms() if _is_metal(a.GetSymbol())}
    n_atoms = mol.GetNumAtoms()
    adj: List[List[int]] = [[] for _ in range(n_atoms)]
    for bond in mol.GetBonds():
        a = int(bond.GetBeginAtomIdx())
        b = int(bond.GetEndAtomIdx())
        if a in metal_atoms or b in metal_atoms:
            continue
        adj[a].append(b)
        adj[b].append(a)
    # BFS from each donor atom up to max_path_len.
    pairs: List[Tuple[int, int]] = []
    donor_set = {idx: i for i, idx in enumerate(donor_idxs)}
    for i, start in enumerate(donor_idxs):
        # BFS up to max_path_len; record any donor reached with a shorter
        # path.
        dist = {start: 0}
        frontier = [start]
        while frontier:
            new_frontier: List[int] = []
            for u in frontier:
                d = dist[u]
                if d >= max_path_len:
                    continue
                for v in adj[u]:
                    if v in dist:
                        continue
                    dist[v] = d + 1
                    new_frontier.append(v)
            frontier = new_frontier
        for atom_idx, d in dist.items():
            if atom_idx == start:
                continue
            if 0 < d <= max_path_len and atom_idx in donor_set:
                j = donor_set[atom_idx]
                if i < j:
                    pairs.append((i, j))
    return sorted(set(pairs))


# ---------------------------------------------------------------------------
# Top-level convenience.
# ---------------------------------------------------------------------------
def flag_active() -> bool:
    """True when the env flag is set to enable polyhedron-vertex enumeration."""
    return os.environ.get(
        "DELFIN_FFFREE_POLYHEDRON_VERTEX_POLYA", "0"
    ).strip().lower() in ("1", "true", "yes", "on")


def enumerate_orbits_for_smiles(
    smiles: str,
    *,
    max_orbits: Optional[int] = None,
) -> Optional[Dict[str, object]]:
    """Top-level: parse ``smiles`` + emit orbit list ready for ETKDG dispatch.

    Returns ``None`` when the SMILES doesn't carry an enumerable metal
    polyhedron (organic, no metal, CN 0/1, unknown polyhedron, RDKit
    parse failure).  On success returns a dict ::

        {
            "polyhedron":          str,                  # canonical geometry string (primary)
            "metal_idx":           int,                  # RDKit atom index of metal
            "donor_atoms":         List[int],            # RDKit atom indices, in order
            "donor_types":         List[str],            # length CN, type label per donor
            "chelate_pairs":       List[Tuple[int,int]], # indices INTO donor_atoms
            "orbits":              List[Tuple[str,...]], # canonical reps for ``polyhedron`` only
            "polyhedron_orbits":   List[{"polyhedron": str, "orbits": List[...]}],
                                                        # multi-polyhedron breakdown (CN3/CN5)
        }

    The caller dispatches one ETKDG embed (with the corresponding CoordMap)
    per (polyhedron, orbit) pair.  Default-OFF: callers must check
    :func:`flag_active` before invoking this — the function itself does NOT
    consult the env-flag, so it remains test-callable without side-effects.

    CN3 / CN5 multi-polyhedron behaviour (2026-06-07 extension)
    -----------------------------------------------------------
    For CN3 and CN5, ``polyhedron_orbits`` enumerates EVERY polyhedron
    registered in :func:`polyhedra_for_cn` (e.g. TBP-5 AND SPY-5 for CN5,
    SP-3 + T-3 + Y-3 for CN3), so callers can dispatch one ETKDG seed per
    (polyhedron, orbit) and cover the textbook geometric-isomer count for
    ambivalent metals (Cu(II) CN5, Au(I) CN3, ...).  The legacy ``orbits``
    field stays populated for the *primary* polyhedron only -- back-compat
    for callers that don't yet consume ``polyhedron_orbits``.

    For other CNs (2, 4, 6+), ``polyhedron_orbits`` carries a single-entry
    list containing exactly the primary polyhedron's orbits -- behaviour is
    bit-exact to the pre-extension code.
    """
    det = detect_from_smiles(smiles)
    if det is None:
        return None
    polyhedron, donor_types, chelate_pairs, metal_idx = det
    # Re-extract donor atom indices in deterministic order (same order used
    # by detect_from_smiles -- by RDKit atom index of the metal's neighbours).
    try:
        from rdkit import Chem
        from delfin._bond_decollapse import _is_metal
    except Exception:
        return None
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    donor_atoms = [
        int(nbr.GetIdx())
        for nbr in mol.GetAtomWithIdx(metal_idx).GetNeighbors()
        if not _is_metal(nbr.GetSymbol())
    ]
    if len(donor_atoms) != len(donor_types):
        return None
    try:
        orbits = enumerate_orbits_with_chelates(
            polyhedron, donor_types, chelate_pairs, max_orbits=max_orbits,
        )
    except (KeyError, ValueError):
        return None
    if not orbits:
        return None

    # Multi-polyhedron orbit breakdown for CN3 / CN5 (ambivalent metals).
    # For other CNs the breakdown carries a single entry = the primary
    # polyhedron, preserving bit-exact pre-extension output.
    cn = len(donor_types)
    polyhedron_orbits: List[Dict[str, object]] = []
    if cn in _MULTI_POLY_CNS:
        for poly in polyhedra_for_cn(cn):
            try:
                poly_orbits = enumerate_orbits_with_chelates(
                    poly, donor_types, chelate_pairs, max_orbits=max_orbits,
                )
            except (KeyError, ValueError):
                continue
            if not poly_orbits:
                continue
            polyhedron_orbits.append({
                "polyhedron": poly,
                "orbits": poly_orbits,
            })
    if not polyhedron_orbits:
        # Single-polyhedron fallback (preserves pre-extension behaviour for
        # every CN not in _MULTI_POLY_CNS, AND for CN3/CN5 cases where the
        # multi-polyhedron pass yielded nothing extra).
        polyhedron_orbits = [{
            "polyhedron": polyhedron,
            "orbits": orbits,
        }]

    return {
        "polyhedron": polyhedron,
        "metal_idx": int(metal_idx),
        "donor_atoms": list(donor_atoms),
        "donor_types": list(donor_types),
        "chelate_pairs": list(chelate_pairs),
        "orbits": orbits,
        "polyhedron_orbits": polyhedron_orbits,
    }


def enumerate_orbits_for_smiles_multi(
    smiles: str,
    *,
    max_orbits: Optional[int] = None,
) -> Optional[List[Dict[str, object]]]:
    """Multi-polyhedron variant of :func:`enumerate_orbits_for_smiles`.

    Returns a LIST of result dicts (one per polyhedron) when the detected
    CN is in :data:`polyhedra._MULTI_POLY_CNS` (3, 4, 5, 8, 9, 10, 11, 12).
    For other CNs returns a single-element list (the legacy behaviour
    wrapped in a list, so the caller's loop body is uniform).

    Each dict has the same shape as the single-polyhedron result:
    ``polyhedron``, ``metal_idx``, ``donor_atoms``, ``donor_types``,
    ``chelate_pairs``, ``orbits``.

    ZURMAA universal fix (2026-06-07): the ETKDG-embed-orbit dispatcher
    can now emit BOTH T-4 and SP-4 frames for a CN4 SMILES like
    ``O=C(O)c1ccc(N(C)C(=S)[Au]2(Nc1cccs1)c1ccccc1)cc1``, instead of
    only the metal-table-default geometry.  Mogul-DG severity ranks the
    polyhedra downstream.

    Returns ``None`` on parse / detection failure (same contract as the
    single-polyhedron function).  Never raises.
    """
    det = detect_from_smiles(smiles)
    if det is None:
        return None
    _, donor_types, chelate_pairs, metal_idx = det
    cn = len(donor_types)
    try:
        from rdkit import Chem
        from delfin._bond_decollapse import _is_metal
    except Exception:
        return None
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    try:
        metal_sym = mol.GetAtomWithIdx(metal_idx).GetSymbol()
    except Exception:
        metal_sym = ""
    donor_atoms = [
        int(nbr.GetIdx())
        for nbr in mol.GetAtomWithIdx(metal_idx).GetNeighbors()
        if not _is_metal(nbr.GetSymbol())
    ]
    if len(donor_atoms) != len(donor_types):
        return None
    polyhedra_list = _all_polyhedra(cn, metal_sym)
    if not polyhedra_list:
        return None
    out: List[Dict[str, object]] = []
    for polyhedron in polyhedra_list:
        try:
            orbits = enumerate_orbits_with_chelates(
                polyhedron, donor_types, chelate_pairs, max_orbits=max_orbits,
            )
        except (KeyError, ValueError):
            continue
        if not orbits:
            continue
        out.append({
            "polyhedron": polyhedron,
            "metal_idx": int(metal_idx),
            "donor_atoms": list(donor_atoms),
            "donor_types": list(donor_types),
            "chelate_pairs": list(chelate_pairs),
            "orbits": orbits,
        })
    if not out:
        return None
    return out


def vertex_coords(polyhedron: str, metal: str = "Cu") -> Optional["np.ndarray"]:
    """Ideal vertex coordinates for ``polyhedron``, scaled to ~ Cu-O distance.

    Returns the unit-radius vertex array from
    :mod:`delfin.fffree.polyhedra.ref_vectors`, scaled by a typical M-D bond
    length so the caller can pass directly into an RDKit ``coordMap``.

    The scale is intentionally coarse: ETKDG will refine all internal
    distances, and we only need the donor atoms placed in the correct
    cis / trans pattern.  Scaling uses :func:`md_distance` with a generic
    ``"O"`` donor as the reference (≈ 1.9 Å for first-row TM); the absolute
    scale doesn't affect orbit topology.

    Y-3 fallback: when ``polyhedron`` is ``"Y-3 Y-shape"``, fall back to the
    T-3 vertex array — Y-3 shares the same proper-rotation group (C2) and is
    a continuous deformation of T-3 around the metal; the ETKDG / relaxation
    phase refines the angular geometry from the placed skeleton.  This keeps
    :mod:`polyhedra` as the single source of truth for vertex tables.
    """
    try:
        import numpy as np
        from delfin.fffree import polyhedra as _PLY
    except Exception:
        return None
    try:
        shape = _canonical_shape(polyhedron)
    except KeyError:
        return None
    try:
        V = _PLY.ref_vectors(shape)
    except KeyError:
        # Y-3 graceful fallback: reuse T-3 placement skeleton (same group,
        # same orbit structure; angular relaxation is downstream).
        if shape == "Y-3 Y-shape":
            try:
                V = _PLY.ref_vectors("T-3 T-shape")
            except KeyError:
                return None
        else:
            return None
    try:
        # Reference M-D ≈ 1.9 Å (Cu-O typical); used uniformly across all
        # vertices because orbit topology doesn't depend on exact distance.
        d = _PLY.md_distance(metal, "O")
    except Exception:
        d = 1.9
    return np.asarray(V, dtype=float) * float(d)


# ---------------------------------------------------------------------------
# Public exports
# ---------------------------------------------------------------------------
__all__ = [
    "enumerate_orbits",
    "enumerate_orbits_with_chelates",
    "enumerate_orbits_for_smiles",
    "enumerate_orbits_for_smiles_multi",
    "detect_from_smiles",
    "cis_edges_for",
    "cn_for",
    "flag_active",
    "vertex_coords",
]
