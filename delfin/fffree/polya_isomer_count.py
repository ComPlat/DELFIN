#!/usr/bin/env python3
"""polya_isomer_count.py — EXACT stereoisomer count + enumeration
for a coordination polyhedron given a donor-label multiset, via Burnside's lemma
over the proper rotation group of the polyhedron.

This is the PROVABLE completeness reference (the denominator for Layer-2 coverage):
the generator must produce exactly this many distinct geometric isomers, no more,
no fewer.  Counts enantiomers as distinct (proper-rotation group), the standard
"all stereoisomers" convention; achiral merging is a separate post-step.

Self-tested against textbook values (octahedral MA4B2=2 cis/trans, MA3B3=2 fac/mer,
MA2B2C2=6, MABCDEF=30; square-planar MA2B2=2; tetrahedral MABCD=2).
"""
from __future__ import annotations
import itertools
from collections import Counter
from typing import Dict, List, Optional, Tuple


def _close_group(generators: List[Tuple[int, ...]], n: int) -> List[Tuple[int, ...]]:
    """Close a set of vertex permutations into the full group they generate."""
    identity = tuple(range(n))
    group = {identity}
    frontier = [identity]
    gens = [tuple(g) for g in generators]
    while frontier:
        a = frontier.pop()
        for g in gens:
            b = tuple(g[a[i]] for i in range(n))   # compose
            if b not in group:
                group.add(b); frontier.append(b)
    return sorted(group)


def _linear_group():
    # Phase G: CN2 L-2 linear: 2 vertices 180° apart on the metal axis.
    # Proper rotation group C2 (order 2): identity + C2 swap of 2 trans donors.
    # Cu(I)/Ag(I)/Au(I)/Hg(II) linear coordination.
    C2 = (1, 0)
    return _close_group([C2], 2)


def _trigonal_planar_group():
    # CN3 SP-3 trigonal planar: 3 vertices in a plane at 120°.  Proper rotation group
    # D3 (order 6): C3 about perp axis + 3 C2 through each vertex (and midpoint of
    # opposite edge in the plane).  Iter-32c (User 2026-05-28 ADUMOD).
    C3 = (1, 2, 0)
    C2 = (0, 2, 1)            # C2 through vertex 0
    return _close_group([C3, C2], 3)


def _tshape_group():
    # CN3 T-3 T-shape: 2 trans + 1 cis (e.g., d8 Pt CN3).  Proper rotation group C2
    # (order 2): C2 swaps the two trans vertices, fixes the cis one.
    C2 = (1, 0, 2)
    return _close_group([C2], 3)


def _octahedron_group():
    # vertices 0=+x 1=-x 2=+y 3=-y 4=+z 5=-z
    Rz = (2, 3, 1, 0, 4, 5)   # +x->+y,+y->-x,-x->-y,-y->+x
    Rx = (0, 1, 4, 5, 3, 2)   # +y->+z,+z->-y,-y->-z,-z->+y
    return _close_group([Rz, Rx], 6)


def _tetrahedron_group():
    # 4 vertices of a tetrahedron; proper rotation group T (order 12)
    # generators: a 3-fold (fix vertex 0, cycle 1,2,3) and a 2-fold (swap pairs)
    C3 = (0, 2, 3, 1)
    C2 = (1, 0, 3, 2)
    return _close_group([C3, C2], 4)


def _square_group():
    # 4 vertices around a square (square-planar); proper group D4 (order 8)
    C4 = (1, 2, 3, 0)        # rotation about the perpendicular axis
    C2p = (0, 3, 2, 1)       # in-plane 2-fold (flips the square)
    return _close_group([C4, C2p], 4)


def _tbp_group():
    # CN5 trigonal bipyramid: 0,1 axial; 2,3,4 equatorial; proper group D3 (order 6)
    C3 = (0, 1, 3, 4, 2)     # rotate equatorial, axials fixed
    C2 = (1, 0, 2, 4, 3)     # swap axials + flip two equatorial
    return _close_group([C3, C2], 5)


def _spy_group():
    # CN5 square pyramid: 0 apical; 1,2,3,4 basal; proper group C4 (order 4)
    C4 = (0, 2, 3, 4, 1)
    return _close_group([C4], 5)


def _trigonal_prism_group():
    # CN6 TPR (matches polyhedra CN6 index order): 0-2 top triangle at z=+0.7 (0,120,240
    # deg); 3-5 bottom triangle at z=-0.7 (eclipsed: 0,120,240 deg).  Proper rotation
    # group of an eclipsed trigonal prism is D3 (order 6): C3 about z + three C2 through
    # rectangular-face midpoints.  Iter-31 (User 2026-05-28: early-TM Mo/W CN6 prefer
    # TPR over OC — coverage gap previously missed).
    C3 = (1, 2, 0, 4, 5, 3)           # C3 about z: top triangle 0->1->2->0, bottom 3->4->5->3
    C2 = (3, 5, 4, 0, 2, 1)           # C2 about x-axis (through midpoints of 0-3, 1-5, 2-4 edges)
    return _close_group([C3, C2], 6)


def _pentagonal_bipyramid_group():
    # CN7 PB (matches polyhedra CN7 index order): 0,1 axial (+z,-z); 2-6 equatorial
    # pentagon (0,72,...,288 deg). Proper rotation group D5 (order 10).
    C5 = (0, 1, 3, 4, 5, 6, 2)        # C5 about z: equator 2->3->4->5->6->2, axials fixed
    C2 = (1, 0, 2, 6, 5, 4, 3)        # C2 through vertex 2: swap axials, fix v2, 3<->6 4<->5
    return _close_group([C5, C2], 7)


def _square_antiprism_group():
    # CN8 SQAP (matches polyhedra CN8 index order): 0-3 top square (0,90,180,270 deg),
    # 4-7 bottom square (45,135,225,315 deg). Proper rotation group D4 (order 8).
    C4 = (1, 2, 3, 0, 5, 6, 7, 4)     # C4 about z: each square cycles
    C2 = (4, 7, 6, 5, 0, 3, 2, 1)     # C2 (axis at 22.5 deg) swaps top<->bottom
    return _close_group([C4, C2], 8)


def _tricapped_trigonal_prism_group():
    # CN9 TTP (matches polyhedra CN9 index order): 0-2 top triangle, 3-5 bottom triangle
    # (eclipsed prism), 6-8 caps (60,180,300 deg). Proper rotation group D3 (order 6).
    C3 = (1, 2, 0, 4, 5, 3, 7, 8, 6)  # C3 about z: triangles + caps cycle
    C2 = (4, 3, 5, 1, 0, 2, 6, 8, 7)  # C2 through cap6: swap top<->bottom, fix cap6, 7<->8
    return _close_group([C3, C2], 9)


def _sandwich_10_group():
    """SANDWICH-10 (D5d staggered ferrocene): vertices 0-4 top Cp ring, 5-9 bottom
    Cp ring (staggered 36°).  Proper rotation group D5 (order 10):
      C5 about z cycles top 0-1-2-3-4 and bottom 5-6-7-8-9 simultaneously.
      C2 perpendicular to z swaps top<->bottom (and reverses each ring).
    """
    C5 = (1, 2, 3, 4, 0, 6, 7, 8, 9, 5)
    # C2 through the midpoint of vertex 0 and vertex 5: swap top<->bottom while
    # keeping the antiprismatic pairing.  Vertex 0 (top, azimuth 0°) maps to
    # vertex 5 (bottom, azimuth 36°) and so on around the antiprism.
    C2 = (5, 9, 8, 7, 6, 0, 4, 3, 2, 1)
    return _close_group([C5, C2], 10)


def _piano_stool_8_group():
    """PIANO-STOOL-8 (C3v) — η⁵-Cp + 3 σ-donors.  The strict proper rotation group
    is the identity (order 1): C3 (about z) cycles the σ-tripod (5,6,7) and
    arrests 1/3-rotations of the pentagonal Cp ring (0-4) which is NOT a
    symmetry of the pentagon → only identity preserves the WHOLE polyhedron.

    Choosing identity-only is the safest default: it slightly over-counts
    vs the true geometric isomer count (different σ-donor arrangements that
    are tripod-C3-equivalent stay separate), but the downstream 3D-RMSD
    dedup gate prunes the duplicates.  Tightening to C3 is a future iter.
    """
    return [tuple(range(8))]


def _half_sandwich_9_group():
    """HALF-SANDWICH-9 (C3v) — η⁶-arene + 3 σ-donors.  C3 is a real symmetry
    (the hexagonal arene IS C3v compatible: rotation by 120° cycles vertices
    0->2->4 and 1->3->5).
    Group: C3 (order 3).
      C3 about z: arene 0-2-4 cycle, 1-3-5 cycle, tripod 6-7-8 cycle.
    """
    C3 = (2, 3, 4, 5, 0, 1, 7, 8, 6)
    return _close_group([C3], 9)


_GROUPS = {
    "linear": (_linear_group(), 2),                  # Phase G CN2
    "trigonal_planar": (_trigonal_planar_group(), 3),
    "tshape": (_tshape_group(), 3),
    "octahedron": (_octahedron_group(), 6),
    "square_planar": (_square_group(), 4),
    "tetrahedron": (_tetrahedron_group(), 4),
    "trigonal_bipyramid": (_tbp_group(), 5),
    "square_pyramid": (_spy_group(), 4 + 1),
    "trigonal_prism": (_trigonal_prism_group(), 6),
    "pentagonal_bipyramid": (_pentagonal_bipyramid_group(), 7),
    "square_antiprism": (_square_antiprism_group(), 8),
    "tricapped_trigonal_prism": (_tricapped_trigonal_prism_group(), 9),
    # Task #44 / Mission A7 (2026-06-05): sandwich + piano-stool + half-sandwich.
    "sandwich_10": (_sandwich_10_group(), 10),
    "piano_stool_8": (_piano_stool_8_group(), 8),
    "half_sandwich_9": (_half_sandwich_9_group(), 9),
}


# Phase C (Task #64, 2026-06-04): lazily merge f-block CN8-12 groups from
# delfin.fffree.f_block_polyhedra when they are referenced.  Done lazily to
# avoid a circular import (f_block_polyhedra → polyhedra → polya_isomer_count
# → f_block_polyhedra would loop).
def _fblock_group_for(key: str) -> Optional[Tuple[Tuple[Tuple[int, ...], ...], int]]:
    """Return ``(group, n)`` for a Phase-C f-block polyhedron key or None."""
    try:
        from delfin.fffree import f_block_polyhedra as _FBP
    except ImportError:
        return None
    # Map Pólya group keys to f-block canonical names.
    mapping = {
        "sap_8": "SAP-8 square antiprism",
        "dd_8": "DD-8 dodecahedron",
        "tdd_8": "TDD-8 triangular dodecahedron",
        "ttp_9_fblock": "TTP-9 tricapped trigonal prism",  # alias avoids clash
        "csap_9": "CSAP-9 capped square antiprism",
        "bicap_10": "BICAP-10 bicapped square antiprism",
        "ih_12": "IH-12 icosahedron",
    }
    canonical = mapping.get(key)
    if canonical is None:
        return None
    res = _FBP.fblock_group(canonical)
    if res is None:
        return None
    return res


def _get_group(geometry_key: str) -> Tuple[List[Tuple[int, ...]], int]:
    """Look up a Pólya group, falling back to the f-block table."""
    if geometry_key in _GROUPS:
        return _GROUPS[geometry_key]
    fb = _fblock_group_for(geometry_key)
    if fb is None:
        raise KeyError(geometry_key)
    return fb


def _multiset_from_spec(spec: Dict[str, int], n: int) -> List[str]:
    labels = []
    for lab, cnt in spec.items():
        labels += [lab] * cnt
    if len(labels) != n:
        raise ValueError(f"donor count {len(labels)} != CN {n} for {spec}")
    return labels


def count_isomers(geometry: str, donor_spec: Dict[str, int]) -> int:
    """Burnside: (1/|G|) Σ_g |colorings fixed by g|.  A coloring (assignment of the
    donor multiset to vertices) is fixed by g iff every cycle of g is monochromatic."""
    group, n = _get_group(geometry)
    labels = _multiset_from_spec(donor_spec, n)
    mult = Counter(labels)
    total = 0
    for g in group:
        # cycle decomposition
        seen = [False] * n
        cycles = []
        for i in range(n):
            if seen[i]:
                continue
            c = []; j = i
            while not seen[j]:
                seen[j] = True; c.append(j); j = g[j]
            cycles.append(len(c))
        # count multiset assignments constant on each cycle = multinomial over cycles
        # only feasible if we can partition the multiset into blocks of the cycle sizes
        total += _fixed_colorings(cycles, mult)
    return total // len(group)


def _fixed_colorings(cycle_sizes: List[int], mult: Counter) -> int:
    """# ways to assign labels (with given multiplicities) so each cycle is one label."""
    from functools import lru_cache
    sizes = tuple(sorted(cycle_sizes))
    labels = tuple(sorted(mult.items()))

    @lru_cache(maxsize=None)
    def rec(ci: int, remaining: Tuple[Tuple[str, int], ...]) -> int:
        if ci == len(sizes):
            return 1 if all(c == 0 for _, c in remaining) else 0
        sz = sizes[ci]; tot = 0
        rem = list(remaining)
        for k in range(len(rem)):
            lab, c = rem[k]
            if c >= sz:
                rem2 = list(rem); rem2[k] = (lab, c - sz)
                tot += rec(ci + 1, tuple(rem2))
        return tot
    return rec(0, labels)


def enumerate_isomers(geometry: str, donor_spec: Dict[str, int]) -> List[Tuple[str, ...]]:
    """Distinct vertex-colorings up to the rotation group (one representative each)."""
    group, n = _get_group(geometry)
    labels = _multiset_from_spec(donor_spec, n)
    seen = set(); reps = []
    for perm in set(itertools.permutations(labels)):
        orbit = min(tuple(perm[g[i]] for i in range(n)) for g in group)
        if orbit not in seen:
            seen.add(orbit); reps.append(orbit)
    return reps


_ANTIPODE_FULL = {
    "octahedron": {0: 1, 1: 0, 2: 3, 3: 2, 4: 5, 5: 4},
    "square_planar": {0: 2, 1: 3, 2: 0, 3: 1},
}

# Max ideal vertex-vertex angle (deg) for a pair to host a bidentate chelate ("cis"
# edge).  Larger separations are trans-like and not chelate-feasible.  100deg:
#  - reproduces octahedron/square_planar EXACTLY (their only non-90deg pairs are 180deg),
#  - admits the ~90deg chelate edges of TBP-5 (axial-equatorial) and SPY-5 (apical-basal
#    + adjacent-basal), while excluding TBP equatorial-equatorial (120deg) and SPY basal
#    diagonals (164deg) — those are not real chelate positions, and admitting them would
#    feed self-gate-rejected configs to the all-or-nothing chelate gate (one bad isomer
#    bails the whole complex to legacy).
CHELATE_CIS_MAX_DEG = 115.0   # Phase G6 (2026-05-31): raised from 100 → 115
# to admit tetrahedron T-4 chelates (all T-4 angles = 109.5°, previously
# rejected). Universal fix: works for OC-6 (90°), SP-4 (90°), TBP-5 (90/120°),
# T-4 (109.5°). Octahedron unchanged (90° still well under 115°).

_GEOM_KEY_TO_SHAPE = {
    "octahedron": "OC-6 octahedron",
    "square_planar": "SP-4 square planar",
    "tetrahedron": "T-4 tetrahedron",
    "trigonal_bipyramid": "TBP-5 trigonal bipyramid",
    "square_pyramid": "SPY-5 square pyramid",
    "pentagonal_bipyramid": "PB-7 pentagonal bipyramid",
    "square_antiprism": "SQAP-8 square antiprism",
    "tricapped_trigonal_prism": "TTP-9 tricapped trigonal prism",
    # Phase C f-block CN8-12 (Task #64).
    "sap_8": "SAP-8 square antiprism",
    "dd_8": "DD-8 dodecahedron",
    "tdd_8": "TDD-8 triangular dodecahedron",
    "ttp_9_fblock": "TTP-9 tricapped trigonal prism",
    "csap_9": "CSAP-9 capped square antiprism",
    "bicap_10": "BICAP-10 bicapped square antiprism",
    "ih_12": "IH-12 icosahedron",
    # Task #44 / Mission A7 (2026-06-05): sandwich + piano-stool + half-sandwich.
    "sandwich_10": "SANDWICH-10 bis-eta5-Cp",
    "piano_stool_8": "PIANO-STOOL-8 eta5-Cp+L3",
    "half_sandwich_9": "HALF-SANDWICH-9 eta6+L3",
}

# POLYA-COVERAGE-FIX-v1 (env-gated additions to _GEOM_KEY_TO_SHAPE).
# Consulted ONLY when ``DELFIN_FFFREE_POLYA_COVERAGE_FIX_v1=1``; default-OFF
# preserves HEAD bcf56f8's KeyError-into-_ANTIPODE_FULL path so existing
# ``enumerate_chelate_configs("trigonal_prism", ...)`` callers (e.g.
# converter_backend._enumerate_geometry's TPR-6 branch, line ~979) keep
# returning ``None`` exactly as today.  With the flag flipped, TPR-6 gets a
# real cis-edge map derived from polyhedra.ref_vectors and chelate enumeration
# stops being silently dropped — closing the 26 ``bis-bidentate-NN+NN``
# 0%-coverage cases identified in b00f9a0 forensik.
_GEOM_KEY_TO_SHAPE_FIX_V1: Dict[str, str] = {
    "trigonal_prism": "TPR-6 trigonal prism",
}


def _chelate_cis_edges(geometry: str, n: int):
    """Vertex pairs a bidentate chelate can span (cis edges), derived geometrically from
    the ideal polyhedron: pairs whose ideal angular separation is below
    CHELATE_CIS_MAX_DEG.  Single source of truth = the SAME vertex set the assembler
    places into (delfin.fffree.polyhedra), so enumeration and placement stay consistent
    for every geometry (incl. TBP-5/SPY-5).  Byte-identical to the old antipode table for
    octahedron/square_planar.  Falls back to the antipode table if no reference vectors."""
    import os as _os
    shape = _GEOM_KEY_TO_SHAPE.get(geometry)
    # POLYA-COVERAGE-FIX-v1: extend the shape lookup with v1 additions (e.g.
    # TPR-6) ONLY when the env-flag is on.  Default OFF -> HEAD bcf56f8 behaviour
    # (KeyError into _ANTIPODE_FULL fallback) is preserved bit-exact.
    if shape is None and _os.environ.get("DELFIN_FFFREE_POLYA_COVERAGE_FIX_v1",
                                          "") not in ("", "0"):
        shape = _GEOM_KEY_TO_SHAPE_FIX_V1.get(geometry)
    if shape is not None:
        try:
            import math
            import numpy as np
            from delfin.fffree import polyhedra as _PLY
            V = _PLY.ref_vectors(shape)
            cos_max = math.cos(math.radians(CHELATE_CIS_MAX_DEG))
            return [(i, j) for i in range(n) for j in range(i + 1, n)
                    if float(np.clip(V[i] @ V[j], -1.0, 1.0)) > cos_max]
        except Exception:
            pass
        # Task #44 / Mission A7 fallback: sandwich/piano-stool/half-sandwich
        # polyhedra are NOT in the legacy ``polyhedra.ref_vectors`` catalogue
        # — they live in the standalone module.  Use the standalone vertex
        # set for the cis-edge geometry derivation.
        try:
            import math
            import numpy as np
            from delfin.fffree import sandwich_piano_polyhedra as _SP
            if _SP.is_sandwich_geometry(shape):
                V = _SP.ref_vectors_sandwich(shape)
                cos_max = math.cos(math.radians(CHELATE_CIS_MAX_DEG))
                return [(i, j) for i in range(n) for j in range(i + 1, n)
                        if float(np.clip(V[i] @ V[j], -1.0, 1.0)) > cos_max]
        except Exception:
            pass
    if geometry not in _ANTIPODE_FULL:
        # No antipode table for this geometry (Task #44 sandwich/piano-stool/
        # half-sandwich cases reach here when the standalone-module path also
        # didn't yield a result, e.g. for higher-denticity hapto rings whose
        # placement does not flow through cis-edges anyway).  Return an empty
        # cis-edge list: bidentate-arm enumeration won't fire, hapto-ring
        # ligands are placed by the assembler's hapto-π branch.
        return []
    anti = _ANTIPODE_FULL[geometry]
    return [(i, j) for i in range(n) for j in range(i + 1, n) if anti[i] != j]


def enumerate_chelate_configs(geometry: str, ligand_specs):
    """Universal isomer enumeration for a mix of chelating (bidentate) + monodentate
    ligands.  ligand_specs: one dict per ligand instance with keys ``type``
    (identity label), ``denticity`` (1 or 2), ``asym`` (bool: True if a bidentate's
    two donor arms are distinguishable).  Returns distinct configs; each maps
    vertex_index -> (ligand_instance_index, arm_index), one per orbit under the
    polyhedron's proper rotation group.  Bidentate ligands occupy cis-edges."""
    group, n = _get_group(geometry)
    cis_edges = _chelate_cis_edges(geometry, n)
    order = sorted(range(len(ligand_specs)),
                   key=lambda k: -ligand_specs[k]["denticity"])
    seen = set(); out = []

    def canon_key(assign):
        best = None
        for g in group:
            units = {}
            for v, (li, arm) in assign.items():
                units.setdefault(li, []).append((g[v], arm))
            key_units = []
            for li, vs in units.items():
                spec = ligand_specs[li]
                if spec["denticity"] == 1:
                    key_units.append((spec["type"], (vs[0][0],)))
                elif spec.get("asym"):
                    ordered = tuple(v for v, a in sorted(vs, key=lambda x: x[1]))
                    key_units.append((spec["type"] + "*", ordered))
                else:
                    key_units.append((spec["type"], tuple(sorted(v for v, a in vs))))
            k = tuple(sorted(key_units))
            if best is None or k < best:
                best = k
        return best

    def place(assign):
        placed = {li for (li, _) in assign.values()}
        rem = [k for k in order if k not in placed]
        if not rem:
            key = canon_key(assign)
            if key not in seen:
                seen.add(key); out.append(dict(assign))
            return
        k = rem[0]; spec = ligand_specs[k]
        dent = spec["denticity"]
        if dent == 1:
            for v in [v for v in range(n) if v not in assign]:
                a = dict(assign); a[v] = (k, 0); place(a)
        elif dent == 2:
            for (v1, v2) in cis_edges:
                if v1 in assign or v2 in assign:
                    continue
                a = dict(assign); a[v1] = (k, 0); a[v2] = (k, 1); place(a)
                if spec.get("asym"):
                    b = dict(assign); b[v1] = (k, 1); b[v2] = (k, 0); place(b)
        else:
            # tridentate+ (kappa>=3): occupy any d-vertex subset.  Assembly seats the
            # ligand on the matching mer/fac arrangement (metallacycle embed + best-
            # permutation Kabsch) and the self-gate prunes geometrically infeasible
            # subsets, so the combinatorial enumeration need not know mer-vs-fac.
            free = [v for v in range(n) if v not in assign]
            for combo in itertools.combinations(free, dent):
                a = dict(assign)
                for arm, v in enumerate(combo):
                    a[v] = (k, arm)
                place(a)

    place({})
    return out


def count_chelate_isomers(geometry: str, n_chelate: int) -> int:
    """Count distinct ways to place n_chelate identical symmetric bidentate
    ligands on cis-edges (+ identical monodentate on the rest), up to the proper
    rotation group.  Chelates occupy EDGES (cis vertex pairs), not vertices, so
    this is edge-combinatorics (counts enantiomers as distinct: e.g. octahedral
    en2X2 -> trans + cis-Δ + cis-Λ = 3).

    Default (env-OFF): historical behaviour — only ``octahedron`` supported;
    every other geometry raises ``NotImplementedError`` for bit-exact parity
    with HEAD bcf56f8.

    Env DELFIN_FFFREE_POLYA_COVERAGE_FIX_v1=1 (or the ``universal=True`` kwarg):
    falls back to ``count_chelate_isomers_universal`` for any geometry that has
    a Pólya group + geometric cis-edge map registered, including SP-4, T-4,
    TBP-5, SPY-5, TPR-6, PBP-7, SAP/DD-8, TTP-9 and the f-block CN8-12 set.
    The octahedron branch keeps its original edge-combinatorial algorithm
    byte-identical when env is unset.
    """
    import os as _os
    universal = bool(_os.environ.get("DELFIN_FFFREE_POLYA_COVERAGE_FIX_v1", "")
                     not in ("", "0"))
    group, n = _get_group(geometry)
    if geometry == "octahedron":
        antipode = {0: 1, 1: 0, 2: 3, 3: 2, 4: 5, 5: 4}
    elif universal:
        return count_chelate_isomers_universal(
            geometry, n_chelate, n_monodentate_fillers=n - 2 * n_chelate,
            asymmetric=False,
        )
    else:
        raise NotImplementedError(geometry)
    cis = [(i, j) for i in range(n) for j in range(i + 1, n) if antipode[i] != j]
    seen = set(); count = 0
    # choose n_chelate disjoint cis-edges
    for edges in itertools.combinations(cis, n_chelate):
        used = [v for e in edges for v in e]
        if len(set(used)) != 2 * n_chelate:
            continue
        config = frozenset(frozenset(e) for e in edges)
        # canonical orbit under the group
        orbit_min = None
        for g in group:
            mapped = frozenset(frozenset((g[a], g[b])) for a, b in edges)
            key = tuple(sorted(tuple(sorted(e)) for e in mapped))
            if orbit_min is None or key < orbit_min:
                orbit_min = key
        if orbit_min not in seen:
            seen.add(orbit_min); count += 1
    return count


# =====================================================================
# POLYA-COVERAGE-FIX-v1 (env DELFIN_FFFREE_POLYA_COVERAGE_FIX_v1=1)
# ---------------------------------------------------------------------
# Forensik (b00f9a0 voll-pool isocoverage, 2026-06-04, n=10647) showed that
# 1746 SMILES (16.4%) reach 0% coverage with theoretical >= 2 isomers — a
# pure builder gap. The biggest single contributor is *bidentate chelate
# enumeration for non-octahedral polyhedra*, where ``count_chelate_isomers``
# raised ``NotImplementedError`` and the downstream pipeline silently
# defaulted to "1 isomer" for SP-4 / T-4 / TBP-5 / SPY-5 / TPR-6 / PBP-7 /
# SAP-8 / DD-8 / TTP-9.  The universal helper below treats chelate counting
# as a generalisation of the existing ``enumerate_chelate_configs`` machinery
# (which already drives the FF-free builder via converter_backend / grip_ensemble)
# so the *theoretical reference* and the *builder* now use the same algorithm.
#
# Determinism: pure-Python, no RNG, no hash-dependent dict iteration on the
# critical path — input ``ligand_specs`` is sorted by stable key before the
# orbit canonicalisation pass.  Bit-exact across runs with PYTHONHASHSEED=0
# and idempotent across machines (verified by ``tests/test_polya_coverage_fix_v1.py``).
# =====================================================================

def count_chelate_isomers_universal(
    geometry: str,
    n_chelate: int,
    n_monodentate_fillers: int = 0,
    asymmetric: bool = False,
    monodentate_labels: Optional[Tuple[str, ...]] = None,
) -> int:
    """Universal chelate-isomer counter — generalisation of the octahedron
    edge-combinatorial form above to *every* polyhedron registered in
    ``_GROUPS`` (or the lazy f-block fallback).

    Parameters
    ----------
    geometry
        Polyhedron key, e.g. ``"square_planar"``, ``"tetrahedron"``,
        ``"trigonal_bipyramid"``, ``"square_pyramid"``, ``"octahedron"``,
        ``"trigonal_prism"``, ``"pentagonal_bipyramid"``,
        ``"square_antiprism"``, ``"tricapped_trigonal_prism"`` or any
        f-block CN8-12 key (``"sap_8"``, ``"dd_8"`` ...).
    n_chelate
        Number of identical bidentate chelates (each occupies one cis-edge).
    n_monodentate_fillers
        Remaining vertices to be filled with monodentate ligands.  Defaults
        to ``0`` (homochelate, ``n_chelate * 2 == n``).
    asymmetric
        ``True`` if the two arms of each bidentate are distinguishable
        (e.g. pyridyl-amine bipy-fragment).  Doubles the orbit count where
        the arm-swap symmetry would have collapsed.
    monodentate_labels
        Optional tuple of length ``n_monodentate_fillers`` giving distinct
        labels for the monodentate ligands.  When ``None`` (default) all
        monodentates are treated as identical ("X").

    Returns
    -------
    int
        Number of geometric isomers (orbits under the polyhedron's proper
        rotation group), counting enantiomers as distinct — consistent with
        the rest of this module.

    Notes
    -----
    * The implementation is a thin wrapper around ``enumerate_chelate_configs``
      (which is the same engine the FF-free builder uses for placement),
      ensuring the theoretical count and the builder's enumeration draw from
      one single source of truth.
    * For the octahedron, ``count_chelate_isomers`` retains its original
      edge-combinatorial algorithm and is *byte-identical* — this helper
      is reachable from there only when the env-flag opts into the
      universal path.
    """
    group, n = _get_group(geometry)
    if n_chelate < 0 or n_monodentate_fillers < 0:
        raise ValueError("counts must be non-negative")
    if 2 * n_chelate + n_monodentate_fillers != n:
        raise ValueError(
            f"2*n_chelate + n_monodentate_fillers ({2*n_chelate + n_monodentate_fillers})"
            f" != polyhedron CN ({n}) for geometry={geometry!r}"
        )

    # Build ligand_specs in a deterministic order: chelates first (by id), then
    # monodentates (by index/label).  We use identical "type"='B' for symmetric
    # bidentates so all are interchangeable under the orbit canonicalisation,
    # and labels 'X0','X1',... for the monodentates when none are supplied.
    specs: List[Dict] = []
    for k in range(n_chelate):
        specs.append({"type": "B", "denticity": 2, "asym": bool(asymmetric)})
    if monodentate_labels is None:
        # All monodentates identical.
        for k in range(n_monodentate_fillers):
            specs.append({"type": "X", "denticity": 1})
    else:
        if len(monodentate_labels) != n_monodentate_fillers:
            raise ValueError("monodentate_labels length mismatch")
        for lab in monodentate_labels:
            specs.append({"type": str(lab), "denticity": 1})

    configs = enumerate_chelate_configs(geometry, specs)
    return len(configs)


def list_supported_chelate_geometries() -> List[str]:
    """Return the geometry keys for which ``count_chelate_isomers_universal``
    has a registered Pólya group + cis-edge map.  Useful as a sanity check
    and as documentation of what the v1 fix unlocks vs HEAD bcf56f8."""
    out = list(_GROUPS.keys())
    # f-block keys are added lazily via _fblock_group_for; advertise them too.
    for k in ("sap_8", "dd_8", "tdd_8", "ttp_9_fblock",
              "csap_9", "bicap_10", "ih_12"):
        if k not in out and _fblock_group_for(k) is not None:
            out.append(k)
    return sorted(out)


if __name__ == "__main__":
    tests = [
        ("octahedron", {"A": 6}, 1),
        ("octahedron", {"A": 5, "B": 1}, 1),
        ("octahedron", {"A": 4, "B": 2}, 2),     # cis/trans
        ("octahedron", {"A": 3, "B": 3}, 2),     # fac/mer
        ("octahedron", {"A": 2, "B": 2, "C": 2}, 6),
        ("octahedron", {"A": 1, "B": 1, "C": 1, "D": 1, "E": 1, "F": 1}, 30),
        ("square_planar", {"A": 2, "B": 2}, 2),  # cis/trans
        ("square_planar", {"A": 1, "B": 1, "C": 1, "D": 1}, 3),
        ("tetrahedron", {"A": 2, "B": 2}, 1),
        ("tetrahedron", {"A": 1, "B": 1, "C": 1, "D": 1}, 2),  # enantiomers
    ]
    print(f"{'geometry':<20}{'donors':<26}{'Burnside':>9}{'enum':>6}{'expect':>8}  ok")
    allok = True
    for geom, spec, exp in tests:
        c = count_isomers(geom, spec); e = len(enumerate_isomers(geom, spec))
        ok = (c == exp == e); allok &= ok
        print(f"{geom:<20}{str(spec):<26}{c:>9}{e:>6}{exp:>8}  {'OK' if ok else 'FAIL'}")
    print("\nALL TEXTBOOK TESTS PASS" if allok else "\nSOME TESTS FAILED")
