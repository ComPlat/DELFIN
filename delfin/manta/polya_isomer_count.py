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
import os
from collections import Counter
from typing import Dict, List, Tuple


def _delfin_env_int(name: str, default: int) -> int:
    """Local copy of ``smiles_converter._delfin_env_int`` (avoid import cycle)."""
    try:
        return int(os.environ.get(name, str(default)))
    except Exception:
        return default


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
    # CN2 L-2 linear: 2 antipodal vertices (180°).  Proper rotation group C2 (order 2):
    # identity + the C2 about an axis perpendicular to the molecular axis, which swaps
    # the two vertices.  iter-32f (DELFIN_FFFREE_CN_EXTEND).
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


def _trigonal_pyramidal_group():
    # CN3 TPY-3 trigonal pyramidal: 3 donors at the base of a pyramid with the metal
    # at the apex above the donor plane (vacant-tetrahedron / NH3 lone-pair geometry).
    # Point group C3v; its PROPER rotation subgroup is C3 (order 3): the C3 about the
    # pyramid axis (0->1->2->0).  Unlike the PLANAR SP-3 (D3), the pyramid has NO
    # in-plane C2 (the apex breaks it) and no horizontal mirror as a proper rotation,
    # so only the 3-fold survives.  iter-32g (User 2026-06-19 ATENET).
    C3 = (1, 2, 0)
    return _close_group([C3], 3)


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


_GROUPS = {
    "linear": (_linear_group(), 2),
    "trigonal_planar": (_trigonal_planar_group(), 3),
    "tshape": (_tshape_group(), 3),
    "trigonal_pyramidal": (_trigonal_pyramidal_group(), 3),
    "octahedron": (_octahedron_group(), 6),
    "square_planar": (_square_group(), 4),
    "tetrahedron": (_tetrahedron_group(), 4),
    "trigonal_bipyramid": (_tbp_group(), 5),
    "square_pyramid": (_spy_group(), 4 + 1),
    "trigonal_prism": (_trigonal_prism_group(), 6),
    "pentagonal_bipyramid": (_pentagonal_bipyramid_group(), 7),
    "square_antiprism": (_square_antiprism_group(), 8),
    "tricapped_trigonal_prism": (_tricapped_trigonal_prism_group(), 9),
}


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
    group, n = _GROUPS[geometry]
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
    group, n = _GROUPS[geometry]
    labels = _multiset_from_spec(donor_spec, n)
    seen = set(); reps = []
    # Deterministic enumeration order INDEPENDENT of PYTHONHASHSEED
    # (DELFIN_FFFREE_DETERMINISTIC_ENUM, default OFF).  The unique permutations were
    # iterated straight out of a ``set`` whose ordering is salted-hash dependent, so
    # the emitted isomer order — and the downstream fac/mer label numbering — varied
    # across processes UNLESS the caller pinned PYTHONHASHSEED=0 in the environment
    # (the shipped pools do).  Iterating in canonical sorted order makes the manifold
    # reproducible regardless of the hash seed.  OFF reproduces the legacy
    # set-iteration order (byte-identical under a pinned PYTHONHASHSEED=0); the SET of
    # representatives is identical either way (only the order/labelling changes).
    _perms = set(itertools.permutations(labels))
    if os.environ.get("DELFIN_FFFREE_DETERMINISTIC_ENUM", "0") == "1":
        _perms = sorted(_perms)
    for perm in _perms:
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
# NOTE: the 100deg default EXCLUDES the tetrahedron's only pair-angle (109.47deg), so a
#  Td bis-chelate yields 0 configs.  Stream-B Fix 1 (DELFIN_FFFREE_TET_CHELATE=1) raises
#  the ceiling to 111deg to admit the 109.47deg Td chelate edges.  Default-OFF keeps the
#  100deg ceiling -> byte-identical (OC-6/SP-4 only have 90/120/180deg pairs, none of which
#  fall in the (100, 111] band, so they are unchanged with the flag either way).
CHELATE_CIS_MAX_DEG = 100.0
CHELATE_CIS_MAX_DEG_TET = 111.0

_GEOM_KEY_TO_SHAPE = {
    "linear": "L-2 linear",
    "octahedron": "OC-6 octahedron",
    "square_planar": "SP-4 square planar",
    "tetrahedron": "T-4 tetrahedron",
    "trigonal_bipyramid": "TBP-5 trigonal bipyramid",
    "square_pyramid": "SPY-5 square pyramid",
    "pentagonal_bipyramid": "PB-7 pentagonal bipyramid",
    "square_antiprism": "SQAP-8 square antiprism",
    "tricapped_trigonal_prism": "TTP-9 tricapped trigonal prism",
}


def _chelate_cis_edges(geometry: str, n: int):
    """Vertex pairs a bidentate chelate can span (cis edges), derived geometrically from
    the ideal polyhedron: pairs whose ideal angular separation is below
    CHELATE_CIS_MAX_DEG.  Single source of truth = the SAME vertex set the assembler
    places into (delfin.manta.polyhedra), so enumeration and placement stay consistent
    for every geometry (incl. TBP-5/SPY-5).  Byte-identical to the old antipode table for
    octahedron/square_planar.  Falls back to the antipode table if no reference vectors."""
    # Stream-B Fix 1: gate the cis-edge ceiling so the tetrahedron's only pair-angle
    # (109.47deg) becomes chelate-feasible when DELFIN_FFFREE_TET_CHELATE=1.  Default-OFF
    # keeps the 100deg ceiling -> byte-identical (no geometry has a pair in (100, 111]
    # except the Td 109.47deg, which is exactly what we want to admit when ON).
    ceiling = (CHELATE_CIS_MAX_DEG_TET
               if _delfin_env_int("DELFIN_FFFREE_TET_CHELATE", 0)
               else CHELATE_CIS_MAX_DEG)
    shape = _GEOM_KEY_TO_SHAPE.get(geometry)
    if shape is not None:
        try:
            import math
            import numpy as np
            from delfin.manta import polyhedra as _PLY
            V = _PLY.ref_vectors(shape)
            cos_max = math.cos(math.radians(ceiling))
            return [(i, j) for i in range(n) for j in range(i + 1, n)
                    if float(np.clip(V[i] @ V[j], -1.0, 1.0)) > cos_max]
        except Exception:
            pass
    anti = _ANTIPODE_FULL[geometry]
    return [(i, j) for i in range(n) for j in range(i + 1, n) if anti[i] != j]


def _meridional_triples(geometry: str, n: int):
    """Vertex triples on which a RIGID PLANAR tridentate (terpy / pincer) can seat
    MERIDIONALLY: 3 vertices that are coplanar through the metal with one antipodal
    (~180deg) outer pair and the central vertex at ~90deg to both — the meridian arc
    a flat conjugated tridentate physically requires.  A facial triple (3 mutually
    ~90deg vertices, no antipodal pair) is EXCLUDED: a rigid flat ligand cannot fold
    to it.  Derived GEOMETRICALLY from the same ideal vertex set the assembler places
    into (delfin.manta.polyhedra), so it is correct for any CN6 polyhedron (and
    returns [] for geometries with no meridional triple).  Deterministic.

    A meridional triple is detected as: the largest pairwise angle is ~180deg (a
    near-antipodal outer pair) AND the three vertices are coplanar through the metal
    (their normal is well-defined and the middle vertex lies in the outer plane)."""
    import math
    import itertools
    import numpy as np
    shape = _GEOM_KEY_TO_SHAPE.get(geometry)
    if shape is None:
        return []
    try:
        from delfin.manta import polyhedra as _PLY
        V = _PLY.ref_vectors(shape)
    except Exception:
        return []
    # Outer-pair floor.  Default 170deg cleanly separates a true antipodal (180deg)
    # outer pair from facial caps (which top out ~90-110deg).  SPY-5 has NO 180deg
    # pair: its widest trans pair is the trans-BASAL arc (~163.9deg, central = the
    # apical vertex, coplanar through the metal) — the meridian a flat tridentate
    # binds in a square-pyramidal field.  Admit it ONLY under the CN5 rigid-planar
    # flag (DELFIN_FFFREE_PLANAR_MER_CN5) so OC-6/TBP-5 stay byte-identical when off.
    amax_floor = 170.0
    if (shape == "SPY-5 square pyramid"
            and os.environ.get("DELFIN_FFFREE_PLANAR_MER_CN5", "0") == "1"):
        amax_floor = 155.0
    out = []
    for tri in itertools.combinations(range(n), 3):
        vv = [V[i] for i in tri]
        angs = []
        for a in range(3):
            for b in range(a + 1, 3):
                c = float(np.clip(vv[a] @ vv[b], -1.0, 1.0))
                angs.append(math.degrees(math.acos(c)))
        amax = max(angs)
        # meridional: one near-antipodal outer pair; facial caps have NO such pair.
        if amax < amax_floor:
            continue
        # coplanarity through the metal: the 3 unit vectors span a plane (rank-2), i.e.
        # the smallest singular value is ~0 (they pass through the origin = the metal).
        try:
            _, sv, _ = np.linalg.svd(np.array(vv, float))
        except Exception:
            continue
        if sv[0] > 1e-9 and (sv[2] / sv[0]) < 0.05:
            out.append(tri)
    return out


def _equatorial_squares(geometry: str, n: int):
    """Vertex 4-sets on which a RIGID PLANAR tetradentate (porphyrin / phthalocyanine /
    planar Schiff-base macrocycle) can seat as an EQUATORIAL SQUARE: 4 vertices
    coplanar through the metal forming a square (two ~180deg antipodal diagonals +
    four ~90deg sides) — the equatorial plane of OC-6 (axial sites left free for
    axial ligands) or the full SP-4 square.  A flat conjugated macrocycle physically
    requires its 4 donors coplanar with the metal; any non-coplanar 4-subset (a
    folded / cis-bent arrangement) is EXCLUDED.  Derived GEOMETRICALLY from the same
    ideal vertex set the assembler places into (delfin.manta.polyhedra), so it is
    correct for any CN>=4 polyhedron (and returns [] for geometries with no
    equatorial square, e.g. TBP-5 / SPY-5).  Deterministic."""
    import math
    import itertools
    import numpy as np
    shape = _GEOM_KEY_TO_SHAPE.get(geometry)
    if shape is None:
        return []
    try:
        from delfin.manta import polyhedra as _PLY
        V = _PLY.ref_vectors(shape)
    except Exception:
        return []
    out = []
    for quad in itertools.combinations(range(n), 4):
        vv = [V[i] for i in quad]
        # coplanar through the metal: the 4 unit vectors span a plane (rank-2), i.e.
        # the smallest singular value is ~0 (they pass through the origin = the metal).
        try:
            _, sv, _ = np.linalg.svd(np.array(vv, float))
        except Exception:
            continue
        if not (sv[0] > 1e-9 and (sv[2] / sv[0]) < 0.05):
            continue
        # square signature: of the 6 pairwise angles, exactly TWO are ~180deg (the
        # diagonals) and the remaining FOUR are ~90deg (the sides).  A coplanar but
        # non-square 4-set (e.g. trapezoidal) is rejected.
        angs = []
        for a in range(4):
            for b in range(a + 1, 4):
                c = float(np.clip(vv[a] @ vv[b], -1.0, 1.0))
                angs.append(math.degrees(math.acos(c)))
        n180 = sum(1 for x in angs if x >= 170.0)
        n90 = sum(1 for x in angs if 75.0 <= x <= 105.0)
        if n180 == 2 and n90 == 4:
            out.append(quad)
    return out


def enumerate_chelate_configs(geometry: str, ligand_specs):
    """Universal isomer enumeration for a mix of chelating (bidentate) + monodentate
    ligands.  ligand_specs: one dict per ligand instance with keys ``type``
    (identity label), ``denticity`` (1 or 2), ``asym`` (bool: True if a bidentate's
    two donor arms are distinguishable).  Returns distinct configs; each maps
    vertex_index -> (ligand_instance_index, arm_index), one per orbit under the
    polyhedron's proper rotation group.  Bidentate ligands occupy cis-edges."""
    group, n = _GROUPS[geometry]
    cis_edges = _chelate_cis_edges(geometry, n)
    order = sorted(range(len(ligand_specs)),
                   key=lambda k: -ligand_specs[k]["denticity"])
    seen = set(); out = []
    # Geometry-aware meridional restriction (DELFIN_FFFREE_PLANAR_MER): a RIGID PLANAR
    # tridentate (terpy / pincer, spec["rigid_planar"] True) may occupy ONLY meridional
    # vertex triples (coplanar through the metal, outer-outer ~180deg) — never a facial
    # cap.  Computed once; empty if any rigid-planar spec is present but the geometry
    # has no meridional triple (the caller then yields no config -> legacy fallback).
    _need_mer = any(s.get("rigid_planar") and s["denticity"] == 3 for s in ligand_specs)
    _mer_triples = _meridional_triples(geometry, n) if _need_mer else None
    # Geometry-aware EQUATORIAL-SQUARE restriction (DELFIN_FFFREE_KAPPA4): a RIGID
    # PLANAR tetradentate (porphyrin / phthalocyanine, spec["rigid_planar"] True) may
    # occupy ONLY a coplanar-through-the-metal square (the OC-6 equatorial plane /
    # SP-4 square), leaving any axial vertices free — never a folded 4-subset.  Empty
    # if a rigid-planar κ4 spec is present but the geometry has no equatorial square
    # (the caller then yields no config -> legacy fallback).
    _need_sq = any(s.get("rigid_planar") and s["denticity"] == 4 for s in ligand_specs)
    _eq_squares = _equatorial_squares(geometry, n) if _need_sq else None

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
        elif dent == 3 and spec.get("rigid_planar") and _mer_triples is not None:
            # RIGID PLANAR tridentate (terpy / pincer): occupy ONLY meridional vertex
            # triples (a flat conjugated tridentate cannot fold to a facial cap).  Seat
            # the CENTRAL donor (canonical arm index 1) on the meridian's central
            # (~90deg) vertex and the two OUTER donors (arms 0,2) on the antipodal
            # (~180deg) outer pair, so the placed donors are coplanar with the metal.
            import numpy as _np
            from delfin.manta import polyhedra as _PLY2
            try:
                _V = _PLY2.ref_vectors(_GEOM_KEY_TO_SHAPE[geometry])
            except Exception:
                _V = None
            for tri in _mer_triples:
                if any(v in assign for v in tri):
                    continue
                # identify the central vertex = the one NOT in the ~180deg outer pair
                central = tri[0]
                if _V is not None:
                    worst = (-1.0, (tri[0], tri[1]))
                    for ia in range(3):
                        for ib in range(ia + 1, 3):
                            ang = -float(_np.clip(_V[tri[ia]] @ _V[tri[ib]], -1.0, 1.0))
                            if ang > worst[0]:
                                worst = (ang, (tri[ia], tri[ib]))
                    outer = set(worst[1])
                    central = next(v for v in tri if v not in outer)
                outer_v = [v for v in tri if v != central]
                # arm 1 = central donor -> central vertex; arms 0,2 = outer donors ->
                # outer vertices (both orientations enumerated for asym backbones).
                a = dict(assign)
                a[outer_v[0]] = (k, 0); a[central] = (k, 1); a[outer_v[1]] = (k, 2)
                place(a)
                if spec.get("asym"):
                    b = dict(assign)
                    b[outer_v[1]] = (k, 0); b[central] = (k, 1); b[outer_v[0]] = (k, 2)
                    place(b)
        elif dent == 4 and spec.get("rigid_planar") and _eq_squares is not None:
            # RIGID PLANAR tetradentate (porphyrin / phthalocyanine): occupy ONLY an
            # EQUATORIAL SQUARE (4 vertices coplanar through the metal), leaving any
            # axial vertices free for axial ligands.  The macrocycle is a closed
            # 4-cycle arm0-arm1-arm2-arm3-arm0, so consecutive donor arms must sit on
            # ADJACENT (~90deg) square vertices, never across the diagonal.  Enumerate
            # the square's dihedral arrangements; canon_key dedups to the orbit(s).
            import numpy as _np
            from delfin.manta import polyhedra as _PLY3
            try:
                _Vq = _PLY3.ref_vectors(_GEOM_KEY_TO_SHAPE[geometry])
            except Exception:
                _Vq = None
            for sq in _eq_squares:
                if _Vq is None or any(v in assign for v in sq):
                    continue
                # cyclic ring order [w0, adjacent, diagonal, adjacent] so that adjacent
                # ring positions are the ~90deg (square-edge) neighbours.
                w0 = sq[0]
                dpart = min((v for v in sq if v != w0),
                            key=lambda v: float(_Vq[w0] @ _Vq[v]))
                adj = [v for v in sq if v != w0 and v != dpart]
                ring = [w0, adj[0], dpart, adj[1]]
                # dihedral group of the square (4 rotations x 2 reflections); canon dedups.
                for refl in (False, True):
                    for rot in range(4):
                        a = dict(assign)
                        for arm in range(4):
                            pos = (rot - arm) % 4 if refl else (rot + arm) % 4
                            a[ring[pos]] = (k, arm)
                        place(a)
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
    en2X2 -> trans + cis-Δ + cis-Λ = 3)."""
    group, n = _GROUPS[geometry]
    if geometry == "octahedron":
        antipode = {0: 1, 1: 0, 2: 3, 3: 2, 4: 5, 5: 4}
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
