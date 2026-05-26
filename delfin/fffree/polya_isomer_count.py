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
from typing import Dict, List, Tuple


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


_GROUPS = {
    "octahedron": (_octahedron_group(), 6),
    "square_planar": (_square_group(), 4),
    "tetrahedron": (_tetrahedron_group(), 4),
    "trigonal_bipyramid": (_tbp_group(), 5),
    "square_pyramid": (_spy_group(), 4 + 1),
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
CHELATE_CIS_MAX_DEG = 100.0

_GEOM_KEY_TO_SHAPE = {
    "octahedron": "OC-6 octahedron",
    "square_planar": "SP-4 square planar",
    "tetrahedron": "T-4 tetrahedron",
    "trigonal_bipyramid": "TBP-5 trigonal bipyramid",
    "square_pyramid": "SPY-5 square pyramid",
}


def _chelate_cis_edges(geometry: str, n: int):
    """Vertex pairs a bidentate chelate can span (cis edges), derived geometrically from
    the ideal polyhedron: pairs whose ideal angular separation is below
    CHELATE_CIS_MAX_DEG.  Single source of truth = the SAME vertex set the assembler
    places into (delfin.fffree.polyhedra), so enumeration and placement stay consistent
    for every geometry (incl. TBP-5/SPY-5).  Byte-identical to the old antipode table for
    octahedron/square_planar.  Falls back to the antipode table if no reference vectors."""
    shape = _GEOM_KEY_TO_SHAPE.get(geometry)
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
    anti = _ANTIPODE_FULL[geometry]
    return [(i, j) for i in range(n) for j in range(i + 1, n) if anti[i] != j]


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
        if spec["denticity"] == 1:
            for v in [v for v in range(n) if v not in assign]:
                a = dict(assign); a[v] = (k, 0); place(a)
        else:
            for (v1, v2) in cis_edges:
                if v1 in assign or v2 in assign:
                    continue
                a = dict(assign); a[v1] = (k, 0); a[v2] = (k, 1); place(a)
                if spec.get("asym"):
                    b = dict(assign); b[v1] = (k, 1); b[v2] = (k, 0); place(b)

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
