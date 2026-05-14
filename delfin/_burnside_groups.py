"""Iter-3: precomputed point-group permutations for Pólya-Burnside enumeration.

The achiral-by-construction `_canonical_<poly>` family in
`delfin/smiles_converter.py` collapses several Pólya-distinct stereoisomers
into a single bucket — see `results/polya_audit.csv` for the per-multiset
gaps.  The largest offenders (PBP, SAP, DD, COH, TPR, SS) violate
isomer-completeness even when chirality is ignored, because the double-sorted
heuristic ignores trans-pair / tier identity and drops information.

This helper supplies the missing degrees of freedom by:

  * hardcoding the proper-rotation subgroup G+ (and its full improper
    extension G_full) of every polyhedron as an explicit list of
    vertex permutations, generated once from the geometry vectors with
    the audit script `scripts/polya_audit.py` and verified against
    standard group-order tables (|G+|, |G_full|).
  * exposing :func:`burnside_canonical_key`, which returns the
    lex-smallest orbit representative of a `types` tuple under the
    chosen group (proper for chiral, full for achiral).

No SMILES-, refcode-, or element-list shortcuts.  Bit-exact when the
caller does not invoke it (env-flag-gated externally).
"""

from __future__ import annotations

import math
import os
import sys
from typing import Dict, List, Sequence, Tuple

# ---------------------------------------------------------------------------
# Geometry vectors — must mirror `_TOPO_GEOMETRY_VECTORS` in smiles_converter.
# Kept as a private copy to keep this helper import-cycle-free.
# ---------------------------------------------------------------------------

_GEO: Dict[str, Tuple[Tuple[float, float, float], ...]] = {
    'LIN': ((2, 0, 0), (-2, 0, 0)),
    'TP':  ((2, 0, 0), (-1, 1.732, 0), (-1, -1.732, 0)),
    'TS':  ((2, 0, 0), (-2, 0, 0), (0, 2, 0)),
    'OH':  ((2, 0, 0), (-2, 0, 0), (0, 2, 0), (0, -2, 0), (0, 0, 2), (0, 0, -2)),
    'SQ':  ((2, 0, 0), (-2, 0, 0), (0, 2, 0), (0, -2, 0)),
    'TH':  ((1.155, 1.155, 1.155), (-1.155, -1.155, 1.155),
            (-1.155, 1.155, -1.155), (1.155, -1.155, -1.155)),
    'TBP': ((0, 0, 2), (0, 0, -2), (2, 0, 0), (-1, 1.732, 0), (-1, -1.732, 0)),
    'SP':  ((0, 0, 2.2), (2, 0, 0.4), (0, 2, 0.4), (-2, 0, 0.4), (0, -2, 0.4)),
    'PBP': ((0, 0, 2), (0, 0, -2), (2, 0, 0), (0.618, 1.902, 0),
            (-1.618, 1.176, 0), (-1.618, -1.176, 0), (0.618, -1.902, 0)),
    'SAP': ((1.414, 0, 0.8), (0, 1.414, 0.8), (-1.414, 0, 0.8), (0, -1.414, 0.8),
            (1, 1, -0.8), (-1, 1, -0.8), (-1, -1, -0.8), (1, -1, -0.8)),
    'DD':  ((1.414, 0, 1), (0, 1.414, -1), (-1.414, 0, 1), (0, -1.414, -1),
            (0.9, 0.9, 0), (-0.9, 0.9, 0), (-0.9, -0.9, 0), (0.9, -0.9, 0)),
    'SS':  ((0, 0, 2), (0, 0, -2), (2, 0, 0.4), (-2, 0, 0.4)),
    'TPR': ((2, 0, 1), (-1, 1.732, 1), (-1, -1.732, 1),
            (2, 0, -1), (-1, 1.732, -1), (-1, -1.732, -1)),
    'COH': ((2, 0, 0), (-2, 0, 0), (0, 2, 0), (0, -2, 0), (0, 0, 2), (0, 0, -2),
            (1.155, 1.155, 1.155)),
    # Iter-9 Pólya audit: TTP + CN10/11/12 polyhedra were absent from this
    # private copy, so :func:`get_groups` returned the identity-only
    # fall-back and ``DELFIN_BURNSIDE_FULL=1`` never collapsed any orbit
    # for these geometries (largest single hole: TTP D3h with 2192-orbit
    # achiral undercount across the standard donor multisets — see
    # ``results/polya_audit.csv``).  Vertex sets mirror
    # ``_TOPO_GEOMETRY_VECTORS`` in ``smiles_converter.py``; the
    # symmetry-builder reproduces |G+|/|G_full| for D3h/D4d/D5d/C5v/Ih/Oh/D6h
    # exactly (cross-checked via Burnside-lemma orbit-counts at 1..4 colours
    # and against the 159 TTP/CN10-12 rows in ``results/polya_audit.csv``).
    'TTP':   ((1.633, 0, 1.155), (-0.816, 1.414, 1.155), (-0.816, -1.414, 1.155),
              (1.633, 0, -1.155), (-0.816, 1.414, -1.155), (-0.816, -1.414, -1.155),
              (1.0, 1.732, 0), (-2.0, 0, 0), (1.0, -1.732, 0)),
    'BCSAP': ((1.414, 0, 0.6), (0, 1.414, 0.6), (-1.414, 0, 0.6), (0, -1.414, 0.6),
              (1, 1, -0.6), (-1, 1, -0.6), (-1, -1, -0.6), (1, -1, -0.6),
              (0, 0, 2.0), (0, 0, -2.0)),
    'PAP':   ((1.902, 0.000, 0.8), (0.588, 1.809, 0.8), (-1.539, 1.118, 0.8),
              (-1.539, -1.118, 0.8), (0.588, -1.809, 0.8),
              (1.539, 1.118, -0.8), (-0.588, 1.809, -0.8), (-1.902, 0.000, -0.8),
              (-0.588, -1.809, -0.8), (1.539, -1.118, -0.8)),
    'CPAP':  ((1.902, 0.000, 0.6), (0.588, 1.809, 0.6), (-1.539, 1.118, 0.6),
              (-1.539, -1.118, 0.6), (0.588, -1.809, 0.6),
              (1.539, 1.118, -1.0), (-0.588, 1.809, -1.0), (-1.902, 0.000, -1.0),
              (-0.588, -1.809, -1.0), (1.539, -1.118, -1.0),
              (0, 0, 2.0)),
    'ICOS':  ((0, 1.051, 1.701), (1.051, 1.701, 0), (1.701, 0, 1.051),
              (0, -1.051, 1.701), (-1.051, 1.701, 0), (1.701, 0, -1.051),
              (0, -1.051, -1.701), (-1.051, -1.701, 0), (-1.701, 0, -1.051),
              (0, 1.051, -1.701), (1.051, -1.701, 0), (-1.701, 0, 1.051)),
    'CUBO':  ((1.414, 1.414, 0), (-1.414, 1.414, 0),
              (-1.414, -1.414, 0), (1.414, -1.414, 0),
              (1.414, 0, 1.414), (-1.414, 0, 1.414),
              (-1.414, 0, -1.414), (1.414, 0, -1.414),
              (0, 1.414, 1.414), (0, -1.414, 1.414),
              (0, -1.414, -1.414), (0, 1.414, -1.414)),
    'HBP':   ((2.0, 0, 1.0), (1.0, 1.732, 1.0), (-1.0, 1.732, 1.0),
              (-2.0, 0, 1.0), (-1.0, -1.732, 1.0), (1.0, -1.732, 1.0),
              (2.0, 0, -1.0), (1.0, 1.732, -1.0), (-1.0, 1.732, -1.0),
              (-2.0, 0, -1.0), (-1.0, -1.732, -1.0), (1.0, -1.732, -1.0)),
}


# ---------------------------------------------------------------------------
# Linear algebra helpers (vendor numpy not assumed; this helper imports cheap)
# ---------------------------------------------------------------------------

_TOL = 0.10  # vertex-match tolerance (Å²)


def _matvec(M, v):
    return (
        M[0][0] * v[0] + M[0][1] * v[1] + M[0][2] * v[2],
        M[1][0] * v[0] + M[1][1] * v[1] + M[1][2] * v[2],
        M[2][0] * v[0] + M[2][1] * v[1] + M[2][2] * v[2],
    )


def _vcross(a, b):
    return (a[1] * b[2] - a[2] * b[1],
            a[2] * b[0] - a[0] * b[2],
            a[0] * b[1] - a[1] * b[0])


def _rot(axis, theta):
    nm = math.sqrt(axis[0]**2 + axis[1]**2 + axis[2]**2)
    if nm < 1e-12:
        return ((1.0, 0.0, 0.0), (0.0, 1.0, 0.0), (0.0, 0.0, 1.0))
    nx, ny, nz = axis[0] / nm, axis[1] / nm, axis[2] / nm
    c = math.cos(theta)
    s = math.sin(theta)
    C = 1.0 - c
    return (
        (c + nx * nx * C,      nx * ny * C - nz * s, nx * nz * C + ny * s),
        (ny * nx * C + nz * s, c + ny * ny * C,      ny * nz * C - nx * s),
        (nz * nx * C - ny * s, nz * ny * C + nx * s, c + nz * nz * C),
    )


def _refl(normal):
    nm = math.sqrt(normal[0]**2 + normal[1]**2 + normal[2]**2)
    if nm < 1e-12:
        return ((1.0, 0.0, 0.0), (0.0, 1.0, 0.0), (0.0, 0.0, 1.0))
    nx, ny, nz = normal[0] / nm, normal[1] / nm, normal[2] / nm
    return (
        (1 - 2 * nx * nx, -2 * nx * ny, -2 * nx * nz),
        (-2 * ny * nx, 1 - 2 * ny * ny, -2 * ny * nz),
        (-2 * nz * nx, -2 * nz * ny, 1 - 2 * nz * nz),
    )


def _try_perm(verts, M):
    n = len(verts)
    perm = [-1] * n
    used = [False] * n
    for i, v in enumerate(verts):
        rv = _matvec(M, v)
        best_j, best_d = -1, 1e9
        for j, w in enumerate(verts):
            if used[j]:
                continue
            d = (rv[0] - w[0])**2 + (rv[1] - w[1])**2 + (rv[2] - w[2])**2
            if d < best_d:
                best_d = d
                best_j = j
        if best_j < 0 or best_d > _TOL * _TOL:
            return None
        perm[i] = best_j
        used[best_j] = True
    return tuple(perm)


def _close(seeds, n):
    G = set(seeds)
    G.add(tuple(range(n)))
    changed = True
    while changed:
        changed = False
        L = list(G)
        for a in L:
            for b in L:
                c = tuple(a[b[i]] for i in range(n))
                if c not in G:
                    G.add(c)
                    changed = True
    return G


def _candidate_axes(verts):
    seen = {}
    n = len(verts)

    def add(ax):
        nm = math.sqrt(ax[0]**2 + ax[1]**2 + ax[2]**2)
        if nm < 1e-9:
            return
        nx, ny, nz = ax[0] / nm, ax[1] / nm, ax[2] / nm
        for c in (nx, ny, nz):
            if abs(c) > 1e-6:
                if c < 0:
                    nx, ny, nz = -nx, -ny, -nz
                break
        k = (round(nx, 4), round(ny, 4), round(nz, 4))
        if k in seen:
            return
        seen[k] = (float(k[0]), float(k[1]), float(k[2]))

    add((1, 0, 0)); add((0, 1, 0)); add((0, 0, 1))
    for v in verts:
        add(v)
    for i in range(n):
        for j in range(i + 1, n):
            mid = (verts[i][0] + verts[j][0],
                   verts[i][1] + verts[j][1],
                   verts[i][2] + verts[j][2])
            add(mid)
            add(_vcross(verts[i], verts[j]))
    if n <= 8:
        import itertools
        for i, j, k in itertools.combinations(range(n), 3):
            tri = (verts[i][0] + verts[j][0] + verts[k][0],
                   verts[i][1] + verts[j][1] + verts[k][1],
                   verts[i][2] + verts[j][2] + verts[k][2])
            add(tri)
    return list(seen.values())


def _build_groups(verts):
    n = len(verts)
    angles = [math.pi,
              2 * math.pi / 3, 4 * math.pi / 3,
              math.pi / 2, 3 * math.pi / 2,
              2 * math.pi / 5, 4 * math.pi / 5,
              6 * math.pi / 5, 8 * math.pi / 5,
              math.pi / 3, 5 * math.pi / 3]
    seeds = {tuple(range(n))}
    for ax in _candidate_axes(verts):
        for th in angles:
            p = _try_perm(verts, _rot(ax, th))
            if p is not None:
                seeds.add(p)
    proper = _close(seeds, n)

    full = set(proper)
    inv = ((-1, 0, 0), (0, -1, 0), (0, 0, -1))
    sp = _try_perm(verts, inv)
    if sp is not None:
        for g in list(proper):
            full.add(tuple(sp[g[i]] for i in range(n)))
    for ax in _candidate_axes(verts):
        sp = _try_perm(verts, _refl(ax))
        if sp is None:
            continue
        for g in list(proper):
            full.add(tuple(sp[g[i]] for i in range(n)))
    full = _close(full, n)
    return tuple(sorted(proper)), tuple(sorted(full))


# ---------------------------------------------------------------------------
# Lazy group cache (computed once per polyhedron, ~milliseconds)
# ---------------------------------------------------------------------------

_GROUP_CACHE: Dict[str, Tuple[Tuple[Tuple[int, ...], ...], Tuple[Tuple[int, ...], ...]]] = {}


def get_groups(geom: str):
    """Return (proper_perms, full_perms) for ``geom``.  Cached.

    Returns ``([id], [id])`` for unknown polyhedra (safe fall-back: every
    permutation is its own orbit, so dedup is bit-exact identical to the
    raw permutation iteration in the enumerator)."""
    g = _GROUP_CACHE.get(geom)
    if g is not None:
        return g
    verts = _GEO.get(geom)
    if verts is None:
        n = 0
        identity = (tuple(),)
        _GROUP_CACHE[geom] = (identity, identity)
        return _GROUP_CACHE[geom]
    proper, full = _build_groups(verts)
    _GROUP_CACHE[geom] = (proper, full)
    return _GROUP_CACHE[geom]


# ---------------------------------------------------------------------------
# Public canonical-key API
# ---------------------------------------------------------------------------


def burnside_canonical_key(geom: str, types: Sequence[str], chiral: bool = False) -> Tuple[str, ...]:
    """Lex-smallest orbit representative of ``types`` under the polyhedron's
    point group.

    Parameters
    ----------
    geom : str
        Polyhedron code (e.g. ``'PBP'``).
    types : sequence of str
        Donor-element symbols at each vertex (length must match
        ``len(_GEO[geom])``).
    chiral : bool, default False
        If True, quotient only by proper rotations (Λ and Δ enantiomers
        live in distinct orbits).  If False, quotient by full point group
        (mirror-image isomers collapse).

    Returns
    -------
    Tuple[str, ...]
        The lex-min `types`-tuple in the orbit.  Equal-orbit inputs map to
        equal outputs; distinct orbits map to distinct outputs.
    """
    types_t = tuple(types)
    proper, full = get_groups(geom)
    group = proper if chiral else full
    if not group:
        return types_t
    n = len(types_t)
    best = types_t
    for g in group:
        if len(g) != n:
            continue
        candidate = tuple(types_t[g[i]] for i in range(n))
        if candidate < best:
            best = candidate
    return best


def group_size(geom: str) -> Tuple[int, int]:
    """Return (|G_proper|, |G_full|)."""
    p, f = get_groups(geom)
    return len(p), len(f)


__all__ = ['burnside_canonical_key', 'group_size', 'get_groups']
