"""Core-frozen clash-relief relaxation for conformer candidates.

The Welle-5o conformer generator rotates backbone torsions RIGIDLY.  In a packed
coordination sphere a rigid rotation routinely drives the moving fragment into gross
overlap with the rest of the complex — measured worst contacts of 0.06-0.5 A — so the
downstream clash gate (correctly) rejects 100% of candidates and only the base frame
survives (the under-yield root, 2026-06-26).  This module relieves those clashes with
a short geometry-only spring relaxation that FREEZES the coordination core (metal +
every first-shell heavy donor), pulls 1,2-bonded pairs toward their ideal length, and
pushes non-bonded pairs (excluding 1,2/1,3 neighbours, exactly as the clash gate
exempts them) apart to their clash floor.  Geometry-only (covalent-radii springs, no
force field / RDKit sanitisation) so it is robust for any metal and fully
deterministic.  UFF is unfit for RANKING TM complexes but ordinary spring relaxation
is perfectly fine for RELIEVING a clash; GFN-FF still does the final ranking.

OPT-IN: the caller gates on DELFIN_CONF_RELAX_CANDIDATES; with the flag off this
module is never imported.
"""
from __future__ import annotations

from typing import Dict, List, Sequence, Set, Tuple

import numpy as np

from delfin._bond_decollapse import _ideal_bond


def core_indices(graph: Dict) -> Set[int]:
    """Metal atoms + their first-shell heavy neighbours (donors) — the atoms that
    must stay pinned so the backbone re-folds around a fixed polyhedron."""
    n = graph.get("n_atoms", 0)
    is_metal = graph["is_metal"]
    atomic_nums = graph["atomic_nums"]
    neighbours = graph["neighbours"]
    core: Set[int] = set()
    for m in range(n):
        if not is_metal[m]:
            continue
        core.add(m)
        for d in neighbours[m]:
            if atomic_nums[d] != 1:
                core.add(d)
    return core


def relax_clashes(
    symbols: Sequence[str],
    coords: Sequence[Tuple[float, float, float]],
    graph: Dict,
    frozen: Set[int],
    cut_hh: float,
    cut_xh: float,
    cut_xx: float,
    n_iter: int = 150,
    max_step: float = 0.20,
) -> List[Tuple[float, float, float]]:
    """Return clash-relieved coordinates (core atoms in *frozen* never move).

    Deterministic: fixed iteration count, no RNG.  If anything is malformed the
    input coordinates are returned unchanged."""
    try:
        P = np.asarray(coords, dtype=float)
        n = len(symbols)
        if P.shape != (n, 3) or n < 3:
            return list(coords)

        # 1,2 bonds + ideal lengths
        bonded: Set[Tuple[int, int]] = set()
        for (a, b, _o, _ar, _r) in graph.get("bonds", []):
            bonded.add((min(a, b), max(a, b)))
        adj: Dict[int, List[int]] = {i: [] for i in range(n)}
        for (i, j) in bonded:
            adj[i].append(j)
            adj[j].append(i)
        # 1,3 exclusion (same exemption the clash gate uses)
        excl: Set[Tuple[int, int]] = set(bonded)
        for k in range(n):
            nb = adj[k]
            for x in range(len(nb)):
                for y in range(x + 1, len(nb)):
                    excl.add((min(nb[x], nb[y]), max(nb[x], nb[y])))

        def floor(i: int, j: int) -> float:
            si, sj = symbols[i], symbols[j]
            if si == "H" and sj == "H":
                return cut_hh
            if si == "H" or sj == "H":
                return cut_xh
            return cut_xx

        # ---- precompute the per-iteration-invariant index arrays ONCE (vectorized;
        # a per-pair Python loop over the O(n^2) non-bonded set x 150 passes is far
        # too slow for the ~500 candidates per isomer) ----
        if bonded:
            _bi = np.fromiter((i for (i, j) in bonded), dtype=np.intp, count=len(bonded))
            _bj = np.fromiter((j for (i, j) in bonded), dtype=np.intp, count=len(bonded))
            # Restrain each 1,2 bond toward the SHORTER of its native length and the
            # covalent-radii ideal — keep bonds intact while the non-bonded repulsion
            # relieves the clash.  (Pure native-length restraint over-stiffens and
            # starves the clash relief; pure generic-ideal can stretch a real bond.)
            _native = np.sqrt(((P[_bj] - P[_bi]) ** 2).sum(1))
            _ideal = np.array([_ideal_bond(symbols[i], symbols[j]) for (i, j) in bonded],
                              dtype=float)
            _btgt = np.minimum(_native, _ideal)
        else:
            _bi = _bj = np.empty(0, dtype=np.intp); _btgt = np.empty(0, dtype=float)
        bscatter = np.empty(2 * _bi.shape[0], dtype=np.intp)
        bscatter[0::2] = _bi; bscatter[1::2] = _bj

        _ri: List[int] = []; _rj: List[int] = []; _rf: List[float] = []
        for i in range(n):
            for j in range(i + 1, n):
                if (i, j) in excl:
                    continue
                _ri.append(i); _rj.append(j); _rf.append(floor(i, j))
        ri = np.asarray(_ri, dtype=np.intp); rj = np.asarray(_rj, dtype=np.intp)
        rfloor = np.asarray(_rf, dtype=float)
        rscatter = np.empty(2 * ri.shape[0], dtype=np.intp)
        rscatter[0::2] = ri; rscatter[1::2] = rj

        free = np.array([i not in frozen for i in range(n)], dtype=bool)

        for _ in range(n_iter):
            F = np.zeros((n, 3), dtype=float)
            moved = False
            # bonded springs -> ideal length (vectorized)
            if _bi.shape[0]:
                d = P[_bj] - P[_bi]
                r = np.sqrt((d * d).sum(1))
                act = (r > 1e-6) & (np.abs(r - _btgt) > 0.02)
                if act.any():
                    coef = np.where(act, 0.5 * (r - _btgt), 0.0)
                    f = (coef / np.where(r > 1e-6, r, 1.0))[:, None] * d
                    vals = np.empty((2 * _bi.shape[0], 3)); vals[0::2] = f; vals[1::2] = -f
                    np.add.at(F, bscatter, vals)
                    moved = True
            # non-bonded repulsion -> push apart to floor (vectorized)
            if ri.shape[0]:
                d = P[rj] - P[ri]
                r = np.sqrt((d * d).sum(1))
                act = (r < rfloor) & (r > 1e-9)
                if act.any():
                    coef = np.where(act, 0.5 * (r - rfloor), 0.0)
                    f = (coef / np.where(r > 1e-9, r, 1.0))[:, None] * d
                    vals = np.empty((2 * ri.shape[0], 3)); vals[0::2] = f; vals[1::2] = -f
                    np.add.at(F, rscatter, vals)
                    moved = True
            if not moved:
                break
            sn = np.sqrt((F * F).sum(1))
            scale = np.where(sn > max_step, max_step / np.where(sn > 1e-12, sn, 1.0), 1.0)
            P = P + (F * scale[:, None]) * free[:, None]
            if not np.all(np.isfinite(P)):
                return list(coords)

        return [tuple(float(x) for x in row) for row in P]
    except Exception:
        return list(coords)
