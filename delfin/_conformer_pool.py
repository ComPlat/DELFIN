"""Universal per-isomer conformer-pool generator (Welle-5o).

Goal
----
For every emitted XYZ frame (one prescribed coordination isomer's
3D realisation), generate **K diverse conformers** spanning the relevant
conformational space, so that downstream xTB / DFT local-opt finds the
*global* minimum without needing a CREST / GOAT global conformer search.

This is the *CREST-obsolescence enabler* per the user direktive 2026-05-18:

    "gesamte konformationsraum der koordinationsisomere muss
     abgedeckt sein damit relevante konformere erhalten bleiben"

Without per-isomer conformer-pool the downstream local-opt may converge
to a non-global minimum.  With it, each pool member is a different basin
of attraction; one of them will be the global-minimum basin.

Architecture (5 sampling layers, then quality filter + top-K diverse)
---------------------------------------------------------------------
1. **Torsion-grid** — staggered rotamer states around every rotatable
   single bond (re-uses :func:`delfin._rotamer_diversity.identify_rotamer_dofs`).
2. **Ring-pucker** — for each non-aromatic sp3 ring of size 5 / 6, two
   pucker modes (chair-like inversion + twist-flip) on the ring atom
   pair displaced perpendicular to the mean ring plane.
3. **Chelate-twist** — for each ring that contains the metal and at
   least one sp3 backbone atom of size ≤ 7, two δ-vs-λ twist conformers
   obtained by sign-flipping the backbone dihedral of the longest sp3
   sequence in the ring.
4. **Macrocycle modes** — for each ring of size ≥ 7 (universal
   macrocycle detection), three concerted out-of-plane modes (saddle,
   ruffled, dome) implemented as low-order Fourier displacement of the
   ring atoms along the ring normal.
5. **Quality filter** — every candidate must preserve M–D bond
   distances (±MD-tol), keep the heavy-atom bond multiset (topology
   hash), avoid any non-bonded atom overlap < 1.5 Å (heavy-heavy) or
   < 1.7 Å (X–H).
6. **Top-K diverse selection** — rank by UFF energy + clash penalty
   (same scorer as Welle-5l T6), then greedily pick the cheapest
   candidate whose pairwise RMSD vs. all already-selected pool members
   is ≥ ``diversity_rmsd_min`` (default 0.5 Å).  Base XYZ is always
   selected first.

Universal-fundamental compliance (per :doc:`feedback_universal_fundamental_doctrine`)
------------------------------------------------------------------------------------
* Every trigger uses **graph features only** — no SMILES patterns, no
  element allowlists, no refcode prefixes.
* Ring membership: ``GetRingInfo().AtomRings()`` from RDKit (or OB).
* Macrocycle: ring size ≥ 7 atoms — universal across all elements.
* Chelate-twist: ring contains a metal atom + at least one sp3 atom +
  ring size ≤ 7 (graph property — works for Fe(en)3, Co(diamine)3,
  Cu(acac)2, etc.).
* M–D invariant via :func:`delfin._rotamer_diversity._coord_bond_invariant_holds`.

Wire-in
-------
:mod:`delfin.smiles_converter` calls :func:`apply_if_enabled` per emitted
isomer.  Default-OFF byte-identical when
``DELFIN_5O_CONFORMER_POOL=0`` (master flag unset).

Env-flags
---------
``DELFIN_5O_CONFORMER_POOL``     (default ``0``) — master switch
``DELFIN_5O_K_TARGET``           (default ``10``) — pool size (K)
``DELFIN_5O_DIVERSITY_RMSD_MIN`` (default ``0.5``) — Å, pairwise RMSD floor
``DELFIN_5O_K_TORSION_STATES``   (default ``3``) — torsion grid states / DOF
``DELFIN_5O_MAX_DOFS``           (default ``6``) — torsion DOF cap
``DELFIN_5O_TORSION_GRID_CAP``   (default ``64``) — torsion combinations cap
``DELFIN_5O_RING_PUCKER_AMPL``   (default ``0.35``) — Å, pucker displacement
``DELFIN_5O_MACROCYCLE_AMPL``    (default ``0.40``) — Å, macrocycle mode amplitude
``DELFIN_5O_CHELATE_TWIST_DEG``  (default ``30``) — chelate δ/λ twist degrees
``DELFIN_5O_MD_TOL``             (default ``0.05``) — Å, M–D invariant tolerance
``DELFIN_5O_CLASH_HH``           (default ``1.7``) — Å, H–H clash cutoff
``DELFIN_5O_CLASH_XH``           (default ``1.7``) — Å, X–H clash cutoff
``DELFIN_5O_CLASH_XX``           (default ``1.5``) — Å, heavy-heavy clash cutoff

Welle-5p-A hard-gate (default ON when ``DELFIN_5P_A_TOPOLOGY_HARDGATE=1``)
Every layer's candidate is additionally passed through
:func:`delfin._topology_hash.topology_preserved` before the energy
filter.  This rejects chelate-ring twists that flip an amine-H toward
the metal or perturb M–D edges beyond ``DELFIN_5P_A_MD_TOL``.
"""

from __future__ import annotations

import math
import os
from typing import Dict, List, Optional, Sequence, Tuple

from delfin.common.logging import get_logger

# Reuse helpers from the T6 rotamer-diversity module to keep the M–D
# invariant, topology-hash and DOF-detection logic single-sourced.
from delfin import _rotamer_diversity as _rot

logger = get_logger(__name__)


# ---------------------------------------------------------------------------
# Environment configuration helpers (same style as _rotamer_diversity)
# ---------------------------------------------------------------------------


def _env_bool(name: str, default: bool = False) -> bool:
    raw = os.environ.get(name)
    if raw is None:
        return default
    return raw.strip().lower() in ("1", "true", "yes", "on")


def _env_int(name: str, default: int, lo: int = 1, hi: int = 4096) -> int:
    raw = os.environ.get(name)
    if raw is None:
        return default
    try:
        return max(lo, min(hi, int(raw)))
    except (TypeError, ValueError):
        return default


def _env_float(name: str, default: float, lo: float = 0.0, hi: float = 10.0) -> float:
    raw = os.environ.get(name)
    if raw is None:
        return default
    try:
        return max(lo, min(hi, float(raw)))
    except (TypeError, ValueError):
        return default


def _is_enabled() -> bool:
    # Either the legacy 5o master switch OR the conformer-COVERAGE master switch
    # (CONFORMER_REACHABILITY_2026_06_17) activates the per-isomer pool.
    return _env_bool("DELFIN_5O_CONFORMER_POOL", False) or _coverage_enabled()


def _coverage_enabled() -> bool:
    """Master switch for the gas-phase backbone-CONFORMER-coverage offensive
    (CONFORMER_REACHABILITY_2026_06_17).  Default ``0`` -> byte-identical: the
    pool path is untouched (legacy 5o behaviour, off by default).  When ``1``
    the pool runs in coverage mode: tolerant XYZ parse, backbone-heavy-moving
    torsion DOFs (bug-2), dense Cremer-Pople ring-pucker templates (bug-3),
    denser torsion grid, and a hard absolute core-freeze guard (+-0.05 A)."""
    return _env_bool("DELFIN_FFFREE_CONFORMER_COVERAGE", False)


# ---------------------------------------------------------------------------
# Geometry helpers (pure Python — no numpy dependency, mirroring _rotamer_diversity)
# ---------------------------------------------------------------------------


Coord = Tuple[float, float, float]


def _vec_sub(a: Coord, b: Coord) -> Coord:
    return (a[0] - b[0], a[1] - b[1], a[2] - b[2])


def _vec_add(a: Coord, b: Coord) -> Coord:
    return (a[0] + b[0], a[1] + b[1], a[2] + b[2])


def _vec_scale(a: Coord, s: float) -> Coord:
    return (a[0] * s, a[1] * s, a[2] * s)


def _vec_dot(a: Coord, b: Coord) -> float:
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]


def _vec_cross(a: Coord, b: Coord) -> Coord:
    return (
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    )


def _vec_norm(a: Coord) -> float:
    return math.sqrt(a[0] * a[0] + a[1] * a[1] + a[2] * a[2])


def _vec_unit(a: Coord) -> Coord:
    n = _vec_norm(a)
    if n < 1e-12:
        return (0.0, 0.0, 0.0)
    return (a[0] / n, a[1] / n, a[2] / n)


def _centroid(coords: Sequence[Coord], indices: Sequence[int]) -> Coord:
    if not indices:
        return (0.0, 0.0, 0.0)
    sx = sy = sz = 0.0
    for i in indices:
        x, y, z = coords[i]
        sx += x
        sy += y
        sz += z
    n = float(len(indices))
    return (sx / n, sy / n, sz / n)


def _ring_plane_normal(
    coords: Sequence[Coord], ring_atoms: Sequence[int]
) -> Tuple[Coord, Coord]:
    """Return ``(centroid, unit_normal)`` of the best-fit plane through ring atoms.

    Uses a cheap covariance-eigenvector approximation (3×3 symmetric power
    iteration on the *minor* axis — the eigenvector with smallest variance is
    the plane normal).  Avoids any external linear-algebra dependency.
    """
    c = _centroid(coords, ring_atoms)
    # Build 3×3 covariance
    cov = [[0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0]]
    for i in ring_atoms:
        x, y, z = coords[i]
        dx = x - c[0]
        dy = y - c[1]
        dz = z - c[2]
        cov[0][0] += dx * dx
        cov[0][1] += dx * dy
        cov[0][2] += dx * dz
        cov[1][1] += dy * dy
        cov[1][2] += dy * dz
        cov[2][2] += dz * dz
    cov[1][0] = cov[0][1]
    cov[2][0] = cov[0][2]
    cov[2][1] = cov[1][2]

    # Power iteration on the *inverse* problem: smallest eigenvector.
    # We approximate by power-iterating on (trace*I - cov) which has the
    # smallest eigenvalue of cov as its largest.
    trace = cov[0][0] + cov[1][1] + cov[2][2]
    a = [[trace - cov[i][j] if i == j else -cov[i][j] for j in range(3)] for i in range(3)]
    v: Coord = (1.0, 1.0, 1.0)
    for _ in range(30):
        nv = (
            a[0][0] * v[0] + a[0][1] * v[1] + a[0][2] * v[2],
            a[1][0] * v[0] + a[1][1] * v[1] + a[1][2] * v[2],
            a[2][0] * v[0] + a[2][1] * v[1] + a[2][2] * v[2],
        )
        n = _vec_norm(nv)
        if n < 1e-12:
            break
        v = (nv[0] / n, nv[1] / n, nv[2] / n)
    return c, _vec_unit(v)


# ---------------------------------------------------------------------------
# Rotation around an arbitrary axis (re-export Rodrigues from _rotamer_diversity)
# ---------------------------------------------------------------------------


def _rotate_atoms(
    coords: List[Coord],
    origin: Coord,
    axis_dir: Coord,
    angle_rad: float,
    atom_indices: Sequence[int],
) -> List[Coord]:
    return _rot._rodrigues_rotate(coords, origin, axis_dir, angle_rad, atom_indices)


# ---------------------------------------------------------------------------
# Layer 1 — Torsion-grid (reuses T6 rotamer DOFs)
# ---------------------------------------------------------------------------


def _layer1_torsion_candidates(
    symbols: Sequence[str],
    base_coords: List[Coord],
    graph: Dict,
    n_states: int,
    max_dofs: int,
    grid_cap: int,
    coverage: bool = False,
) -> List[Tuple[List[Coord], str]]:
    """Sample staggered torsion conformers around every rotatable DOF.

    Returns a list of ``(coords, mode_tag)`` candidates *excluding* the
    identity combination (which is the input ``base_coords``).

    ``coverage=True`` selects backbone-heavy-moving torsions only (drops the
    methyl/terminal rotors that produce 0 backbone diversity) — bug-2 fix.
    """
    dofs = _rot.identify_rotamer_dofs(graph, max_dofs=max_dofs, coverage=coverage)
    if not dofs:
        return []

    step_rad = 2.0 * math.pi / float(max(2, n_states))
    out: List[Tuple[List[Coord], str]] = []
    for combo in _rot._grid_iter(len(dofs), n_states, grid_cap):
        if all(s == 0 for s in combo):
            continue
        coords = list(base_coords)
        for dof_idx, state in enumerate(combo):
            if state == 0:
                continue
            dof = dofs[dof_idx]
            anchor = dof["anchor"]
            pivot = dof["pivot"]
            ox, oy, oz = coords[anchor]
            px, py, pz = coords[pivot]
            axis = (px - ox, py - oy, pz - oz)
            coords = _rotate_atoms(
                coords, (ox, oy, oz), axis, step_rad * state, dof["rotating"]
            )
        out.append((coords, f"torsion-{'.'.join(str(s) for s in combo)}"))
    return out


# ---------------------------------------------------------------------------
# Layer 2 — Ring-pucker (universal — every non-aromatic, non-metal ring of size 5..6)
# ---------------------------------------------------------------------------


def _rings_from_graph(graph: Dict) -> List[List[int]]:
    """Return SSSR-like ring atom lists from the OB-perceived bonds.

    We collect the bond-level ``is_ring`` flags from
    :func:`_rotamer_diversity._graph_from_ob` and do a simple
    cycle-finding via DFS on the ring-bond subgraph.  Pure Python — no
    extra dependencies.
    """
    ring_bonds = []
    for (a, b, _order, _arom, is_ring) in graph.get("bonds", []):
        if is_ring:
            ring_bonds.append((a, b))
    if not ring_bonds:
        return []

    # Build ring-subgraph adjacency
    adj: Dict[int, List[int]] = {}
    for (a, b) in ring_bonds:
        adj.setdefault(a, []).append(b)
        adj.setdefault(b, []).append(a)

    # Find all simple cycles up to length 12 — universal across mono- /
    # bi- / tri- / tetra-dentate chelates and macrocyclic ligands.
    rings: List[List[int]] = []
    seen_ring_keys: set = set()
    nodes = sorted(adj.keys())

    def _dfs(start: int, cur: int, path: List[int], depth: int):
        if depth > 16:  # cap ring size at 16 — macrocycle ceiling
            return
        for nbr in adj.get(cur, []):
            if nbr == start and len(path) >= 3:
                # Found cycle of length len(path)
                key = tuple(sorted(path))
                if key not in seen_ring_keys:
                    seen_ring_keys.add(key)
                    rings.append(list(path))
                continue
            if nbr in path:
                continue
            path.append(nbr)
            _dfs(start, nbr, path, depth + 1)
            path.pop()

    for n in nodes:
        _dfs(n, n, [n], 0)

    # Deduplicate: keep only smallest unique rings (SSSR-like)
    # We sort by size and only keep rings whose vertices aren't fully
    # covered by a smaller already-accepted ring.
    rings.sort(key=lambda r: (len(r), r))
    accepted: List[List[int]] = []
    for r in rings:
        rset = set(r)
        # only suppress when an EXISTING smaller ring covers all atoms
        if any(set(a).issuperset(rset) for a in accepted if len(a) <= len(r)):
            continue
        accepted.append(r)
    return accepted


def _layer2_ring_pucker_candidates(
    symbols: Sequence[str],
    base_coords: List[Coord],
    graph: Dict,
    rings: Sequence[Sequence[int]],
    amplitude: float,
) -> List[Tuple[List[Coord], str]]:
    """For each saturated non-aromatic 5- or 6-ring, generate two pucker modes.

    Universal: triggers on any ring with size in [4, 7] that does NOT
    contain any aromatic bond and does NOT contain a metal atom.  For
    each such ring we displace two opposite ring atoms by ±amplitude
    along the ring normal — this is the chair-vs-twist mode pair used
    in standard cyclohexane conformational analysis.
    """
    atomic_nums = graph["atomic_nums"]
    is_metal = graph["is_metal"]
    bonds = graph.get("bonds", [])
    aromatic_bond_set = set()
    bonded_h: Dict[int, List[int]] = {}
    for (a, b, _o, arom, _ring) in bonds:
        if arom:
            aromatic_bond_set.add((min(a, b), max(a, b)))
        # Iter-18 rigid-H drag (per-atom forensik 2026-05-18): without
        # this, pucker displacement of ring heavy atoms leaves bonded H
        # behind → CH3 umbrella collapse 54.78% pool-wide (vs 123a130
        # 21.98% baseline) — 100-500× violations per frame.
        if symbols[a] == "H" and symbols[b] != "H":
            bonded_h.setdefault(b, []).append(a)
        elif symbols[b] == "H" and symbols[a] != "H":
            bonded_h.setdefault(a, []).append(b)

    out: List[Tuple[List[Coord], str]] = []
    for ring in rings:
        size = len(ring)
        if size < 4 or size > 7:
            continue
        if any(is_metal[i] for i in ring):
            continue  # chelate rings handled in Layer 3
        # Must be saturated (no aromatic bond between any consecutive
        # ring atoms).  We check both ring-edge bonds.
        is_aromatic = False
        for k in range(size):
            a = ring[k]
            b = ring[(k + 1) % size]
            if (min(a, b), max(a, b)) in aromatic_bond_set:
                is_aromatic = True
                break
        if is_aromatic:
            continue
        # Ring atoms must be majority sp3-like (heavy atomic number ≥ 6
        # for C and analogues).  Skip pure hetero or single-element rings
        # that aren't relevant (e.g. all-N).
        n_sp3_like = sum(1 for i in ring if atomic_nums[i] in (6, 7, 8, 15, 16))
        if n_sp3_like < max(3, size - 1):
            continue

        c, normal = _ring_plane_normal(base_coords, ring)
        if _vec_norm(normal) < 1e-9:
            continue

        # Two pucker modes: ±amplitude on opposite atoms.  Pick the pair
        # at indices (0, size//2) — the chemically meaningful pair for
        # chair-vs-boat and twist-vs-half-chair sampling.
        i0 = ring[0]
        i1 = ring[size // 2]

        for sign_a, sign_b, tag in (
            (+1.0, -1.0, "chair"),
            (-1.0, +1.0, "boat"),
            (+1.0, +1.0, "dome"),
        ):
            coords = list(base_coords)
            disp_a = _vec_scale(normal, amplitude * sign_a)
            disp_b = _vec_scale(normal, amplitude * sign_b)
            coords[i0] = _vec_add(coords[i0], disp_a)
            for h_idx in bonded_h.get(i0, ()):
                coords[h_idx] = _vec_add(coords[h_idx], disp_a)
            coords[i1] = _vec_add(coords[i1], disp_b)
            for h_idx in bonded_h.get(i1, ()):
                coords[h_idx] = _vec_add(coords[h_idx], disp_b)
            out.append((coords, f"pucker-r{size}-{tag}"))
    return out


# ---------------------------------------------------------------------------
# Layer 3 — Chelate-twist (universal — rings containing the metal + sp3 backbone)
# ---------------------------------------------------------------------------


def _layer3_chelate_twist_candidates(
    symbols: Sequence[str],
    base_coords: List[Coord],
    graph: Dict,
    rings: Sequence[Sequence[int]],
    twist_deg: float,
) -> List[Tuple[List[Coord], str]]:
    """Generate δ / λ twist conformers for chelate rings.

    Universal trigger:

    * ring contains exactly one metal atom (otherwise it's not a chelate
      in the standard sense),
    * ring size ≤ 7 (5- and 6-membered chelates predominantly, 7-rings
      for diamine / diether / acac-extended cases),
    * ring contains at least one non-aromatic single-bond pair of
      adjacent atoms (the dihedral candidates).

    For each qualifying ring we generate two conformers: a δ-twist
    (positive rotation around the *non-metal* axis of two opposite
    backbone atoms) and a λ-twist (negative rotation).  The rotation
    axis is the line between the two backbone atoms that are *not*
    bonded directly to the metal.
    """
    atomic_nums = graph["atomic_nums"]
    is_metal = graph["is_metal"]
    bonds = graph.get("bonds", [])
    aromatic_bond_set = set()
    for (a, b, _o, arom, _r) in bonds:
        if arom:
            aromatic_bond_set.add((min(a, b), max(a, b)))

    twist_rad = math.radians(twist_deg)
    out: List[Tuple[List[Coord], str]] = []

    for ring in rings:
        size = len(ring)
        if size < 4 or size > 7:
            continue
        metals_in_ring = [i for i in ring if is_metal[i]]
        if len(metals_in_ring) != 1:
            continue  # not a single-metal chelate ring
        m_idx = metals_in_ring[0]
        # Backbone atoms = ring atoms not bonded directly to the metal
        backbone = []
        for i in ring:
            if i == m_idx:
                continue
            # is i adjacent to m_idx in the ring?
            pos_in_ring = ring.index(i)
            m_pos = ring.index(m_idx)
            ring_n = size
            is_adj_to_m = (
                pos_in_ring == (m_pos + 1) % ring_n
                or pos_in_ring == (m_pos - 1) % ring_n
            )
            if not is_adj_to_m:
                backbone.append(i)
        if len(backbone) < 2:
            continue
        # Pick the two backbone atoms farthest from the metal in the ring
        # (= the "outer" pair for the chelate twist).
        # Use ring-position-distance to metal.
        m_pos = ring.index(m_idx)
        def _ring_dist(a_idx: int) -> int:
            ap = ring.index(a_idx)
            d = (ap - m_pos) % size
            return min(d, size - d)
        backbone.sort(key=_ring_dist, reverse=True)
        a_idx = backbone[0]
        b_idx = backbone[1]

        # Edge between a_idx and b_idx must be single, non-aromatic — but
        # they're allowed to be non-adjacent (in a 6-ring chelate the
        # backbone has 4 atoms; we pick the central two).
        pair = (min(a_idx, b_idx), max(a_idx, b_idx))
        if pair in aromatic_bond_set:
            continue

        axis_origin = base_coords[a_idx]
        axis_dir = _vec_sub(base_coords[b_idx], base_coords[a_idx])
        if _vec_norm(axis_dir) < 1e-6:
            continue

        # Atoms to rotate = one side of the ring (BFS in the molecular
        # graph from a_idx without crossing the bond a_idx–b_idx and
        # without crossing the metal — keep the metal/donor side fixed).
        neighbours = graph["neighbours"]

        def _bfs_side(start: int, blocked: int, banned_metal: bool) -> List[int]:
            visited = {start, blocked}
            if banned_metal:
                visited.add(m_idx)
            out_list = [start]
            stack = [start]
            while stack:
                cur = stack.pop()
                for nbr in neighbours[cur]:
                    if nbr in visited:
                        continue
                    visited.add(nbr)
                    out_list.append(nbr)
                    stack.append(nbr)
            return out_list

        side_atoms = _bfs_side(a_idx, b_idx, banned_metal=True)
        if any(is_metal[i] for i in side_atoms):
            # We hit a metal — abort to keep M–D invariant safe
            continue
        if len(side_atoms) < 2:
            continue

        for sign, tag in ((+1.0, "delta"), (-1.0, "lambda")):
            coords = list(base_coords)
            coords = _rotate_atoms(
                coords, axis_origin, axis_dir, sign * twist_rad, side_atoms
            )
            out.append((coords, f"chelate-{tag}-r{size}-m{m_idx}"))
    return out


# ---------------------------------------------------------------------------
# Layer 4 — Macrocycle modes (universal — ring size ≥ 7 atoms)
# ---------------------------------------------------------------------------


def _layer4_macrocycle_candidates(
    symbols: Sequence[str],
    base_coords: List[Coord],
    graph: Dict,
    rings: Sequence[Sequence[int]],
    amplitude: float,
) -> List[Tuple[List[Coord], str]]:
    """Sample macrocyclic deformation modes (saddle, ruffled, dome).

    Universal trigger: any ring of size ≥ 7 — this covers porphyrin
    24-atom inner ring (when reduced to the 16-atom core or considered
    as the full 24-membered macrocycle), corrin, calix-arenes, and
    any large ligand backbone.  Excludes pure-aromatic rings only if
    every ring bond is aromatic — partial-aromatic macrocycles like
    porphyrin remain valid since the deformation modes are physically
    real (porphyrin saddle / ruffled / dome are textbook).

    The deformation is a low-order Fourier mode along the ring normal:

    * saddle  — :math:`\\cos(2\\theta_k)` displacement per ring atom k
    * ruffled — :math:`\\cos(4\\theta_k)`
    * dome    — uniform offset of an *inner* subset (mean of ring) along
      the normal (zero-mode + radial gradient)

    For metal-containing macrocycles (porphyrin Fe / Co / Mn) the metal
    is treated as a normal ring atom — it deforms with the ring,
    preserving M–N bonds since neighbours follow.  M–D invariant guard
    enforces the bond-length tolerance downstream.
    """
    out: List[Tuple[List[Coord], str]] = []
    bonds = graph.get("bonds", [])
    aromatic_bond_set = set()
    for (a, b, _o, arom, _r) in bonds:
        if arom:
            aromatic_bond_set.add((min(a, b), max(a, b)))

    for ring in rings:
        size = len(ring)
        if size < 7 or size > 30:
            continue
        # We allow partial-aromatic macrocycles (porphyrin) — the
        # deformation modes are physically defined as collective
        # rigid-displacement of every ring atom along the ring normal,
        # which preserves bond lengths to 2nd order.
        # But pure-aromatic rings (e.g. 18-annulene) we skip to avoid
        # breaking sp2 planarity.
        n_aromatic_bonds = 0
        for k in range(size):
            a = ring[k]
            b = ring[(k + 1) % size]
            if (min(a, b), max(a, b)) in aromatic_bond_set:
                n_aromatic_bonds += 1
        if n_aromatic_bonds >= size:  # fully aromatic — skip
            continue

        c, normal = _ring_plane_normal(base_coords, ring)
        if _vec_norm(normal) < 1e-9:
            continue

        # Establish angular ordering θ_k around the centroid using the
        # projection of (atom – centroid) onto a tangent basis.  Pick a
        # reference radial vector from centroid to the first ring atom,
        # build orthonormal in-plane basis (e1, e2).
        r0_vec = _vec_sub(base_coords[ring[0]], c)
        r0_in_plane = _vec_sub(r0_vec, _vec_scale(normal, _vec_dot(r0_vec, normal)))
        e1 = _vec_unit(r0_in_plane)
        if _vec_norm(e1) < 1e-9:
            continue
        e2 = _vec_cross(normal, e1)
        e2 = _vec_unit(e2)
        if _vec_norm(e2) < 1e-9:
            continue

        angles: List[float] = []
        for atom_idx in ring:
            r = _vec_sub(base_coords[atom_idx], c)
            r_in_plane = _vec_sub(r, _vec_scale(normal, _vec_dot(r, normal)))
            x = _vec_dot(r_in_plane, e1)
            y = _vec_dot(r_in_plane, e2)
            angles.append(math.atan2(y, x))

        for n_fold, tag in ((2, "saddle"), (4, "ruffled"), (0, "dome")):
            coords = list(base_coords)
            for k, atom_idx in enumerate(ring):
                if n_fold == 0:
                    # Dome: only displace inner atoms (every other atom)
                    # to avoid pure translation — universal heuristic.
                    if k % 2 != 0:
                        continue
                    disp_amp = amplitude
                else:
                    disp_amp = amplitude * math.cos(n_fold * angles[k])
                # Drag bonded hydrogens with the ring atom (rigid-H,
                # universal — fixes the Welle-5l "B5 v2 H stretched"
                # bug per project_rigid_h_tracking_b5_bug).
                disp = _vec_scale(normal, disp_amp)
                coords[atom_idx] = _vec_add(coords[atom_idx], disp)
                for nbr in graph["neighbours"][atom_idx]:
                    if graph["atomic_nums"][nbr] == 1:
                        coords[nbr] = _vec_add(coords[nbr], disp)
            out.append((coords, f"macrocycle-{tag}-r{size}"))
    return out


# ---------------------------------------------------------------------------
# Quality filter — clash + M–D invariant + topology hash
# ---------------------------------------------------------------------------


def _core_frozen(
    graph: Dict,
    base_coords: Sequence[Coord],
    new_coords: Sequence[Coord],
    tol: float,
) -> bool:
    """Hard core-invariant guard (#82/#100): the metal AND every donor atom
    (first-shell heavy neighbour of any metal) must stay at their NATIVE
    absolute positions within *tol* Å.  Stricter than the M-D bond-length
    invariant (which only checks distances) — this pins the coordination core
    in space so the backbone re-folds around a fixed polyhedron.

    Returns True iff the absolute displacement of every core atom is ≤ tol.
    """
    n = graph.get("n_atoms", 0)
    if n == 0 or n != len(new_coords):
        return False
    atomic_nums = graph["atomic_nums"]
    neighbours = graph["neighbours"]
    is_metal = graph["is_metal"]
    core: set = set()
    for m_idx in range(n):
        if not is_metal[m_idx]:
            continue
        core.add(m_idx)
        for d_idx in neighbours[m_idx]:
            if atomic_nums[d_idx] != 1:
                core.add(d_idx)
    tol2 = tol * tol
    for i in core:
        dx = new_coords[i][0] - base_coords[i][0]
        dy = new_coords[i][1] - base_coords[i][1]
        dz = new_coords[i][2] - base_coords[i][2]
        if dx * dx + dy * dy + dz * dz > tol2:
            return False
    return True


def _no_clash(
    symbols: Sequence[str],
    coords: Sequence[Coord],
    graph: Dict,
    cut_hh: float,
    cut_xh: float,
    cut_xx: float,
) -> bool:
    """Return True iff no non-bonded atom pair is closer than its cutoff.

    Bonded 1,2 and 1,3 pairs are exempt — they are VSEPR-determined
    geometric neighbours and would always trigger.
    """
    n = len(symbols)
    bonded: set = set()
    for (a, b, _o, _arom, _r) in graph.get("bonds", []):
        bonded.add((min(a, b), max(a, b)))
    # 1,3 pairs
    atom_bonds: Dict[int, List[int]] = {i: [] for i in range(n)}
    for (i, j) in bonded:
        atom_bonds[i].append(j)
        atom_bonds[j].append(i)
    one_three: set = set()
    for mid in range(n):
        nbrs = atom_bonds[mid]
        for ii in range(len(nbrs)):
            for jj in range(ii + 1, len(nbrs)):
                a, b = nbrs[ii], nbrs[jj]
                one_three.add((min(a, b), max(a, b)))

    for i in range(n):
        for j in range(i + 1, n):
            pair = (i, j)
            if pair in bonded or pair in one_three:
                continue
            dx = coords[i][0] - coords[j][0]
            dy = coords[i][1] - coords[j][1]
            dz = coords[i][2] - coords[j][2]
            d2 = dx * dx + dy * dy + dz * dz
            si = symbols[i]
            sj = symbols[j]
            if si == "H" and sj == "H":
                cut = cut_hh
            elif si == "H" or sj == "H":
                cut = cut_xh
            else:
                cut = cut_xx
            if d2 < cut * cut:
                return False
    return True


# ---------------------------------------------------------------------------
# Pairwise RMSD (heavy atoms only, no alignment — coords are already aligned
# because they all derive from the same base XYZ by local displacements)
# ---------------------------------------------------------------------------


def _heavy_rmsd(
    symbols: Sequence[str],
    a: Sequence[Coord],
    b: Sequence[Coord],
) -> float:
    n = len(symbols)
    heavy_indices = [i for i in range(n) if symbols[i] != "H"]
    if not heavy_indices:
        heavy_indices = list(range(n))
    s = 0.0
    for i in heavy_indices:
        dx = a[i][0] - b[i][0]
        dy = a[i][1] - b[i][1]
        dz = a[i][2] - b[i][2]
        s += dx * dx + dy * dy + dz * dz
    return math.sqrt(s / float(len(heavy_indices)))


# ---------------------------------------------------------------------------
# Top-K diverse selection
# ---------------------------------------------------------------------------


def _select_top_k_diverse(
    candidates: List[Tuple[float, str, str, List[Coord]]],
    base_xyz: str,
    base_coords: List[Coord],
    symbols: Sequence[str],
    k: int,
    rmsd_min: float,
) -> List[Tuple[str, str]]:
    """Greedy top-K diverse selection.

    The base XYZ is always position 0 (with energy ranking among
    candidates).  We then iterate candidates by energy ascending and
    accept each that has RMSD ≥ ``rmsd_min`` vs. every already-selected
    pool member.  Stops when k pool members are accepted.

    Returns ``[(xyz_string, mode_tag), ...]`` with length ≤ k.
    """
    if k < 1:
        return [(base_xyz, "base")]
    candidates.sort(key=lambda c: c[0])  # by energy ascending
    selected: List[Tuple[str, str, List[Coord]]] = [(base_xyz, "base", base_coords)]
    for energy, xyz_str, tag, coords in candidates:
        if len(selected) >= k:
            break
        if xyz_str == base_xyz:
            continue
        accept = True
        for (_xs, _tg, sc) in selected:
            if _heavy_rmsd(symbols, coords, sc) < rmsd_min:
                accept = False
                break
        if accept:
            selected.append((xyz_str, tag, coords))
    return [(xs, tg) for (xs, tg, _c) in selected]


# ---------------------------------------------------------------------------
# Main entry — build the conformer pool for a single XYZ
# ---------------------------------------------------------------------------


def generate_conformer_pool(
    xyz: str,
    k_target: int = 10,
    diversity_rmsd_min: float = 0.5,
    n_torsion_states: int = 3,
    max_dofs: int = 6,
    torsion_grid_cap: int = 64,
    ring_pucker_ampl: float = 0.35,
    macrocycle_ampl: float = 0.40,
    chelate_twist_deg: float = 30.0,
    md_tol: float = 0.05,
    clash_hh: float = 1.7,
    clash_xh: float = 1.7,
    clash_xx: float = 1.5,
    coverage: bool = False,
    core_freeze_tol: float = 0.05,
) -> List[Tuple[str, str]]:
    """Return ``[(xyz, mode_tag), ...]`` of length ≤ ``k_target``.

    ``(base_xyz, "base")`` is always position 0 if ``k_target ≥ 1``.

    If the molecule has no relevant DOFs (e.g. ``[Pt(NH3)4]`` — rigid
    tetraamminecation) the function returns ``[(base_xyz, "base")]``
    only.  This satisfies the "rigid → K=1" contract.

    The function never crashes — on any internal failure it falls back
    to ``[(base_xyz, "base")]``.
    """
    if k_target < 1:
        return [(xyz, "base")]
    # Build the graph + parse coords once
    ob_mol = _rot._build_ob_mol_from_xyz(xyz)
    if ob_mol is None:
        return [(xyz, "base")]
    graph = _rot._graph_from_ob(ob_mol)
    if not graph:
        return [(xyz, "base")]
    try:
        symbols, base_coords_t = _rot._parse_delfin_xyz(xyz)
    except Exception:
        return [(xyz, "base")]
    if len(symbols) != graph["n_atoms"]:
        return [(xyz, "base")]
    base_coords: List[Coord] = [tuple(c) for c in base_coords_t]
    base_topo = _rot._topology_hash(graph)

    # Detect rings (universal — RDKit-free)
    rings = _rings_from_graph(graph)

    # ---------- Layered sampling ----------
    raw_candidates: List[Tuple[List[Coord], str]] = []
    try:
        raw_candidates.extend(
            _layer1_torsion_candidates(
                symbols,
                base_coords,
                graph,
                n_torsion_states,
                max_dofs,
                torsion_grid_cap,
                coverage=coverage,
            )
        )
    except Exception as exc:
        logger.debug("5o Layer-1 (torsion) failed: %s", exc)
    try:
        # Welle-5p-B: chemistry-template ring-conformer generator
        # supersedes the naive 2-atom Layer-2 pucker when its env-flag
        # is set.  Default-OFF byte-identical fall-through.  Coverage mode
        # (bug-3, dense Cremer-Pople-style pucker) force-enables it so the
        # ring-bound backbone (azamacrocycle/chelate-backbone/saturated rings
        # where the rotatable-bond detector finds 0 backbone torsions) is
        # actually re-folded — the chair/boat/twist/envelope + macrocycle
        # saddle/ruffle/dome templates are the local-minima the torsion grid
        # cannot reach (ILUFOK/NUMBIH/SEHXOU/KADYIA-class).
        use_5p_b_templates = coverage or _env_bool(
            "DELFIN_5P_B_RING_TEMPLATES", False
        )
        if use_5p_b_templates:
            from delfin import _ring_conformer_templates as _ring_tpl

            k_per_ring = _env_int(
                "DELFIN_5P_B_K_PER_RING", 6 if coverage else 3, lo=1, hi=20
            )
            max_var = _env_int(
                "DELFIN_5P_B_MAX_RING_VARIANTS", 24 if coverage else 6,
                lo=1, hi=64
            )
            amp_frac = _env_float(
                "DELFIN_5P_B_AMPLITUDE_FRACTION",
                0.50,
                lo=0.0,
                hi=2.0,
            )
            raw_candidates.extend(
                _ring_tpl.layer2_ring_template_candidates(
                    symbols,
                    base_coords,
                    graph,
                    rings,
                    amp_frac,
                    k_per_ring,
                    max_var,
                )
            )
        else:
            raw_candidates.extend(
                _layer2_ring_pucker_candidates(
                    symbols, base_coords, graph, rings, ring_pucker_ampl
                )
            )
    except Exception as exc:
        logger.debug("5o Layer-2 (ring pucker) failed: %s", exc)
    try:
        raw_candidates.extend(
            _layer3_chelate_twist_candidates(
                symbols, base_coords, graph, rings, chelate_twist_deg
            )
        )
    except Exception as exc:
        logger.debug("5o Layer-3 (chelate twist) failed: %s", exc)
    try:
        raw_candidates.extend(
            _layer4_macrocycle_candidates(
                symbols, base_coords, graph, rings, macrocycle_ampl
            )
        )
    except Exception as exc:
        logger.debug("5o Layer-4 (macrocycle) failed: %s", exc)

    if not raw_candidates:
        return [(xyz, "base")]

    # Welle-5p-A: precompute hard-gate parameters (read env vars once)
    # Default reverted 1 -> 0 on 2026-05-18 (Iter-17) — voll-pool b5defcd
    # showed pool-wide sigma -2783 isomere over-rejection.  Per-class
    # adaptive thresholds deferred to Iter-18.  Env-flag opt-in preserved.
    _hardgate_on = _env_bool("DELFIN_5P_A_TOPOLOGY_HARDGATE", False)
    if _hardgate_on:
        _amine_min_deg = _env_float(
            "DELFIN_5P_A_AMINE_H_MIN_DEG", 60.0, lo=0.0, hi=180.0
        )
        _hm_min = _env_float("DELFIN_5P_A_CLASH_HM_MIN", 2.30, lo=0.0, hi=10.0)
        _bond_tol_th = _env_float("DELFIN_5P_A_BOND_TOL", 1.30, lo=0.5, hi=3.0)
        try:
            from delfin import _topology_hash as _th  # local import
        except Exception:
            _th = None  # type: ignore
            _hardgate_on = False
    else:
        _th = None  # type: ignore

    # ---------- Quality filter ----------
    # Under-yield fix (DELFIN_CONF_RELAX_CANDIDATES, default OFF): a rigidly-rotated
    # torsion candidate clashes grossly in the packed coordination sphere (measured
    # worst contacts 0.06-0.5 A), so the clash gate below rejects ~100% of them and
    # only the base survives.  Relieve the clash first with a core-frozen spring
    # relaxation (metal + every donor pinned), then gate the RELAXED frame.  Geometry
    # -only + deterministic; byte-identical when off (coords pass through untouched).
    _relax_on = os.environ.get("DELFIN_CONF_RELAX_CANDIDATES", "0") == "1"
    _cr = None
    _relax_frozen = None
    if _relax_on:
        try:
            from delfin.fffree import _conf_relax as _cr
            _relax_frozen = _cr.core_indices(graph)
        except Exception:
            _relax_on = False

    scored_candidates: List[Tuple[float, str, str, List[Coord]]] = []
    for (coords, tag) in raw_candidates:
        try:
            if _relax_on and _relax_frozen is not None and _cr is not None:
                coords = _cr.relax_clashes(
                    symbols, coords, graph, _relax_frozen,
                    clash_hh, clash_xh, clash_xx,
                )
            if not _rot._coord_bond_invariant_holds(
                graph, base_coords, coords, tol=md_tol
            ):
                continue
            # Coverage hard core-freeze guard (#82/#100): metal + every donor
            # pinned to native absolute positions within core_freeze_tol.  Any
            # frame whose coordination core drifts is REJECTED (the backbone
            # must re-fold around a FIXED polyhedron, not move the core).
            if coverage and not _core_frozen(
                graph, base_coords, coords, core_freeze_tol
            ):
                continue
            if not _no_clash(
                symbols, coords, graph, clash_hh, clash_xh, clash_xx
            ):
                continue
            cand_xyz = _rot._format_delfin_xyz(symbols, coords)
            cand_mol = _rot._build_ob_mol_from_xyz(cand_xyz)
            if cand_mol is None:
                continue
            cand_graph = _rot._graph_from_ob(cand_mol)
            if _rot._topology_hash(cand_graph) != base_topo:
                continue
            # Welle-5p-A universal hard-gate — single source of truth
            if _hardgate_on and _th is not None:
                _gate_res = _th.topology_preserved(
                    symbols,
                    base_coords,
                    coords,
                    md_tol=md_tol,
                    amine_h_min_deg=_amine_min_deg,
                    h_to_metal_min_dist=_hm_min,
                    bond_tol=_bond_tol_th,
                )
                if not _gate_res.passed:
                    logger.debug(
                        "5p-A pool reject (%s): %s",
                        tag,
                        _gate_res.violations[:3],
                    )
                    continue
            energy = _rot._evaluate_xyz(cand_xyz)
            if energy is None:
                continue
            if not math.isfinite(energy):
                continue
            scored_candidates.append((energy, cand_xyz, tag, coords))
        except Exception as exc:
            logger.debug("5o candidate skipped: %s", exc)
            continue

    if not scored_candidates:
        return [(xyz, "base")]

    # ---------- GFN-FF re-rank (DELFIN_CONF_GFNFF_RANK, default OFF) ----------
    # The UFF scoring above is only a cheap pre-filter: on transition-metal complexes
    # OpenBabel-UFF energies are unusable for ranking (measured 170-300 kcal/mol
    # spreads + 1e16 blow-ups where the physical conformer window is a few kcal/mol).
    # When enabled, re-score the lowest-UFF shortlist with GFN-FF (xtb CLI,
    # license-clean, parametrised for metals; measured 8 kcal/mol spread on a pool
    # UFF placed at 301) and rank by that.  Byte-identical when off (block skipped).
    if os.environ.get("DELFIN_CONF_GFNFF_RANK", "0") == "1":
        try:
            from delfin.fffree import _gfnff_rank as _gff
            _top_m = _env_int(
                "DELFIN_CONF_GFNFF_TOPM", max(3 * k_target, 12), lo=1, hi=256
            )
            _chg = _env_int("DELFIN_CONF_GFNFF_CHARGE", 0, lo=-12, hi=12)
            _reranked = _gff.rerank(scored_candidates, _top_m, charge=_chg)
            if _reranked:
                scored_candidates = _reranked
        except Exception as exc:
            logger.debug("GFN-FF rerank skipped: %s", exc)

    # ---------- Top-K diverse selection ----------
    pool = _select_top_k_diverse(
        scored_candidates,
        xyz,
        base_coords,
        symbols,
        k_target,
        diversity_rmsd_min,
    )
    return pool


# ---------------------------------------------------------------------------
# Public wire-in helper used by smiles_converter
# ---------------------------------------------------------------------------


def apply_if_enabled(xyz: str) -> List[Tuple[str, str]]:
    """Wire-in entry: returns ``[(xyz, mode_tag), ...]`` when the env-flag is set.

    Default-OFF semantics: if ``DELFIN_5O_CONFORMER_POOL`` is unset,
    returns ``[(xyz, "base")]`` unchanged.  The wire-in caller iterates
    the returned list and appends each extra frame as a labelled isomer
    suffixed with the mode-tag.
    """
    if not _is_enabled():
        return [(xyz, "base")]
    coverage = _coverage_enabled()
    # Coverage mode raises the density defaults (bug-3): more pool members,
    # finer torsion grid, more DOFs/combinations, tighter dedup so genuine
    # backbone folds survive.  Legacy 5o defaults preserved when coverage off.
    k = _env_int(
        "DELFIN_FFFREE_COV_K_TARGET" if coverage else "DELFIN_5O_K_TARGET",
        24 if coverage else 10, lo=1, hi=128,
    )
    rmsd_min = _env_float(
        "DELFIN_FFFREE_COV_DIVERSITY_RMSD_MIN" if coverage
        else "DELFIN_5O_DIVERSITY_RMSD_MIN",
        0.30 if coverage else 0.5, lo=0.0, hi=5.0,
    )
    n_states = _env_int(
        "DELFIN_FFFREE_COV_TORSION_STATES" if coverage
        else "DELFIN_5O_K_TORSION_STATES",
        6 if coverage else 3, lo=2, hi=12,
    )
    max_dofs = _env_int(
        "DELFIN_FFFREE_COV_MAX_DOFS" if coverage else "DELFIN_5O_MAX_DOFS",
        8 if coverage else 6, lo=1, hi=32,
    )
    grid_cap = _env_int(
        "DELFIN_FFFREE_COV_TORSION_GRID_CAP" if coverage
        else "DELFIN_5O_TORSION_GRID_CAP",
        512 if coverage else 64, lo=2, hi=4096,
    )
    pucker_ampl = _env_float("DELFIN_5O_RING_PUCKER_AMPL", 0.35, lo=0.0, hi=2.0)
    macro_ampl = _env_float("DELFIN_5O_MACROCYCLE_AMPL", 0.40, lo=0.0, hi=2.0)
    twist_deg = _env_float("DELFIN_5O_CHELATE_TWIST_DEG", 30.0, lo=0.0, hi=180.0)
    md_tol = _env_float("DELFIN_5O_MD_TOL", 0.05, lo=0.0, hi=2.0)
    core_freeze_tol = _env_float(
        "DELFIN_FFFREE_COV_CORE_FREEZE_TOL", 0.05, lo=0.0, hi=2.0
    )
    # Coverage mode relaxes the non-bonded clash floors toward the physical
    # fusion limit: crystal backbone folds (esp. crowded biaryl / o-terphenyl
    # like CIYROT) legitimately adopt sterically TIGHT contacts that the
    # conservative 1.5/1.7 A selection floors reject — and COVERAGE keeps the
    # best-of-ensemble MIN, so admitting a tight-but-valid fold can only help,
    # never hurt (a downstream local-opt relaxes residual strain).  Floors stay
    # ABOVE the bond-fusion regime (heavy-heavy 1.35, X-H 1.35, H-H 1.45) so no
    # atoms-on-top-of-each-other frame is ever emitted.
    clash_hh = _env_float(
        "DELFIN_5O_CLASH_HH", 1.45 if coverage else 1.7, lo=0.5, hi=3.0
    )
    clash_xh = _env_float(
        "DELFIN_5O_CLASH_XH", 1.35 if coverage else 1.7, lo=0.5, hi=3.0
    )
    clash_xx = _env_float(
        "DELFIN_5O_CLASH_XX", 1.35 if coverage else 1.5, lo=0.5, hi=3.0
    )
    try:
        return generate_conformer_pool(
            xyz,
            k_target=k,
            diversity_rmsd_min=rmsd_min,
            n_torsion_states=n_states,
            max_dofs=max_dofs,
            torsion_grid_cap=grid_cap,
            ring_pucker_ampl=pucker_ampl,
            macrocycle_ampl=macro_ampl,
            chelate_twist_deg=twist_deg,
            md_tol=md_tol,
            clash_hh=clash_hh,
            clash_xh=clash_xh,
            clash_xx=clash_xx,
            coverage=coverage,
            core_freeze_tol=core_freeze_tol,
        )
    except Exception as exc:  # pragma: no cover - safety net
        logger.debug("Welle-5o conformer-pool failed: %s", exc)
        return [(xyz, "base")]


# ---------------------------------------------------------------------------
# Diagnostic helper — useful for tests / future Pólya-completeness audits
# ---------------------------------------------------------------------------


def count_dofs(xyz: str) -> Dict[str, int]:
    """Return per-layer DOF counts for the given XYZ.

    Pure analytic — no candidate generation, no scoring.  Useful for
    validation reports listing torsion / chelate / macrocycle DOFs.
    Returns
    ``{"torsion": int, "ring_pucker": int, "chelate_twist": int, "macrocycle": int}``.
    """
    out = {"torsion": 0, "ring_pucker": 0, "chelate_twist": 0, "macrocycle": 0}
    ob_mol = _rot._build_ob_mol_from_xyz(xyz)
    if ob_mol is None:
        return out
    graph = _rot._graph_from_ob(ob_mol)
    if not graph:
        return out
    out["torsion"] = len(_rot.identify_rotamer_dofs(graph, max_dofs=999))
    rings = _rings_from_graph(graph)
    is_metal = graph["is_metal"]
    bonds = graph.get("bonds", [])
    aromatic_bond_set = set()
    for (a, b, _o, arom, _r) in bonds:
        if arom:
            aromatic_bond_set.add((min(a, b), max(a, b)))
    for ring in rings:
        size = len(ring)
        if size < 4 or size > 7:
            if size >= 7 and size <= 30:
                # Layer 4 candidate
                n_arom = sum(
                    1
                    for k in range(size)
                    if (min(ring[k], ring[(k + 1) % size]),
                        max(ring[k], ring[(k + 1) % size])) in aromatic_bond_set
                )
                if n_arom < size:
                    out["macrocycle"] += 1
            continue
        metals_in_ring = [i for i in ring if is_metal[i]]
        if metals_in_ring:
            if len(metals_in_ring) == 1 and size <= 7:
                out["chelate_twist"] += 1
        else:
            # Layer 2 candidate gate
            is_aromatic = False
            for k in range(size):
                a = ring[k]
                b = ring[(k + 1) % size]
                if (min(a, b), max(a, b)) in aromatic_bond_set:
                    is_aromatic = True
                    break
            if not is_aromatic:
                out["ring_pucker"] += 1
    return out


__all__ = [
    "generate_conformer_pool",
    "apply_if_enabled",
    "count_dofs",
]
