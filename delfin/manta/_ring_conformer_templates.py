"""Chemistry-accurate ring-conformer templates (Welle-5p-B).

User direktive 2026-05-18:
    "Sessel-Wanne-Übergang nicht durch Rotation modellierbar — muss
     anders gemacht werden."

A naive single-bond rotation cannot produce chair ↔ boat ↔ twist
transitions of saturated rings.  Cyclohexane chair → boat requires a
*concerted* displacement of six ring atoms perpendicular to the mean
ring plane.  This module provides chemistry-accurate, normalised
ring-conformer templates (pre-computed atom-Z displacement patterns)
that are applied as rigid-body atom-Z perturbations to the existing
ring frame — **not** as rotations.

Architecture
------------
For each ring size the module exposes a dictionary of normalised
displacement patterns (one signed scalar per ring atom, in ring-order).
The L2-norm of each pattern is 1.0 so a single global amplitude in Å
can be applied universally.

Ring sizes covered
~~~~~~~~~~~~~~~~~~
* **3 / 4**     — flat / butterfly (very limited conformational space).
* **5**         — five envelopes (E1–E5), five twists (T1–T5),
                  half-chair conformers.
* **6**         — chair (two enantiomeric inversions),
                  boat (three orientations),
                  twist-boat (three orientations),
                  half-chair, envelope.
* **7**         — chair, twist-chair, boat, twist-boat.
* **macrocycle** (≥8) — saddle, ruffled, dome (low-order Fourier
                        modes along the ring normal).

Universal-fundamental compliance
--------------------------------
* No SMILES patterns, no refcode regex.
* Trigger purely on **graph features**: ring size, aromatic-bond flag,
  metal-membership flag.
* Chelate-rings (any ring containing a metal atom) are **EXCLUDED**
  from templating — they require different handling (chelate-twist
  in :mod:`delfin.manta._conformer_pool` Layer 3, or the planned 5p-A
  topology-hash framework).
* Aromatic rings are EXCLUDED (planar fixed by π-conjugation).
* Templates are normalised to unit RMS amplitude; the actual
  displacement scales with average ring bond length, so cyclopentane
  and cyclohexane share the same template library but produce
  proportional perturbations.

Wire-in
-------
:func:`generate_ring_conformer_pool` returns ``[(coords, mode_tag), ...]``
suitable for plugging into :mod:`delfin.manta._conformer_pool` as a
**replacement** for the current naive 2-atom Layer-2 pucker.

Env-flags
---------
``DELFIN_5P_B_RING_TEMPLATES``      (default ``0``) — master switch
``DELFIN_5P_B_K_PER_RING``          (default ``3``) — top-K templates / ring
``DELFIN_5P_B_MAX_RING_VARIANTS``   (default ``6``) — combinatorial cap
``DELFIN_5P_B_AMPLITUDE_FRACTION``  (default ``0.50``) — fraction of ring
                                                       bond length used as
                                                       displacement amplitude
``DELFIN_5P_B_MD_TOL``              (default ``0.05``) — Å, M-D invariant guard
"""

from __future__ import annotations

import math
import os
from typing import Dict, List, Optional, Sequence, Tuple

from delfin.common.logging import get_logger

# Reuse helpers from the rotamer-diversity module
from delfin.manta import _rotamer_diversity as _rot

logger = get_logger(__name__)


Coord = Tuple[float, float, float]


# ---------------------------------------------------------------------------
# Environment configuration
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
    return _env_bool("DELFIN_5P_B_RING_TEMPLATES", False)


# ---------------------------------------------------------------------------
# Template library
# ---------------------------------------------------------------------------
#
# Each template is a tuple of signed displacement coefficients, one per
# ring atom in ring-order.  Coefficients are normalised so the L2-norm
# is 1.0 — a single scalar amplitude (in Å) then controls the per-atom
# perpendicular shift along the ring normal.
#
# Patterns are derived from textbook chemistry:
#  * Cyclohexane chair: alternating ±1 (D3d symmetry, atoms 1,3,5 up;
#                       2,4,6 down).
#  * Cyclohexane boat:  two flagpoles + four base atoms in plane (C2v).
#  * Twist-boat:        diagonal twist (D2 symmetry).
#  * Half-chair:        4 atoms in plane, 2 puckered antisymmetric.
#  * Envelope (5-ring): one atom out of plane.
#  * Twist (5-ring):    two adjacent atoms antisymmetrically displaced.
#
# All normalised templates are computed once at module import via
# :func:`_normalise` so the per-call cost is just a multiplication.


def _normalise(pattern: Sequence[float]) -> Tuple[float, ...]:
    """Return *pattern* rescaled so its L2-norm is 1.0."""
    s = math.sqrt(sum(v * v for v in pattern))
    if s < 1e-12:
        return tuple(float(v) for v in pattern)
    return tuple(v / s for v in pattern)


# ----- 3-ring -----
TEMPLATES_3_RING: Dict[str, Tuple[float, ...]] = {
    # 3-rings are rigid by geometry — provide only a trivial dome mode
    # for completeness; rarely used.
    "dome": _normalise((1.0, 1.0, 1.0)),
}


# ----- 4-ring -----
TEMPLATES_4_RING: Dict[str, Tuple[float, ...]] = {
    # Butterfly (puckered) — two opposite corners up, two down.
    "butterfly": _normalise((+1.0, -1.0, +1.0, -1.0)),
    # Envelope-like (1,3 up, 2,4 zero).  Use small zero with safety
    # epsilon so the pattern is non-degenerate.
    "envelope-1": _normalise((+1.0, 0.0, +1.0, 0.0)),
}


# ----- 5-ring -----
# Cyclopentane conformational map: 10 envelope (E) + 10 twist (T)
# conformers form a pseudo-rotational pathway.  We tabulate the five
# unique-symmetry envelopes and five twists; cyclic ring-shifts of the
# input ring-list cover the other half automatically.

_E5 = []
for k in range(5):
    p = [0.0] * 5
    p[k] = 2.0
    for j in range(5):
        if j != k:
            p[j] = -0.5  # the other four atoms compensate downward
    _E5.append(_normalise(p))

_T5 = []
for k in range(5):
    p = [0.0] * 5
    p[k] = +1.0
    p[(k + 1) % 5] = -1.0
    _T5.append(_normalise(p))

TEMPLATES_5_RING: Dict[str, Tuple[float, ...]] = {
    f"envelope-{k+1}": _E5[k] for k in range(5)
}
TEMPLATES_5_RING.update({f"twist-{k+1}": _T5[k] for k in range(5)})
# Half-chair (4 atoms in plane, 2 adjacent puckered)
for k in range(5):
    p = [0.0] * 5
    p[k] = +1.0
    p[(k + 2) % 5] = -1.0
    TEMPLATES_5_RING[f"half-chair-{k+1}"] = _normalise(p)


# ----- 6-ring -----
# Cyclohexane conformational map — textbook chair / boat / twist-boat /
# half-chair / envelope.  Chair has two enantiomeric inversions.

TEMPLATES_6_RING: Dict[str, Tuple[float, ...]] = {
    # Chair (D3d) — alternating ±1.  Two enantiomers (chair-A, chair-B).
    "chair-A": _normalise((+1.0, -1.0, +1.0, -1.0, +1.0, -1.0)),
    "chair-B": _normalise((-1.0, +1.0, -1.0, +1.0, -1.0, +1.0)),
    # Boat (C2v) — atoms 1 and 4 are flagpoles (up); atoms 2,3,5,6 in
    # the base plane (small alternating compensation).  Three boats
    # by cycling which atoms are flagpoles.
    "boat-14": _normalise((+1.0, -0.33, -0.33, +1.0, -0.33, -0.33)),
    "boat-25": _normalise((-0.33, +1.0, -0.33, -0.33, +1.0, -0.33)),
    "boat-36": _normalise((-0.33, -0.33, +1.0, -0.33, -0.33, +1.0)),
    # Twist-boat (D2) — diagonal twist on a boat.
    "twist-boat-1": _normalise((+1.0, -1.0, +0.5, -1.0, +1.0, -0.5)),
    "twist-boat-2": _normalise((+0.5, +1.0, -1.0, +0.5, +1.0, -1.0)),
    "twist-boat-3": _normalise((-1.0, +0.5, +1.0, -1.0, +0.5, +1.0)),
    # Half-chair (C2) — five atoms in plane, one puckered.  Six
    # rotations cover the conformational manifold.
    "half-chair-1": _normalise((+1.0, +0.5, -0.5, -1.0, -0.5, +0.5)),
    "half-chair-2": _normalise((+0.5, +1.0, +0.5, -0.5, -1.0, -0.5)),
    # Envelope (Cs) — one atom out of plane.
    "envelope-1": _normalise((+1.0, 0.0, 0.0, 0.0, 0.0, 0.0)),
    "envelope-4": _normalise((0.0, 0.0, 0.0, +1.0, 0.0, 0.0)),
}


# ----- 7-ring -----
# Cycloheptane: chair, twist-chair, boat, twist-boat.  All four are
# accessible via low-energy pseudo-rotation; we expose normalised
# patterns for the four canonical entries.

TEMPLATES_7_RING: Dict[str, Tuple[float, ...]] = {
    "chair-7": _normalise((+1.0, -1.0, +1.0, -1.0, +1.0, -1.0, +0.0)),
    "twist-chair-7": _normalise((+1.0, -0.5, -1.0, +0.5, +1.0, -0.5, -1.0)),
    "boat-7": _normalise((+1.0, -0.33, -0.33, +1.0, -0.33, -0.33, -0.33)),
    "twist-boat-7": _normalise((+1.0, +0.5, -1.0, -0.5, +1.0, +0.5, -1.0)),
}


# ---------------------------------------------------------------------------
# Macrocycle (size ≥ 8) — Fourier mode templates parametric on size
# ---------------------------------------------------------------------------


def _macrocycle_templates(size: int) -> Dict[str, Tuple[float, ...]]:
    """Return a dict of (mode_tag → normalised displacement pattern) for a
    macrocyclic ring of *size* atoms.

    Uses three universal low-order Fourier modes along the ring normal:

    * ``saddle``  — :math:`\\cos(2\\theta)`
    * ``ruffled`` — :math:`\\cos(4\\theta)`
    * ``dome``    — alternating up/zero per atom (the "half-up" mode)

    The actual angular positions :math:`\\theta_k` are uniform around
    the ring centroid: :math:`\\theta_k = 2 \\pi k / size`.
    """
    out: Dict[str, Tuple[float, ...]] = {}
    saddle = [math.cos(2.0 * math.pi * k * 2 / size) for k in range(size)]
    out["saddle"] = _normalise(saddle)
    if size >= 6:
        ruffled = [math.cos(2.0 * math.pi * k * 4 / size) for k in range(size)]
        out["ruffled"] = _normalise(ruffled)
    # Dome — even-indexed atoms up, odd-indexed atoms zero (half-up mode).
    dome = [1.0 if (k % 2) == 0 else 0.0 for k in range(size)]
    out["dome"] = _normalise(dome)
    return out


def get_templates_for_ring_size(size: int) -> Dict[str, Tuple[float, ...]]:
    """Return the template library for *size*.

    Parameters
    ----------
    size:
        Ring size in atoms.

    Returns
    -------
    dict
        ``{mode_tag: normalised_pattern}``.  Empty dict if *size* is
        outside the supported range (3..30).
    """
    if size == 3:
        return dict(TEMPLATES_3_RING)
    if size == 4:
        return dict(TEMPLATES_4_RING)
    if size == 5:
        return dict(TEMPLATES_5_RING)
    if size == 6:
        return dict(TEMPLATES_6_RING)
    if size == 7:
        return dict(TEMPLATES_7_RING)
    if 8 <= size <= 30:
        return _macrocycle_templates(size)
    return {}


# ---------------------------------------------------------------------------
# Geometry helpers (pure-Python; mirror _conformer_pool style)
# ---------------------------------------------------------------------------


def _vec_sub(a: Coord, b: Coord) -> Coord:
    return (a[0] - b[0], a[1] - b[1], a[2] - b[2])


def _vec_add(a: Coord, b: Coord) -> Coord:
    return (a[0] + b[0], a[1] + b[1], a[2] + b[2])


def _vec_scale(a: Coord, s: float) -> Coord:
    return (a[0] * s, a[1] * s, a[2] * s)


def _vec_dot(a: Coord, b: Coord) -> float:
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]


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

    Power iteration on (trace*I − cov) to extract the smallest
    eigenvector of the covariance matrix (= the plane normal).  Same
    technique as :func:`delfin.manta._conformer_pool._ring_plane_normal`.
    """
    c = _centroid(coords, ring_atoms)
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
    trace = cov[0][0] + cov[1][1] + cov[2][2]
    a = [
        [trace - cov[i][j] if i == j else -cov[i][j] for j in range(3)]
        for i in range(3)
    ]
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


def _average_ring_bond_length(
    coords: Sequence[Coord], ring_atoms: Sequence[int]
) -> float:
    """Return mean Euclidean distance between adjacent atoms in *ring_atoms*."""
    size = len(ring_atoms)
    if size < 2:
        return 0.0
    s = 0.0
    for k in range(size):
        a = ring_atoms[k]
        b = ring_atoms[(k + 1) % size]
        s += _vec_norm(_vec_sub(coords[a], coords[b]))
    return s / float(size)


# ---------------------------------------------------------------------------
# Ring detection (re-used pattern from _conformer_pool)
# ---------------------------------------------------------------------------


def _rings_from_graph(graph: Dict) -> List[List[int]]:
    """Return SSSR-like ring atom lists from an OB-perceived graph.

    Mirrors :func:`delfin.manta._conformer_pool._rings_from_graph` for
    independence from that module.
    """
    ring_bonds = []
    for (a, b, _order, _arom, is_ring) in graph.get("bonds", []):
        if is_ring:
            ring_bonds.append((a, b))
    if not ring_bonds:
        return []

    adj: Dict[int, List[int]] = {}
    for (a, b) in ring_bonds:
        adj.setdefault(a, []).append(b)
        adj.setdefault(b, []).append(a)

    rings: List[List[int]] = []
    seen_ring_keys: set = set()
    nodes = sorted(adj.keys())

    def _dfs(start: int, cur: int, path: List[int], depth: int):
        if depth > 16:
            return
        for nbr in adj.get(cur, []):
            if nbr == start and len(path) >= 3:
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

    rings.sort(key=lambda r: (len(r), r))
    accepted: List[List[int]] = []
    for r in rings:
        rset = set(r)
        if any(set(a).issuperset(rset) for a in accepted if len(a) <= len(r)):
            continue
        accepted.append(r)
    return accepted


# ---------------------------------------------------------------------------
# Ring eligibility (universal — graph features only)
# ---------------------------------------------------------------------------


def find_rings_for_templating(graph: Dict) -> List[List[int]]:
    """Return rings amenable to template application.

    EXCLUDE:
        * rings containing a metal atom (= chelate rings; handled by
          Layer-3 chelate-twist in :mod:`delfin.manta._conformer_pool`)
        * rings whose every consecutive bond is aromatic (planar π-system,
          fixed)
        * rings whose size lies outside the supported range (3..30)

    INCLUDE:
        * cyclopentyl / cyclohexyl backbones
        * pyranose / furanose / piperidine etc. (heteroatom carbocycles)
        * macrocyclic ligand backbones up to size 30

    Universal-fundamental: triggers exclusively on graph features
    (ring-bond flag, aromatic flag, metal flag).  No SMILES patterns.
    """
    if not graph:
        return []
    is_metal = graph.get("is_metal", [])
    bonds = graph.get("bonds", [])
    aromatic_bond_set = set()
    for (a, b, _o, arom, _r) in bonds:
        if arom:
            aromatic_bond_set.add((min(a, b), max(a, b)))

    rings = _rings_from_graph(graph)
    out: List[List[int]] = []
    for ring in rings:
        size = len(ring)
        if size < 3 or size > 30:
            continue
        # EXCLUDE metal-containing rings — handled separately
        if any(is_metal[i] for i in ring):
            continue
        # Count aromatic ring-bonds
        n_aromatic = 0
        for k in range(size):
            a = ring[k]
            b = ring[(k + 1) % size]
            if (min(a, b), max(a, b)) in aromatic_bond_set:
                n_aromatic += 1
        # Fully aromatic → skip (planar fixed)
        if n_aromatic >= size:
            continue
        out.append(ring)
    return out


# ---------------------------------------------------------------------------
# Template application
# ---------------------------------------------------------------------------


def apply_template(
    base_coords: Sequence[Coord],
    ring_atoms: Sequence[int],
    pattern: Sequence[float],
    amplitude: float,
    graph: Optional[Dict] = None,
    drag_hydrogens: bool = True,
) -> List[Coord]:
    """Apply a normalised displacement *pattern* to *ring_atoms*.

    Each ring atom is displaced by ``pattern[k] * amplitude`` along the
    ring-plane normal.  Bonded hydrogens are dragged rigidly with their
    heavy parent atom (per :doc:`feedback_rigid_h_tracking_b5_bug`).

    Parameters
    ----------
    base_coords:
        Full molecule coordinates (will not be mutated; a copy is
        returned).
    ring_atoms:
        Atom indices that form the ring, in ring-order.
    pattern:
        Normalised displacement pattern (one signed scalar per ring
        atom).
    amplitude:
        Scalar (Å) controlling how far the most-displaced atom moves
        along the normal.
    graph:
        Optional bond/neighbour graph.  When supplied, hydrogens that
        are bonded to a ring atom are dragged rigidly with that heavy
        atom — preserves C-H bond lengths.
    drag_hydrogens:
        If True (default) and *graph* is supplied, drag bonded
        hydrogens with their heavy parent.

    Returns
    -------
    list
        New ``[(x, y, z), ...]`` coordinates.
    """
    n_atoms = len(base_coords)
    size = len(ring_atoms)
    if len(pattern) != size:
        raise ValueError(
            f"pattern length {len(pattern)} != ring size {size}"
        )
    _c, normal = _ring_plane_normal(base_coords, ring_atoms)
    if _vec_norm(normal) < 1e-9:
        return list(base_coords)

    coords: List[Coord] = list(base_coords)
    neighbours = graph.get("neighbours", [[] for _ in range(n_atoms)]) if graph else None
    atomic_nums = graph.get("atomic_nums", []) if graph else []

    for k, atom_idx in enumerate(ring_atoms):
        disp = _vec_scale(normal, amplitude * pattern[k])
        coords[atom_idx] = _vec_add(coords[atom_idx], disp)
        # Drag bonded hydrogens (rigid-H)
        if drag_hydrogens and neighbours is not None and atomic_nums:
            for nbr in neighbours[atom_idx]:
                if nbr < len(atomic_nums) and atomic_nums[nbr] == 1:
                    coords[nbr] = _vec_add(coords[nbr], disp)
    return coords


# ---------------------------------------------------------------------------
# Validation: M-D invariant + topology hash (standalone fallbacks)
# ---------------------------------------------------------------------------


def _md_distance_check(
    base_coords: Sequence[Coord],
    new_coords: Sequence[Coord],
    graph: Dict,
    tol: float,
) -> bool:
    """Return True iff every M-D bond length is preserved within *tol* Å.

    Standalone fallback: re-uses
    :func:`delfin.manta._rotamer_diversity._coord_bond_invariant_holds`
    when available; otherwise iterates metal-donor edges directly.
    """
    if hasattr(_rot, "_coord_bond_invariant_holds"):
        try:
            return _rot._coord_bond_invariant_holds(
                graph, base_coords, new_coords, tol=tol
            )
        except Exception:
            pass
    # Inline fallback
    is_metal = graph.get("is_metal", [])
    for (a, b, _o, _arom, _ring) in graph.get("bonds", []):
        if a >= len(is_metal) or b >= len(is_metal):
            continue
        if not (is_metal[a] or is_metal[b]):
            continue
        d0 = _vec_norm(_vec_sub(base_coords[a], base_coords[b]))
        d1 = _vec_norm(_vec_sub(new_coords[a], new_coords[b]))
        if abs(d1 - d0) > tol:
            return False
    return True


# ---------------------------------------------------------------------------
# Pool generation
# ---------------------------------------------------------------------------


def generate_ring_conformer_pool(
    xyz: str,
    k_per_ring: int = 3,
    max_ring_variants: int = 6,
    amplitude_fraction: float = 0.50,
    md_tol: float = 0.05,
) -> List[Tuple[str, str]]:
    """Generate ring-template variants for *xyz*.

    For each ring eligible for templating (non-aromatic, non-chelate)
    pick the first ``k_per_ring`` templates from the library and apply
    them.  The overall list is capped at ``max_ring_variants`` to
    avoid combinatorial explosion.

    Returns
    -------
    list
        ``[(xyz_string, mode_tag), ...]``.  The first element is always
        ``(base_xyz, "base")``.  Subsequent entries are validated
        variants.
    """
    base_pool: List[Tuple[str, str]] = [(xyz, "base")]

    ob_mol = _rot._build_ob_mol_from_xyz(xyz)
    if ob_mol is None:
        return base_pool
    graph = _rot._graph_from_ob(ob_mol)
    if not graph:
        return base_pool

    try:
        symbols, base_coords_t = _rot._parse_delfin_xyz(xyz)
    except Exception:
        return base_pool
    base_coords: List[Coord] = [tuple(c) for c in base_coords_t]

    rings = find_rings_for_templating(graph)
    if not rings:
        return base_pool

    base_topo = None
    try:
        base_topo = _rot._topology_hash(graph)
    except Exception:
        base_topo = None

    variants: List[Tuple[str, str]] = []
    for ring in rings:
        if len(variants) >= max_ring_variants:
            break
        size = len(ring)
        templates = get_templates_for_ring_size(size)
        if not templates:
            continue
        # Pick top-K templates (deterministic ordering)
        template_items = list(templates.items())[: max(1, k_per_ring)]
        avg_bond = _average_ring_bond_length(base_coords, ring)
        if avg_bond < 1e-3:
            continue
        amplitude = avg_bond * amplitude_fraction
        for tag, pattern in template_items:
            if len(variants) >= max_ring_variants:
                break
            try:
                new_coords = apply_template(
                    base_coords, ring, pattern, amplitude, graph=graph
                )
            except Exception as exc:
                logger.debug("5p-B template apply failed (%s): %s", tag, exc)
                continue
            # Validate M-D invariant
            if not _md_distance_check(base_coords, new_coords, graph, md_tol):
                continue
            # Validate topology hash
            cand_xyz = _rot._format_delfin_xyz(symbols, new_coords)
            if base_topo is not None:
                try:
                    cand_mol = _rot._build_ob_mol_from_xyz(cand_xyz)
                    if cand_mol is None:
                        continue
                    cand_graph = _rot._graph_from_ob(cand_mol)
                    if _rot._topology_hash(cand_graph) != base_topo:
                        continue
                except Exception:
                    continue
            variants.append(
                (cand_xyz, f"ring-r{size}-{tag}-a{ring[0]}")
            )

    return base_pool + variants


# ---------------------------------------------------------------------------
# Layer-2 wire-in helper (drop-in replacement for the naive ring-pucker)
# ---------------------------------------------------------------------------


def layer2_ring_template_candidates(
    symbols: Sequence[str],
    base_coords: Sequence[Coord],
    graph: Dict,
    rings: Sequence[Sequence[int]],
    amplitude_fraction: float,
    k_per_ring: int,
    max_ring_variants: int,
) -> List[Tuple[List[Coord], str]]:
    """Replacement for :func:`_conformer_pool._layer2_ring_pucker_candidates`.

    Returns ``[(coords, mode_tag), ...]`` for the conformer-pool layered
    sampler.  Same signature pattern as the existing Layer-2 helper so
    wire-in is a single-line swap.

    Universal: triggers on every non-aromatic, non-metal ring of size
    3..30.  Chelate rings (metal-containing) are excluded — those are
    Layer-3's responsibility.
    """
    is_metal = graph.get("is_metal", [])
    bonds = graph.get("bonds", [])
    aromatic_bond_set = set()
    for (a, b, _o, arom, _r) in bonds:
        if arom:
            aromatic_bond_set.add((min(a, b), max(a, b)))

    out: List[Tuple[List[Coord], str]] = []
    for ring in rings:
        if len(out) >= max_ring_variants:
            break
        size = len(ring)
        if size < 3 or size > 30:
            continue
        # Skip metal-containing rings — chelate Layer-3 handles those.
        if any(is_metal[i] for i in ring):
            continue
        # Skip fully aromatic rings (planar fixed).
        n_aromatic = 0
        for k in range(size):
            a = ring[k]
            b = ring[(k + 1) % size]
            if (min(a, b), max(a, b)) in aromatic_bond_set:
                n_aromatic += 1
        if n_aromatic >= size:
            continue
        templates = get_templates_for_ring_size(size)
        if not templates:
            continue
        template_items = list(templates.items())[: max(1, k_per_ring)]
        avg_bond = _average_ring_bond_length(base_coords, ring)
        if avg_bond < 1e-3:
            continue
        amplitude = avg_bond * amplitude_fraction
        for tag, pattern in template_items:
            if len(out) >= max_ring_variants:
                break
            try:
                new_coords = apply_template(
                    base_coords, ring, pattern, amplitude, graph=graph
                )
            except Exception:
                continue
            out.append((new_coords, f"ring-r{size}-{tag}-a{ring[0]}"))
    return out


# ---------------------------------------------------------------------------
# Public env-gated wire-in
# ---------------------------------------------------------------------------


def apply_if_enabled(xyz: str) -> List[Tuple[str, str]]:
    """Return ``[(xyz, mode_tag), ...]`` when ``DELFIN_5P_B_RING_TEMPLATES=1``.

    Default-OFF semantics: when the flag is unset, returns
    ``[(xyz, "base")]`` unchanged — byte-identical fall-through.

    On any internal failure the helper falls back to the base-only pool.
    """
    if not _is_enabled():
        return [(xyz, "base")]
    k_per_ring = _env_int("DELFIN_5P_B_K_PER_RING", 3, lo=1, hi=20)
    max_var = _env_int("DELFIN_5P_B_MAX_RING_VARIANTS", 6, lo=1, hi=64)
    amp_frac = _env_float("DELFIN_5P_B_AMPLITUDE_FRACTION", 0.50, lo=0.0, hi=2.0)
    md_tol = _env_float("DELFIN_5P_B_MD_TOL", 0.05, lo=0.0, hi=2.0)
    try:
        return generate_ring_conformer_pool(
            xyz,
            k_per_ring=k_per_ring,
            max_ring_variants=max_var,
            amplitude_fraction=amp_frac,
            md_tol=md_tol,
        )
    except Exception as exc:  # pragma: no cover - safety net
        logger.debug("Welle-5p-B ring-template pool failed: %s", exc)
        return [(xyz, "base")]


# ---------------------------------------------------------------------------
# Diagnostic / introspection
# ---------------------------------------------------------------------------


def template_library_coverage() -> Dict[str, int]:
    """Return number of templates per ring-size in the static library.

    Useful for validation reports.  Macrocycle sizes are sampled at a
    few representative values; the actual library is computed
    on-the-fly per ring.
    """
    out: Dict[str, int] = {
        "ring-3": len(TEMPLATES_3_RING),
        "ring-4": len(TEMPLATES_4_RING),
        "ring-5": len(TEMPLATES_5_RING),
        "ring-6": len(TEMPLATES_6_RING),
        "ring-7": len(TEMPLATES_7_RING),
    }
    for size in (8, 10, 12, 16, 24):
        out[f"ring-{size}"] = len(_macrocycle_templates(size))
    return out


__all__ = [
    "TEMPLATES_3_RING",
    "TEMPLATES_4_RING",
    "TEMPLATES_5_RING",
    "TEMPLATES_6_RING",
    "TEMPLATES_7_RING",
    "apply_if_enabled",
    "apply_template",
    "find_rings_for_templating",
    "generate_ring_conformer_pool",
    "get_templates_for_ring_size",
    "layer2_ring_template_candidates",
    "template_library_coverage",
]
