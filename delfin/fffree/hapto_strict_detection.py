"""delfin.fffree.hapto_strict_detection — Strict hapto-π eligibility filter.

Fundamental architecture extension (2026-06-07): the existing hapto-π
detection (:func:`delfin.fffree.coordination_mode_enum._detect_hapto_modes`)
walks the molecular graph with SMARTS patterns matching any aromatic / open-π
fragment and emits η²/η⁴/η⁶ candidates for each match.  This is correct
chemistry for STAND-ALONE arene/Cp/COT/butadiene SMILES (the "what hapto
modes can this ligand support in principle?" question) but produces
**false-positive** hapto orbits when the SMILES describes a full
organometallic complex in which most aromatic rings are σ-bond
substituents (phenyl groups on a chelate backbone, pyridyl pendant to a
σ-N donor, …) that the metal is NOT hapto-bonded to.

Example (YUHRUP-class, V3 voll-pool 2026-06-06)
-----------------------------------------------
SMILES (Re σ-S/σ-N chelate with phenyl substituent on the chelate)::

    CCN(CC)C1=NC(C2=CC=CC=C2)=[N+]2N=C(N3CCCCCC3)[S][Re]2([O-])([Cl])[S]1

Re is σ-bound (S, N₂, O⁻, Cl, S — CN5).  The phenyl ring ``c1ccccc1`` is
a pure substituent on the chelate backbone, NOT a hapto donor.  The legacy
``_detect_hapto_modes`` nevertheless emits ``eta2-arene``, ``eta4-arene``
and ``eta6-arene`` for that phenyl because the SMARTS ``[c]1[c][c][c][c][c]1``
matches it.  The orbit set then contains binding-mode mixtures (κ⁵-σ ⊕ η⁶-arene)
that have no geometric realisation and waste enumeration budget downstream.

Universal graph topology criterion
----------------------------------
A ring (or open-π chain) is a **legitimate hapto candidate** if and only if
at least one of its atoms is within ``D_max`` bonds (in the metal-removed
graph) of an atom that IS directly bonded to a metal in the original
graph.  ``D_max`` defaults to 1, which is the textbook criterion (one of
the ring atoms must itself be the metal's neighbour to be hapto-bonded);
the bound is exposed via the env-flag ``DELFIN_FFFREE_HAPTO_STRICT_MAX_DIST``
so μ-η-bridging arenes (rare; ring sits between two metals) can be admitted
with ``D_max = 2`` if required.

Why universal & graph-only
--------------------------
* No SMILES patterns — only graph BFS on the molecular adjacency.
* No element heuristics — the metal set is the standard
  :func:`delfin._bond_decollapse._is_metal` predicate.
* No SMILES-specific shortcuts — the rule is the same whether the ring is
  a phenyl, pyridyl, naphthyl, thiophene, … any ring whose atoms are too
  far from a metal in the graph cannot be hapto.
* Mathematically: hapto-bond = π-electron donation requires the ring's
  π-cloud to overlap with the metal's empty d-orbital.  Overlap requires
  proximity in real space, which (modulo conformer) is bounded by graph
  distance.  Distance > 1 bond from a metal-coordinated atom means the
  ring is in the second coordination sphere or further — never hapto.

Determinism & env-gate
----------------------
* Pure-Python; sorted containers; no random sampling.
* Default OFF (byte-identical) gated by
  ``DELFIN_FFFREE_HAPTO_STRICT_DETECTION=1``.
* When OFF, callers MUST treat the helper as identity (always returning
  True for legitimacy).  This module's public entry-point
  :func:`is_legitimate_hapto_candidate` already does so via its
  ``strict`` argument (default = False).
"""
from __future__ import annotations

import os
from typing import Iterable, List, Optional, Sequence, Set


# ---------------------------------------------------------------------------
# Env flags
# ---------------------------------------------------------------------------
def strict_detection_enabled() -> bool:
    """Live-evaluated env-gate for the strict hapto-eligibility filter.

    Default OFF (byte-identical).  Returns True when
    ``DELFIN_FFFREE_HAPTO_STRICT_DETECTION`` is set to ``1``/``true``/
    ``yes``/``on`` (case-insensitive).
    """
    return os.environ.get(
        "DELFIN_FFFREE_HAPTO_STRICT_DETECTION", "0"
    ).strip().lower() in ("1", "true", "yes", "on")


def _max_dist() -> int:
    """Max graph distance (in bonds, metal-removed graph) from a
    metal-coordinated atom to a ring atom for the ring to be considered
    a hapto candidate.

    Default = 1 (textbook hapto: at least one ring atom IS a metal neighbour).
    Env-overridable via ``DELFIN_FFFREE_HAPTO_STRICT_MAX_DIST``.  A value of
    2 admits μ-η-bridging arenes; 0 would require every ring atom to be a
    direct neighbour (too strict for asymmetric η-modes).
    """
    raw = os.environ.get("DELFIN_FFFREE_HAPTO_STRICT_MAX_DIST", "").strip()
    if raw:
        try:
            v = int(raw)
            if v >= 0:
                return v
        except (TypeError, ValueError):
            pass
    return 1


# ---------------------------------------------------------------------------
# Graph helpers (RDKit-mol agnostic — accept either an RDKit mol or a
# (n_atoms, bonds, metal_idxs) triple so callers don't have to materialise
# the full molecule object).
# ---------------------------------------------------------------------------
def _is_metal_safe(symbol: str) -> bool:
    try:
        from delfin._bond_decollapse import _is_metal as _im
    except Exception:
        return False
    try:
        return bool(_im(symbol))
    except Exception:
        return False


def _metal_atoms_from_mol(mol) -> List[int]:
    """Return sorted list of RDKit atom indices that are metals."""
    out: List[int] = []
    try:
        for a in mol.GetAtoms():
            try:
                if _is_metal_safe(a.GetSymbol()):
                    out.append(int(a.GetIdx()))
            except Exception:
                continue
    except Exception:
        return []
    return sorted(out)


def _metal_coord_atoms(mol, metal_idxs: Sequence[int]) -> Set[int]:
    """Return the set of (non-metal) atom indices that are directly bonded
    to ANY metal in ``metal_idxs`` (the "first coordination sphere").

    Pure graph; deterministic.  Returns an empty set if mol has no metals
    or no metal-X bonds.
    """
    out: Set[int] = set()
    if not metal_idxs:
        return out
    metal_set = {int(m) for m in metal_idxs}
    try:
        for m in metal_set:
            try:
                atom = mol.GetAtomWithIdx(int(m))
            except Exception:
                continue
            for nbr in atom.GetNeighbors():
                try:
                    ni = int(nbr.GetIdx())
                except Exception:
                    continue
                if ni in metal_set:
                    continue           # M-M bond, not a donor
                out.add(ni)
    except Exception:
        return out
    return out


def _bfs_distance(
    n_atoms: int,
    adj: Sequence[Sequence[int]],
    sources: Iterable[int],
    max_dist: int,
) -> List[int]:
    """Multi-source BFS up to ``max_dist`` bonds.

    Returns a list ``dist`` of length ``n_atoms`` where ``dist[i]`` is the
    minimum bond-graph distance from any source to atom ``i``, or a large
    sentinel (= n_atoms + 1) if unreached within ``max_dist``.

    The graph is ``adj`` (already metal-removed by the caller, so paths
    through the metal are excluded).  Pure Python, deterministic.
    """
    INF = n_atoms + 1
    dist = [INF] * n_atoms
    frontier: List[int] = []
    for s in sources:
        s = int(s)
        if 0 <= s < n_atoms and dist[s] > 0:
            dist[s] = 0
            frontier.append(s)
    d = 0
    while frontier and d < max_dist:
        new_frontier: List[int] = []
        d += 1
        for u in frontier:
            for v in adj[u]:
                if dist[v] > d:
                    dist[v] = d
                    new_frontier.append(v)
        frontier = new_frontier
    return dist


def _build_metal_removed_adj(mol, metal_idxs: Sequence[int]) -> List[List[int]]:
    """Return a deterministic adjacency list (sorted neighbour lists) for
    ``mol`` with all metal atoms removed (their bonds skipped).
    """
    n = int(mol.GetNumAtoms())
    metal_set = {int(m) for m in metal_idxs}
    adj: List[List[int]] = [[] for _ in range(n)]
    for bond in mol.GetBonds():
        try:
            a = int(bond.GetBeginAtomIdx())
            b = int(bond.GetEndAtomIdx())
        except Exception:
            continue
        if a in metal_set or b in metal_set:
            continue
        adj[a].append(b)
        adj[b].append(a)
    for i in range(n):
        adj[i] = sorted(set(adj[i]))
    return adj


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------
def is_legitimate_hapto_candidate(
    mol,
    ring_atoms: Sequence[int],
    *,
    metal_idxs: Optional[Sequence[int]] = None,
    max_dist: Optional[int] = None,
    strict: Optional[bool] = None,
) -> bool:
    """Return True iff ``ring_atoms`` is a plausible hapto donor for ``mol``.

    Strict graph-topology criterion (when ``strict=True``):

    1. At least one ring atom is directly bonded to a metal in ``mol``
       (the ring is in the first coordination sphere = textbook hapto), OR
    2. At least one ring atom is within ``max_dist`` bonds of an atom that
       is directly bonded to a metal, in the metal-removed graph (admits
       μ-η-bridging arenes when ``max_dist`` ≥ 2).

    Parameters
    ----------
    mol : RDKit Mol
        The molecular graph to test against.  Must support
        :meth:`GetAtoms`, :meth:`GetBonds`, :meth:`GetAtomWithIdx`,
        :meth:`GetNumAtoms`.
    ring_atoms : sequence of int
        RDKit atom indices that constitute the candidate hapto ring /
        open-π chain.
    metal_idxs : sequence of int, optional
        Pre-computed metal atom indices.  When omitted, computed via
        :func:`delfin._bond_decollapse._is_metal`.
    max_dist : int, optional
        Override the default max graph distance (defaults to the env flag
        ``DELFIN_FFFREE_HAPTO_STRICT_MAX_DIST`` or 1).
    strict : bool, optional
        When ``False`` (default), the helper returns True unconditionally
        (byte-identical legacy behaviour).  When ``True``, applies the
        topology filter.  When ``None``, consults
        :func:`strict_detection_enabled`.

    Returns
    -------
    bool
        True  → ring is a legitimate hapto candidate (or filter disabled),
        False → ring is too far from the metal to be hapto-bonded; the
                caller should SKIP enumerating η-modes for this ring.

    Universal: graph-only, no SMILES-specific shortcuts.  When the
    molecule contains no metal (e.g. a stand-alone ligand SMILES like
    ``c1ccccc1``) the filter returns True (we cannot tell from a metal-
    free graph which rings would be hapto-bound to a hypothetical metal,
    so we preserve the question "what hapto modes can this ligand
    support in principle?").

    Deterministic; no side-effects.
    """
    if strict is None:
        strict = strict_detection_enabled()
    if not strict:
        return True
    if mol is None or ring_atoms is None or len(ring_atoms) == 0:
        return True
    try:
        n_atoms = int(mol.GetNumAtoms())
    except Exception:
        return True
    if metal_idxs is None:
        metal_idxs = _metal_atoms_from_mol(mol)
    if not metal_idxs:
        # No metal context — preserve legacy behaviour (the caller is
        # asking "can this ligand support hapto modes?", which is a yes
        # for any aromatic ring).
        return True
    # 1. Direct metal-neighbour atoms (first coordination sphere).
    coord = _metal_coord_atoms(mol, metal_idxs)
    if not coord:
        # Metal is in the graph but has no detectable neighbours
        # (e.g. metal-only fragment) -- can't be hapto.
        return False
    # Fast path: any ring atom IS a metal neighbour ⇒ definitely hapto-able.
    ring_set: Set[int] = {int(i) for i in ring_atoms}
    if ring_set & coord:
        return True
    if max_dist is None:
        max_dist = _max_dist()
    if max_dist <= 0:
        # The fast path already returned True for ring-atom ∈ coord.
        # max_dist = 0 means we require ring atoms themselves to be coord
        # atoms.  Since they aren't, fail.
        return False
    # 2. BFS from the coordination atoms in the metal-removed graph;
    # any ring atom reached within ``max_dist`` bonds qualifies.
    adj = _build_metal_removed_adj(mol, metal_idxs)
    dist = _bfs_distance(n_atoms, adj, sources=coord, max_dist=max_dist)
    for r in ring_set:
        if 0 <= r < n_atoms and dist[r] <= max_dist:
            return True
    return False


__all__ = [
    "strict_detection_enabled",
    "is_legitimate_hapto_candidate",
]
