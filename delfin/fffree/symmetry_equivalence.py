"""delfin.fffree.symmetry_equivalence — graph-automorphism bond/angle classes.

Per memory ``project_degeneracy_master_doctrine_2026_06_04``: detect
equivalence classes via graph automorphism (networkx), then enforce that
symmetry-equivalent bonds share a parameter — i.e. their lengths should
converge (variance penalty in the GRIP loss).

For LUHMOT: the aromatic 6-ring carries six C-C bonds, all symmetry-equivalent.
The 1.130 Å outlier should be penalised by the within-class variance term.

This module gives two services:

1. ``find_equivalence_classes(syms, bonds)`` -> list of lists of bond tuples
   (each list is an automorphism orbit).
2. ``variance_penalty(P, classes, alpha=10.0)`` -> scalar + gradient
   contribution for the GRIP loss aggregator.

Determinism: GraphMatcher iteration order over automorphisms is non-deterministic
in general; we therefore use a lex-sorted orbit-canonicalisation pass so the
returned classes are stable across runs.  No RNG anywhere.

Env-gate: ``DELFIN_FFFREE_SYMMETRY_EQUIVALENCE=1`` (default OFF, byte-identical
when unset).  Auto-enabled under PURE_TRACK3.
"""
from __future__ import annotations

import os
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

import numpy as np


def is_enabled() -> bool:
    if os.environ.get("DELFIN_FFFREE_SYMMETRY_EQUIVALENCE", "0") == "1":
        return True
    if os.environ.get("DELFIN_FFFREE_PURE_TRACK3", "0") == "1":
        return True
    return False


def _build_graph(syms: Sequence[str], bonds: Iterable[Tuple[int, int]]):
    """Build an undirected networkx Graph with atom-symbol node attributes.

    Aromatic flag (when known) is folded into the node label so aromatic
    sp2 carbons are NOT equated with sp3 carbons.  Here we accept a plain
    bond list; aromaticity gets added in via a separate attribute if the
    caller passes ``aromatic_atoms``.
    """
    import networkx as nx
    G = nx.Graph()
    for i, s in enumerate(syms):
        G.add_node(int(i), elem=str(s))
    for (a, b) in bonds:
        G.add_edge(int(a), int(b))
    return G


def find_equivalence_classes(
    syms: Sequence[str],
    bonds: Iterable[Tuple[int, int]],
    aromatic_atoms: Optional[Sequence[bool]] = None,
    bond_attrs: Optional[Dict[Tuple[int, int], Dict[str, object]]] = None,
) -> List[List[Tuple[int, int]]]:
    """Return bond-orbits under the graph's automorphism group.

    Parameters
    ----------
    syms
        Atom symbols.
    bonds
        Iterable of unordered (a, b) pairs (atom indices).
    aromatic_atoms
        Optional aromaticity flag per atom — folded into the node label so
        an aromatic C is NOT equated with an sp3 C.  Highly recommended.
    bond_attrs
        Optional per-bond attributes (e.g. {'aromatic': True/False}) used as
        an edge label.  Bonds with mismatching attrs are NOT equated.

    Returns
    -------
    A list of bond-classes; each class is a lex-sorted list of (a, b) tuples.
    Singleton classes (bonds with no symmetry partner) are also returned.
    The outer list itself is lex-sorted by the canonical (smallest) bond
    of each class — fully deterministic.

    Implementation
    --------------
    We compute the automorphism group of the *node-labelled* graph using
    networkx's GraphMatcher with the identity isomorphism (G iso G), then
    apply each automorphism to every bond and accumulate orbits.
    """
    import networkx as nx
    from networkx.algorithms.isomorphism import GraphMatcher

    bond_list = sorted({tuple(sorted((int(a), int(b)))) for (a, b) in bonds})
    if not bond_list:
        return []

    G = _build_graph(syms, bond_list)
    if aromatic_atoms is not None:
        for i, fl in enumerate(aromatic_atoms):
            if i in G.nodes:
                G.nodes[i]["aromatic"] = bool(fl)
    if bond_attrs:
        for (a, b), attrs in bond_attrs.items():
            ab = tuple(sorted((int(a), int(b))))
            if G.has_edge(ab[0], ab[1]):
                for k, v in attrs.items():
                    G.edges[ab[0], ab[1]][k] = v

    # Node match: same element AND same aromatic flag (if provided).
    def nm(n1, n2):
        if n1.get("elem") != n2.get("elem"):
            return False
        if "aromatic" in n1 and "aromatic" in n2:
            if n1["aromatic"] != n2["aromatic"]:
                return False
        return True

    def em(e1, e2):
        # If aromatic attr present on either, must match
        if "aromatic" in e1 or "aromatic" in e2:
            return e1.get("aromatic", False) == e2.get("aromatic", False)
        return True

    GM = GraphMatcher(G, G, node_match=nm, edge_match=em)
    # Cap the search: aut group can be huge; orbit closure converges fast,
    # so we union edges progressively and bail when a configurable cap of
    # mappings has been processed.
    max_aut = int(os.environ.get("DELFIN_FFFREE_SYMM_MAX_AUT", "256"))

    # Initialise: each bond is its own class.
    parent = {b: b for b in bond_list}

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    def union(a, b):
        ra, rb = find(a), find(b)
        if ra == rb:
            return
        if ra <= rb:
            parent[rb] = ra
        else:
            parent[ra] = rb

    count = 0
    try:
        for mapping in GM.isomorphisms_iter():
            for (a, b) in bond_list:
                ma, mb = mapping[a], mapping[b]
                mab = tuple(sorted((ma, mb)))
                if mab in parent:
                    union((a, b), mab)
            count += 1
            if count >= max_aut:
                break
    except Exception:
        # Defensive: if the matcher chokes, return identity classes.
        pass

    orbits: Dict[Tuple[int, int], List[Tuple[int, int]]] = {}
    for b in bond_list:
        root = find(b)
        orbits.setdefault(root, []).append(b)
    # Sort each orbit + sort the outer list canonically.
    classes = []
    for root in sorted(orbits.keys()):
        classes.append(sorted(orbits[root]))
    classes.sort(key=lambda lst: lst[0])
    return classes


def class_means_and_stds(
    P: np.ndarray, classes: List[List[Tuple[int, int]]]
) -> List[Tuple[float, float]]:
    """For each class, compute (mean, std) of the constituent bond lengths.

    Useful for forensik / metrics output (e.g. show how a per-bond pull has
    collapsed within-class variance).
    """
    out: List[Tuple[float, float]] = []
    for cls in classes:
        ds = [float(np.linalg.norm(P[a] - P[b])) for (a, b) in cls]
        if not ds:
            out.append((0.0, 0.0))
            continue
        mu = float(np.mean(ds))
        sd = float(np.std(ds))
        out.append((mu, sd))
    return out


def variance_penalty(
    P: np.ndarray,
    classes: List[List[Tuple[int, int]]],
    alpha: float = 10.0,
) -> Tuple[float, np.ndarray]:
    """Compute Σ_class α · n_class · Var(class_lengths), with the analytic
    gradient on positions.

    Var(class) = (1/n) Σ_k (d_k - μ)² where d_k = |P[a]-P[b]| and μ is the
    class mean.  Gradient on each endpoint:

        d Var / d d_k = (2/n)(d_k - μ)(1 - 1/n) ≈ (2/n)(d_k - μ)  for large n
        d d_k / d P[a] = (P[a] - P[b]) / d_k

    Singleton classes contribute zero (no variance penalty when only one
    bond).
    """
    P = np.asarray(P, dtype=np.float64)
    grad = np.zeros_like(P)
    total = 0.0
    for cls in classes:
        n = len(cls)
        if n < 2:
            continue
        ds = []
        for (a, b) in cls:
            v = P[a] - P[b]
            d = float(np.linalg.norm(v))
            ds.append(d)
        ds = np.array(ds, dtype=np.float64)
        mu = float(np.mean(ds))
        var = float(np.mean((ds - mu) ** 2))
        total += alpha * n * var

        # Gradient: ∂/∂P[a] [α n Var] = α · ∂Σ_k (d_k - μ)² / ∂P[a]
        #         = 2α Σ_k (d_k - μ) ∂d_k/∂P[a] · (1 - 1/n)·[k==i]+ correction
        # For determinism we use the closed form:
        #   ∂Var/∂d_i = (2/n)(d_i - μ)(1 - 1/n)
        for k, (a, b) in enumerate(cls):
            d_k = ds[k]
            if d_k < 1e-9:
                continue
            coef = alpha * 2.0 * (d_k - mu) * (1.0 - 1.0 / n)
            unit = (P[a] - P[b]) / d_k
            grad[a] += coef * unit
            grad[b] -= coef * unit

    return total, grad


__all__ = [
    "is_enabled",
    "find_equivalence_classes",
    "class_means_and_stds",
    "variance_penalty",
]
