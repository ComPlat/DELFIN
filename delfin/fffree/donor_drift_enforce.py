"""delfin.fffree.donor_drift_enforce — universal post-construction donor-shell
distance enforcement (Bug class #1: donor_drift, 38.9% of V3 voll-pool).

Diagnosis
=========
Auto-Diagnostic on the V3 voll-pool (commit ``c1e0fde``, 6627 files):
2578 / 6627 (38.9 %) structures have a declared donor sitting at 2.4-3.0 Å from
the metal — drifted OUT of the first coordination shell.  ~92 % of those cases
originate from the ``native`` constructive path: the donor is placed at the
correct M-D distance during ``_orient_chelate_to_vertices`` / monodentate
placement, but a downstream corrector (sp2-flatten, GRIP polish, amide-VSEPR,
oxoanion-VSEPR, refine clash-relief, …) shifts it outwards until ``r(M-D)`` is
above the bonding range.  No existing module pulls it back.

The complementary detection module :mod:`construction_sanity` already labels
this mode as ``md_drifted`` and **rejects** the entire build when the drift
exceeds 1.20 × ``md_target``.  This module is the *constructive* counterpart:
instead of rejecting, it projects the drifted donor back onto the polyhedron
vertex direction at the ideal ``md_target`` and translates the bonded ligand
subtree by the same delta, preserving ligand-internal geometry.

Universal contract
==================
* No per-class branches.  Pure geometry: polyhedron vertex direction + scalar
  M-D distance.  Reuses :func:`polyhedra.ref_vectors` and
  :func:`polyhedra.md_distance` — the single sources of truth.
* For every declared σ-donor we hard-enforce
  ``r(M-D) <= md_target * MAX_FACTOR`` (default ``MAX_FACTOR = 1.15``,
  configurable via ``DELFIN_FFFREE_DONOR_DRIFT_MAX_FACTOR``).
* If the donor has a known polyhedron-vertex index we project it back along
  the vertex unit direction; otherwise we project along its CURRENT
  M-to-donor unit direction (orientation-preserving fallback).
* The donor moves by ``delta = d_new - d_now``.  The same ``delta`` is added
  to every atom in the bonded subtree rooted at the donor, blocked at the
  metal and at every OTHER donor (so other coordination sites stay pinned and
  ligand-internal geometry is preserved by translation only).
* Hapto-π donors (denticity-1 site spread over a ring of carbons) are skipped
  by default — they are governed by the M-centroid distance, not per-atom
  M-D, and have their own corrector (:mod:`hapto_honest_construction`).
* The metal at ``metal_idx=0`` is NEVER moved.

Env flags
=========
``DELFIN_FFFREE_DONOR_DRIFT_ENFORCE=1``
    Activate the post-construction drift repair.  Default OFF =
    byte-identical to HEAD when unset.
``DELFIN_FFFREE_DONOR_DRIFT_MAX_FACTOR=<float>``
    Maximum allowed drift ratio ``r(M-D) / md_target``.  Default ``1.15``.

Determinism
===========
The enforcement walks donors in sorted (ascending global index) order and the
subtree BFS uses sorted neighbour iteration, so the result is byte-identical
across runs / machines for a given input (matches the rest of the FF-free
pipeline).
"""
from __future__ import annotations

import os
from typing import Iterable, List, Mapping, Optional, Sequence, Set, Tuple

import numpy as np

# ---------------------------------------------------------------------------
# Public env-gate helpers.
# ---------------------------------------------------------------------------


def enforce_active() -> bool:
    """``True`` when ``DELFIN_FFFREE_DONOR_DRIFT_ENFORCE=1``."""
    return os.environ.get("DELFIN_FFFREE_DONOR_DRIFT_ENFORCE", "0") == "1"


def max_factor() -> float:
    """Maximum allowed ``r(M-D) / md_target`` before drift repair fires.

    Configurable via ``DELFIN_FFFREE_DONOR_DRIFT_MAX_FACTOR`` (default 1.15).
    Clamped to ``[1.01, 2.0]`` to keep the repair physically meaningful.
    """
    try:
        v = float(os.environ.get("DELFIN_FFFREE_DONOR_DRIFT_MAX_FACTOR", "1.15"))
    except (TypeError, ValueError):
        v = 1.15
    if v < 1.01:
        v = 1.01
    elif v > 2.0:
        v = 2.0
    return v


# ---------------------------------------------------------------------------
# Internal helpers.
# ---------------------------------------------------------------------------


def _bfs_subtree(
    root: int,
    bonds: Sequence[Tuple[int, int]],
    blocked: Set[int],
    n_atoms: int,
) -> List[int]:
    """Return the set of atom indices in the bonded subtree rooted at ``root``,
    blocked at every index in ``blocked`` (which always contains the metal and
    the other donors).  ``root`` itself IS included.  Deterministic — neighbour
    iteration is performed in sorted order, BFS preserves first-discovery
    order.

    Bonds are treated as undirected.  Self-loops and out-of-range indices are
    ignored.  When the graph is empty (e.g. monatomic ligand), the result is
    simply ``[root]``.
    """
    if not (0 <= root < n_atoms):
        return [root]
    adj: List[List[int]] = [[] for _ in range(n_atoms)]
    for b in bonds:
        try:
            i, j = int(b[0]), int(b[1])
        except (TypeError, ValueError, IndexError):
            continue
        if i == j:
            continue
        if not (0 <= i < n_atoms) or not (0 <= j < n_atoms):
            continue
        adj[i].append(j)
        adj[j].append(i)
    for k in range(n_atoms):
        # Deterministic: sort+dedup so the BFS order is reproducible
        # cross-machine / cross-RDKit-version.
        adj[k] = sorted(set(adj[k]))

    visited: Set[int] = {root}
    order: List[int] = [root]
    frontier: List[int] = [root]
    blocked_set = set(int(x) for x in blocked)
    blocked_set.discard(root)
    while frontier:
        nxt: List[int] = []
        for u in frontier:
            for v in adj[u]:
                if v in visited or v in blocked_set:
                    continue
                visited.add(v)
                order.append(v)
                nxt.append(v)
        # Deterministic BFS layer expansion: walk frontier in sorted order so
        # the resulting subtree-list is independent of the input bond order.
        frontier = sorted(nxt)
    return order


def _vertex_unit(
    geometry: str,
    vertex_idx: int,
) -> Optional[np.ndarray]:
    """Return the unit vector of polyhedron vertex ``vertex_idx`` for
    ``geometry``, or ``None`` if the geometry / vertex is unknown."""
    try:
        from delfin.fffree.polyhedra import ref_vectors as _ref
    except ImportError:
        return None
    try:
        V = _ref(str(geometry))
    except Exception:
        return None
    try:
        v = np.asarray(V[int(vertex_idx)], dtype=float)
    except (IndexError, TypeError, ValueError):
        return None
    n = float(np.linalg.norm(v))
    if not np.isfinite(n) or n <= 1e-9:
        return None
    return v / n


def _md_target(metal_sym: str, donor_sym: str) -> float:
    """Wrapper for :func:`polyhedra.md_distance` with safe fallback."""
    try:
        from delfin.fffree.polyhedra import md_distance as _md
    except ImportError:
        return 2.0
    try:
        v = float(_md(str(metal_sym), str(donor_sym)))
    except Exception:
        return 2.0
    if not np.isfinite(v) or v <= 0.5:
        return 2.0
    return v


# ---------------------------------------------------------------------------
# Public API.
# ---------------------------------------------------------------------------


def enforce_donor_shell(
    syms: Sequence[str],
    P: np.ndarray,
    metal_idx: int,
    donor_idxs: Sequence[int],
    bonds: Sequence[Tuple[int, int]],
    *,
    geometry: str = "",
    donor_to_vertex: Optional[Mapping[int, int]] = None,
    hapto_donors: Optional[Iterable[int]] = None,
    max_factor_override: Optional[float] = None,
) -> Tuple[np.ndarray, List[dict]]:
    """Post-construction donor-shell distance enforcement.

    For every σ-donor in ``donor_idxs`` (``hapto_donors`` are skipped), check
    ``r(M-D) <= max_factor * md_target``.  If exceeded, project the donor back
    to ``M + v_ideal * md_target`` where ``v_ideal`` is the polyhedron vertex
    unit direction (or the current ``(D-M)/||D-M||`` direction when no vertex
    map is provided).  Apply the same translation delta to every atom in the
    donor's bonded subtree (blocked at the metal and the other donors) so
    ligand-internal geometry is preserved.

    Parameters
    ----------
    syms
        Element symbols, length ``n``.
    P
        Coordinates ``(n, 3)``.
    metal_idx
        Global index of the central metal (NEVER moved).
    donor_idxs
        Global indices of every declared coordination donor.
    bonds
        Iterable of ``(i, j)`` pairs giving the heavy bond connectivity
        (metal-donor bonds are silently filtered out by the blocked-set).
    geometry
        Polyhedron name (e.g. ``"OC-6 octahedron"``).  Used together with
        ``donor_to_vertex`` to pick the ideal vertex direction.
    donor_to_vertex
        Optional mapping ``donor_global_idx -> polyhedron_vertex_idx``.
        When provided AND the geometry is known, the donor is projected to
        the IDEAL vertex direction (orientation correction + distance fix);
        otherwise we keep the current direction (distance fix only).
    hapto_donors
        Optional iterable of donor global indices that belong to a hapto-π
        site (η³-η⁸).  Skipped by the enforcement because their M-D target is
        per-atom, not per-site (covered by :mod:`hapto_honest_construction`).
    max_factor_override
        Optional ratio override (default = :func:`max_factor`).

    Returns
    -------
    (P_new, repairs)
        ``P_new`` is a NEW array — the input is not mutated.  ``repairs`` is
        a list of ``dict`` entries (one per repaired donor) carrying
        ``donor_idx``, ``r_before``, ``r_after``, ``md_target``, ``delta`` and
        ``mode`` (``"vertex"`` or ``"current_dir"``).
    """
    P_in = np.asarray(P, dtype=float)
    n = P_in.shape[0]
    out = P_in.copy()
    repairs: List[dict] = []

    if n == 0 or len(donor_idxs) == 0:
        return out, repairs
    try:
        m = int(metal_idx)
    except (TypeError, ValueError):
        return out, repairs
    if not (0 <= m < n):
        return out, repairs

    cap = float(max_factor_override) if max_factor_override is not None else max_factor()
    M = out[m].copy()

    skip_donors = set(int(x) for x in (hapto_donors or []) if 0 <= int(x) < n)
    blocked = set(int(d) for d in donor_idxs if 0 <= int(d) < n)
    blocked.add(m)

    # Walk donors deterministically (sorted by ascending global index).
    for d_raw in sorted(set(int(d) for d in donor_idxs)):
        if not (0 <= d_raw < n) or d_raw == m or d_raw in skip_donors:
            continue
        donor_sym = str(syms[d_raw]) if d_raw < len(syms) else ""
        metal_sym = str(syms[m]) if m < len(syms) else ""
        md_tgt = _md_target(metal_sym, donor_sym)
        vec = out[d_raw] - M
        r_now = float(np.linalg.norm(vec))
        if r_now <= 1e-9:
            # Degenerate (collapsed) donor — leave to other modules.
            continue
        if r_now <= cap * md_tgt:
            continue

        # Direction: polyhedron vertex when available, else current dir.
        mode = "current_dir"
        u: Optional[np.ndarray] = None
        if donor_to_vertex is not None and geometry:
            v_idx = donor_to_vertex.get(int(d_raw))
            if v_idx is not None:
                u = _vertex_unit(geometry, int(v_idx))
                if u is not None:
                    mode = "vertex"
        if u is None:
            u = vec / r_now

        d_new = M + u * md_tgt
        delta = d_new - out[d_raw]
        if not np.all(np.isfinite(delta)):
            continue

        # BFS subtree rooted at the donor, blocked at metal + other donors.
        # The donor itself is included in the result.
        # We exclude OTHER donors from the blocked set ONLY for this donor's
        # own root traversal (the donor must move); other donors stay pinned.
        local_blocked = set(blocked)
        local_blocked.discard(d_raw)
        subtree = _bfs_subtree(d_raw, bonds, local_blocked, n)

        # Translate the entire subtree by ``delta`` so ligand-internal
        # geometry is preserved (pure rigid translation -- no internal
        # distortion).  The metal (and other donors) are NOT in the subtree
        # because they were blocked.
        for k in subtree:
            out[k] = out[k] + delta

        r_after = float(np.linalg.norm(out[d_raw] - M))
        repairs.append({
            "donor_idx": int(d_raw),
            "donor_sym": donor_sym,
            "md_target": float(md_tgt),
            "r_before": r_now,
            "r_after": r_after,
            "delta": [float(delta[0]), float(delta[1]), float(delta[2])],
            "mode": mode,
            "subtree_size": len(subtree),
        })

    return out, repairs


# ---------------------------------------------------------------------------
# Convenience: detect drift without repairing (useful for tests / diagnostics).
# ---------------------------------------------------------------------------


def find_drifted_donors(
    syms: Sequence[str],
    P: np.ndarray,
    metal_idx: int,
    donor_idxs: Sequence[int],
    *,
    max_factor_override: Optional[float] = None,
) -> List[dict]:
    """List declared donors whose ``r(M-D) > max_factor * md_target``.

    Pure read-only — coordinates are not modified.  Returns a list of dicts
    with ``donor_idx``, ``donor_sym``, ``md_target``, ``r`` and ``ratio``.
    """
    P_in = np.asarray(P, dtype=float)
    n = P_in.shape[0]
    out: List[dict] = []
    if n == 0:
        return out
    try:
        m = int(metal_idx)
    except (TypeError, ValueError):
        return out
    if not (0 <= m < n):
        return out
    cap = float(max_factor_override) if max_factor_override is not None else max_factor()
    M = P_in[m]
    metal_sym = str(syms[m]) if m < len(syms) else ""
    for d_raw in sorted(set(int(d) for d in donor_idxs)):
        if not (0 <= d_raw < n) or d_raw == m:
            continue
        donor_sym = str(syms[d_raw]) if d_raw < len(syms) else ""
        md_tgt = _md_target(metal_sym, donor_sym)
        r = float(np.linalg.norm(P_in[d_raw] - M))
        if md_tgt <= 0:
            continue
        ratio = r / md_tgt
        if ratio > cap:
            out.append({
                "donor_idx": int(d_raw),
                "donor_sym": donor_sym,
                "md_target": float(md_tgt),
                "r": r,
                "ratio": ratio,
            })
    return out
