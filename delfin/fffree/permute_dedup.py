"""permute_dedup.py — UNIVERSAL permutation-invariant duplicate removal (FF-free).

ROOT-CAUSE fix for the #1 user-reported pool defect: "es gibt sehr viele gleiche
Strukturen im Pool — Ununterscheidbarkeit von Atomen — Duplikate müssen raus".

The defect
----------
The pre-existing ensemble dedup (``conformer_complete.apply_to_ensemble`` /
``assemble_complex._complex_rmsd`` / ``_conformer_pool``) compares two frames with
a heavy-atom **Kabsch RMSD that uses the FIXED atom order** of the two coordinate
sets.  Atom indices are an arbitrary labelling artefact, NOT chemistry.  Two
structures that are IDENTICAL up to relabeling of indistinguishable atoms (a pair
of identical ligands swapped, a methyl rotated onto a symmetry-equivalent
position, any molecular-graph automorphism) therefore receive a LARGE
fixed-order RMSD and are BOTH kept — genuine duplicates survive.

The fix (this module)
---------------------
A **permutation-invariant** duplicate test: two frames A, B are duplicates iff

    min over automorphism permutations π of  heavy-RMSD(A, π(B))  <  threshold

where π ranges over the automorphisms of the molecular graph (permutations that
(a) preserve element identity, (b) preserve adjacency, and (c) map every atom
only WITHIN its own graph-symmetry orbit).  Only indistinguishable atoms can be
swapped, so the test can NEVER make two genuinely-different structures look equal:
the geometric RMSD floor still has to be cleared, and distinct geometric isomers /
conformers / stereoisomers do not coincide under a graph relabelling.

This runs as a STRICTER dedup LAYER at the SAME ensemble chokepoint where the
fixed-order RMSD dedup runs (the public ``smiles_to_xyz_isomers`` boundary, over
the final native+legacy+conformer-completeness union).  The first member of each
permutation-equivalence cluster is kept; deterministic input order is preserved.

Never-collapse guard (most important property)
----------------------------------------------
* Atoms are permuted ONLY within their graph-symmetry orbit (same element + same
  iterated-degree refinement color), so a permutation is always a true graph
  automorphism — it can never relabel an atom onto a chemically-different atom.
* The duplicate decision is the GEOMETRIC RMSD under the best automorphism; two
  distinct geometric isomers (all-cis vs all-trans octahedral, mer vs fac, …) or
  two distinct conformers have geometries that differ by MORE than the threshold
  even after the best symmetry relabelling, so they are KEPT.
* Chirality: the automorphism search is purely topological, but the test that
  decides duplicate-or-not is the proper-rotation (reflection-FORBIDDEN) Kabsch
  RMSD.  Two mirror-image (enantiomeric) frames cannot be superposed by a proper
  rotation, so their RMSD stays high and they are KEPT.  ``use_chirality`` (on by
  default) additionally seeds each atom's refinement color with a local geometric
  parity, so a stereocenter's neighbours are never treated as interchangeable —
  belt-and-braces against ever merging a true stereoisomer.
* On ANY ambiguity / error / cap overflow the module falls back to the IDENTITY
  permutation (= the existing fixed-order behaviour) and never raises, so it can
  only ever remove TRUE duplicates, never distinct structures.

Design constraints
------------------
* UNIVERSAL / graph-only.  The automorphism group is derived purely from the
  perceived connectivity (the same Open Babel perception the rest of the pipeline
  uses for DOF / topology) plus element identity.  NEVER special-cased on a
  SMILES string, refcode or element.
* FF-FREE.  Pure graph theory + Kabsch geometry.  No force field.
* DETERMINISTIC.  Fixed refinement ordering, fixed backtracking order, no RNG /
  time / hash dependence.  ``PYTHONHASHSEED`` independent (no set / dict iteration
  order leaks into output).
* env-gated, default-OFF, BYTE-IDENTICAL when off.  Master flag
  ``DELFIN_FFFREE_PERMUTE_DEDUP`` (default ``0``).  Unset / "0" -> the public
  ``dedup_ensemble`` is the identity (returns the input list untouched).
* NEVER raises.  Any failure returns the input unchanged.

Env-flags (read once per call)
------------------------------
``DELFIN_FFFREE_PERMUTE_DEDUP``        (default ``0``) — master switch.
``DELFIN_FFFREE_PERMUTE_DEDUP_RMSD``   (default ``0.5``) — Angstrom heavy-RMSD
                                         duplicate threshold (matches the existing
                                         ``DELFIN_FFFREE_CONF_RMSD`` default).
``DELFIN_FFFREE_PERMUTE_DEDUP_MAXPERM``(default ``4096``) — cap on enumerated
                                         automorphisms; above it, fall back to the
                                         identity permutation (= current behaviour)
                                         and log.
``DELFIN_FFFREE_PERMUTE_DEDUP_CHIRAL`` (default ``1``) — seed refinement colors
                                         with local geometric parity so atoms
                                         around a stereocenter are never swapped.

Spec lives in this docstring (re-buildable from it).
"""

from __future__ import annotations

import math
import os
from typing import Dict, List, Optional, Sequence, Tuple

from delfin.common.logging import get_logger

logger = get_logger(__name__)


# ---------------------------------------------------------------------------
# Environment configuration helpers.
# ---------------------------------------------------------------------------


def _env_bool(name: str, default: bool = False) -> bool:
    raw = os.environ.get(name)
    if raw is None:
        return default
    return raw.strip().lower() in ("1", "true", "yes", "on")


def _env_int(name: str, default: int, lo: int = 1, hi: int = 1 << 24) -> int:
    raw = os.environ.get(name)
    if raw is None:
        return default
    try:
        return max(lo, min(hi, int(float(raw.strip()))))
    except Exception:
        return default


def _env_float(name: str, default: float, lo: float = 0.0, hi: float = 1e6) -> float:
    raw = os.environ.get(name)
    if raw is None:
        return default
    try:
        return max(lo, min(hi, float(raw.strip())))
    except Exception:
        return default


def _is_enabled() -> bool:
    return os.environ.get("DELFIN_FFFREE_PERMUTE_DEDUP", "0") == "1"


# ---------------------------------------------------------------------------
# Heavy-atom Kabsch RMSD under a permutation (FF-free, proper rotation only).
# ---------------------------------------------------------------------------


def _kabsch_rmsd_perm(
    heavy_idx: Sequence[int],
    A,
    B,
    perm: Sequence[int],
) -> float:
    """Heavy-atom RMSD of A vs the permuted B after a proper-rotation Kabsch
    superposition, using the heavy atoms only.

    ``A`` and ``B`` are full (N,3) coordinate arrays; ``perm`` maps every atom
    index i -> perm[i] (B atom that plays the role of A atom i).  Reflections are
    FORBIDDEN (determinant-corrected) so enantiomeric frames never align.

    Returns a large number on degeneracy / error.  Deterministic; never raises.
    """
    try:
        import numpy as np
        P = A[list(heavy_idx)]
        Q = B[[perm[i] for i in heavy_idx]]
        if P.shape != Q.shape or P.shape[0] < 1:
            return 1e9
        Pc = P - P.mean(0)
        Qc = Q - Q.mean(0)
        if P.shape[0] >= 3:
            H = Pc.T @ Qc
            U, _S, Vt = np.linalg.svd(H)
            d = np.sign(np.linalg.det(Vt.T @ U.T))
            D = np.diag([1.0, 1.0, d])
            R = Vt.T @ D @ U.T
            diff = (R @ Pc.T).T - Qc
        else:
            diff = Pc - Qc
        return float(math.sqrt((diff * diff).sum() / P.shape[0]))
    except Exception:
        return 1e9


# ---------------------------------------------------------------------------
# Graph-symmetry orbits via iterated-degree (color) refinement.
# ---------------------------------------------------------------------------


def _initial_colors(
    atomic_nums: Sequence[int],
    neighbours: Sequence[Sequence[int]],
    coords=None,
    use_chirality: bool = True,
) -> List[int]:
    """Initial atom colors: (element, degree[, local geometric parity]).

    The optional geometric-parity seed splits the two faces of any atom whose
    neighbourhood is chiral in 3D, so the refinement (and hence the automorphism
    search) never treats the neighbours of a stereocenter as interchangeable.
    Deterministic; parity is quantised by SIGN only (robust to tiny coordinate
    noise; a flat / achiral center contributes 0).
    """
    n = len(atomic_nums)
    base: List[Tuple] = []
    for i in range(n):
        parity = 0
        if use_chirality and coords is not None:
            parity = _local_parity(i, neighbours[i], coords)
        base.append((int(atomic_nums[i]), len(neighbours[i]), parity))
    # map distinct tuples -> small dense ints in sorted (deterministic) order
    order = {key: c for c, key in enumerate(sorted(set(base)))}
    return [order[k] for k in base]


def _local_parity(i: int, nbrs: Sequence[int], coords) -> int:
    """Sign of the signed volume of the first four distinct neighbour directions
    of atom *i* (a CIP-free, label-independent 3D chirality seed).  0 for fewer
    than 4 neighbours or a (near-)planar / degenerate center.

    Uses the SIGN only, so it is robust and deterministic.  The neighbour order
    is the perceived-graph order; for a TRUE automorphism the swap of two
    symmetry-equivalent neighbours flips the sign, which correctly forbids the
    swap when the two faces are geometrically distinct (a stereocenter).
    """
    try:
        import numpy as np
        if len(nbrs) < 4:
            return 0
        o = np.asarray(coords[i], float)
        v = [np.asarray(coords[j], float) - o for j in nbrs[:4]]
        vol = float(np.dot(np.cross(v[1] - v[0], v[2] - v[0]), v[3] - v[0]))
        if vol > 1e-6:
            return 1
        if vol < -1e-6:
            return -1
        return 0
    except Exception:
        return 0


def _refine_colors(
    colors: List[int],
    neighbours: Sequence[Sequence[int]],
) -> List[int]:
    """1-WL color refinement to a fixed point: a node's new color is a function
    of its own color and the MULTISET of its neighbours' colors.  Deterministic
    (sorted signatures -> dense relabel in sorted order).  Returns stable
    orbit-like color classes (a sound over-approximation of automorphism orbits;
    the backtracking search below enforces exact automorphism).
    """
    n = len(colors)
    cur = list(colors)
    for _ in range(n + 1):
        sig = []
        for i in range(n):
            sig.append((cur[i], tuple(sorted(cur[j] for j in neighbours[i]))))
        order = {key: c for c, key in enumerate(sorted(set(sig)))}
        new = [order[s] for s in sig]
        if new == cur:
            break
        cur = new
    return cur


# ---------------------------------------------------------------------------
# Automorphism enumeration (deterministic backtracking, capped).
# ---------------------------------------------------------------------------


def _enumerate_automorphisms(
    atomic_nums: Sequence[int],
    neighbours: Sequence[Sequence[int]],
    coords=None,
    use_chirality: bool = True,
    max_perm: int = 4096,
) -> Optional[List[List[int]]]:
    """Enumerate graph automorphisms (including the identity) as a list of
    permutations ``perm`` with ``perm[i]`` = image of atom i.

    A permutation is an automorphism iff it preserves the refined color of every
    atom AND maps the neighbour SET of every atom onto the neighbour set of its
    image.  Candidate images for each atom are restricted to its color class, so
    only graph-indistinguishable atoms are ever swapped.

    Deterministic backtracking with the most-constrained-atom-first variable
    order.  Returns ``None`` (= "fall back to identity only") if the count would
    exceed ``max_perm`` (e.g. a fullerene-like or huge-symmetry case), so the
    pass degrades gracefully to the existing fixed-order behaviour.  Never raises.
    """
    try:
        n = len(atomic_nums)
        if n == 0:
            return [[]]
        nbr_sets = [set(neighbours[i]) for i in range(n)]
        colors = _initial_colors(atomic_nums, neighbours, coords, use_chirality)
        colors = _refine_colors(colors, neighbours)

        # candidate image set per atom = same refined color (deterministic order)
        by_color: Dict[int, List[int]] = {}
        for i in range(n):
            by_color.setdefault(colors[i], []).append(i)
        candidates = [sorted(by_color[colors[i]]) for i in range(n)]

        # If every atom is in a singleton color class the only automorphism is
        # the identity (common for low-symmetry complexes) -> cheap exit.
        if all(len(candidates[i]) == 1 for i in range(n)):
            return [list(range(n))]

        # variable order: fewest candidates first (then lower index) for a tight,
        # deterministic search tree.
        var_order = sorted(range(n), key=lambda i: (len(candidates[i]), i))

        results: List[List[int]] = []
        perm = [-1] * n
        used = [False] * n
        overflow = [False]

        def backtrack(k: int) -> None:
            if overflow[0]:
                return
            if k == n:
                results.append(list(perm))
                if len(results) > max_perm:
                    overflow[0] = True
                return
            i = var_order[k]
            for j in candidates[i]:
                if used[j]:
                    continue
                # color must match (guaranteed by candidates) and the mapping must
                # respect adjacency for every ALREADY-assigned neighbour of i.
                ok = True
                for nb in neighbours[i]:
                    pj = perm[nb]
                    if pj == -1:
                        continue
                    # i~nb  =>  j~pj   (and the reverse is enforced symmetrically
                    # when nb is processed); also forbid spurious edges.
                    if (pj not in nbr_sets[j]):
                        ok = False
                        break
                if not ok:
                    continue
                # Additional pruning: j's already-assigned neighbours among the
                # images must be neighbours-of-i preimages (forbid extra edges).
                # Enforced implicitly by the symmetric check above across atoms.
                perm[i] = j
                used[j] = True
                backtrack(k + 1)
                perm[i] = -1
                used[j] = False
                if overflow[0]:
                    return

        backtrack(0)
        if overflow[0]:
            logger.debug("permute-dedup automorphism cap %d exceeded -> identity", max_perm)
            return None
        if not results:
            return [list(range(n))]
        # ensure identity present and deterministic order
        identity = list(range(n))
        results.sort()
        if identity not in results:
            results.insert(0, identity)
        return results
    except Exception:
        return None


# ---------------------------------------------------------------------------
# Per-structure automorphism cache (the graph is shared by all frames of one
# structure: same atom order, same connectivity).
# ---------------------------------------------------------------------------


def _automorphisms_for_xyz(
    xyz: str,
    use_chirality: bool,
    max_perm: int,
) -> Tuple[Optional[List[str]], Optional[List[List[int]]], Optional[List[int]]]:
    """Build (symbols, automorphism-permutation-list, heavy-index-list) for a
    structure's representative frame.  Returns (None, None, None) on failure
    (caller then falls back to identity / fixed-order behaviour).
    """
    try:
        from delfin._rotamer_diversity import (
            _build_ob_mol_from_xyz, _graph_from_ob, _parse_delfin_xyz,
        )
        symbols, coords = _parse_delfin_xyz(xyz)
        ob_mol = _build_ob_mol_from_xyz(xyz)
        if ob_mol is None:
            return symbols, [list(range(len(symbols)))], _heavy(symbols)
        graph = _graph_from_ob(ob_mol)
        if not graph or int(graph.get("n_atoms", 0)) != len(symbols):
            return symbols, [list(range(len(symbols)))], _heavy(symbols)
        autos = _enumerate_automorphisms(
            graph["atomic_nums"], graph["neighbours"],
            coords=coords, use_chirality=use_chirality, max_perm=max_perm,
        )
        if autos is None:
            # cap overflow / error -> identity only (current behaviour)
            autos = [list(range(len(symbols)))]
        return symbols, autos, _heavy(symbols)
    except Exception:
        return None, None, None


def _heavy(symbols: Sequence[str]) -> List[int]:
    heavy = [k for k in range(len(symbols)) if symbols[k] != "H"]
    if len(heavy) < 3:
        return list(range(len(symbols)))
    return heavy


def _coords_array(xyz: str):
    """Parse one frame -> (N,3) numpy array, or None."""
    try:
        import numpy as np
        from delfin._rotamer_diversity import _parse_delfin_xyz
        _syms, co = _parse_delfin_xyz(xyz)
        return np.asarray(co, float)
    except Exception:
        return None


def is_permutation_duplicate(
    xyz_a: str,
    xyz_b: str,
    rmsd: float = 0.5,
    use_chirality: bool = True,
    max_perm: int = 4096,
) -> bool:
    """Public single-pair test: True iff frame B is a permutation-duplicate of A
    (min over graph automorphisms of heavy-RMSD(A, π(B)) < ``rmsd``).

    Deterministic; never raises (returns False on any error so distinct frames
    are never wrongly merged)."""
    try:
        symbols, autos, heavy = _automorphisms_for_xyz(xyz_a, use_chirality, max_perm)
        if symbols is None or autos is None:
            return False
        A = _coords_array(xyz_a)
        B = _coords_array(xyz_b)
        if A is None or B is None or A.shape != B.shape:
            return False
        if len(symbols) != A.shape[0]:
            return False
        best = 1e9
        for perm in autos:
            r = _kabsch_rmsd_perm(heavy, A, B, perm)
            if r < best:
                best = r
            if best < rmsd:
                return True
        return best < rmsd
    except Exception:
        return False


# ---------------------------------------------------------------------------
# Public ensemble-level dedup pass.
# ---------------------------------------------------------------------------


def dedup_ensemble(isomers):
    """Permutation-invariant duplicate removal over the emitted ``(xyz, label[, …])``
    ensemble of ONE structure.

    Identity (returns *isomers* untouched) unless
    ``DELFIN_FFFREE_PERMUTE_DEDUP=1`` -> output BYTE-IDENTICAL when off.

    Contract:
      * Frames are grouped by atom-count signature (different N -> different
        structure -> never compared / never merged).
      * Within a group, a candidate frame is dropped iff it is a permutation-
        duplicate (min-over-automorphisms heavy-Kabsch-RMSD < threshold) of an
        ALREADY-KEPT frame.  The automorphism group is taken from the KEPT frame
        being compared against (cached per kept frame), NOT from a single group
        representative: perceived connectivity is geometry-dependent (a bent
        chelate / close contact can be perceived as bonded in one conformer but
        not another), so different frames legitimately expose different symmetry,
        and the kept frame's own group is the correct relabelling set for that
        comparison.
      * The FIRST member of every cluster is kept; deterministic input order is
        preserved among kept frames.
      * Distinct geometric isomers / conformers / stereoisomers have geometries
        that differ by more than the threshold even under the best symmetry
        relabelling and are therefore KEPT (never-collapse guard).
      * On any error / cap overflow the pass falls back to identity-permutation
        comparison (= existing fixed-order behaviour) and never raises.

    Deterministic; never raises (falls back to *isomers* on error)."""
    if not isomers or not _is_enabled():
        return isomers
    try:
        rmsd = _env_float("DELFIN_FFFREE_PERMUTE_DEDUP_RMSD", 0.5, lo=0.0, hi=10.0)
        max_perm = _env_int("DELFIN_FFFREE_PERMUTE_DEDUP_MAXPERM", 4096, lo=1, hi=1 << 22)
        use_chirality = _env_bool("DELFIN_FFFREE_PERMUTE_DEDUP_CHIRAL", True)
    except Exception:
        return isomers

    try:
        # normalise items to (xyz, label, rest) while preserving the original tuple
        items = []
        for item in isomers:
            if isinstance(item, (tuple, list)) and len(item) >= 2:
                items.append((item[0], item[1], tuple(item[2:]), item))
            else:
                items.append((item, "", (), item))

        # cache of (symbols, autos, heavy) keyed by atom-count group rep index
        kept_flags = [True] * len(items)
        # group by parsed atom-count signature (cheap, deterministic)
        groups: Dict[int, List[int]] = {}
        coords_cache: Dict[int, object] = {}
        for idx, (xyz, _lbl, _rest, _orig) in enumerate(items):
            C = _coords_array(xyz)
            coords_cache[idx] = C
            key = -1 if C is None else int(C.shape[0])
            groups.setdefault(key, []).append(idx)

        import numpy as np  # noqa: F401 — ensures numpy path is active

        removed = 0
        auto_cache: Dict[int, Tuple[object, object]] = {}  # kept-frame idx -> (autos, heavy)
        for key in sorted(groups.keys()):
            members = groups[key]
            if key < 0 or len(members) < 2:
                continue  # unparseable group or singleton -> nothing to dedup
            kept_in_group: List[int] = []
            for idx in members:
                B = coords_cache[idx]
                if B is None or B.shape[0] != key:
                    kept_in_group.append(idx)  # keep odd ones, never drop blindly
                    continue
                is_dup = False
                for kidx in kept_in_group:
                    A = coords_cache[kidx]
                    if A is None or A.shape != B.shape:
                        continue
                    # automorphisms from the KEPT frame A (its own perceived graph),
                    # cached per kept frame.  Perceived connectivity is geometry-
                    # dependent, so the kept frame's own symmetry group is the
                    # correct relabelling set for the A-vs-B comparison.
                    if kidx not in auto_cache:
                        _s, autos, heavy = _automorphisms_for_xyz(
                            items[kidx][0], use_chirality, max_perm)
                        if autos is None:
                            autos = [list(range(key))]
                            heavy = list(range(key))
                        auto_cache[kidx] = (autos, heavy)
                    autos, heavy = auto_cache[kidx]
                    best = 1e9
                    for perm in autos:
                        r = _kabsch_rmsd_perm(heavy, A, B, perm)
                        if r < best:
                            best = r
                        if best < rmsd:
                            break
                    if best < rmsd:
                        is_dup = True
                        break
                if is_dup:
                    kept_flags[idx] = False
                    removed += 1
                else:
                    kept_in_group.append(idx)

        if removed:
            logger.debug("permute-dedup removed %d permutation-duplicate frame(s)", removed)
        out = [items[i][3] for i in range(len(items)) if kept_flags[i]]
        return out or isomers
    except Exception as exc:  # pragma: no cover — safety net
        logger.debug("permute-dedup ensemble pass failed: %s", exc)
        return isomers
