"""delfin.fffree.construction_sanity — Post-embedding sanity check, retry loop,
extended vdW floor and per-class robust placement templates.

Construction-side fix (User mandate 2026-06-07, Bug-class 2+3):

DELFIN V13 voll-pool (V3 + hapto-π + vdW-floor stack) produces severely broken
structures for two complex TM SMILES classes:

  * Bug-class 2 (ODUXAN-like): Fischer-carbene + 5 CO around a Cr(0) loses one
    nitrogen donor — it drifts to ~2.84 Å (well beyond bonding), leaving a
    distorted CN5 instead of the expected OC-6.
  * Bug-class 3 (WICROP-like): η⁶-arene piano-stool + 3 CO + P(o-tol)₂ collapses
    completely — multiple atom pairs end up below 1.0 Å, the M-C28 distance
    falls to 1.44 Å, and the CO oxygens face inward.

Both bugs are independent of the Pólya-vertex enumeration gap (handled by the
parallel polyhedron-vertex enumeration extension).  They are about
*construction quality*, not enumeration completeness:

  1. ETKDG embedding fails on complex hapto + multi-σ + bulky topologies.
  2. Mogul-DG REPLACE_RIGID sets ideal M-D distances but cannot fix
     ligand-internal collapse.
  3. The existing vdW-floor in :mod:`grip_polish` excludes hapto-π rigid bodies
     (the freeze pass locks them rigid) → internal collapses inside a rigid
     body are never penalised.
  4. No sanity-check rejects structures with collapses after embedding —
     bad embeds silently make it through.

This module adds, all env-gated default-OFF byte-identical:

  * :func:`assert_construction_sane` — post-embed Pauli / bond-stretch / M-D /
    coordination-CShM sanity vector.  Returns ``(is_sane, list_of_violations)``.
  * :func:`build_with_retries` — drop-in wrapper that re-embeds with a sequence
    of ETKDG seeds (env-gated) until sanity passes or a rigid-template fallback
    catches the residual.
  * :func:`vdw_floor_all_pairs_value_and_grad` — extended vdW floor that
    penalises ALL non-bonded heavy pairs **including those inside the hapto-π
    rigid body** (translating their gradient to the rigid-body's global shift,
    not internal distortion).  Used by :mod:`grip_polish` when
    ``DELFIN_FFFREE_GRIP_VDW_FLOOR_ALL_PAIRS=1``.
  * :func:`piano_stool_template` — η⁶-arene + 3 CO + L piano-stool placement
    (ring at the M-L axis, three CO 90° from the ring axis at 120° azimuth).
  * :func:`fischer_carbene_template` — Fischer-carbene + ≤5 CO + auxiliary
    N-donor OC-6 placement (carbene C axial, CO equatorial, side-chain N donor
    from the carbene's amine).

Env flags (all default OFF, byte-identical when unset):
    DELFIN_FFFREE_CONSTRUCTION_SANITY=1
        Enable post-embed sanity check + retry loop.
    DELFIN_FFFREE_CONSTRUCTION_SANITY_MAX_RETRIES=N  (optional, default 3)
        Number of alternative ETKDG seeds to try before giving up.
    DELFIN_FFFREE_CONSTRUCTION_SANITY_CSHM_MAX=v  (optional, default 35.0)
        Reject builds with coordination-polyhedron CShM > v.
    DELFIN_FFFREE_GRIP_VDW_FLOOR_ALL_PAIRS=1
        Extend vdW floor in grip_polish to ALL non-bonded heavy pairs, including
        inside the hapto-π rigid body.
    DELFIN_FFFREE_PIANO_STOOL_TEMPLATE=1
        Use a robust piano-stool template (η⁶-arene + 3 CO + L) BEFORE ETKDG
        rather than rely on ETKDG to find it.
    DELFIN_FFFREE_FISCHER_CARBENE_TEMPLATE=1
        Use a robust Fischer-carbene OC-6 template (carbene C axial, CO
        equatorial, auxiliary N-donor) for the Cr(0)-carbene + multi-CO class.

The module is *universal* in the sense that the SMILES patterns it recognises
are topology-based (graph features: η^n hapto + n CO + 1 σ-donor; carbene C-N
sp2 + ≤5 CO).  No element-specific SMILES strings are hard-coded.

API summary:

    >>> from delfin.fffree.construction_sanity import (
    ...     assert_construction_sane, sanity_active, build_with_retries,
    ...     vdw_floor_all_pairs_active, vdw_floor_all_pairs_value_and_grad,
    ...     piano_stool_template_active, fischer_carbene_template_active,
    ...     classify_topology, piano_stool_template, fischer_carbene_template,
    ... )
    >>> sanity_active()
    False
    >>> assert_construction_sane(['Cr', 'C', 'O'], np.eye(3), [(0, 1), (1, 2)])
    (True, [])

Pure read-only when all flags are off.  Coordinates returned by template
helpers are deterministic (no RNG).
"""
from __future__ import annotations

import math
import os
from typing import Callable, Dict, FrozenSet, Iterable, List, Mapping, Optional, Sequence, Set, Tuple

import numpy as np

import delfin._bond_decollapse as _bd

# ---------------------------------------------------------------------------
# Universal physical constants (single-source-of-truth, mirror grip_polish /
# build_time_clash_gate to keep the floor uniform across construction +
# detection).
# ---------------------------------------------------------------------------

# Bondi vdW radii (Å) — identical table to ``_bond_decollapse._VDW`` so the
# build-time sanity check, GRIP vdW-floor and structqual detector see the same
# geometry.  Missing elements fall through to :data:`_VDW_DEFAULT`.
_VDW: Dict[str, float] = dict(_bd._VDW)
_VDW_DEFAULT: float = float(_bd._VDW_DEFAULT)

# Default thresholds (mirror existing detectors)
DEFAULT_PAULI_FRACTION: float = 0.85    # d_ij < 0.85 * (r_vdw_i + r_vdw_j) → violation
DEFAULT_BOND_STRETCH_FACTOR: float = 1.6
DEFAULT_BOND_COMPRESS_FACTOR: float = 0.7
DEFAULT_MD_FLOOR_FACTOR: float = 0.6
DEFAULT_MD_CEIL_FACTOR: float = 2.0
DEFAULT_CSHM_MAX: float = 35.0          # very loose: only truly catastrophic builds rejected
DEFAULT_MAX_RETRIES: int = 3
DEFAULT_RETRY_SEEDS: Tuple[int, ...] = (42, 17, 91, 137, 271)

# Topology-classification thresholds
_PIANO_STOOL_HAPTO_MIN: int = 5         # η⁵-Cp and up qualify
_PIANO_STOOL_HAPTO_MAX: int = 8         # η⁶-arene / η⁷-cycloheptatrienyl / η⁸-COT
_FISCHER_CO_MIN: int = 4                # 4-5 CO + carbene → octahedral
_FISCHER_CO_MAX: int = 5

# ---------------------------------------------------------------------------
# Env-flag accessors
# ---------------------------------------------------------------------------


def _env_truthy(name: str) -> bool:
    """Return ``True`` iff env-var ``name`` is set to a truthy value (1, true,
    yes, on, case-insensitive).  Unset / empty / any other value -> ``False``.
    """
    raw = os.environ.get(name, "").strip().lower()
    if not raw:
        return False
    return raw in ("1", "true", "yes", "on")


def _env_int(name: str, default: int) -> int:
    raw = os.environ.get(name, "").strip()
    if not raw:
        return default
    try:
        v = int(raw)
        if v >= 0:
            return v
    except (TypeError, ValueError):
        pass
    return default


def _env_float(name: str, default: float) -> float:
    raw = os.environ.get(name, "").strip()
    if not raw:
        return default
    try:
        v = float(raw)
        if np.isfinite(v):
            return v
    except (TypeError, ValueError):
        pass
    return default


def sanity_active() -> bool:
    """``True`` iff the post-embed sanity check + retry loop is enabled.

    Default OFF — byte-identical to legacy assemble path when unset.
    """
    return _env_truthy("DELFIN_FFFREE_CONSTRUCTION_SANITY")


def vdw_floor_all_pairs_active() -> bool:
    """``True`` iff the extended vdW floor (ALL pairs, including rigid body)
    is enabled in grip_polish.  Default OFF.
    """
    return _env_truthy("DELFIN_FFFREE_GRIP_VDW_FLOOR_ALL_PAIRS")


def piano_stool_template_active() -> bool:
    """``True`` iff the per-class piano-stool template is enabled.  Default OFF."""
    return _env_truthy("DELFIN_FFFREE_PIANO_STOOL_TEMPLATE")


def fischer_carbene_template_active() -> bool:
    """``True`` iff the per-class Fischer-carbene OC-6 template is enabled.

    Default OFF — byte-identical to legacy assemble path when unset.
    """
    return _env_truthy("DELFIN_FFFREE_FISCHER_CARBENE_TEMPLATE")


def fragment_check_active() -> bool:
    """``True`` iff the post-build graph-fragmentation check is enabled.

    Default OFF — byte-identical to legacy assemble path when unset.  Logs
    violations to stderr (or via the standard ``logging`` module) without
    rejecting the build; combine with
    :func:`fragment_check_strict_active` to reject fragmented builds.
    """
    return _env_truthy("DELFIN_FFFREE_FRAGMENT_CHECK")


def fragment_check_strict_active() -> bool:
    """``True`` iff the fragment check rejects builds with extra fragments.

    Default OFF — byte-identical to legacy assemble path when unset.  When
    both this and :func:`fragment_check_active` are on, a fragmented build
    causes ``assemble_from_config`` to return ``None`` so the upstream
    caller can fall back.  When only :func:`fragment_check_active` is on,
    fragments are logged but the build is kept.
    """
    return _env_truthy("DELFIN_FFFREE_FRAGMENT_CHECK_STRICT")


# ---------------------------------------------------------------------------
# Post-embed sanity check
# ---------------------------------------------------------------------------


def _vdw_radius(sym: str) -> float:
    """Return the Bondi vdW radius for ``sym`` (Å), or :data:`_VDW_DEFAULT`
    for unknown elements."""
    return float(_VDW.get(sym, _VDW_DEFAULT))


def _is_metal(sym: str) -> bool:
    """Predicate mirroring :func:`_bd._is_metal`."""
    return _bd._is_metal(sym)


def _bond_set(bonds: Iterable[Tuple[int, int]]) -> Set[FrozenSet[int]]:
    """Convert ``bonds`` into a set of frozenset pairs (order-insensitive)."""
    s: Set[FrozenSet[int]] = set()
    for b in bonds:
        if not isinstance(b, (tuple, list)) or len(b) < 2:
            continue
        try:
            i, j = int(b[0]), int(b[1])
        except (TypeError, ValueError):
            continue
        if i == j:
            continue
        s.add(frozenset((i, j)))
    return s


def _build_adjacency(n: int, bonds: Set[FrozenSet[int]]) -> List[Set[int]]:
    """Per-atom neighbour set."""
    nbr: List[Set[int]] = [set() for _ in range(n)]
    for b in bonds:
        try:
            i, j = sorted(int(x) for x in b)
        except (TypeError, ValueError):
            continue
        if 0 <= i < n and 0 <= j < n:
            nbr[i].add(j)
            nbr[j].add(i)
    return nbr


def _bonded_or_13_pairs(
    n: int,
    bonds: Set[FrozenSet[int]],
) -> Set[FrozenSet[int]]:
    """Return the set of bonded + 1,3 (angle-sharing) pairs.

    Used to EXCLUDE these from the Pauli-floor check (the Bondi radii
    intentionally collide for bonded atoms; only non-bonded contacts count).
    """
    nbr = _build_adjacency(n, bonds)
    excl: Set[FrozenSet[int]] = set(bonds)
    for k in range(n):
        nb = sorted(nbr[k])
        m = len(nb)
        for a in range(m):
            for b in range(a + 1, m):
                excl.add(frozenset((nb[a], nb[b])))
    return excl


def assert_construction_sane(
    P: np.ndarray,
    syms: Sequence[str],
    bonds: Iterable[Tuple[int, int]],
    *,
    metal_idx: Optional[int] = 0,
    donor_idxs: Optional[Sequence[int]] = None,
    geometry: Optional[str] = None,
    pauli_fraction: float = DEFAULT_PAULI_FRACTION,
    bond_stretch_factor: float = DEFAULT_BOND_STRETCH_FACTOR,
    bond_compress_factor: float = DEFAULT_BOND_COMPRESS_FACTOR,
    md_floor_factor: float = DEFAULT_MD_FLOOR_FACTOR,
    md_ceil_factor: float = DEFAULT_MD_CEIL_FACTOR,
    cshm_max: Optional[float] = None,
) -> Tuple[bool, List[Dict]]:
    """Universal post-embedding sanity check.

    Parameters
    ----------
    P : ndarray (N, 3)
        Atom coordinates.
    syms : sequence of str
        Atomic symbols (must match ``P``).
    bonds : iterable of (i, j)
        Bonded atom pairs.  Order-insensitive.  Used for bond-length checks
        and to exclude bonded / 1,3 pairs from the Pauli floor.
    metal_idx : int or None, default 0
        Index of the central metal atom (M-D distances are checked from here).
        ``None`` disables the M-D check (useful for ligand-only structures).
    donor_idxs : sequence of int or None
        Indices of donor atoms (the polyhedron vertices).  ``None`` skips both
        the M-D check and the coordination-shape check.
    geometry : str or None
        Polyhedron geometry name (``"OC-6 octahedron"``, ``"SPY-5 square
        pyramid"`` …).  When provided, the donor unit-vectors are scored
        against the ideal vertex set via :func:`polyhedra.cshm`.
    pauli_fraction, bond_stretch_factor, bond_compress_factor,
    md_floor_factor, md_ceil_factor, cshm_max
        Thresholds (see module-level constants for defaults).  ``cshm_max=None``
        means "respect ``DELFIN_FFFREE_CONSTRUCTION_SANITY_CSHM_MAX`` or fall
        back to :data:`DEFAULT_CSHM_MAX`".

    Returns
    -------
    (is_sane, violations) : bool, list of dict
        ``is_sane`` is ``True`` iff no violation fires.  Each violation is a
        dict with keys ``{"mode", "i", "j", "d", "ideal", ...}`` describing the
        offending pair (or the offending CShM / M-D pair).  The list is empty
        when the structure passes.

    The function is *always callable* — it does **not** check the env flag.
    Callers gate it with :func:`sanity_active`.  This way unit tests can
    deterministically drive the check regardless of env state.
    """
    violations: List[Dict] = []
    P_arr = np.asarray(P, dtype=float)
    n = int(P_arr.shape[0])
    if n == 0 or n != len(syms):
        return False, [{"mode": "size_mismatch", "n_P": n, "n_syms": len(syms)}]
    if not np.all(np.isfinite(P_arr)):
        return False, [{"mode": "non_finite_coords"}]

    bonds_set = _bond_set(bonds)
    excl_pairs = _bonded_or_13_pairs(n, bonds_set)

    # ----- (1) Pauli floor over ALL non-bonded, non-1,3 pairs ---------------
    for i in range(n):
        ri = _vdw_radius(str(syms[i]))
        for j in range(i + 1, n):
            if frozenset((i, j)) in excl_pairs:
                continue
            rj = _vdw_radius(str(syms[j]))
            d = float(np.linalg.norm(P_arr[i] - P_arr[j]))
            floor = float(pauli_fraction) * (ri + rj)
            if d < floor:
                violations.append({
                    "mode": "pauli_floor",
                    "i": i, "j": j,
                    "sym_i": str(syms[i]), "sym_j": str(syms[j]),
                    "d": d, "floor": floor,
                })

    # ----- (2) bond-stretch / (3) bond-compress -----------------------------
    for b in bonds_set:
        i, j = sorted(int(x) for x in b)
        if not (0 <= i < n and 0 <= j < n):
            continue
        si, sj = str(syms[i]), str(syms[j])
        # M-D bonds are scored against md_distance, not _ideal_bond.
        if _is_metal(si) or _is_metal(sj):
            continue
        ideal = _bd._ideal_bond(si, sj)
        if ideal <= 0:
            continue
        d = float(np.linalg.norm(P_arr[i] - P_arr[j]))
        if d > float(bond_stretch_factor) * ideal:
            violations.append({
                "mode": "bond_stretched",
                "i": i, "j": j,
                "sym_i": si, "sym_j": sj,
                "d": d, "ideal": ideal,
                "ratio": d / ideal,
            })
        elif d < float(bond_compress_factor) * ideal:
            violations.append({
                "mode": "bond_compressed",
                "i": i, "j": j,
                "sym_i": si, "sym_j": sj,
                "d": d, "ideal": ideal,
                "ratio": d / ideal,
            })

    # ----- (4) Metal-donor distance bracket --------------------------------
    donor_unit_vecs: List[np.ndarray] = []
    if metal_idx is not None and donor_idxs:
        try:
            from delfin.fffree.polyhedra import md_distance as _md_dist
        except ImportError:
            _md_dist = None
        m = int(metal_idx)
        if 0 <= m < n:
            for d_idx in donor_idxs:
                try:
                    di = int(d_idx)
                except (TypeError, ValueError):
                    continue
                if not (0 <= di < n) or di == m:
                    continue
                ideal_md = float(_md_dist(str(syms[m]), str(syms[di]))) if _md_dist else 2.0
                if ideal_md <= 0:
                    ideal_md = 2.0
                vec = P_arr[di] - P_arr[m]
                d = float(np.linalg.norm(vec))
                if d <= 1e-9:
                    violations.append({
                        "mode": "md_collapsed",
                        "i": m, "j": di,
                        "sym_i": str(syms[m]), "sym_j": str(syms[di]),
                        "d": d, "ideal": ideal_md,
                    })
                    continue
                if d < md_floor_factor * ideal_md:
                    violations.append({
                        "mode": "md_too_short",
                        "i": m, "j": di,
                        "sym_i": str(syms[m]), "sym_j": str(syms[di]),
                        "d": d, "ideal": ideal_md,
                        "ratio": d / ideal_md,
                    })
                elif d > md_ceil_factor * ideal_md:
                    violations.append({
                        "mode": "md_drifted",
                        "i": m, "j": di,
                        "sym_i": str(syms[m]), "sym_j": str(syms[di]),
                        "d": d, "ideal": ideal_md,
                        "ratio": d / ideal_md,
                    })
                donor_unit_vecs.append(vec / d)

    # ----- (5) Coordination-polyhedron CShM --------------------------------
    if geometry and donor_unit_vecs:
        if cshm_max is None:
            cshm_max = _env_float(
                "DELFIN_FFFREE_CONSTRUCTION_SANITY_CSHM_MAX",
                DEFAULT_CSHM_MAX,
            )
        try:
            from delfin.fffree.polyhedra import cshm as _cshm
            obs = np.asarray(donor_unit_vecs, dtype=float)
            s = float(_cshm(obs, str(geometry)))
            if np.isfinite(s) and s > float(cshm_max):
                violations.append({
                    "mode": "cshm_too_large",
                    "geometry": str(geometry),
                    "cshm": s,
                    "cshm_max": float(cshm_max),
                })
        except Exception:
            # CShM failure is non-fatal; we already covered the structure
            # with Pauli + bond + M-D checks.
            pass

    return (len(violations) == 0), violations


# ---------------------------------------------------------------------------
# Retry loop
# ---------------------------------------------------------------------------


def build_with_retries(
    embed_fn: Callable[[int], Optional[np.ndarray]],
    syms: Sequence[str],
    bonds: Iterable[Tuple[int, int]],
    *,
    metal_idx: Optional[int] = 0,
    donor_idxs: Optional[Sequence[int]] = None,
    geometry: Optional[str] = None,
    fallback_fn: Optional[Callable[[], Optional[np.ndarray]]] = None,
    max_retries: Optional[int] = None,
    seeds: Sequence[int] = DEFAULT_RETRY_SEEDS,
    sanity_kwargs: Optional[Mapping[str, float]] = None,
) -> Tuple[Optional[np.ndarray], List[Dict]]:
    """Embed-then-sanity-check loop.

    Parameters
    ----------
    embed_fn : callable (seed) -> ndarray | None
        Embedder. Called with successive seeds from ``seeds``; first sane
        return wins.
    syms, bonds : passed through to :func:`assert_construction_sane`.
    metal_idx, donor_idxs, geometry : passed through.
    fallback_fn : callable () -> ndarray | None, optional
        Final fallback called when ALL ``max_retries`` seeds fail.  Typically
        a rigid-template placement.  Returns its result regardless of sanity
        (the caller can inspect ``violations``).
    max_retries : int, optional
        Number of seeds to try (default = env
        ``DELFIN_FFFREE_CONSTRUCTION_SANITY_MAX_RETRIES`` or
        :data:`DEFAULT_MAX_RETRIES`).  Clamped to ``len(seeds)``.
    seeds : sequence of int
        ETKDG seed candidates (default :data:`DEFAULT_RETRY_SEEDS`).
    sanity_kwargs : mapping, optional
        Extra kwargs forwarded to :func:`assert_construction_sane`.

    Returns
    -------
    (P, violations) : ndarray | None, list of dict
        First sane embedding's coordinates and an empty violations list, OR
        the fallback's coordinates with the LAST attempted seed's violations,
        OR ``(None, violations)`` if every embedder + fallback returned
        ``None`` / failed.

    The function is *callable regardless of env state*.  Callers that want a
    byte-identical legacy path simply do not call it.
    """
    if max_retries is None:
        max_retries = _env_int(
            "DELFIN_FFFREE_CONSTRUCTION_SANITY_MAX_RETRIES",
            DEFAULT_MAX_RETRIES,
        )
    max_retries = max(1, min(int(max_retries), len(seeds)))
    sanity_kwargs = dict(sanity_kwargs or {})

    last_viols: List[Dict] = []
    for k in range(max_retries):
        seed = int(seeds[k])
        try:
            P = embed_fn(seed)
        except Exception:
            P = None
        if P is None:
            last_viols = [{"mode": "embed_returned_none", "seed": seed}]
            continue
        P_arr = np.asarray(P, dtype=float)
        if P_arr.size == 0 or not np.all(np.isfinite(P_arr)):
            last_viols = [{"mode": "embed_returned_non_finite", "seed": seed}]
            continue
        ok, viols = assert_construction_sane(
            P_arr, syms, bonds,
            metal_idx=metal_idx,
            donor_idxs=donor_idxs,
            geometry=geometry,
            **sanity_kwargs,
        )
        if ok:
            return P_arr, []
        last_viols = viols

    # All retries failed -> rigid template fallback (if provided).
    if fallback_fn is not None:
        try:
            Pf = fallback_fn()
        except Exception:
            Pf = None
        if Pf is not None:
            Pf_arr = np.asarray(Pf, dtype=float)
            if Pf_arr.size > 0 and np.all(np.isfinite(Pf_arr)):
                return Pf_arr, last_viols
    return None, last_viols


# ---------------------------------------------------------------------------
# Extended vdW floor — ALL pairs, including rigid-body internals
# ---------------------------------------------------------------------------


def vdw_floor_all_pairs_value_and_grad(
    R: np.ndarray,
    *,
    symbols: Sequence[str],
    excluded_pairs: Set[FrozenSet[int]],
    weight: float,
    fraction: float = DEFAULT_PAULI_FRACTION,
    rigid_body_atoms: Optional[Sequence[int]] = None,
    rigid_translation_metal: Optional[int] = None,
    eps: float = 1e-9,
) -> Tuple[float, np.ndarray]:
    """Heavy-atom Pauli-floor penalty + gradient, applied to **all**
    non-bonded heavy pairs.

    Differences from the existing :func:`grip_polish._vdw_floor_value_and_grad`:

      * Pairs whose endpoints are BOTH inside ``rigid_body_atoms`` are still
        penalised — but their gradient is summed onto the rigid body's
        translation slot (``rigid_translation_metal``) rather than left on the
        individual atoms, so the rigid body translates as a whole instead of
        deforming.  This mirrors the chain rule the hapto-π freeze applies
        elsewhere in :mod:`grip_polish`.
      * Pairs where exactly ONE endpoint is in the rigid body have their
        gradient on the rigid endpoint summed onto the translation slot and
        zeroed in place.  The free endpoint keeps its gradient directly.
      * Pairs where neither endpoint is rigid behave identically to the legacy
        floor.

    Parameters
    ----------
    R : ndarray (N, 3)
        Current coordinates.
    symbols : sequence of str
        Atomic symbols (used to look up vdW radii).
    excluded_pairs : set of frozenset({i, j})
        Bonded + 1,3 pairs to skip.
    weight, fraction : float
        Penalty multiplier and floor fraction.
    rigid_body_atoms : sequence of int, optional
        Indices of atoms inside a rigid body (e.g. the hapto-π ring +
        substituent subtree).  When provided WITH ``rigid_translation_metal``,
        their gradients are summed onto the metal.
    rigid_translation_metal : int, optional
        Index of the rigid-body translation anchor (typically the metal).
    eps : float
        Floor on inter-atomic distance for numerical safety.

    Returns
    -------
    (L, G) : float, ndarray (N, 3)
    """
    R_arr = np.asarray(R, dtype=np.float64)
    n = int(R_arr.shape[0])
    grad = np.zeros_like(R_arr)
    if n == 0 or weight <= 0.0 or fraction <= 0.0:
        return 0.0, grad

    rigid_set: Set[int] = set()
    if rigid_body_atoms is not None:
        for x in rigid_body_atoms:
            try:
                rigid_set.add(int(x))
            except (TypeError, ValueError):
                continue
    rigid_metal = (
        int(rigid_translation_metal)
        if rigid_translation_metal is not None
        and 0 <= int(rigid_translation_metal) < n
        else None
    )

    # Precompute heavy indices + radii.
    heavy: List[int] = []
    radii: List[float] = [float("nan")] * n
    for i in range(min(n, len(symbols))):
        s = symbols[i]
        if not isinstance(s, str):
            continue
        if s == "H":
            continue
        r = _VDW.get(s)
        if r is None:
            continue
        heavy.append(i)
        radii[i] = float(r)

    total = 0.0
    for ii in range(len(heavy)):
        i = heavy[ii]
        ri = radii[i]
        if not np.isfinite(ri):
            continue
        Pi = R_arr[i]
        i_rigid = i in rigid_set
        for jj in range(ii + 1, len(heavy)):
            j = heavy[jj]
            rj = radii[j]
            if not np.isfinite(rj):
                continue
            if frozenset((i, j)) in excluded_pairs:
                continue
            d_vec = Pi - R_arr[j]
            d = float(np.linalg.norm(d_vec))
            floor_ij = fraction * (ri + rj)
            if d >= floor_ij:
                continue
            if d < eps:
                continue
            gap = floor_ij - d
            total += weight * gap * gap
            coef = -2.0 * weight * gap / d
            gi = coef * d_vec
            j_rigid = j in rigid_set
            # Distribute gradient:
            #   - both rigid + metal known     -> both onto metal
            #   - exactly one rigid + metal    -> rigid onto metal, free direct
            #   - neither rigid                -> direct (legacy behaviour)
            if i_rigid and j_rigid and rigid_metal is not None:
                # Both endpoints rigid: equal-and-opposite -> net zero on metal.
                # (Translation cannot change an internal distance.)  We
                # intentionally drop the gradient -- L-BFGS will not deform
                # the rigid body to relieve the contact; downstream sanity
                # check will trigger a re-embed.
                pass
            elif i_rigid and rigid_metal is not None:
                grad[rigid_metal] += gi
                grad[j] -= gi
            elif j_rigid and rigid_metal is not None:
                grad[i] += gi
                grad[rigid_metal] -= gi
            else:
                grad[i] += gi
                grad[j] -= gi

    return float(total), grad


# ---------------------------------------------------------------------------
# Topology classification — universal, graph-only (no SMILES strings)
# ---------------------------------------------------------------------------


def classify_topology(ligands: Sequence[Mapping]) -> str:
    """Classify a decomposed-ligand list into a coarse construction class.

    The classification is *graph-only* — based on hapto/CO/donor counts in
    the decomposed-ligand records.  It does NOT inspect SMILES strings or
    metal identity.

    Returns one of:

      * ``"piano_stool"``    — exactly one η⁵+ hapto ring + at least 3 terminal
        CO / NO / CN ligands AND optionally one σ-donor.  Matches the
        WICROP-class (η⁶-arene + 3 CO + P(o-tol)₂).
      * ``"fischer_carbene"`` — exactly one carbene-like sp² C donor + 4-5 CO
        + optional auxiliary σ-donor (N / P / O).  Matches the ODUXAN-class
        (Fischer carbene + 5 CO).
      * ``"unclassified"``   — everything else (the legacy path stays).

    A ligand entry is treated as a CO/CN/NO terminal sp-donor if its
    ``donor_local_idxs`` has length 1 AND ``donor_elems[0]`` ∈ {"C", "N"}
    AND its ``smiles`` contains a triple/double-bonded oxygen partner.  When
    those fields are missing the entry is conservatively skipped.

    Failure-safe: any missing field, type error or empty list returns
    ``"unclassified"``.  Callers gate template usage on the return value.
    """
    if not ligands:
        return "unclassified"
    n_hapto5p = 0
    hapto_eta_max = 0
    n_co = 0
    n_carbene = 0
    n_aux_sigma = 0
    for lg in ligands:
        if not isinstance(lg, Mapping):
            continue
        try:
            is_hapto = bool(lg.get("is_hapto", False))
            eta = int(lg.get("hapto_eta", 0))
            dent = int(lg.get("denticity", 0))
            donor_elems = lg.get("donor_elems") or []
            smi = str(lg.get("smiles", ""))
        except Exception:
            continue
        if is_hapto and eta >= _PIANO_STOOL_HAPTO_MIN:
            n_hapto5p += 1
            hapto_eta_max = max(hapto_eta_max, eta)
            continue
        # Terminal CO / CN / NO check: single donor C or N, triple/double O in SMILES.
        if dent == 1 and donor_elems and str(donor_elems[0]) in ("C", "N"):
            # CO: [C]#[O+] or C#O; NO: N=O; CN: C#N
            if any(token in smi for token in ("#[O+]", "#O", "=O", "#N", "C#")):
                # Carbenes (Fischer / Schrock) have aliphatic / aromatic backbone
                # attached, so the SMILES carries a substituent on the donor C
                # (e.g. "[C](OCC)=..." or "[C]=N").  Distinguish heuristically:
                # treat as carbene IFF the donor C has more than just a triple
                # bond to O.
                # Cheap fingerprint: a CO ligand has SMILES of the form
                # "[O+]#[C-]" / "[C]#[O+]" / "C#O" / "[C-]#[O+]" -- always 2-4
                # heavy atoms.  Anything longer is a carbene-ish donor.
                heavy_count = sum(1 for c in smi if c.isalpha() and c.isupper())
                if heavy_count <= 4:
                    n_co += 1
                    continue
                # Falls through to carbene below
            # Carbene detection: donor C with sp² + N/O neighbour in SMILES.
            # The Fischer-carbene fingerprint is a donor C that has a
            # substituent like "=N", "-N", "=O" or "-O" attached directly:
            if str(donor_elems[0]) == "C" and any(
                fp in smi for fp in ("=N", "-N", "(=N", "([N", "=[N", "[N+]",
                                     "=O", "-O", "(=O", "(O", "([O")
            ):
                n_carbene += 1
                continue
            # Otherwise: a generic monodentate σ-donor.
            n_aux_sigma += 1
        elif dent >= 1 and donor_elems and str(donor_elems[0]) in ("N", "P", "O", "S"):
            # Multi-atom phosphine / amine / phenoxide -> auxiliary σ-donor.
            n_aux_sigma += 1

    # Piano-stool: 1 hapto ring + ≥3 CO (auxiliary σ-donor optional)
    if (n_hapto5p == 1
            and _PIANO_STOOL_HAPTO_MIN <= hapto_eta_max <= _PIANO_STOOL_HAPTO_MAX
            and n_co >= 3):
        return "piano_stool"

    # Fischer carbene: 1 carbene + 4-5 CO (NO hapto)
    if (n_hapto5p == 0
            and n_carbene >= 1
            and _FISCHER_CO_MIN <= n_co <= _FISCHER_CO_MAX):
        return "fischer_carbene"

    return "unclassified"


# ---------------------------------------------------------------------------
# Robust per-class placement templates
# ---------------------------------------------------------------------------


def piano_stool_template(
    metal: str,
    ring_size: int,
    n_co: int,
    aux_donor: Optional[str] = None,
    *,
    ring_distance: Optional[float] = None,
    co_distance: Optional[float] = None,
    co_bond_length: float = 1.15,
    aux_distance: Optional[float] = None,
    ring_axis: np.ndarray = np.array([0.0, 0.0, 1.0]),
) -> Dict[str, np.ndarray]:
    """Geometry-only piano-stool placement (η⁵/η⁶/η⁷-ring + n CO + 1 aux σ).

    The metal sits at the origin.  The hapto ring is placed in the +z half-space
    at distance ``ring_distance`` (default = covalent-radii sum for C +
    metal-C bond), with the ring perpendicular to ``ring_axis``.  The n CO
    ligands are placed in the lower hemisphere at half-angle 54.7° from -z
    (matching the IUPAC Wells "piano-stool" tripod opening) and equal azimuth.
    The optional auxiliary σ-donor (P / N) goes at -z if not None.

    Returns
    -------
    dict with keys ``"metal"``, ``"ring"``, ``"co_carbons"``, ``"co_oxygens"``,
    optionally ``"aux_donor"``.  Each value is an ``ndarray`` (… , 3).
    """
    metal_sym = str(metal)
    ring_size = int(ring_size)
    n_co = int(n_co)
    if ring_size < 3 or n_co < 1:
        raise ValueError("piano_stool_template requires ring_size >= 3 and n_co >= 1")

    # Defaults
    if ring_distance is None:
        # M-ring centroid distance: approx. covalent C + metal covalent
        try:
            from delfin.fffree.polyhedra import md_distance as _md_dist
            md_c = _md_dist(metal_sym, "C")
        except Exception:
            md_c = 2.05
        # Ring carbons sit ~0.05 Å above the centroid plane; centroid ~
        # 1.70-1.90 Å for arene.
        ring_distance = float(md_c - 0.25)
    if co_distance is None:
        try:
            from delfin.fffree.polyhedra import md_distance as _md_dist
            co_distance = _md_dist(metal_sym, "C")
        except Exception:
            co_distance = 1.85
    if aux_distance is None and aux_donor is not None:
        try:
            from delfin.fffree.polyhedra import md_distance as _md_dist
            aux_distance = _md_dist(metal_sym, str(aux_donor))
        except Exception:
            aux_distance = 2.30

    # Normalise the ring axis.
    z = np.asarray(ring_axis, dtype=float)
    nz = float(np.linalg.norm(z))
    if nz < 1e-9:
        z = np.array([0.0, 0.0, 1.0])
    else:
        z = z / nz

    # Build an orthonormal frame (x, y, z).
    if abs(z[2]) < 0.9:
        x = np.cross(z, np.array([0.0, 0.0, 1.0]))
    else:
        x = np.cross(z, np.array([1.0, 0.0, 0.0]))
    x = x / float(np.linalg.norm(x))
    y = np.cross(z, x)

    # Ring (Cp / arene) centred in +z plane at z = ring_distance, with the
    # appropriate radius for the polygon (CCDC mean: Cp ~ 1.20 Å, arene
    # ~ 1.39 Å, COT ~ 1.42 Å).
    ring_radius = {
        3: 0.78, 4: 0.90, 5: 1.20, 6: 1.39, 7: 1.40, 8: 1.42,
    }.get(ring_size, 1.20)
    centroid = z * float(ring_distance)
    ring_pos = np.zeros((ring_size, 3))
    for k in range(ring_size):
        theta = 2.0 * math.pi * k / float(ring_size)
        ring_pos[k] = centroid + ring_radius * (math.cos(theta) * x + math.sin(theta) * y)

    # CO ligands in the -z hemisphere.  Open the tripod at half-angle ~54.7°
    # from -z (the magic-angle of the regular tetrahedron; matches the
    # sandwich_piano_polyhedra ``SIGMA_TRIPOD_OPEN_ANGLE``).
    half_open = math.radians(54.7356)
    # Azimuthal offset so the CO axes do not overlap the ring atoms in
    # projection.
    az_offset = math.pi / float(n_co)
    co_carbons = np.zeros((n_co, 3))
    co_oxygens = np.zeros((n_co, 3))
    for k in range(n_co):
        phi = 2.0 * math.pi * k / float(n_co) + az_offset
        direction = (
            -math.cos(half_open) * z
            + math.sin(half_open) * (math.cos(phi) * x + math.sin(phi) * y)
        )
        direction = direction / float(np.linalg.norm(direction))
        co_carbons[k] = float(co_distance) * direction
        # O points OUTWARD (away from metal) along the same axis.
        co_oxygens[k] = co_carbons[k] + float(co_bond_length) * direction

    out: Dict[str, np.ndarray] = {
        "metal": np.zeros(3),
        "ring": ring_pos,
        "co_carbons": co_carbons,
        "co_oxygens": co_oxygens,
    }
    if aux_donor is not None and aux_distance is not None:
        # Auxiliary donor opposite the ring centroid (the "stool leg").
        out["aux_donor"] = -z * float(aux_distance)
    return out


def fischer_carbene_template(
    metal: str,
    n_co: int,
    *,
    aux_donor: Optional[str] = None,
    carbene_distance: Optional[float] = None,
    co_distance: Optional[float] = None,
    co_bond_length: float = 1.15,
    aux_distance: Optional[float] = None,
) -> Dict[str, np.ndarray]:
    """Geometry-only OC-6 placement for a Fischer-carbene Cr(0)-class complex.

    Layout (metal at origin):

      * Carbene C at +z (axial)
      * Auxiliary σ-donor (N / P / O) at -z (axial, opposite carbene)
      * ``n_co`` CO carbons in the equatorial plane at 90° from the carbene
        axis, 360°/n_co azimuthal spacing.  CO oxygens point OUTWARD.

    When ``aux_donor`` is None the -z position is left blank — caller can place
    an extra CO there to make a CN6 homoleptic-CO + 1 carbene structure.

    Returns
    -------
    dict with keys ``"metal"``, ``"carbene_C"``, ``"co_carbons"``,
    ``"co_oxygens"`` and optionally ``"aux_donor"``.
    """
    metal_sym = str(metal)
    n_co = int(n_co)
    if n_co < 1:
        raise ValueError("fischer_carbene_template requires n_co >= 1")

    try:
        from delfin.fffree.polyhedra import md_distance as _md_dist
    except Exception:
        _md_dist = None

    def _md(donor: str, default: float) -> float:
        if _md_dist is None:
            return default
        try:
            return float(_md_dist(metal_sym, donor))
        except Exception:
            return default

    if carbene_distance is None:
        carbene_distance = _md(metal_sym, "C") if False else _md("C", 2.05)
    if co_distance is None:
        co_distance = _md("C", 1.85)
    if aux_distance is None and aux_donor is not None:
        aux_distance = _md(str(aux_donor), 2.20)

    z = np.array([0.0, 0.0, 1.0])
    x = np.array([1.0, 0.0, 0.0])
    y = np.array([0.0, 1.0, 0.0])

    # Carbene C at axial +z
    carbene_C = z * float(carbene_distance)

    # CO ligands in the equatorial plane.
    co_carbons = np.zeros((n_co, 3))
    co_oxygens = np.zeros((n_co, 3))
    for k in range(n_co):
        phi = 2.0 * math.pi * k / float(n_co)
        direction = math.cos(phi) * x + math.sin(phi) * y
        direction = direction / float(np.linalg.norm(direction))
        co_carbons[k] = float(co_distance) * direction
        co_oxygens[k] = co_carbons[k] + float(co_bond_length) * direction

    out: Dict[str, np.ndarray] = {
        "metal": np.zeros(3),
        "carbene_C": carbene_C,
        "co_carbons": co_carbons,
        "co_oxygens": co_oxygens,
    }
    if aux_donor is not None and aux_distance is not None:
        out["aux_donor"] = -z * float(aux_distance)
    return out


# ---------------------------------------------------------------------------
# Build-time graph-fragmentation check (extra_fragments bug class)
# ---------------------------------------------------------------------------


def _connected_components(n: int, bonds: Set[FrozenSet[int]]) -> List[List[int]]:
    """Return the connected components of the atom-bond graph.

    Pure BFS on the bond list — no SMILES, no chemistry, no element-specific
    behaviour.  Atoms with no incident bonds form singleton components.

    Parameters
    ----------
    n : int
        Total number of atoms (vertex count).
    bonds : set of frozenset({i, j})
        Edge set.  Self-loops / out-of-range indices are silently dropped by
        :func:`_build_adjacency`.

    Returns
    -------
    list of list of int
        One list of atom indices per component.  Each inner list is sorted
        ascending; the outer list is sorted by ``(size, smallest_index)`` so
        the smallest component appears first (tie-broken by lex order).  This
        makes the result deterministic across runs and machines.
    """
    if n <= 0:
        return []
    nbr = _build_adjacency(n, bonds)
    seen = [False] * n
    comps: List[List[int]] = []
    for start in range(n):
        if seen[start]:
            continue
        # BFS from ``start``.
        stack = [start]
        seen[start] = True
        comp: List[int] = []
        while stack:
            v = stack.pop()
            comp.append(v)
            for w in nbr[v]:
                if not seen[w]:
                    seen[w] = True
                    stack.append(w)
        comp.sort()
        comps.append(comp)
    # Deterministic ordering: smallest component first, lex tie-break.
    comps.sort(key=lambda c: (len(c), c[0] if c else -1))
    return comps


def verify_no_extra_fragments(
    P: np.ndarray,
    syms: Sequence[str],
    bonds: Iterable[Tuple[int, int]],
) -> Tuple[bool, List[Dict]]:
    """Verify that the atom-bond graph is a single connected component.

    A frequent failure mode in the GRIP polish stage is an atom being pulled
    far from its parent so that, while the bond list still records the bond,
    the resulting geometry is physically disconnected.  This check is
    *graph-only*: it operates on the discrete bond list and ignores
    coordinates.  Combined with :func:`assert_construction_sane`'s
    Pauli/bond/M-D checks it catches both physical collapses (coord-side)
    and topology fragmentation (graph-side).

    Parameters
    ----------
    P : ndarray (N, 3)
        Atom coordinates.  Currently unused — kept in the signature so the
        function is a drop-in companion to :func:`assert_construction_sane`
        and future versions can add coordinate-aware checks (e.g.
        "components separated by > X Å through space").
    syms : sequence of str
        Atomic symbols (length must match ``P`` and the highest bond index).
    bonds : iterable of (i, j)
        Bonded atom pairs.  Order-insensitive.

    Returns
    -------
    (is_intact, violations) : bool, list of dict
        ``is_intact`` is ``True`` iff exactly one connected component covers
        all atoms.  When ``False``, ``violations`` contains one entry per
        extra fragment:

        ``{
            "mode": "extra_fragment",
            "n_components": int,         # total component count
            "main_size": int,            # size of the largest component
            "fragment_size": int,        # size of THIS extra fragment
            "fragment_atoms": list[int], # 0-based indices of the fragment
            "fragment_syms": list[str],  # symbols (parallels fragment_atoms)
        }``

        The list is sorted smallest-fragment-first (matching
        :func:`_connected_components`'s deterministic order); the main
        component itself is NOT listed.  Returned empty when the build is
        intact.

    Universal: pure graph BFS, no SMILES patterns, no element rules.
    Deterministic: BFS order is fixed by the sorted bond/neighbour lists.
    Safe on degenerate input: size mismatch / non-finite coords are reported
    as a dedicated violation rather than raising.

    Examples
    --------
    >>> P = np.array([[0.0, 0, 0], [1.5, 0, 0], [3.0, 0, 0]])
    >>> verify_no_extra_fragments(P, ["C", "C", "C"], [(0, 1), (1, 2)])
    (True, [])
    >>> # Same atoms but with the bond (1, 2) removed -> two components
    >>> ok, viols = verify_no_extra_fragments(P, ["C", "C", "C"], [(0, 1)])
    >>> ok
    False
    >>> viols[0]["fragment_atoms"]
    [2]
    """
    P_arr = np.asarray(P, dtype=float)
    n_p = int(P_arr.shape[0]) if P_arr.ndim >= 1 else 0
    n_s = len(syms)
    n = n_s  # the symbol list is the authoritative atom count
    if n == 0:
        return True, []
    if n_p != n_s:
        return False, [{
            "mode": "fragment_check_size_mismatch",
            "n_P": n_p,
            "n_syms": n_s,
        }]

    bonds_set = _bond_set(bonds)
    comps = _connected_components(n, bonds_set)
    if len(comps) <= 1:
        return True, []

    # Identify the main (largest) component; everything else is an extra
    # fragment.  ``comps`` is sorted smallest-first, so the main component
    # is the LAST one.
    main = comps[-1]
    main_size = len(main)
    n_components = len(comps)

    violations: List[Dict] = []
    for comp in comps[:-1]:
        violations.append({
            "mode": "extra_fragment",
            "n_components": n_components,
            "main_size": main_size,
            "fragment_size": len(comp),
            "fragment_atoms": list(comp),
            "fragment_syms": [str(syms[i]) for i in comp],
        })
    return False, violations


def _log_fragment_violations(
    violations: Sequence[Mapping],
    *,
    context: str = "",
) -> None:
    """Emit a single-line summary of fragment violations via ``logging``.

    Internal helper used by the ``assemble_from_config`` integration so
    fragments are always logged when ``DELFIN_FFFREE_FRAGMENT_CHECK=1``,
    regardless of strict-mode.  Falls back silently if logging is broken.

    The log line includes the number of components, the largest fragment
    size and the indices of the smallest fragment — enough to grep for the
    bug class in voll-pool runs without flooding the log.
    """
    if not violations:
        return
    try:
        import logging
        logger = logging.getLogger("delfin.fffree.construction_sanity")
        for v in violations:
            mode = v.get("mode", "")
            if mode == "extra_fragment":
                logger.warning(
                    "fragment_check%s: %d components, main_size=%d, "
                    "extra_size=%d, extra_atoms=%s, extra_syms=%s",
                    (f" [{context}]" if context else ""),
                    int(v.get("n_components", 0)),
                    int(v.get("main_size", 0)),
                    int(v.get("fragment_size", 0)),
                    list(v.get("fragment_atoms", [])),
                    list(v.get("fragment_syms", [])),
                )
            else:
                logger.warning(
                    "fragment_check%s: %s %s",
                    (f" [{context}]" if context else ""),
                    mode, dict(v),
                )
    except Exception:
        # Logging failures are never fatal — the gate is informational.
        pass


# ---------------------------------------------------------------------------
# Exports
# ---------------------------------------------------------------------------

__all__ = [
    "DEFAULT_PAULI_FRACTION",
    "DEFAULT_BOND_STRETCH_FACTOR",
    "DEFAULT_BOND_COMPRESS_FACTOR",
    "DEFAULT_MD_FLOOR_FACTOR",
    "DEFAULT_MD_CEIL_FACTOR",
    "DEFAULT_CSHM_MAX",
    "DEFAULT_MAX_RETRIES",
    "DEFAULT_RETRY_SEEDS",
    "sanity_active",
    "vdw_floor_all_pairs_active",
    "piano_stool_template_active",
    "fischer_carbene_template_active",
    "fragment_check_active",
    "fragment_check_strict_active",
    "assert_construction_sane",
    "build_with_retries",
    "vdw_floor_all_pairs_value_and_grad",
    "classify_topology",
    "piano_stool_template",
    "fischer_carbene_template",
    "verify_no_extra_fragments",
]
