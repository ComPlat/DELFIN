"""delfin.fffree.inter_ligand_clash_gate — pre-GRIP-polish inter-ligand
clash filter for combinatorial rotamer enumeration.

Spec: ``project_grip_combinatorial_clash_gate_2026_06_02.md``.

Problem
-------
Naive combinatorial assembly of per-ligand rotamer pools blows up as
``(states_per_bond)**(n_rotors * n_ligands)`` -- 43M attempts for 4 ligands
with 3 rotors each at 6 states.  ~99 % of those configs have inter-ligand
clashes and are unusable; running them all through GRIP polish wastes
compute.

Solution
--------
A lightweight pre-polish "clash gate" that:

  1. counts inter-ligand vdW-floor violations using the ligand-subgraph
     map from :func:`grip_ensemble.identify_ligand_subgraphs`;
  2. rejects any candidate with ``count > clash_threshold`` (default 5);
  3. caps the total number of yielded candidates at ``max_total`` (default
     200);
  4. yields valid candidates in deterministic (input-product) order.

The gate is purely additive: when disabled, the caller iterates the raw
``itertools.product`` and produces the same set of candidates as before.

Universal: depends only on the molecular graph + Cartesian coordinates +
vdW radii table.  Deterministic: stable iteration order; quick distance
checks instead of full-CShM.  No SMILES patterns.

Env-flags (all default OFF -> byte-identical to HEAD when unset):

  - ``DELFIN_FFFREE_PRE_POLISH_CLASH_GATE``      -> "1" enables the gate
                                                    (default "1" only when
                                                    ``ENUMERATE_ROTAMERS=1``,
                                                    otherwise honored as-is).
  - ``DELFIN_FFFREE_PRE_POLISH_CLASH_MAX``       -> int, max clashes
                                                    allowed (default 5).
  - ``DELFIN_FFFREE_PRE_POLISH_CLASH_THRESHOLD`` -> float, vdW fraction
                                                    (default 0.85).
  - ``DELFIN_FFFREE_ENUMERATION_MAX_TOTAL``      -> int, hard cap on
                                                    yielded candidates
                                                    (default 200).
"""
from __future__ import annotations

import os
from typing import Callable, Iterable, Iterator, List, Optional, Sequence, Tuple, Set

import numpy as np


# ------------------------------------------------------------------
# Defaults
# ------------------------------------------------------------------

DEFAULT_CLASH_VDW_FRACTION = 0.85   # match grip_constraints / count_inter_ligand_clashes
DEFAULT_CLASH_THRESHOLD = 5          # max acceptable inter-ligand clash count
DEFAULT_MAX_TOTAL = 200              # hard cap on yielded candidates


# ------------------------------------------------------------------
# Env helpers
# ------------------------------------------------------------------


def gate_enabled(enumerate_rotamers_on: bool = False) -> bool:
    """Return ``True`` when the pre-polish clash gate is active.

    Precedence:

      1. ``DELFIN_FFFREE_PRE_POLISH_CLASH_GATE`` explicit "0"/"1".
      2. Auto-on when the rotamer-enumeration env-flag is set (the gate's
         primary use case is to bound the combinatorial product).
      3. Default OFF otherwise (byte-identical to HEAD).
    """
    raw = os.environ.get("DELFIN_FFFREE_PRE_POLISH_CLASH_GATE")
    if raw is not None:
        return str(raw).strip() in ("1", "true", "True", "on")
    return bool(enumerate_rotamers_on)


def _env_int(name: str, default: int) -> int:
    try:
        return int(os.environ.get(name, str(default)))
    except (TypeError, ValueError):
        return default


def _env_float(name: str, default: float) -> float:
    try:
        return float(os.environ.get(name, str(default)))
    except (TypeError, ValueError):
        return default


def env_clash_threshold(default: int = DEFAULT_CLASH_THRESHOLD) -> int:
    return _env_int("DELFIN_FFFREE_PRE_POLISH_CLASH_MAX", default)


def env_clash_vdw_fraction(default: float = DEFAULT_CLASH_VDW_FRACTION) -> float:
    return _env_float("DELFIN_FFFREE_PRE_POLISH_CLASH_THRESHOLD", default)


def env_max_total(default: int = DEFAULT_MAX_TOTAL) -> int:
    return _env_int("DELFIN_FFFREE_ENUMERATION_MAX_TOTAL", default)


# ------------------------------------------------------------------
# Quick distance-based clash counter (subgraph-aware)
# ------------------------------------------------------------------


def count_inter_ligand_clashes_quick(
    P: np.ndarray,
    syms: Sequence[str],
    ligand_subgraphs: Sequence[Iterable[int]],
    metal_idx: Optional[int] = None,
    threshold: float = DEFAULT_CLASH_VDW_FRACTION,
    vdw_radii: Optional[dict] = None,
) -> int:
    """Count atom pairs ``(i, j)`` from DIFFERENT ligand subgraphs whose
    distance is below ``threshold * (r_vdW_i + r_vdW_j)``.

    This is a stripped-down sibling of
    :func:`grip_ensemble.count_inter_ligand_clashes` that does NOT require
    an RDKit ``mol`` -- it consumes pre-computed subgraph membership +
    raw coordinates + element symbols.  Optimised for tight inner loops
    (no bond-graph lookup, no per-atom mol inspection).

    Parameters
    ----------
    P : (N, 3) ndarray
        Cartesian coordinates of the assembled complex.
    syms : sequence of str
        Element symbol per atom, aligned with ``P``.
    ligand_subgraphs : sequence of iterable[int]
        One iterable per ligand giving the atom indices belonging to that
        ligand.  Atoms NOT covered by any subgraph (e.g. the metal centre)
        are ignored.
    metal_idx : int, optional
        Atom index of the metal centre; ignored when counting clashes.
    threshold : float, default 0.85
        Multiplier of the vdW radius sum below which a pair counts as a
        clash.  Matches the grip_ensemble convention.
    vdw_radii : dict, optional
        Override the default vdW-radius table (``{symbol: r_Å}``).  When
        ``None`` a small canonical table is used.

    Returns
    -------
    int
        Number of inter-ligand clashing pairs.

    Universal: graph-agnostic; identical inputs -> identical output.
    """
    P = np.asarray(P, dtype=float)
    if P.ndim == 1:
        if P.size % 3 != 0:
            raise ValueError("P must reshape to (N, 3)")
        P = P.reshape(-1, 3)
    n = P.shape[0]
    if vdw_radii is None:
        vdw_radii = _DEFAULT_VDW_RADII
    metal_set: Set[int] = {int(metal_idx)} if metal_idx is not None else set()

    # Build per-atom subgraph id (atoms in no subgraph and the metal -> -1).
    sg_of: List[int] = [-1] * n
    for k, comp in enumerate(ligand_subgraphs):
        for a in comp:
            ai = int(a)
            if 0 <= ai < n and ai not in metal_set:
                sg_of[ai] = k

    # Pre-extract radii (NaN -> skip).
    radii = np.full(n, np.nan, dtype=float)
    for i, s in enumerate(syms):
        if i in metal_set:
            continue
        r = vdw_radii.get(str(s))
        if r is not None:
            radii[i] = float(r)

    # H-inclusion check: when DELFIN_FFFREE_HH_CLASH_INCLUDE=0 (default),
    # include H exactly as before -- H IS already in _DEFAULT_VDW_RADII
    # (r_H=1.20) so inter-ligand H-H pairs already contribute.  This
    # branch is retained for explicitness and future H-only-skip toggling.
    count = 0
    thr = float(threshold)
    for i in range(n):
        si = sg_of[i]
        if si < 0:
            continue
        ri = radii[i]
        if not np.isfinite(ri):
            continue
        Pi = P[i]
        for j in range(i + 1, n):
            sj = sg_of[j]
            if sj < 0 or sj == si:
                continue
            rj = radii[j]
            if not np.isfinite(rj):
                continue
            d = float(np.linalg.norm(Pi - P[j]))
            if d < thr * (ri + rj):
                count += 1
    return count


def count_inter_ligand_hh_clashes_quick(
    P: np.ndarray,
    syms: Sequence[str],
    ligand_subgraphs: Sequence[Iterable[int]],
    metal_idx: Optional[int] = None,
    threshold: float = DEFAULT_CLASH_VDW_FRACTION,
) -> int:
    """Count INTER-ligand H-H pairs below ``threshold × (r_H + r_H)``.

    Companion helper to :func:`count_inter_ligand_clashes_quick` used
    by the post-build forensics to attribute "how many of the inter-
    ligand clashes are H-H" without re-running the full pair scan.

    Returns 0 when the env-flag ``DELFIN_FFFREE_HH_CLASH_INCLUDE`` is
    unset (default OFF -> byte-identical with HEAD callers that don't
    import this function).
    """
    if os.environ.get("DELFIN_FFFREE_HH_CLASH_INCLUDE", "0").strip().lower() not in ("1", "true", "on", "yes"):
        return 0
    P = np.asarray(P, dtype=float)
    if P.ndim == 1:
        if P.size % 3 != 0:
            raise ValueError("P must reshape to (N, 3)")
        P = P.reshape(-1, 3)
    n = P.shape[0]
    metal_set: Set[int] = {int(metal_idx)} if metal_idx is not None else set()

    sg_of: List[int] = [-1] * n
    for k, comp in enumerate(ligand_subgraphs):
        for a in comp:
            ai = int(a)
            if 0 <= ai < n and ai not in metal_set:
                sg_of[ai] = k

    r_H = 1.20
    d_min = float(threshold) * (r_H + r_H)
    count = 0
    for i in range(n):
        if sg_of[i] < 0 or str(syms[i]) != "H":
            continue
        Pi = P[i]
        for j in range(i + 1, n):
            if sg_of[j] < 0 or sg_of[j] == sg_of[i]:
                continue
            if str(syms[j]) != "H":
                continue
            d = float(np.linalg.norm(Pi - P[j]))
            if d < d_min:
                count += 1
    return count


# Minimal vdW table copied from grip_ensemble.DEFAULT_VDW_RADII (Bondi-style).
# Kept local so this module has no hard dep on grip_ensemble at import time.
_DEFAULT_VDW_RADII = {
    "H": 1.20, "He": 1.40,
    "Li": 1.82, "Be": 1.53, "B": 1.92, "C": 1.70, "N": 1.55, "O": 1.52,
    "F": 1.47, "Ne": 1.54,
    "Na": 2.27, "Mg": 1.73, "Al": 1.84, "Si": 2.10, "P": 1.80, "S": 1.80,
    "Cl": 1.75, "Ar": 1.88,
    "K": 2.75, "Ca": 2.31,
    "Sc": 2.11, "Ti": 2.00, "V": 2.00, "Cr": 2.00, "Mn": 2.00, "Fe": 2.00,
    "Co": 2.00, "Ni": 1.63, "Cu": 1.40, "Zn": 1.39, "Ga": 1.87, "Ge": 2.11,
    "As": 1.85, "Se": 1.90, "Br": 1.85, "Kr": 2.02,
    "Rb": 3.03, "Sr": 2.49,
    "Y": 2.10, "Zr": 2.05, "Nb": 2.00, "Mo": 2.00, "Tc": 2.00, "Ru": 2.00,
    "Rh": 2.00, "Pd": 1.63, "Ag": 1.72, "Cd": 1.58, "In": 1.93, "Sn": 2.17,
    "Sb": 2.06, "Te": 2.06, "I": 1.98, "Xe": 2.16,
    "Cs": 3.43, "Ba": 2.68,
    "La": 2.40, "Ce": 2.35, "Pr": 2.39, "Nd": 2.29, "Pm": 2.36, "Sm": 2.29,
    "Eu": 2.33, "Gd": 2.37, "Tb": 2.21, "Dy": 2.29, "Ho": 2.16, "Er": 2.35,
    "Tm": 2.27, "Yb": 2.42, "Lu": 2.21,
    "Hf": 2.05, "Ta": 2.00, "W": 2.00, "Re": 2.00, "Os": 2.00, "Ir": 2.00,
    "Pt": 1.75, "Au": 1.66, "Hg": 1.55, "Tl": 1.96, "Pb": 2.02, "Bi": 2.07,
    "Po": 1.97, "At": 2.02, "Rn": 2.20,
}


# ------------------------------------------------------------------
# Generator wrapper: combinatorial assembly with early reject
# ------------------------------------------------------------------


def enumerate_with_clash_gate(
    combos_iter: Iterable,
    assembler: Callable,
    *,
    max_total: Optional[int] = None,
    clash_threshold: Optional[int] = None,
    vdw_fraction: Optional[float] = None,
) -> Iterator:
    """Iterate ``combos_iter`` of assembly inputs, build each via
    ``assembler(combo)`` and yield only those whose inter-ligand clash
    count is ``<= clash_threshold``.

    Parameters
    ----------
    combos_iter : iterable
        Yields opaque ``combo`` values that ``assembler`` accepts.  Typical
        usage: ``itertools.product(*rotamer_lists_per_ligand)``.
    assembler : callable
        ``assembler(combo) -> (P, syms, ligand_subgraphs, metal_idx)`` or
        ``None`` on failure.  ``metal_idx`` may be ``None``.
    max_total : int, optional
        Cap on yielded combinations.  Default from env (200).
    clash_threshold : int, optional
        Max inter-ligand clash count to accept.  Default from env (5).
    vdw_fraction : float, optional
        vdW-floor multiplier passed to
        :func:`count_inter_ligand_clashes_quick`.  Default from env (0.85).

    Yields
    ------
    The same tuple that ``assembler`` returns for each accepted candidate.

    Determinism: order matches ``combos_iter``; the gate only filters,
    never reorders.
    """
    cap = int(max_total) if max_total is not None else env_max_total()
    thr_clash = int(clash_threshold) if clash_threshold is not None else env_clash_threshold()
    thr_vdw = float(vdw_fraction) if vdw_fraction is not None else env_clash_vdw_fraction()
    n_yield = 0
    for combo in combos_iter:
        if n_yield >= cap:
            return
        try:
            assembled = assembler(combo)
        except Exception:
            continue
        if assembled is None:
            continue
        try:
            P, syms, subgraphs, metal_idx = assembled
        except Exception:
            continue
        try:
            n_clash = count_inter_ligand_clashes_quick(
                P, syms, subgraphs, metal_idx=metal_idx,
                threshold=thr_vdw,
            )
        except Exception:
            continue
        if n_clash > thr_clash:
            continue
        yield assembled
        n_yield += 1


# ------------------------------------------------------------------
# Diagnostics helper: count without yielding (for forensics)
# ------------------------------------------------------------------


def gate_stats(
    combos_iter: Iterable,
    assembler: Callable,
    *,
    max_attempts: int = 1000,
    clash_threshold: Optional[int] = None,
    vdw_fraction: Optional[float] = None,
) -> Tuple[int, int, int]:
    """Return ``(attempts, assembled_ok, kept)`` -- pure counting for
    forensic A/B comparisons.  Does NOT yield candidates.
    """
    thr_clash = int(clash_threshold) if clash_threshold is not None else env_clash_threshold()
    thr_vdw = float(vdw_fraction) if vdw_fraction is not None else env_clash_vdw_fraction()
    attempts = 0
    assembled_ok = 0
    kept = 0
    for combo in combos_iter:
        if attempts >= max_attempts:
            break
        attempts += 1
        try:
            assembled = assembler(combo)
        except Exception:
            continue
        if assembled is None:
            continue
        assembled_ok += 1
        try:
            P, syms, subgraphs, metal_idx = assembled
            n_clash = count_inter_ligand_clashes_quick(
                P, syms, subgraphs, metal_idx=metal_idx,
                threshold=thr_vdw,
            )
            if n_clash <= thr_clash:
                kept += 1
        except Exception:
            continue
    return attempts, assembled_ok, kept
