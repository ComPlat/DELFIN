"""sp³-H degenerate-umbrella healing pass.

Companion to :mod:`delfin.fffree.sp3_h_umbrella` — packaged separately so the
build-time hook (umbrella) and the GRIP pre-polish hook (heal) can be
toggled independently via distinct env flags.

Healing semantics
=================
* Acts as a defensive *post-hoc* pass after the structure has been built
  and partially relaxed.
* Only triggers on *degenerate* sp³ centres (H-X-H angles outside the
  band ``(95°, 125°)`` or near-collinear, |cos| > 0.9).  Non-degenerate
  centres are left untouched.
* Frozen-atom contract: metal indices and donor indices supplied by the
  caller are NEVER moved.  The healing only modifies H positions; the
  centre and its heavy neighbours are also untouched.
* Accept-if-better gate: returns the healed array only when the maximum
  H-X-H deviation strictly improves across all rewritten centres; on
  regression the original ``coords`` is returned verbatim.

Env flag
========
``DELFIN_FFFREE_SP3_H_HEAL`` (default OFF, byte-identical when unset).

Public API
==========
* :func:`heal_active() -> bool`
* :func:`detect_degenerate_sp3_h(coords, mol, ...) -> list[Sp3Center]`
* :func:`heal_degenerate_sp3_h(coords, mol, ...) -> np.ndarray`

GRIP integration
================
Designed to be invoked from
:func:`delfin.fffree.grip_polish.grip_polish` as a pre-polish stage,
parallel to ``_run_pre_polish_topology_healing`` and
``_run_pre_polish_grip_healing``.  See the wiring snippet in
``grip_polish._run_pre_polish_sp3_h_heal`` (added in the same commit as
this module).
"""
from __future__ import annotations

import os
from typing import Dict, FrozenSet, List, Optional, Sequence, Set, Tuple

import numpy as np

from .sp3_h_umbrella import (
    Sp3Center,
    Sp3Diagnostic,
    check_umbrella_geometry,
    detect_sp3_centers,
    enforce_umbrella,
)

# ---------------------------------------------------------------------------
# Env flag
# ---------------------------------------------------------------------------

_HEAL_ENV: str = "DELFIN_FFFREE_SP3_H_HEAL"


def heal_active() -> bool:
    """``True`` iff sp³-H heal is enabled (default OFF byte-identical)."""
    raw = os.environ.get(_HEAL_ENV, "").strip().lower()
    return raw in ("1", "true", "yes", "on")


# ---------------------------------------------------------------------------
# Detector
# ---------------------------------------------------------------------------


def detect_degenerate_sp3_h(
    coords: np.ndarray,
    mol=None,
    *,
    syms: Optional[Sequence[str]] = None,
    bond_topology: Optional[Sequence[Sequence[int]]] = None,
    frozen_atoms: Optional[Set[int]] = None,
) -> List[Sp3Center]:
    """Return the sp³ centres whose H umbrella is degenerate.

    "Degenerate" = any H-X-H angle outside ``(95°, 125°)`` OR any H-H
    pair with |cos(angle)| > 0.9 (near-collinear).  See
    :func:`check_umbrella_geometry`.

    ``frozen_atoms`` (optional) — centres whose ``center_idx`` belongs to
    this set are skipped entirely; centres whose ``h_indices`` overlap
    the frozen set are kept (the heal preserves frozen atoms anyway).
    """
    coords = np.asarray(coords, dtype=np.float64)
    if syms is None:
        if mol is None:
            return []
        try:
            syms = [str(a.GetSymbol()) for a in mol.GetAtoms()]
        except Exception:
            return []
    if len(syms) != coords.shape[0]:
        return []
    centres = detect_sp3_centers(coords, syms, bond_topology=bond_topology,
                                 mol=mol)
    frozen = frozen_atoms or set()
    out: List[Sp3Center] = []
    for cen in centres:
        if cen.center_idx in frozen:
            continue
        diag = check_umbrella_geometry(cen, coords)
        if diag.flag_degenerate:
            out.append(cen)
    return out


# ---------------------------------------------------------------------------
# Heal
# ---------------------------------------------------------------------------


def heal_degenerate_sp3_h(
    coords: np.ndarray,
    mol=None,
    *,
    syms: Optional[Sequence[str]] = None,
    bond_topology: Optional[Sequence[Sequence[int]]] = None,
    frozen_atoms: Optional[Set[int]] = None,
    preserve_bond_lengths: bool = True,
    accept_if_better: bool = True,
) -> Tuple[np.ndarray, Dict[str, int]]:
    """Re-orient H atoms of all degenerate sp³ centres to a Td umbrella.

    Parameters
    ----------
    coords
        ``(N, 3)`` Cartesian positions; treated as the working copy.
    mol
        Optional RDKit ``Mol``; when provided, hybridization-aware
        detection is used and ``syms`` is derived from ``mol``.
    syms
        Explicit symbol list (overrides ``mol``-derived symbols when both
        are supplied).
    bond_topology
        Optional adjacency override.
    frozen_atoms
        Atoms in this set are never moved.  Centres whose ``center_idx``
        is frozen are not processed at all.  ``h_indices`` overlapping
        ``frozen_atoms`` (rare — H is usually free) are also preserved.
    preserve_bond_lengths
        Forward to :func:`enforce_umbrella`.
    accept_if_better
        When True (default), the healed array is only returned when the
        maximum H-X-H deviation strictly decreases (over the union of
        rewritten centres).  On regression, ``coords`` is returned
        unchanged.

    Returns
    -------
    (new_coords, report)
        ``report`` keys: ``"degenerate_detected"``,
        ``"centres_healed"``, ``"hs_moved"``,
        ``"max_hh_dev_before_deg"``, ``"max_hh_dev_after_deg"``,
        ``"accepted"``.
    """
    coords = np.asarray(coords, dtype=np.float64).copy()
    report = {
        "degenerate_detected": 0,
        "centres_healed": 0,
        "hs_moved": 0,
        "max_hh_dev_before_deg": 0.0,
        "max_hh_dev_after_deg": 0.0,
        "accepted": False,
    }
    if coords.shape[0] == 0:
        return coords, report
    if syms is None and mol is not None:
        try:
            syms = [str(a.GetSymbol()) for a in mol.GetAtoms()]
        except Exception:
            return coords, report
    if syms is None or len(syms) != coords.shape[0]:
        return coords, report
    centres = detect_degenerate_sp3_h(
        coords, mol=mol, syms=syms, bond_topology=bond_topology,
        frozen_atoms=frozen_atoms,
    )
    report["degenerate_detected"] = len(centres)
    if not centres:
        return coords, report
    frozen = frozen_atoms or set()
    # Snapshot original
    original = coords.copy()
    max_before = 0.0
    for cen in centres:
        diag = check_umbrella_geometry(cen, original)
        if diag.max_hh_dev_deg > max_before:
            max_before = diag.max_hh_dev_deg
    # Apply
    working = original.copy()
    hs_moved = 0
    centres_healed = 0
    for cen in centres:
        new_arr = enforce_umbrella(
            working, cen, preserve_bond_lengths=preserve_bond_lengths,
        )
        # Honour frozen contract on H atoms (rare case)
        moved_here = 0
        for h_idx in cen.h_indices:
            if h_idx in frozen:
                new_arr[h_idx] = working[h_idx]
                continue
            if not np.allclose(new_arr[h_idx], working[h_idx], atol=1e-6):
                moved_here += 1
        if moved_here > 0:
            centres_healed += 1
            hs_moved += moved_here
        working = new_arr
    # Recompute post deviation
    max_after = 0.0
    for cen in centres:
        d2 = check_umbrella_geometry(cen, working)
        if d2.max_hh_dev_deg > max_after:
            max_after = d2.max_hh_dev_deg
    report["max_hh_dev_before_deg"] = float(max_before)
    report["max_hh_dev_after_deg"] = float(max_after)
    report["centres_healed"] = centres_healed
    report["hs_moved"] = hs_moved
    if accept_if_better:
        if max_after < max_before - 1e-6:
            report["accepted"] = True
            return working, report
        # Regression or no improvement → return original
        return original, report
    report["accepted"] = True
    return working, report


__all__ = [
    "heal_active",
    "detect_degenerate_sp3_h",
    "heal_degenerate_sp3_h",
]
