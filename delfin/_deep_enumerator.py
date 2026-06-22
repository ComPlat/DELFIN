"""Welle-5m-X: deep coordination-isomer enumerator (e6761e4 mechanism port).

Background
----------

Commit ``e6761e4`` (Mar 23 2026, ``Add extended smiles_converter with
universal hapto coordination support``) emits more named-isomer candidates
than HEAD on hetero coordination spheres.  X10-YIRQIC (Re CN=6 with five
distinct donor classes: Br, C-carbene, C-CO, P-PPh3, N-pyridyl) returns
4 distinct fingerprints in e6761e4 vs 2 in HEAD — an empirical gap of
2× on the universal trigger regime ``(CN >= 5 AND distinct_donor_classes
>= 3)``.

Forensic root cause
-------------------

``smiles_to_xyz_isomers`` HEAD adds the post-Iter-6 ``_classify_one_conf``
filter cascade that converts e6761e4's permissive strict-pass into a
stricter one.  Specifically, HEAD's strict-pass rejects (returns
``None``) on:

1. ``_md_distance_in_tolerance`` -- pre-UFF M-D gate (env-flag).
2. ``_has_atom_clash(min_dist=0.3)``                        SHARED.
3. ``_metal_donor_distances_realistic`` ``min_abs_ml=1.70`` NEW vs e6761e4.
4. ``_has_pi_ring_nonplanarity``                            SHARED.

e6761e4 strict pass (file: ``smiles_converter.py:13768-13783`` of
``e6761e4``) only rejected on (2) and (4); (3) is the discriminator.
Items (5)-(8) are *penalties* under both code paths.  The
``_metal_donor_distances_realistic`` reject is geometry-based (1.70 A
absolute floor) and catches partly-collapsed sampling conformers that
e6761e4 still keeps for fingerprint dedup -- those extra survivors are
what yields the larger isomer count.

Universal port mechanism
------------------------

When the trigger predicate fires, demote ``_metal_donor_distances_realistic``
from REJECT to PENALTY (similar to ``_has_unphysical_metal_nonbonded_contact``
in the strict-pass on HEAD) and skip ``_has_pi_ring_nonplanarity`` reject
for non-aromatic-ring complexes (e6761e4 kept it as REJECT for all classes
but pure trigger compatibility is enforced via the donor-class predicate
so this stays universal).

Trigger logic (universal, graph-only)
-------------------------------------

A SMILES qualifies for deeper enumeration iff *both*:

  1. coordination number ``CN >= 5`` at any metal centre.
  2. number of distinct donor classes ``>= 3`` (donor class = atomic
     symbol + Morgan-invariant tuple).

The Morgan invariant ensures chemically-distinct nitrogens (pyridyl vs
amine vs imine) count separately, matching ``_donor_type_map`` semantics
already used downstream.

Both conditions are pure RDKit graph queries -- no SMILES substring,
no refcode pattern, no element allow-list.  Per
``feedback_universal_fundamental_doctrine.md``.

Module shape
------------

Pure RDKit graph + ``os.environ``; no numpy / openbabel imports.
Bit-exact default-OFF behaviour preserved by ``is_deep_enum_enabled``
returning ``False`` when the env-flag is unset.

Public API
~~~~~~~~~~

* ``is_deep_enum_enabled() -> bool`` -- env-flag check.
* ``count_distinct_donor_classes(mol, metal_idx, dtype_map) -> int``.
* ``should_apply_deep_enum(mol, dtype_map) -> bool`` -- top-level
  predicate combining both trigger conditions across all metal centres.

Env-flags
~~~~~~~~~

``DELFIN_5M_X_DEEP_ENUM`` (default ``0``) -- master gate.  When set to
``1`` the universal trigger is consulted; when set to ``2`` the deep
enumeration is forced for every SMILES regardless of trigger
(diagnostic / smoke).

``DELFIN_5M_X_DEEP_ENUM_MIN_CN`` (default ``5``) -- override the CN
threshold.

``DELFIN_5M_X_DEEP_ENUM_MIN_CLASSES`` (default ``3``) -- override the
donor-class threshold.

Tests
~~~~~

Co-located in ``tests/test_welle5m_x_deep_enum.py``:

* ``test_default_off_byte_identical`` -- gate returns False when env-flag unset.
* ``test_trigger_yirqic`` -- YIRQIC normalised SMILES triggers deep enum.
* ``test_trigger_simple_cisplatin_skips`` -- cisplatin (CN=4, 2 donors) does NOT trigger.
* ``test_donor_class_counting`` -- two-N-environment ligand counts as
  2 distinct classes via Morgan invariant.

References
----------

* ``feedback_universal_fundamental_doctrine.md``
* ``feedback_named_isomer_coverage.md``
* ``project_acpicf_planarity_gap.md`` (motivation: more isomers needed
  for mixed-donor CN5+)
* ``feedback_iter12_baustein3_success.md`` -- architecture pattern
  (helper file + single insertion + env-flag gating)
"""

from __future__ import annotations

import os
from typing import Dict, Optional, Tuple

try:
    from rdkit import Chem  # type: ignore  # noqa: F401
    _RDKIT_AVAILABLE = True
except Exception:  # pragma: no cover
    _RDKIT_AVAILABLE = False


# Metal-set mirror -- kept in sync with ``smiles_converter._METAL_SET``.
# Duplicated rather than imported so this module stays import-cycle-safe.
_METAL_SET_FALLBACK = {
    "Li", "Be", "Na", "Mg", "Al", "K", "Ca", "Sc", "Ti", "V", "Cr", "Mn",
    "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge", "Rb", "Sr", "Y", "Zr",
    "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn", "Sb",
    "Cs", "Ba", "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb",
    "Dy", "Ho", "Er", "Tm", "Yb", "Lu", "Hf", "Ta", "W", "Re", "Os",
    "Ir", "Pt", "Au", "Hg", "Tl", "Pb", "Bi", "Po", "Fr", "Ra", "Ac",
    "Th", "Pa", "U", "Np", "Pu",
}


def _env_int(name: str, default: int) -> int:
    """Read an integer env-var with safe fallback."""
    try:
        raw = os.environ.get(name)
        if raw is None or raw == "":
            return int(default)
        return int(raw)
    except (TypeError, ValueError):
        return int(default)


def is_deep_enum_enabled() -> bool:
    """Return True iff the master gate ``DELFIN_5M_X_DEEP_ENUM`` is non-zero.

    Default ``0`` -- bit-exact pre-patch behaviour.  Operator override:
    set to ``1`` for trigger-conditional activation, ``2`` for forced
    activation on every SMILES (diagnostic only).
    """
    return _env_int("DELFIN_5M_X_DEEP_ENUM", 0) != 0


def _force_deep_enum() -> bool:
    """Return True iff the master gate is set to the diagnostic ``2`` value."""
    return _env_int("DELFIN_5M_X_DEEP_ENUM", 0) == 2


def count_distinct_donor_classes(
    mol,
    metal_idx: int,
    dtype_map: Optional[Dict[int, Tuple]] = None,
) -> int:
    """Return the number of distinct donor classes coordinated to ``metal_idx``.

    A donor class is the tuple ``(element, morgan_invariant)`` from the
    caller-supplied ``dtype_map`` -- which mirrors the same Morgan-enriched
    typing used by ``_compute_coordination_fingerprint`` so the trigger
    aligns with downstream dedup keys.

    When ``dtype_map`` is omitted, falls back to element-only classes
    (the conservative direction: fewer classes -> trigger fires LESS
    often, never MORE).

    Universal: queries the RDKit graph only; no element allow-list,
    no SMILES inspection.
    """
    if not _RDKIT_AVAILABLE or mol is None:
        return 0
    try:
        atom = mol.GetAtomWithIdx(int(metal_idx))
    except Exception:
        return 0
    if atom.GetSymbol() not in _METAL_SET_FALLBACK:
        return 0

    classes = set()
    for nbr in atom.GetNeighbors():
        if nbr.GetSymbol() in _METAL_SET_FALLBACK:
            # Skip M-M donors -- they are not coordination-isomer drivers.
            continue
        if dtype_map is not None:
            key = dtype_map.get(nbr.GetIdx())
            if key is None:
                key = (nbr.GetSymbol(), frozenset())
        else:
            key = (nbr.GetSymbol(), frozenset())
        # Normalise to a hashable representation.  ``frozenset`` is
        # already hashable; the (sym, frozenset) tuple ditto.
        try:
            classes.add(key)
        except TypeError:
            classes.add((str(key[0]), tuple(sorted(key[1])) if len(key) > 1 else ()))
    return len(classes)


def coordination_number(mol, metal_idx: int) -> int:
    """Return the coordination number (count of non-metal neighbours)."""
    if not _RDKIT_AVAILABLE or mol is None:
        return 0
    try:
        atom = mol.GetAtomWithIdx(int(metal_idx))
    except Exception:
        return 0
    if atom.GetSymbol() not in _METAL_SET_FALLBACK:
        return 0
    return sum(
        1 for nbr in atom.GetNeighbors()
        if nbr.GetSymbol() not in _METAL_SET_FALLBACK
    )


def should_apply_deep_enum(
    mol,
    dtype_map: Optional[Dict[int, Tuple]] = None,
) -> bool:
    """Top-level predicate: should the deep enumerator activate for ``mol``?

    Returns True iff
    ``DELFIN_5M_X_DEEP_ENUM`` is non-zero AND
    at least one metal centre satisfies
    ``CN >= DELFIN_5M_X_DEEP_ENUM_MIN_CN`` and
    ``distinct_donor_classes >= DELFIN_5M_X_DEEP_ENUM_MIN_CLASSES``.

    When the gate is set to the diagnostic ``2`` value the trigger
    bypasses the predicate and returns True unconditionally (intended
    for full-pool smoke runs to verify universal application).

    Universal: graph features only; no SMILES inspection, no
    element / refcode allow-list.
    """
    if not is_deep_enum_enabled():
        return False
    if not _RDKIT_AVAILABLE or mol is None:
        return False
    if _force_deep_enum():
        return True
    min_cn = max(2, _env_int("DELFIN_5M_X_DEEP_ENUM_MIN_CN", 5))
    min_classes = max(2, _env_int("DELFIN_5M_X_DEEP_ENUM_MIN_CLASSES", 3))
    try:
        for atom in mol.GetAtoms():
            if atom.GetSymbol() not in _METAL_SET_FALLBACK:
                continue
            m_idx = atom.GetIdx()
            cn = coordination_number(mol, m_idx)
            if cn < min_cn:
                continue
            n_classes = count_distinct_donor_classes(mol, m_idx, dtype_map)
            if n_classes >= min_classes:
                return True
    except Exception:
        return False
    return False


def deep_enum_relax_md_realistic_reject(
    mol,
    dtype_map: Optional[Dict[int, Tuple]] = None,
) -> bool:
    """Return True iff ``_metal_donor_distances_realistic`` should be
    DEMOTED to a penalty (instead of REJECT) for ``mol``'s strict pass.

    This is the e6761e4-style permissive strict-pass behaviour: keep
    sampling conformers that fall slightly below the 1.70 A absolute
    M-D floor so their fingerprints can still contribute to the
    coordination-isomer dedup pool.

    Strict-equivalent to ``should_apply_deep_enum`` -- exposed as a
    separate name so future per-feature gating can split the levers.
    """
    return should_apply_deep_enum(mol, dtype_map=dtype_map)


def deep_enum_pi_planarity_softgate(
    mol,
    dtype_map: Optional[Dict[int, Tuple]] = None,
) -> bool:
    """Return True iff ``_has_pi_ring_nonplanarity`` should be DEMOTED
    to a penalty for ``mol``'s strict pass.

    This is a second-level deep-enum knob: e6761e4 still rejected on
    pi-ring nonplanarity, so by default this returns False even when
    the deep-enum master gate is active.  Operator override:
    ``DELFIN_5M_X_DEEP_ENUM_RELAX_PI=1`` enables the pi-ring softgate
    additionally.
    """
    if not should_apply_deep_enum(mol, dtype_map=dtype_map):
        return False
    return _env_int("DELFIN_5M_X_DEEP_ENUM_RELAX_PI", 0) != 0


def deep_enum_unlabeled_rmsd_threshold(
    mol,
    dtype_map: Optional[Dict[int, Tuple]] = None,
    *,
    default_threshold: float = 2.5,
    deep_threshold: float = 0.8,
) -> float:
    """Return the RMSD threshold to use when both candidates carry an
    empty (unclassified) coordination label.

    Background
    ----------

    HEAD's RMSD dedup uses an aggressive 2.5 A threshold whenever two
    isomers share the same base label.  Empty string == empty string
    counts as "same label", so two geometrically distinct but
    unclassified candidates get merged.  On hetero CN >= 5 complexes
    where ``_classify_isomer_label`` returns ``''`` for every candidate
    (no recognised MA*B*C* pattern), this collapses what should be
    multiple distinct isomers into one.

    Mechanism
    ---------

    When the deep-enum trigger fires, return ``deep_threshold`` (0.8 A,
    the same value already used for cross-label dedup) so unlabeled
    candidates are merged only when geometrically near-identical.
    Labeled candidates keep the legacy 2.5 A behaviour -- the caller
    consults this helper only when both labels are empty.

    Returns
    -------

    ``deep_threshold`` iff ``should_apply_deep_enum`` returns True;
    otherwise ``default_threshold`` (bit-exact pre-patch).
    """
    if should_apply_deep_enum(mol, dtype_map=dtype_map):
        return float(deep_threshold)
    return float(default_threshold)


__all__ = [
    "is_deep_enum_enabled",
    "count_distinct_donor_classes",
    "coordination_number",
    "should_apply_deep_enum",
    "deep_enum_relax_md_realistic_reject",
    "deep_enum_pi_planarity_softgate",
    "deep_enum_unlabeled_rmsd_threshold",
]
