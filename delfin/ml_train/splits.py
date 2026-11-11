"""Train/validation/test splitting and a leakage check for delfin/ml_train.

Splits are expressed as a mapping ``split_name -> list[record index]`` over
the full record list, so every split is provably disjoint and covers the
whole dataset with no index in two splits. Four strategies are provided:

* :func:`random_split` -- plain seeded random assignment per record.
* :func:`scaffold_split` -- Bemis-Murcko scaffold grouping (rdkit); a whole
  scaffold goes to one split, which is the standard way to avoid
  near-duplicate-bearing scaffolds leaking between train and test.
* :func:`stratified_split` -- random split *within* target-value bins, so each
  split keeps the target distribution's bins represented.
* :func:`grouped_split` -- assign whole groups (e.g. by metal or by source)
  to one split, so no group appears in two splits.

Every split takes a ``seed`` and is fully reproducible from it; the seed is
also recorded in the split file written by :func:`save_split`.

:func:`check_leakage` verifies that no molecule (canonical InChIKey) and no
scaffold appears in two *different* splits, and reports any overlap it
finds. rdkit is imported lazily (it is a core DELFIN dependency).
"""
from __future__ import annotations

import json
import random
from typing import Any, Callable, Dict, List, Optional

from .datasets import Record, canonical_inchikey

__all__ = [
    "check_leakage",
    "grouped_split",
    "load_split",
    "random_split",
    "save_split",
    "scaffold_split",
    "stratified_split",
    "LeakageReport",
]


# ---------------------------------------------------------------------------
# Shared greedy assignment primitive
# ---------------------------------------------------------------------------


def _assign_groups_to_splits(
    groups: List[tuple], fractions: Dict[str, float], seed: int
) -> Dict[str, List[int]]:
    """Assign whole groups to splits so each split's fill tracks its fraction.

    ``groups`` is a list of ``(key, [record indices])``. Groups are visited
    in a seeded-random order; each group is placed into the split currently
    furthest below its target fraction. Disjoint by construction (each index
    belongs to one group, one split).
    """
    rng = random.Random(seed)
    names = list(fractions)
    if not names:
        return {}
    order = list(range(len(groups)))
    rng.shuffle(order)
    splits: Dict[str, List[int]] = {name: [] for name in names}
    counts = {name: 0 for name in names}

    def fill_ratio(name: str) -> float:
        target = fractions[name]
        return counts[name] / target if target > 0 else float("inf")

    for gi in order:
        _key, indices = groups[gi]
        chosen = min(names, key=fill_ratio)
        splits[chosen].extend(indices)
        counts[chosen] += len(indices)
    return splits


def _assert_fractions(fractions: Dict[str, float]) -> List[str]:
    names = list(fractions)
    if not names:
        raise ValueError("fractions must name at least one split")
    for name in fractions:
        if fractions[name] < 0:
            raise ValueError(f"split {name!r} has negative fraction {fractions[name]}")
    return names


# ---------------------------------------------------------------------------
# Individual strategies
# ---------------------------------------------------------------------------


def random_split(
    n: int, fractions: Dict[str, float], seed: int
) -> Dict[str, List[int]]:
    """Split ``n`` records (indices 0..n-1) randomly, reproducibly by seed."""
    _assert_fractions(fractions)
    groups = [(i, [i]) for i in range(n)]
    return _assign_groups_to_splits(groups, fractions, seed)


def scaffold_split(
    records: List[Record], fractions: Dict[str, float], seed: int
) -> Dict[str, List[int]]:
    """Bemis-Murcko scaffold split: each scaffold goes wholly to one split."""
    _assert_fractions(fractions)
    by_scaffold: Dict[str, List[int]] = {}
    for i, rec in enumerate(records):
        scaf = _murcko_scaffold(rec.smiles) if rec.smiles else f"smiles:{i}"
        by_scaffold.setdefault(scaf, []).append(i)
    groups = sorted(by_scaffold.items(), key=lambda kv: (kv[0]))
    return _assign_groups_to_splits(groups, fractions, seed)


def stratified_split(
    records: List[Record],
    target_name: str,
    n_bins: int,
    fractions: Dict[str, float],
    seed: int,
) -> Dict[str, List[int]]:
    """Stratified split by target-value bins: random split within each bin."""
    _assert_fractions(fractions)
    if n_bins < 1:
        raise ValueError(f"n_bins must be >= 1, got {n_bins}")

    values = [_target_value(records[i], target_name) for i in range(len(records))]
    finite = [v for v in values if v is not None]
    if not finite:
        raise ValueError(f"no record has target {target_name!r}; cannot stratify")
    vmin, vmax = min(finite), max(finite)
    span = (vmax - vmin) if vmax > vmin else 1.0

    bins: Dict[int, List[int]] = {b: [] for b in range(n_bins)}
    unbinned: List[int] = []
    for i, v in enumerate(values):
        if v is None:
            unbinned.append(i)
            continue
        b = min(int((v - vmin) / span * n_bins), n_bins - 1)
        bins[b].append(i)

    combined: Dict[str, List[int]] = {name: [] for name in fractions}
    for b in range(n_bins):
        part = _assign_groups_to_splits(
            [(i, [i]) for i in bins[b]], fractions, seed=b * 1000 + seed
        )
        for name in fractions:
            combined[name].extend(part[name])
    # unstratified leftovers: single shot
    if unbinned:
        part = _assign_groups_to_splits(
            [(i, [i]) for i in unbinned], fractions, seed=seed
        )
        for name in fractions:
            combined[name].extend(part[name])
    return combined


def grouped_split(
    records: List[Record],
    key: Callable[[Record], Any],
    fractions: Dict[str, float],
    seed: int,
) -> Dict[str, List[int]]:
    """Split by a grouping key (e.g. metal or source); a group never splits."""
    _assert_fractions(fractions)
    by_group: Dict[Any, List[int]] = {}
    for i, rec in enumerate(records):
        by_group.setdefault(key(rec), []).append(i)
    groups = sorted(by_group.items(), key=lambda kv: str(kv[0]))
    return _assign_groups_to_splits(groups, fractions, seed)


# ---------------------------------------------------------------------------
# Leakage check
# ---------------------------------------------------------------------------


class LeakageReport:
    """Result of :func:`check_leakage`: overlaps of molecules / scaffolds."""

    def __init__(self) -> None:
        self.inchikey_overlaps: List[tuple] = []  # (a, b, inchikey)
        self.scaffold_overlaps: List[tuple] = []  # (a, b, scaffold)

    @property
    def has_leakage(self) -> bool:
        return bool(self.inchikey_overlaps or self.scaffold_overlaps)


def check_leakage(
    splits: Dict[str, List[int]], records: List[Record]
) -> LeakageReport:
    """Report any molecule (InChIKey) or scaffold shared between two splits.

    ``splits`` maps split-name to record indices; ``records`` is the full
    record list. Records without a parseable SMILES fall back to their own
    identifier as InChIKey (never a false overlap). Only *different* splits
    are compared.
    """
    report = LeakageReport()
    names = list(splits)
    # molecule and scaffold present in each split
    mol_of: List[Dict[int, str]] = [{} for _ in names]
    scaf_of: List[Dict[int, str]] = [{} for _ in names]
    for si, name in enumerate(names):
        for i in splits[name]:
            rec = records[i]
            mol_of[si][i] = canonical_inchikey(rec.smiles) if rec.smiles else rec.identifier
            scaf_of[si][i] = (
                _murcko_scaffold(rec.smiles) if rec.smiles else rec.identifier
            )
    for a in range(len(names)):
        for b in range(a + 1, len(names)):
            a_mol = {v for v in mol_of[a].values() if v}
            b_mol = {v for v in mol_of[b].values() if v}
            for key in a_mol & b_mol:
                report.inchikey_overlaps.append((names[a], names[b], key))
            a_scaf = {v for v in scaf_of[a].values() if v}
            b_scaf = {v for v in scaf_of[b].values() if v}
            for key in a_scaf & b_scaf:
                report.scaffold_overlaps.append((names[a], names[b], key))
    return report


# ---------------------------------------------------------------------------
# Split-file serialisation (records the seed + method)
# ---------------------------------------------------------------------------


def save_split(
    splits: Dict[str, List[int]],
    path: str,
    seed: int,
    method: str,
    extra: Optional[Dict[str, Any]] = None,
) -> None:
    """Write a split definition with its seed and method to JSON."""
    payload = {
        "seed": seed,
        "method": method,
        "splits": {name: list(idx) for name, idx in splits.items()},
    }
    if extra:
        payload["extra"] = extra
    with open(path, "w", encoding="utf-8") as fh:
        json.dump(payload, fh, indent=2, sort_keys=True)


def load_split(path: str) -> Dict[str, Any]:
    """Read a split file written by :func:`save_split`."""
    with open(path, "r", encoding="utf-8") as fh:
        return json.load(fh)


# ---------------------------------------------------------------------------
# private helpers
# ---------------------------------------------------------------------------


def _target_value(rec: Record, name: str) -> Optional[float]:
    for t in rec.targets:
        if t.name == name:
            return t.value
    return None


def _murcko_scaffold(smiles: str) -> str:
    from rdkit import Chem
    from rdkit.Chem.Scaffolds import MurckoScaffold

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return smiles
    scaf = MurckoScaffold.GetScaffoldForMol(mol)
    return MurckoScaffold.MurckoScaffoldSmiles(mol=scaf)
