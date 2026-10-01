"""SMILES parsing, metal and hapto-group detection, complex-class classification and the hapto approximation of the MANTA constructor.

Moved verbatim from delfin/smiles_converter.py (MANTA split, 2026-10); every
statement is the original text, only the imports below are new.
"""
from __future__ import annotations

import os
import re
from typing import Dict, List, Optional, Tuple

from delfin.common.logging import get_logger
from delfin.manta.ml_tables import (
    Chem,
    RDKIT_AVAILABLE,
    _METALS,
    _METAL_SET,
    _prefer_no_sanitize,
)

logger = get_logger("delfin.smiles_converter")


def contains_metal(smiles: str) -> bool:
    """Return True if SMILES likely contains a metal atom."""
    for metal in _METALS:
        if re.search(rf'\[{metal}[+\-\d\]@]', smiles, re.IGNORECASE):
            return True
        if re.search(rf'\[{metal}\]', smiles, re.IGNORECASE):
            return True
    return False


def _hapto_approx_enabled(flag: Optional[bool] = None) -> bool:
    """Return True when the hapto approximation/first-class path is enabled.

    Defaults to ENABLED so hapto-detected SMILES (Cp, arene-π, allyl,
    diene complexes) produce η-labelled output instead of fail-fast.
    The underlying hapto pipeline (introduced fully in `327b4e6` and
    refined through `44fce9e`) is universal and stable; the historical
    fail-fast gate was added in `7411aee` only as transitional safety
    while hapto code matured.

    Disable explicitly via DELFIN_HAPTO_APPROX=0 to revert to fail-fast
    when diagnosing hapto-builder issues.

    Returns True when:
      - explicit flag is True, OR
      - DELFIN_HAPTO_APPROX env var is truthy or unset, OR
      - DELFIN_HAPTO_FIRST_CLASS env var is truthy.
    Returns False only when DELFIN_HAPTO_APPROX is explicitly truthy-false
    (``0``/``false``/``no``/``off``) AND DELFIN_HAPTO_FIRST_CLASS is not set.
    """
    if flag is not None:
        return bool(flag)
    env = os.environ.get("DELFIN_HAPTO_APPROX", "").strip().lower()
    if env in {"0", "false", "no", "off"}:
        env_fc = os.environ.get("DELFIN_HAPTO_FIRST_CLASS", "").strip().lower()
        return env_fc in {"1", "true", "yes", "on"}
    return True


def _hapto_label_for_group(mol, c_indices: List[int]) -> str:
    """Build canonical η{n}-{ring_type} label for a hapto donor block.

    n = number of metal-bound C atoms in the contiguous block.
    ring_type heuristic:
      n == 2  → 'η2-alkene' (or 'η2-π')
      n == 3  → 'η3-allyl'
      n == 4  → 'η4-diene'
      n == 5  → 'η5-Cp' (cyclopentadienyl-like)
      n == 6  → 'η6-arene'
      n == 7  → 'η7-cycloheptatrienyl'
      n == 8  → 'η8-cyclooctatetraenyl'
    Universal — driven only by atom count of the contiguous block,
    not by SMILES patterns.
    """
    n = len(c_indices)
    name_map = {
        2: 'alkene', 3: 'allyl', 4: 'diene',
        5: 'Cp', 6: 'arene', 7: 'cht', 8: 'cot',
    }
    suffix = name_map.get(n, f'C{n}')
    return f'η{n}-{suffix}'


def _find_hapto_groups(mol) -> List[Tuple[int, List[int]]]:
    """Detect likely eta/hapto groups as contiguous metal-bound carbon sets."""
    if not RDKIT_AVAILABLE or mol is None:
        return []

    groups: List[Tuple[int, List[int]]] = []
    for atom in mol.GetAtoms():
        if atom.GetSymbol() not in _METAL_SET:
            continue
        metal_idx = atom.GetIdx()
        all_neighbors = list(atom.GetNeighbors())
        c_neighbors = [n.GetIdx() for n in all_neighbors if n.GetAtomicNum() == 6]
        if len(c_neighbors) < 2:
            continue
        # Skip cluster compounds (borane/carborane cages): if >50% of
        # metal neighbours are B atoms, the C atoms are cage vertices,
        # not part of a cyclopentadienyl-type hapto ligand.
        b_count = sum(1 for n in all_neighbors if n.GetAtomicNum() == 5)
        if b_count > len(all_neighbors) / 2:
            continue

        c_set = set(c_neighbors)
        seen: set = set()
        for start in c_neighbors:
            if start in seen:
                continue
            comp: List[int] = []
            stack = [start]
            seen.add(start)
            while stack:
                cur = stack.pop()
                comp.append(cur)
                cur_atom = mol.GetAtomWithIdx(cur)
                for nbr in cur_atom.GetNeighbors():
                    ni = nbr.GetIdx()
                    if ni not in c_set or ni in seen:
                        continue
                    seen.add(ni)
                    stack.append(ni)

            if len(comp) < 2:
                continue
            # Avoid false positives for ordinary C,C chelation: classify as hapto
            # when the contiguous donor block is ring-like or has >=3 atoms.
            ring_like = any(mol.GetAtomWithIdx(i).IsInRing() for i in comp)
            if len(comp) >= 3 or ring_like:
                groups.append((metal_idx, sorted(comp)))

    return groups


def _probe_hapto_groups_from_smiles(smiles: str) -> List[Tuple[int, List[int]]]:
    """Parse SMILES quickly and report likely hapto groups."""
    if not RDKIT_AVAILABLE or not contains_metal(smiles):
        return []

    mol, _note = mol_from_smiles_rdkit(smiles, allow_metal=True)
    if mol is None:
        try:
            p = Chem.SmilesParserParams()
            p.sanitize = False
            p.removeHs = False
            p.strictParsing = False
            mol = Chem.MolFromSmiles(smiles, p)
        except Exception:
            try:
                mol = Chem.MolFromSmiles(smiles, sanitize=False)
            except Exception:
                mol = None
    return _find_hapto_groups(mol)


# ============================================================================
# Iter-8 FPCFD: per-class function dispatch infrastructure
# ============================================================================
#
# Forensics-Driven Per-Class Champion Function Dispatch (FPCFD) is the Iter-8+
# strategy that replaces the env-flag-stacking approach of Iter-1..7.  Per
# (chemistry_class × function) cell, the best-known champion's function-body
# is forward-ported into HEAD with suffix `_<champion_sha>` and dispatched at
# call-site via `_classify_complex_class()` below.
#
# Classes (5):
#   no_metal      — pure organic / inorganic SMILES with no transition metal
#   sigma         — single metal, all donors are σ-only (Cl, P, N, O, S, ...)
#   hapto         — single metal with at least one η-coordinated ligand (Cp, arene, alkene)
#   multi_sigma   — two or more metals, all σ-only coordination
#   multi_hapto   — two or more metals, at least one metal η-coordinated
#
# Design doctrine (HANDOFF2.md sections 14, 15, 16):
#   - Champion bodies live permanently in HEAD as named functions.
#   - Per call-site: dispatcher selects exactly ONE body per SMILES.
#   - Forward-only: no revert; on regression → forward-only fix-commit.
#   - Per-class Δ-matrix as acceptance gate (no aggregate-trap).
#
# This dispatch infrastructure is *pure infrastructure* (Iter-8.0): adding the
# helper does NOT change behavior — no champion-bodies are ported yet, no
# call-sites dispatch yet.  Subsequent Iter-8.1..8.7 patches port one cell
# each and wire the dispatcher.
def _classify_complex_class(mol) -> str:
    """5-class chemistry classification based on metal count + hapto coordination.

    Returns one of:
      - 'no_metal'    — n_metals == 0
      - 'sigma'       — n_metals == 1, no hapto group
      - 'hapto'       — n_metals == 1, at least one hapto group
      - 'multi_sigma' — n_metals >= 2, no hapto-coordinated metal
      - 'multi_hapto' — n_metals >= 2, at least one hapto-coordinated metal

    Used by Iter-8+ FPCFD per-class function dispatch.  Mirrors
    `quality_framework/scripts/find_topology_loss.py:chemistry_class()` so
    that runtime-class-detection matches detector-class-detection (essential
    for per-class Δ-matrix attribution).

    Defensive: returns 'no_metal' on any exception (mol parsing failure etc).
    """
    if mol is None:
        return "no_metal"
    try:
        n_metals = sum(
            1 for a in mol.GetAtoms() if a.GetSymbol() in _METAL_SET
        )
    except Exception:
        return "no_metal"
    if n_metals == 0:
        return "no_metal"
    try:
        hapto_metal_idx = {m_idx for m_idx, _ in _find_hapto_groups(mol)}
    except Exception:
        hapto_metal_idx = set()
    n_hapto_metals = len(hapto_metal_idx)
    if n_metals == 1:
        return "hapto" if n_hapto_metals else "sigma"
    return "multi_hapto" if n_hapto_metals else "multi_sigma"


def _hapto_failfast_error(hapto_groups: List[Tuple[int, List[int]]]) -> str:
    """Build a concise user-facing error for unsupported hapto coordination."""
    if not hapto_groups:
        return (
            "Hapto (eta) coordination detected. "
            "Enable DELFIN_HAPTO_APPROX=1 for experimental approximation mode."
        )
    max_group = max(len(g[1]) for g in hapto_groups)
    return (
        "Hapto (eta) coordination detected "
        f"({len(hapto_groups)} group(s), max eta~{max_group}). "
        "Standard conversion does not support this reliably. "
        "Set DELFIN_HAPTO_APPROX=1 to enable experimental approximation mode."
    )


def _select_multihapto_anchors(
    mol,
    metal_idx: int,
    metal_groups: List[List[int]],
) -> List[int]:
    """Choose one anchor per hapto group, maximally spread around metal_idx.

    For metals with 2+ groups: pick anchors from opposite ends of the graph
    to avoid two anchors being direct neighbours.
    """
    if len(metal_groups) <= 1:
        # Single group: just pick best anchor
        grp = metal_groups[0]
        grp_set = set(grp)

        def _akey(idx: int) -> Tuple[int, int, int, int]:
            a = mol.GetAtomWithIdx(idx)
            aromatic = 1 if a.GetIsAromatic() else 0
            in_ring = 1 if a.IsInRing() else 0
            c_in_grp = sum(
                1 for n in a.GetNeighbors() if n.GetAtomicNum() == 6 and n.GetIdx() in grp_set
            )
            return (aromatic, in_ring, c_in_grp, -idx)
        return [max(grp, key=_akey)]

    # For each group, compute a "position score" = average graph distance
    # from the metal through the group members.  Pick anchors that are
    # maximally spread by assigning them greedily.
    anchors: List[int] = []
    used_anchors: set = set()
    for grp in metal_groups:
        grp_set = set(grp)

        def _anchor_score(idx: int) -> Tuple[float, int, int, int, int]:
            a = mol.GetAtomWithIdx(idx)
            aromatic = 1 if a.GetIsAromatic() else 0
            in_ring = 1 if a.IsInRing() else 0
            c_in_grp = sum(
                1 for n in a.GetNeighbors() if n.GetAtomicNum() == 6 and n.GetIdx() in grp_set
            )
            # Penalty for being adjacent to an already-chosen anchor
            adj_penalty = 0
            for n in a.GetNeighbors():
                if n.GetIdx() in used_anchors:
                    adj_penalty = -2
                    break
            return (adj_penalty, aromatic, in_ring, c_in_grp, -idx)

        anchor = max(grp, key=_anchor_score)
        anchors.append(anchor)
        used_anchors.add(anchor)

    return anchors


def _apply_hapto_approximation(
    mol,
    hapto_groups: Optional[List[Tuple[int, List[int]]]] = None,
):
    """Experimental eta->anchor approximation for embedding/optimization.

    For each contiguous metal-bound carbon block, keep one representative
    metal-carbon bond as a normal bond and convert the other M-C contacts to
    dative bonds (C->M). For metals with multiple hapto groups, anchors are
    chosen to be maximally spread to improve embedding convergence.
    """
    if not RDKIT_AVAILABLE or mol is None:
        return mol, 0

    groups = hapto_groups if hapto_groups is not None else _find_hapto_groups(mol)
    if not groups:
        return mol, 0

    # Group hapto groups by metal for coordinated anchor selection
    by_metal: Dict[int, List[List[int]]] = {}
    by_metal_order: Dict[int, List[int]] = {}  # track group indices
    for gi, (metal_idx, grp) in enumerate(groups):
        if len(grp) < 2:
            continue
        by_metal.setdefault(metal_idx, []).append(grp)
        by_metal_order.setdefault(metal_idx, []).append(gi)

    # Compute coordinated anchors per metal
    anchor_map: Dict[int, int] = {}  # group_index -> anchor_atom_idx
    for metal_idx, metal_groups in by_metal.items():
        anchors = _select_multihapto_anchors(mol, metal_idx, metal_groups)
        for i, anchor in enumerate(anchors):
            gi = by_metal_order[metal_idx][i]
            anchor_map[gi] = anchor

    rw = Chem.RWMol(mol)
    converted = 0
    for gi, (metal_idx, grp) in enumerate(groups):
        if len(grp) < 2:
            continue
        anchor = anchor_map.get(gi)
        if anchor is None:
            grp_set = set(grp)

            def _anchor_key(idx: int) -> Tuple[int, int, int, int]:
                a = rw.GetAtomWithIdx(idx)
                aromatic = 1 if a.GetIsAromatic() else 0
                in_ring = 1 if a.IsInRing() else 0
                c_in_grp = sum(
                    1 for n in a.GetNeighbors() if n.GetAtomicNum() == 6 and n.GetIdx() in grp_set
                )
                return (aromatic, in_ring, c_in_grp, -idx)
            anchor = max(grp, key=_anchor_key)

        for c_idx in grp:
            c_atom = rw.GetAtomWithIdx(c_idx)
            c_atom.SetNoImplicit(False)
            if c_idx == anchor:
                continue
            bond = rw.GetBondBetweenAtoms(metal_idx, c_idx)
            if bond is None:
                continue
            if bond.GetBondType() == Chem.BondType.DATIVE:
                continue
            rw.RemoveBond(metal_idx, c_idx)
            rw.AddBond(c_idx, metal_idx, Chem.BondType.DATIVE)
            converted += 1

    out = rw.GetMol()
    try:
        out.UpdatePropertyCache(strict=False)
    except Exception:
        pass
    return out, converted


def mol_from_smiles_rdkit(smiles: str, allow_metal: bool = False):
    """Create RDKit Mol, with relaxed sanitizing for metal complexes."""
    try:
        if not _prefer_no_sanitize(smiles):
            mol = Chem.MolFromSmiles(smiles)
            if mol is not None:
                return mol, None
        if not allow_metal:
            return None, "Failed to parse SMILES string"
        mol = Chem.MolFromSmiles(smiles, sanitize=False)
        if mol is None:
            return None, "Failed to parse SMILES (no sanitize)"
        try:
            Chem.SanitizeMol(
                mol,
                sanitizeOps=(
                    Chem.SanitizeFlags.SANITIZE_ALL
                    ^ Chem.SanitizeFlags.SANITIZE_PROPERTIES
                    ^ Chem.SanitizeFlags.SANITIZE_KEKULIZE
                ),
            )
        except Exception:
            pass
        # Update property cache to compute implicit H counts for AddHs
        # Avoid property-cache updates for neutral Ni/Co+[N] SMILES to prevent valence errors
        if not _prefer_no_sanitize(smiles):
            try:
                mol.UpdatePropertyCache(strict=False)
            except Exception:
                pass
        return mol, "partial sanitize"
    except Exception as e:
        return None, f"RDKit error: {e}"
