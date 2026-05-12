"""SMILES → ClassFeatures classifier — Phase 2 of Hybrid-Path.

Standalone classifier for routing decisions. No dependency on smiles_converter
internals (which would create a circular import).

Empirically validated against existing archive class-labels in
quality_framework/results/class_taxonomy.md (~7M classified records).

Per nature_project/15_HYBRID_PATH_FINAL.md.
"""
from __future__ import annotations

import re
from dataclasses import dataclass
from typing import FrozenSet, List, Optional


# Metal blocks (Pyykko classification)
METAL_3D: FrozenSet[str] = frozenset({
    "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn",
})
METAL_4D: FrozenSet[str] = frozenset({
    "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd",
})
METAL_5D: FrozenSet[str] = frozenset({
    "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg",
})
LANTHANIDE: FrozenSet[str] = frozenset({
    "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb",
    "Dy", "Ho", "Er", "Tm", "Yb", "Lu",
})
ACTINIDE: FrozenSet[str] = frozenset({
    "Ac", "Th", "Pa", "U", "Np", "Pu", "Am", "Cm", "Bk",
    "Cf", "Es", "Fm", "Md", "No", "Lr",
})
SBLOCK: FrozenSet[str] = frozenset({
    "Li", "Na", "K", "Rb", "Cs", "Fr",
    "Be", "Mg", "Ca", "Sr", "Ba", "Ra",
})
PBLOCK: FrozenSet[str] = frozenset({
    "Al", "Ga", "In", "Tl",
    "Ge", "Sn", "Pb", "Sb", "Bi",
})

ALL_METALS: FrozenSet[str] = (
    METAL_3D | METAL_4D | METAL_5D | LANTHANIDE | ACTINIDE | SBLOCK | PBLOCK
)


# Coord-class enumeration
COORD_SIGMA = "sigma"
COORD_HAPTO = "hapto"
COORD_MULTI_SIGMA = "multi-sigma"
COORD_MULTI_HAPTO = "multi-hapto"
COORD_NO_METAL = "no_metal"
COORD_BORANE = "borane"  # B-cluster
COORD_CLUSTER = "cluster"  # M-M-bonded ≥3 metals

# Metal-block enumeration
BLOCK_3D = "3d"
BLOCK_4D = "4d"
BLOCK_5D = "5d"
BLOCK_LN = "Ln"
BLOCK_AN = "An"
BLOCK_S = "s"
BLOCK_P = "p"
BLOCK_NONE = "no-metal"
BLOCK_UNKNOWN = "?"


@dataclass(frozen=True)
class ClassFeatures:
    """Class-features for routing decisions.

    Attributes:
        coord_class: One of {sigma, hapto, multi-sigma, multi-hapto, no_metal,
                     borane, cluster}.
        metal_block: One of {3d, 4d, 5d, Ln, An, s, p, no-metal, ?}.
        metals:      Tuple of metal symbols in molecule (empty if no_metal).
        n_metals:    Number of metal atoms in SMILES (>1 → multi-metal coord).
        has_hapto:   At least one η-coordination (heuristic via ring-bonded M).
        has_aromatic: Aromatic ring present (Kekulé or aromatic-flag).
        has_bridge:  μ-bridging ligand (heuristic via M...L...M).
    """
    coord_class: str
    metal_block: str
    metals: tuple[str, ...]
    n_metals: int
    has_hapto: bool
    has_aromatic: bool
    has_bridge: bool


# ----------------------------------------------------------------------
# Helpers
# ----------------------------------------------------------------------

# Match bracketed atoms with optional charge: [Fe-2], [NH3+], [Pt+]
_BRACKET_ATOM_RE = re.compile(r"\[([A-Z][a-z]?)")


def _metal_block_of(sym: str) -> str:
    if sym in METAL_3D: return BLOCK_3D
    if sym in METAL_4D: return BLOCK_4D
    if sym in METAL_5D: return BLOCK_5D
    if sym in LANTHANIDE: return BLOCK_LN
    if sym in ACTINIDE: return BLOCK_AN
    if sym in SBLOCK: return BLOCK_S
    if sym in PBLOCK: return BLOCK_P
    return BLOCK_UNKNOWN


def _detect_metals(smiles: str) -> List[str]:
    """Extract metal symbols from bracketed atoms in SMILES."""
    metals: List[str] = []
    for m in _BRACKET_ATOM_RE.findall(smiles):
        if m in ALL_METALS:
            metals.append(m)
    return metals


def _majority_block(metals: List[str]) -> str:
    """Return the most-common metal-block among given metals."""
    if not metals:
        return BLOCK_NONE
    from collections import Counter
    blocks = [_metal_block_of(m) for m in metals]
    return Counter(blocks).most_common(1)[0][0]


def _heuristic_has_hapto(smiles: str) -> bool:
    """Detect hapto coordination heuristically.

    Hapto = metal bonded to multiple atoms of the same aromatic/conjugated ring.
    Signature in SMILES:
      - Metal-bracket followed by multiple ring-closure digits like [Fe]12345...
      - Or `(=` patterns with ring numbers attached to metal
      - Or metal directly bonded to multiple atoms in same ring
    """
    # Pattern: [<metal>]<digit>{2,} indicates multiple ring closures from metal
    for m in _BRACKET_ATOM_RE.findall(smiles):
        if m not in ALL_METALS:
            continue
    # Stronger heuristic: look for [<metal>][some chars]<digit><digit>+
    # e.g. [Fe]12345, [Cr]([CH]1)([CH]2)...
    pattern = re.compile(r"\[(?:[A-Z][a-z]?)(?:[+\-]\d*)?\](?:\d{2,}|\([^)]+\)\d|\d+\([^)]+\)\d)")
    if pattern.search(smiles):
        return True
    # Or 3+ ring-closure digits immediately after metal bracket
    pat2 = re.compile(r"\[[A-Z][a-z]?(?:[+\-]\d*)?\][0-9]{3,}")
    if pat2.search(smiles):
        return True
    return False


def _heuristic_has_bridge(smiles: str, n_metals: int) -> bool:
    """Detect μ-bridging heuristically.

    Only meaningful for n_metals >= 2.  Signature: same ring closure
    digit appearing in two different metal brackets' regions.
    """
    if n_metals < 2:
        return False
    # crude: count ring-closures, if heavy ring closures + multiple metals →
    # likely bridged or M-M bonded
    metals_brackets = re.findall(r"\[[A-Z][a-z]?[^\]]*\]", smiles)
    metal_brackets_with_digits = [b for b in metals_brackets
                                   if re.search(r"\d{2,}|\d\)\d", b)]
    return len(metal_brackets_with_digits) >= 2


def _heuristic_has_aromatic(smiles: str) -> bool:
    """Detect aromatic features in SMILES."""
    if any(c in smiles for c in "cnops"):
        return True
    if "C1=" in smiles or "c1" in smiles:
        return True
    return False


def _is_borane(smiles: str, metals: List[str]) -> bool:
    """Detect borane/carborane: B atoms forming cluster polyhedron.

    Heuristic: ≥6 B atoms + cage-like ring closures.
    """
    if metals:
        return False
    n_B = len(re.findall(r"\[B[^\]]*\]|B(?![a-zerih])", smiles))
    if n_B < 6:
        return False
    # ring-closure density
    digit_count = sum(1 for c in smiles if c.isdigit())
    return digit_count >= 6


def _is_cluster(smiles: str, n_metals: int) -> bool:
    """Detect M-M-bonded cluster (≥3 metals)."""
    return n_metals >= 3


# ----------------------------------------------------------------------
# Public API
# ----------------------------------------------------------------------

def classify_smiles(smiles: str) -> ClassFeatures:
    """Classify a SMILES string for routing decisions.

    Args:
        smiles: SMILES input string.

    Returns:
        ClassFeatures with coord_class + metal_block + auxiliary flags.

    Guarantees:
        - Pure function: no global state, no I/O.
        - Deterministic: same input → same output.

    Examples:
        >>> features = classify_smiles("[Fe-2]([Cl])([Cl])([Cl])[Cl]")
        >>> features.coord_class == "sigma"
        True
        >>> features.metal_block == "3d"
        True
    """
    metals = _detect_metals(smiles)
    n_metals = len(metals)
    block = _majority_block(metals)
    has_hapto = _heuristic_has_hapto(smiles) if n_metals > 0 else False
    has_aromatic = _heuristic_has_aromatic(smiles)
    has_bridge = _heuristic_has_bridge(smiles, n_metals)

    # Determine coord_class hierarchically
    if _is_borane(smiles, metals):
        coord = COORD_BORANE
    elif _is_cluster(smiles, n_metals):
        coord = COORD_CLUSTER
    elif n_metals == 0:
        coord = COORD_NO_METAL
    elif n_metals == 1:
        coord = COORD_HAPTO if has_hapto else COORD_SIGMA
    else:  # n_metals == 2
        coord = COORD_MULTI_HAPTO if has_hapto else COORD_MULTI_SIGMA

    return ClassFeatures(
        coord_class=coord,
        metal_block=block,
        metals=tuple(metals),
        n_metals=n_metals,
        has_hapto=has_hapto,
        has_aromatic=has_aromatic,
        has_bridge=has_bridge,
    )
