"""delfin.fffree.pseudo_rotation — Universal Berry pseudo-rotation enumeration
for fluxional CN5 trigonal-bipyramidal complexes.

The Berry mechanism interconverts the two distinct trigonal-bipyramidal (TBP)
and square-pyramidal (SPY) configurations of a CN5 metal. For complete
enumeration, all 5 distinct Berry pseudo-rotation TBP↔SPY transitions are
canonical states.

For PF5 (textbook example): 5 distinct TBP arrangements correspond to choosing
which 2 of 5 ligands are axial. Berry rotation cycles through these.

Universal across all CN5 complexes.
FF-free, deterministic, geometric.

Env-gate: DELFIN_FFFREE_BERRY=1 (default OFF, byte-identical when unset).
"""
from __future__ import annotations

import math
import os
from typing import Iterator, List, Sequence, Tuple

import numpy as np


_PURE_TRACK3 = os.environ.get("DELFIN_FFFREE_PURE_TRACK3", "0") == "1"
_BERRY = _PURE_TRACK3 or os.environ.get("DELFIN_FFFREE_BERRY", "0") == "1"


def enumerate_berry_axial_pairs(n_donors: int) -> List[Tuple[int, int]]:
    """For CN5 TBP, enumerate all distinct (axial_a, axial_b) pairs.

    Returns 10 pairs (5 choose 2). Each defines one Berry-distinct configuration.

    Universal across any CN5 complex with N=5 distinguishable donors.
    """
    if not _BERRY or n_donors != 5:
        return [(0, 4)]  # default trivial
    out = []
    for i in range(n_donors):
        for j in range(i + 1, n_donors):
            out.append((i, j))
    return out


def assemble_tbp_from_axial_pair(
    metal_pos: np.ndarray,
    donor_indices: Sequence[int],
    axial_a: int,
    axial_b: int,
    md_distance: float = 2.0,
) -> Dict:
    """Place 5 donors around metal with axial_a/axial_b on the principal axis.

    Returns dict {donor_idx: position} for each donor at TBP positions:
      axial_a: (0, 0, +md)
      axial_b: (0, 0, -md)
      equatorial 3: at 120° in xy plane

    Universal: any 5 donor indices.
    """
    if not _BERRY:
        return {}
    out = {}
    out[donor_indices[axial_a]] = metal_pos + np.array([0, 0, md_distance])
    out[donor_indices[axial_b]] = metal_pos + np.array([0, 0, -md_distance])
    # Equatorial: 3 donors at 0°, 120°, 240° in xy plane
    eq_indices = [i for i in range(5) if i not in (axial_a, axial_b)]
    for k, idx in enumerate(eq_indices):
        angle = 2 * math.pi * k / 3
        out[donor_indices[idx]] = metal_pos + np.array([
            md_distance * math.cos(angle),
            md_distance * math.sin(angle),
            0,
        ])
    return out


def enumerate_berry_tbp_states(
    metal_pos: np.ndarray,
    donor_indices: Sequence[int],
    md_distance: float = 2.0,
) -> Iterator[Tuple[Tuple[int, int], dict]]:
    """Enumerate all Berry-distinct TBP configurations for CN5.

    Yields ((axial_a, axial_b), {donor_idx: position}) for each pair.
    Universal CN5.
    """
    if not _BERRY:
        return
    pairs = enumerate_berry_axial_pairs(len(donor_indices))
    for axial_a, axial_b in pairs:
        config = assemble_tbp_from_axial_pair(
            metal_pos, donor_indices, axial_a, axial_b, md_distance
        )
        yield ((axial_a, axial_b), config)


from typing import Dict


if __name__ == "__main__":
    os.environ["DELFIN_FFFREE_BERRY"] = "1"
    import importlib, sys
    sys.modules.pop("delfin.fffree.pseudo_rotation", None)
    from delfin.fffree.pseudo_rotation import (
        enumerate_berry_axial_pairs, enumerate_berry_tbp_states
    )

    print("=== Berry pseudo-rotation pairs for CN5 ===")
    pairs = enumerate_berry_axial_pairs(5)
    for p in pairs:
        print(f"  axial pair: {p}")
    print(f"Total: {len(pairs)} distinct Berry-TBP configurations")

    # Test assembly
    metal_pos = np.array([0, 0, 0])
    donor_indices = [1, 2, 3, 4, 5]
    print(f"\n=== Sample TBP assembly (axial = 0, 1) ===")
    for pair, config in list(enumerate_berry_tbp_states(metal_pos, donor_indices))[:2]:
        print(f"  pair {pair}:")
        for d_idx, pos in sorted(config.items()):
            print(f"    donor {d_idx}: {pos}")
