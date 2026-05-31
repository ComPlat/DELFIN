"""delfin.fffree.trans_effect — Universal trans-effect/trans-influence detection.

The trans-effect is the kinetic/structural influence a ligand exerts on the
group in the position trans to it. In TMC:
  - Strong trans-effect ligands (CO, CN-, H-, alkyl, NO) labilize the trans
    position (relevant for catalysis kinetics)
  - Trans-influence elongates the trans M-D bond (structural feature)

For each donor, identify its trans-partner. Compute M-D distance variations
based on trans-effect ranking.

Universal across all TMC.
FF-free, deterministic, geometric.

Env-gate: DELFIN_FFFREE_TRANS_EFFECT=1 (default OFF, byte-identical when unset).
Auto-enabled under PURE_TRACK3.
"""
from __future__ import annotations

import math
import os
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np


_PURE_TRACK3 = os.environ.get("DELFIN_FFFREE_PURE_TRACK3", "0") == "1"
_TRANS = _PURE_TRACK3 or os.environ.get("DELFIN_FFFREE_TRANS_EFFECT", "0") == "1"


# Trans-effect ranking (Pauling-Chatt order, ~Pt(II) systems)
# Higher value = stronger trans-effect (more elongation of trans M-D)
TRANS_EFFECT_RANK = {
    "C": 8,     # alkyl/CO/CN strongest
    "H": 7,     # hydride
    "P": 6,     # phosphine
    "Cl": 5, "Br": 5, "I": 5,  # halides moderate
    "N": 4,     # ammine/pyridine
    "O": 3,     # aqua/alkoxide
    "S": 4,     # thio
    "F": 2,     # fluoride weak
}


def find_trans_partner(
    coords: np.ndarray,
    metal_idx: int,
    donor_indices: Sequence[int],
    angle_threshold: float = 150.0,
) -> Dict[int, Optional[int]]:
    """For each donor, find the donor index in the trans position.

    Two donors are trans if D1-M-D2 angle > angle_threshold (e.g., 150°).

    Returns: {donor_idx: trans_partner_idx_or_None}.

    Universal across any coordination geometry.
    """
    if not _TRANS:
        return {}
    metal_pos = coords[int(metal_idx)]
    out = {}
    donors = [int(d) for d in donor_indices]
    for d1 in donors:
        v1 = coords[d1] - metal_pos
        best_partner = None
        best_angle = 0.0
        for d2 in donors:
            if d2 == d1:
                continue
            v2 = coords[d2] - metal_pos
            cos_a = float(np.dot(v1, v2) /
                          (np.linalg.norm(v1) * np.linalg.norm(v2) + 1e-9))
            cos_a = max(-1.0, min(1.0, cos_a))
            angle = math.degrees(math.acos(cos_a))
            if angle > best_angle and angle > angle_threshold:
                best_angle = angle
                best_partner = d2
        out[d1] = best_partner
    return out


def apply_trans_influence_elongation(
    coords: np.ndarray,
    syms: Sequence[str],
    metal_idx: int,
    donor_indices: Sequence[int],
    elongation_factor: float = 1.05,
) -> np.ndarray:
    """Apply trans-influence-driven M-D elongation.

    For each donor pair (D1, D2) that are trans to each other, the donor
    with the higher trans-effect ranking causes its partner to elongate.

    Returns coords with adjusted M-D bond lengths.
    Universal.
    """
    if not _TRANS:
        return coords
    P = coords.copy()
    metal_pos = P[int(metal_idx)]
    partners = find_trans_partner(P, metal_idx, donor_indices)
    for d1, d2 in partners.items():
        if d2 is None:
            continue
        s1 = syms[d1]
        s2 = syms[d2]
        rank1 = TRANS_EFFECT_RANK.get(s1, 3)
        rank2 = TRANS_EFFECT_RANK.get(s2, 3)
        # d1 is donor with higher trans-effect → elongate d2
        if rank1 > rank2:
            v2 = P[d2] - metal_pos
            d_v2 = float(np.linalg.norm(v2))
            if d_v2 > 1e-6:
                P[d2] = metal_pos + v2 * elongation_factor
    return P


if __name__ == "__main__":
    os.environ["DELFIN_FFFREE_TRANS_EFFECT"] = "1"
    import importlib, sys
    sys.modules.pop("delfin.fffree.trans_effect", None)
    from delfin.fffree.trans_effect import (
        find_trans_partner, apply_trans_influence_elongation, TRANS_EFFECT_RANK
    )

    print("=== Trans-effect ranking (Pauling-Chatt) ===")
    for sym, rank in sorted(TRANS_EFFECT_RANK.items(), key=lambda x: -x[1]):
        print(f"  {sym}: {rank}")

    # Test: Pt(NH3)(CO)Cl2 (cis-isomer, CO and Cl are trans to NH3 + Cl resp.)
    coords = np.array([
        [0, 0, 0],         # Pt (idx 0)
        [2.0, 0, 0],       # N (idx 1, NH3)
        [-2.0, 0, 0],      # C (idx 2, CO)
        [0, 2.0, 0],       # Cl (idx 3)
        [0, -2.0, 0],      # Cl (idx 4)
    ])
    syms = ["Pt", "N", "C", "Cl", "Cl"]
    partners = find_trans_partner(coords, 0, [1, 2, 3, 4])
    print(f"\nPt(N)(C)(Cl)2 trans partners:")
    for d1, d2 in partners.items():
        if d2 is not None:
            print(f"  {syms[d1]}(idx={d1}) trans to {syms[d2]}(idx={d2})")

    # Apply trans-influence: C should elongate N (C has higher rank)
    new_coords = apply_trans_influence_elongation(coords, syms, 0, [1, 2, 3, 4])
    print(f"\nM-D distances after trans-influence:")
    for i in range(1, 5):
        d = float(np.linalg.norm(new_coords[i] - new_coords[0]))
        print(f"  {syms[i]}(idx={i}): M-D = {d:.3f} (was 2.000)")
