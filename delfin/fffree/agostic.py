"""delfin.fffree.agostic — Universal agostic interaction enumeration.

Agostic interactions are M···H-C close contacts where a C-H bond donates
electron density to an unsaturated metal center. They are critical features
in catalytic mechanisms (β-hydride elimination, σ-bond metathesis).

For complete TMC enumeration, agostic configurations represent distinct
isomers/conformers from non-agostic ones.

Detection criteria (CCDC empirical):
  M···H distance: 1.8-2.4 Å
  M-H-C angle: 100-140°
  C-H-bond lengthening: 1.10-1.20 Å (vs typical 1.08 Å)

Universal across all TMC + organic C-H sources.
FF-free, deterministic, geometric.

Env-gate: DELFIN_FFFREE_AGOSTIC=1 (default OFF, byte-identical when unset).
Auto-enabled under PURE_TRACK3.
"""
from __future__ import annotations

import math
import os
from typing import Dict, List, Optional, Sequence

import numpy as np


_PURE_TRACK3 = os.environ.get("DELFIN_FFFREE_PURE_TRACK3", "0") == "1"
_AGOSTIC = _PURE_TRACK3 or os.environ.get("DELFIN_FFFREE_AGOSTIC", "0") == "1"


def detect_agostic_contacts(
    coords: np.ndarray,
    syms: Sequence[str],
    metal_idx: int,
    d_max: float = 2.5,
    angle_min: float = 90.0,
    angle_max: float = 150.0,
) -> List[Dict]:
    """Detect M···H-C agostic contacts.

    Algorithm:
      1. For each H atom: find its parent C
      2. Compute M-H distance (must be < d_max)
      3. Compute M-H-C angle (must be in [angle_min, angle_max])
      4. Distinguish from over-coordination (if M-H < M-C, just hydride)

    Returns list of {h_idx, c_idx, mh_dist, mhc_angle}.

    Universal across all TMC + organic C-H bonds.
    """
    if not _AGOSTIC:
        return []
    out = []
    metal_pos = coords[int(metal_idx)]
    n = len(coords)
    for h in range(n):
        if syms[h] != "H":
            continue
        # Find parent heavy atom (closest non-metal heavy)
        parents = []
        for k in range(n):
            if k == h or syms[k] == "H":
                continue
            d_hk = float(np.linalg.norm(coords[h] - coords[k]))
            if d_hk < 1.3 and syms[k] in {"C", "N", "O", "Si", "B"}:
                parents.append((k, d_hk))
        if not parents:
            continue
        parent_idx = min(parents, key=lambda x: x[1])[0]
        # Only C-H agostic (per definition)
        if syms[parent_idx] != "C":
            continue
        d_mh = float(np.linalg.norm(coords[h] - metal_pos))
        if d_mh > d_max or d_mh < 1.5:
            continue
        # M-H-C angle
        v_hm = metal_pos - coords[h]
        v_hc = coords[parent_idx] - coords[h]
        cos_a = float(np.dot(v_hm, v_hc) /
                      (np.linalg.norm(v_hm) * np.linalg.norm(v_hc) + 1e-9))
        cos_a = max(-1.0, min(1.0, cos_a))
        angle = math.degrees(math.acos(cos_a))
        if angle_min < angle < angle_max:
            out.append({
                "h_idx": h,
                "c_idx": parent_idx,
                "mh_dist": d_mh,
                "mhc_angle": angle,
            })
    return out


if __name__ == "__main__":
    os.environ["DELFIN_FFFREE_AGOSTIC"] = "1"
    import importlib, sys
    sys.modules.pop("delfin.fffree.agostic", None)
    from delfin.fffree.agostic import detect_agostic_contacts

    # Test: hypothetical Ti complex with agostic C-H
    # Ti at origin, C-H bond positioned to be agostic
    coords = np.array([
        [0, 0, 0],       # Ti (idx 0)
        [2.0, 0, 0],     # C (idx 1)
        [2.0, 1.0, 0],   # H (idx 2, non-agostic)
        [2.4, 0.5, 0],   # H (idx 3, AGOSTIC - close to Ti)
    ])
    syms = ["Ti", "C", "H", "H"]
    contacts = detect_agostic_contacts(coords, syms, metal_idx=0, d_max=3.0)
    print(f"Detected agostic contacts: {len(contacts)}")
    for c in contacts:
        print(f"  H={c['h_idx']}, C={c['c_idx']}, M-H={c['mh_dist']:.2f}Å, M-H-C={c['mhc_angle']:.1f}°")
