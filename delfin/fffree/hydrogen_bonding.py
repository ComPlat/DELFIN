"""delfin.fffree.hydrogen_bonding — Universal hydrogen-bond detection in TMCs.

Hydrogen bonds D-H···A (D=donor atom, A=acceptor atom) are critical in many
TMC structures (water bridges, NH···Cl, OH···O ligand-ligand interactions).
For complete enumeration: H-bond network topology can distinguish isomers.

Universal SMARTS + geometric detection.
FF-free, deterministic, no atom-specific code.

Env-gate: DELFIN_FFFREE_HBOND=1 (default OFF, byte-identical when unset).
Auto-enabled under PURE_TRACK3.
"""
from __future__ import annotations

import math
import os
from typing import Dict, List, Optional, Sequence

import numpy as np


_PURE_TRACK3 = os.environ.get("DELFIN_FFFREE_PURE_TRACK3", "0") == "1"
_HBOND = _PURE_TRACK3 or os.environ.get("DELFIN_FFFREE_HBOND", "0") == "1"


# Standard H-bond donor / acceptor atom types
_HBOND_DONORS = {"N", "O", "F", "S"}      # heavy atoms that can donate via H
_HBOND_ACCEPTORS = {"N", "O", "F", "Cl", "Br", "S"}


def detect_hbonds(
    coords: np.ndarray,
    syms: Sequence[str],
    d_da_max: float = 3.5,   # heavy-heavy max distance
    d_ha_min: float = 1.5,   # H to acceptor min
    d_ha_max: float = 2.8,   # H to acceptor max
    angle_min: float = 120.0,  # D-H···A min angle
) -> List[Dict]:
    """Detect classical D-H···A hydrogen bonds.

    Algorithm:
      1. For each H atom: identify parent heavy donor D (closest D)
      2. Search for acceptor A with H···A within [d_ha_min, d_ha_max]
      3. Verify D-H···A angle >= angle_min (typically > 120°)
      4. Verify D···A distance < d_da_max

    Returns list of {h_idx, d_idx, a_idx, d_da, d_ha, angle}.

    Universal across all H-bond donor/acceptor combinations.
    """
    if not _HBOND:
        return []
    n = len(coords)
    out = []
    for h in range(n):
        if syms[h] != "H":
            continue
        # Find donor parent
        donor_idx = None
        d_dh = 999.0
        for k in range(n):
            if k == h or syms[k] == "H":
                continue
            if syms[k] not in _HBOND_DONORS:
                continue
            d = float(np.linalg.norm(coords[h] - coords[k]))
            if d < 1.3 and d < d_dh:
                donor_idx = k
                d_dh = d
        if donor_idx is None:
            continue
        # Search acceptors
        for a in range(n):
            if a == h or a == donor_idx or syms[a] == "H":
                continue
            if syms[a] not in _HBOND_ACCEPTORS:
                continue
            d_ha = float(np.linalg.norm(coords[h] - coords[a]))
            if not (d_ha_min < d_ha < d_ha_max):
                continue
            d_da = float(np.linalg.norm(coords[donor_idx] - coords[a]))
            if d_da > d_da_max:
                continue
            # D-H-A angle
            v_hd = coords[donor_idx] - coords[h]
            v_ha = coords[a] - coords[h]
            cos_a = float(np.dot(v_hd, v_ha) /
                          (np.linalg.norm(v_hd) * np.linalg.norm(v_ha) + 1e-9))
            cos_a = max(-1.0, min(1.0, cos_a))
            angle = 180.0 - math.degrees(math.acos(cos_a))
            if angle < angle_min:
                continue
            out.append({
                "h_idx": h, "d_idx": donor_idx, "a_idx": a,
                "d_da": d_da, "d_ha": d_ha, "angle": angle,
                "donor_sym": syms[donor_idx], "acceptor_sym": syms[a],
            })
    return out


if __name__ == "__main__":
    os.environ["DELFIN_FFFREE_HBOND"] = "1"
    import importlib, sys
    sys.modules.pop("delfin.fffree.hydrogen_bonding", None)
    from delfin.fffree.hydrogen_bonding import detect_hbonds

    # Test: water dimer with H-bond
    coords = np.array([
        [0.0, 0.0, 0.0],         # O1 (idx 0)
        [-0.76, -0.59, 0.0],     # H1a
        [0.76, -0.59, 0.0],      # H1b (donor)
        [0.76, -0.59 + 1.8, 0.0],  # O2 (idx 3)
        [-0.76 + 0.76, -1.18 + 1.8, 0],  # H2a
        [0.76 + 0.76, -1.18 + 1.8, 0],   # H2b
    ])
    syms = ["O", "H", "H", "O", "H", "H"]
    hbs = detect_hbonds(coords, syms)
    print(f"Water dimer test: {len(hbs)} H-bonds")
    for hb in hbs:
        print(f"  H={hb['h_idx']}, D={hb['d_idx']}({hb['donor_sym']}), "
              f"A={hb['a_idx']}({hb['acceptor_sym']}), "
              f"D...A={hb['d_da']:.2f}Å, H...A={hb['d_ha']:.2f}Å, "
              f"angle={hb['angle']:.1f}°")
