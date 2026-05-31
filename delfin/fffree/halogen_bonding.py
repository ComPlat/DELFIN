"""delfin.fffree.halogen_bonding — Universal halogen bond (σ-hole) detection.

Halogen bonds R-X···Y (X = Cl/Br/I, Y = N/O/S lone-pair acceptor) are important
non-covalent interactions in TMC + supramolecular chemistry. The σ-hole on the
halogen creates a directional electrophilic site.

Detection criteria (Cavallo 2016):
  X···Y distance: < sum of vdW radii (e.g. <3.5 Å for Br···N)
  R-X···Y angle: > 140° (linear preferred, σ-hole alignment)

Universal across all R-X···Y combinations.
FF-free, deterministic, geometric.

Env-gate: DELFIN_FFFREE_XBOND=1 (default OFF, byte-identical when unset).
Auto-enabled under PURE_TRACK3.
"""
from __future__ import annotations

import math
import os
from typing import Dict, List, Sequence

import numpy as np


_PURE_TRACK3 = os.environ.get("DELFIN_FFFREE_PURE_TRACK3", "0") == "1"
_XBOND = _PURE_TRACK3 or os.environ.get("DELFIN_FFFREE_XBOND", "0") == "1"


# Halogen donor atoms (form σ-hole)
_XBOND_DONORS = {"Cl", "Br", "I"}
# Acceptor atoms (lone-pair available)
_XBOND_ACCEPTORS = {"N", "O", "S", "Se", "F", "Cl", "Br", "I"}

# VdW sum cutoffs (Å) — Bondi radii sum reduced 10%
_VDW_SUMS = {
    ("Cl", "N"): 3.16,  ("Cl", "O"): 3.10,  ("Cl", "S"): 3.45,  ("Cl", "F"): 3.05,
    ("Br", "N"): 3.30,  ("Br", "O"): 3.25,  ("Br", "S"): 3.59,  ("Br", "F"): 3.20,
    ("I",  "N"): 3.50,  ("I",  "O"): 3.45,  ("I",  "S"): 3.75,  ("I",  "F"): 3.40,
}


def detect_halogen_bonds(
    coords: np.ndarray,
    syms: Sequence[str],
    d_max_extra: float = 0.5,
    angle_min: float = 140.0,
) -> List[Dict]:
    """Detect R-X...Y halogen bonds.

    Algorithm:
      1. For each X atom (Cl/Br/I): find parent R (closest non-halogen)
      2. For each Y atom (acceptor): check X-Y distance + R-X...Y angle
      3. Distance: X-Y < vdW_sum + d_max_extra
      4. Angle: R-X-Y > angle_min (σ-hole alignment)

    Returns list of {x_idx, r_idx, y_idx, xy_dist, rxy_angle}.

    Universal across all R-X...Y combinations.
    """
    if not _XBOND:
        return []
    n = len(coords)
    out = []
    for x in range(n):
        if syms[x] not in _XBOND_DONORS:
            continue
        # Find parent R
        r_idx = None
        r_d = 999.0
        for k in range(n):
            if k == x or syms[k] == "H":
                continue
            d = float(np.linalg.norm(coords[x] - coords[k]))
            if d < 2.5 and d < r_d:
                r_idx = k
                r_d = d
        if r_idx is None:
            continue
        # Search acceptors
        for y in range(n):
            if y == x or y == r_idx:
                continue
            if syms[y] not in _XBOND_ACCEPTORS:
                continue
            d_xy = float(np.linalg.norm(coords[x] - coords[y]))
            cutoff = _VDW_SUMS.get((syms[x], syms[y]), 3.5) + d_max_extra
            if d_xy > cutoff or d_xy < 2.0:
                continue
            # R-X-Y angle
            v_xr = coords[r_idx] - coords[x]
            v_xy = coords[y] - coords[x]
            cos_a = float(np.dot(v_xr, v_xy) /
                          (np.linalg.norm(v_xr) * np.linalg.norm(v_xy) + 1e-9))
            cos_a = max(-1.0, min(1.0, cos_a))
            angle = 180.0 - math.degrees(math.acos(cos_a))
            if angle < angle_min:
                continue
            out.append({
                "x_idx": x, "r_idx": r_idx, "y_idx": y,
                "xy_dist": d_xy, "rxy_angle": angle,
                "x_sym": syms[x], "y_sym": syms[y],
            })
    return out


if __name__ == "__main__":
    os.environ["DELFIN_FFFREE_XBOND"] = "1"
    import importlib, sys
    sys.modules.pop("delfin.fffree.halogen_bonding", None)
    from delfin.fffree.halogen_bonding import detect_halogen_bonds

    # Test: simple I-C ... N system
    coords = np.array([
        [0.0, 0.0, 0.0],   # C
        [0.0, 0.0, 2.1],   # I (R-X)
        [0.0, 0.0, 5.2],   # N (acceptor at 3.1 A)
    ])
    syms = ["C", "I", "N"]
    xb = detect_halogen_bonds(coords, syms)
    print(f"Detected halogen bonds: {len(xb)}")
    for x in xb:
        print(f"  {x['x_sym']}(idx={x['x_idx']}) ... {x['y_sym']}(idx={x['y_idx']}) "
              f"d={x['xy_dist']:.2f}Å angle={x['rxy_angle']:.1f}°")
