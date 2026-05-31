"""delfin.fffree.topology_isomers — Universal topology-isomer classifier.

For TMC isomers that are NOT distinguishable from SMILES alone (since SMILES
gives only connectivity, not 3D arrangement), classify:

  fac vs mer    : for κ³-tridentate ligands on CN6 (facial vs meridional)
  cis vs trans  : for κ²-bidentate pairs on CN4/CN6 (90° vs 180° donor-donor)
  Δ vs Λ        : for tris-chelate complexes (Delta vs Lambda helicity)
  axial vs equatorial : for TBP-5 (CN5 trigonal bipyramidal)

These are required for complete enumeration: each topology isomer can be a
distinct DFT-startable structure.

Universal: works on any TMC with bidentate/tridentate ligands.
FF-free, deterministic, geometric.

Env-gate: DELFIN_FFFREE_TOPOLOGY=1 (default OFF, byte-identical when unset).
Auto-enabled under PURE_TRACK3 for complete isomer enumeration.
"""
from __future__ import annotations

import math
import os
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np


_PURE_TRACK3 = os.environ.get("DELFIN_FFFREE_PURE_TRACK3", "0") == "1"
_TOPOLOGY = _PURE_TRACK3 or os.environ.get("DELFIN_FFFREE_TOPOLOGY", "0") == "1"


def classify_fac_mer(donor_positions: Sequence[np.ndarray],
                      metal_position: np.ndarray) -> str:
    """Classify κ³-tridentate as facial (fac) or meridional (mer).

    Algorithm: if all 3 donors lie approximately on a face of the polyhedron
    (i.e., their angular sum to metal is ~270° = 3×90° on a face), it's fac.
    If they're meridional (one in middle, two at 180°), it's mer.

    Returns 'fac' or 'mer'.
    """
    if len(donor_positions) != 3:
        return "unknown"
    d = [np.asarray(p) - np.asarray(metal_position) for p in donor_positions]
    # Compute angles d_i ↔ d_j around metal (D-M-D angles)
    angles = []
    for i in range(3):
        for j in range(i + 1, 3):
            cos_a = float(np.dot(d[i], d[j]) /
                          (np.linalg.norm(d[i]) * np.linalg.norm(d[j])))
            cos_a = max(-1.0, min(1.0, cos_a))
            angles.append(math.degrees(math.acos(cos_a)))
    # fac: all 3 angles ~90°; mer: 2 angles ~90° + 1 angle ~180°
    n_near_180 = sum(1 for a in angles if 150 < a < 210)
    if n_near_180 == 0:
        return "fac"
    if n_near_180 >= 1:
        return "mer"
    return "unknown"


def classify_cis_trans(donor_a: np.ndarray, donor_b: np.ndarray,
                        metal: np.ndarray) -> str:
    """Classify κ²-bidentate as cis (~90°) or trans (~180°)."""
    va = donor_a - metal
    vb = donor_b - metal
    cos_a = float(np.dot(va, vb) / (np.linalg.norm(va) * np.linalg.norm(vb)))
    cos_a = max(-1.0, min(1.0, cos_a))
    angle = math.degrees(math.acos(cos_a))
    if 150 < angle < 210:
        return "trans"
    return "cis"


def classify_delta_lambda(
    chelate_rings_atoms: Sequence[Sequence[int]],
    coords: np.ndarray,
    metal_idx: int,
) -> str:
    """Classify tris-chelate complex chirality as Δ (Delta) or Λ (Lambda).

    Algorithm: compute the sign of the scalar triple product of the
    chelate-ring centroids relative to the metal. Positive = Δ (right-handed
    propeller), negative = Λ (left-handed).

    Universal across [M(N-N)₃], [M(O-O)₃], [M(N-O)₃], etc.
    Requires exactly 3 chelate rings around the metal.
    """
    if len(chelate_rings_atoms) != 3:
        return "achiral"
    metal_pos = coords[int(metal_idx)]
    # Compute chelate-ring centroids (excluding metal)
    centroids = []
    for ring in chelate_rings_atoms:
        ring_no_metal = [i for i in ring if int(i) != int(metal_idx)]
        if not ring_no_metal:
            continue
        c = np.array([coords[int(i)] for i in ring_no_metal]).mean(axis=0)
        centroids.append(c - metal_pos)
    if len(centroids) < 3:
        return "achiral"
    # Scalar triple product (Delta vs Lambda by sign)
    triple = float(np.dot(centroids[0], np.cross(centroids[1], centroids[2])))
    return "Delta" if triple > 0 else "Lambda"


def classify_axial_equatorial(
    donor_position: np.ndarray,
    metal_position: np.ndarray,
    principal_axis: np.ndarray,
) -> str:
    """For CN5 TBP geometry: classify a donor as axial or equatorial.

    Axial: donor along principal axis (angle to axis ~0° or 180°).
    Equatorial: donor in trigonal plane (angle to axis ~90°).
    """
    v = donor_position - metal_position
    v_norm = float(np.linalg.norm(v))
    if v_norm < 1e-6:
        return "unknown"
    cos_a = float(np.dot(v, principal_axis) /
                  (v_norm * float(np.linalg.norm(principal_axis))))
    cos_a = max(-1.0, min(1.0, cos_a))
    angle = math.degrees(math.acos(cos_a))
    if angle < 30 or angle > 150:
        return "axial"
    return "equatorial"


def enumerate_topology_isomers(
    n_chelate_rings: int,
    n_tridentate_groups: int = 0,
    cn: int = 6,
) -> List[Dict]:
    """Enumerate all valid topology isomers for a coordination environment.

    Parameters
    ----------
    n_chelate_rings : number of bidentate chelate rings
    n_tridentate_groups : number of κ³-tridentate groups
    cn : coordination number

    Returns: list of topology dicts:
      {n_topology_isomers, types, labels}
    """
    if not _TOPOLOGY:
        return [{"label": "no-topology-enum", "n": 1}]
    out = []
    # Δ/Λ for tris-chelate
    if n_chelate_rings == 3:
        out.append({"label": "Delta", "type": "chirality"})
        out.append({"label": "Lambda", "type": "chirality"})
    # fac/mer for tridentate on CN6
    if n_tridentate_groups >= 1 and cn == 6:
        out.append({"label": "fac", "type": "tridentate-topology"})
        out.append({"label": "mer", "type": "tridentate-topology"})
    # cis/trans for bidentate-pair
    if n_chelate_rings == 2 and cn in (4, 6):
        out.append({"label": "cis", "type": "bidentate-topology"})
        out.append({"label": "trans", "type": "bidentate-topology"})
    # axial/equatorial-pair for CN5 TBP
    if cn == 5:
        out.append({"label": "axial-axial", "type": "TBP-position"})
        out.append({"label": "axial-equatorial", "type": "TBP-position"})
        out.append({"label": "equatorial-equatorial", "type": "TBP-position"})
    return out


if __name__ == "__main__":
    os.environ["DELFIN_FFFREE_TOPOLOGY"] = "1"
    import importlib, sys
    sys.modules.pop("delfin.fffree.topology_isomers", None)
    from delfin.fffree.topology_isomers import (
        classify_fac_mer, classify_cis_trans, classify_delta_lambda,
        classify_axial_equatorial, enumerate_topology_isomers
    )

    # Test fac/mer
    metal = np.array([0, 0, 0])
    fac_donors = [np.array([1, 0, 0]), np.array([0, 1, 0]), np.array([0, 0, 1])]
    mer_donors = [np.array([1, 0, 0]), np.array([0, 1, 0]), np.array([-1, 0, 0])]
    print(f"fac test donors: {classify_fac_mer(fac_donors, metal)}")
    print(f"mer test donors: {classify_fac_mer(mer_donors, metal)}")

    # Test cis/trans
    print(f"cis bidentate: {classify_cis_trans(np.array([1,0,0]), np.array([0,1,0]), metal)}")
    print(f"trans bidentate: {classify_cis_trans(np.array([1,0,0]), np.array([-1,0,0]), metal)}")

    # Test axial/equatorial
    axis = np.array([0, 0, 1])
    print(f"axial donor: {classify_axial_equatorial(np.array([0,0,2]), metal, axis)}")
    print(f"equatorial donor: {classify_axial_equatorial(np.array([2,0,0]), metal, axis)}")

    # Enumerate
    print("\n=== Topology enumeration cases ===")
    for n_chel, n_trid, cn in [(3, 0, 6), (2, 0, 6), (1, 1, 6), (0, 0, 5)]:
        topos = enumerate_topology_isomers(n_chel, n_trid, cn)
        labels = [t["label"] for t in topos]
        print(f"  n_chel={n_chel}, n_trid={n_trid}, CN={cn}: {labels}")
