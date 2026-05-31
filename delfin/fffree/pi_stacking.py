"""delfin.fffree.pi_stacking — Universal π-π stacking + CH-π interaction detection.

π-π stacking interactions (face-to-face, edge-to-face) are critical in:
  - DNA/RNA structure
  - Metalloporphyrin assemblies
  - Bipyridyl/terpyridyl-bridged complexes
  - Supramolecular hosts

CH-π interactions are weaker but ubiquitous in M-arene complexes + protein
binding sites.

Detection criteria:
  π-π face-to-face: ring centroids < 4.5 Å, normals ∥ (angle < 25°)
  π-π edge-to-face: centroids 4-6 Å, normals ⊥ (angle 60-90°)
  CH-π: H above ring centroid, H...centroid 2.5-3.5 Å

Universal across all aromatic systems.
FF-free, deterministic, geometric.

Env-gate: DELFIN_FFFREE_PI_STACKING=1 (default OFF, byte-identical when unset).
"""
from __future__ import annotations

import math
import os
from typing import Dict, List, Sequence

import numpy as np


_PURE_TRACK3 = os.environ.get("DELFIN_FFFREE_PURE_TRACK3", "0") == "1"
_PI = _PURE_TRACK3 or os.environ.get("DELFIN_FFFREE_PI_STACKING", "0") == "1"


def _ring_centroid_normal(ring_coords: np.ndarray) -> tuple:
    """Return (centroid, normal) of a ring via SVD."""
    centroid = ring_coords.mean(axis=0)
    centered = ring_coords - centroid
    _, _, Vt = np.linalg.svd(centered, full_matrices=False)
    normal = Vt[-1]
    return centroid, normal


def detect_pi_stacking(
    coords: np.ndarray,
    aromatic_rings: List[Sequence[int]],
    face_to_face_max_dist: float = 4.5,
    face_to_face_max_angle: float = 25.0,
    edge_to_face_dist_range: tuple = (3.5, 6.0),
    edge_to_face_angle_range: tuple = (60.0, 90.0),
) -> List[Dict]:
    """Detect π-π stacking between aromatic rings.

    Returns list of {ring_a, ring_b, dist, angle, mode} where mode is
    'face-to-face' or 'edge-to-face'.
    """
    if not _PI:
        return []
    n_rings = len(aromatic_rings)
    out = []
    centroids = []
    normals = []
    for r in aromatic_rings:
        c, nrm = _ring_centroid_normal(np.array([coords[int(i)] for i in r]))
        centroids.append(c)
        normals.append(nrm)
    for i in range(n_rings):
        for j in range(i + 1, n_rings):
            d = float(np.linalg.norm(centroids[i] - centroids[j]))
            cos_a = abs(float(np.dot(normals[i], normals[j])))
            cos_a = max(-1.0, min(1.0, cos_a))
            angle = math.degrees(math.acos(cos_a))
            if d < face_to_face_max_dist and angle < face_to_face_max_angle:
                out.append({"ring_a": i, "ring_b": j, "dist": d, "angle": angle,
                            "mode": "face-to-face"})
            elif (edge_to_face_dist_range[0] < d < edge_to_face_dist_range[1] and
                  edge_to_face_angle_range[0] < angle < edge_to_face_angle_range[1]):
                out.append({"ring_a": i, "ring_b": j, "dist": d, "angle": angle,
                            "mode": "edge-to-face"})
    return out


def detect_ch_pi(
    coords: np.ndarray,
    syms: Sequence[str],
    aromatic_rings: List[Sequence[int]],
    h_centroid_min: float = 2.3,
    h_centroid_max: float = 3.5,
) -> List[Dict]:
    """Detect C-H···π contacts (H above ring centroid)."""
    if not _PI:
        return []
    centroids_normals = []
    for r in aromatic_rings:
        c, nrm = _ring_centroid_normal(np.array([coords[int(i)] for i in r]))
        centroids_normals.append((c, nrm, set(int(i) for i in r)))
    out = []
    n = len(coords)
    for h in range(n):
        if syms[h] != "H":
            continue
        for ridx, (c, nrm, ring_set) in enumerate(centroids_normals):
            if h in ring_set:
                continue
            v = coords[h] - c
            dist = float(np.linalg.norm(v))
            if h_centroid_min < dist < h_centroid_max:
                # H should be above ring plane (cos angle > 0.5 = within ~60° of normal)
                cos_a = abs(float(np.dot(v, nrm) / (dist * float(np.linalg.norm(nrm)) + 1e-9)))
                if cos_a > 0.5:
                    out.append({"h_idx": h, "ring_idx": ridx, "dist": dist,
                                "cos_angle": cos_a})
    return out


if __name__ == "__main__":
    os.environ["DELFIN_FFFREE_PI_STACKING"] = "1"
    import importlib, sys
    sys.modules.pop("delfin.fffree.pi_stacking", None)
    from delfin.fffree.pi_stacking import detect_pi_stacking, detect_ch_pi

    # Test: 2 benzene rings face-to-face at 3.5 Å
    angles = np.array([i * math.pi / 3 for i in range(6)])
    ring_a = np.array([[math.cos(a) * 1.4, math.sin(a) * 1.4, 0] for a in angles])
    ring_b = ring_a + np.array([0, 0, 3.5])
    coords = np.vstack([ring_a, ring_b])
    rings = [[0, 1, 2, 3, 4, 5], [6, 7, 8, 9, 10, 11]]
    stacks = detect_pi_stacking(coords, rings)
    print(f"Face-to-face benzene pair: {len(stacks)} stacking contacts")
    for s in stacks:
        print(f"  rings {s['ring_a']}-{s['ring_b']}: d={s['dist']:.2f}Å "
              f"angle={s['angle']:.1f}° mode={s['mode']}")
