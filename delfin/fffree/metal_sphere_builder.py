#!/usr/bin/env python3
"""metal_sphere_builder.py — Phase 2: metal-FF-free coordination-sphere placement.

The core that replaces OB-UFF for the metal: place donor atoms on COD-ideal
polyhedron vertices (from metric_coord_shape._REFS) at empirical M-D distances —
NO force field on the metal (UFF can't type 4d/5d / has no η-term).  Deterministic
geometric construction.

M-D distance = covalent-radii sum (metal + donor); avoids the _bond_decollapse
_ideal_bond bug (defaults metals to 0.9 A).  Realistic-distortion hook (Phase-0.3
finding: T-4/TBP/SPY have COD median CShM ~2, not 0 -> ideal placement is "too
perfect" and the discriminator would flag it) seeded + deterministic.

Validated: ideal ML6 octahedron -> CShM(OC-6) ~0, all M-D = target, cis 90/trans 180.
"""
from __future__ import annotations
import math
from typing import List, Tuple
import numpy as np
from delfin.fffree import polyhedra as _PH

# self-contained references (delfin.fffree.polyhedra) — no agent_workspace dep
_COV = _PH.COV
_GEOM_BY_CN = _PH.GEOM_BY_CN
md_distance = _PH.md_distance
_ref_vectors = _PH.ref_vectors


def place_sphere(metal: str, donor_elems: List[str], geometry: str,
                 distort: float = 0.0, seed: int = 42
                 ) -> Tuple[List[str], np.ndarray]:
    """Metal at origin + donors on the ideal polyhedron at M-D distances.
    distort>0 adds a seeded radial/angular perturbation (to reach realistic COD
    CShM instead of a too-perfect 0)."""
    ref = _ref_vectors(geometry)
    n = len(donor_elems)
    if n != len(ref):
        raise ValueError(f"{n} donors != {len(ref)} vertices for {geometry}")
    syms = [metal] + list(donor_elems)
    P = [np.zeros(3)]
    rng = np.random.RandomState(seed)
    for i, d in enumerate(donor_elems):
        v = ref[i] / np.linalg.norm(ref[i])
        r = md_distance(metal, d)
        if distort > 0:
            v = v + distort * rng.randn(3)        # angular wobble
            v = v / np.linalg.norm(v)
            r = r * (1.0 + distort * 0.3 * rng.randn())   # radial wobble
        P.append(v * r)
    return syms, np.array(P, float)
