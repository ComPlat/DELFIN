#!/usr/bin/env python3
"""metal_sphere_builder.py — metal-FF-free coordination-sphere placement.

Places donor atoms on ideal polyhedron vertices (delfin.manta.polyhedra) at
empirical metal-donor distances, with NO force field on the metal.  Deterministic
geometric construction.

M-D distance = covalent-radii sum (metal + donor).  Optional seeded distortion
hook produces a realistic (rather than perfectly ideal) coordination sphere.
"""
from __future__ import annotations
import math
from typing import List, Tuple
import numpy as np
from delfin.manta import polyhedra as _PH

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
