"""delfin.fffree.sandwich_piano_polyhedra — Phase C / Task #44 (2026-06-05).

Dedicated coordination polyhedra for **bis-η-sandwich**, **piano-stool**
(η⁵-Cp + 3 σ-donors), and **half-sandwich** (η⁶-arene + 3 σ-donors)
topologies.  These three geometries are NOT representable in the legacy
``delfin.fffree.polyhedra`` catalogue (which only carries point-donor
polyhedra CN2-9): a hapto ring contributes a **set of n vertices** lying
on a parallel circle around the M-centroid axis, NOT a single donor
vertex.  The mismatch is the root cause of the ferrocene-class M-Cp(centroid)
distance collapse observed in pool 2792332 (median 0.99 Å vs ideal 1.65 Å
for Fe and 1.313 Å for arenes vs ideal 1.62-1.69 Å — see Mission A7
forensik, 73 % more H-anomalies / structure on the Cp/arene subset than
on the whole pool).

Scope (Phase C, Task #44)
-------------------------

* CN10 ``SANDWICH-10``     bis-η⁵-Cp (D5d staggered / D5h eclipsed)
                            — Cp₂Fe (ferrocene), Cp*₂Ru, Cp₂Co⁺ ...
* CN8  ``PIANO-STOOL-8``   η⁵-Cp + 3 σ-donors                      (C3v)
                            — CpRe(CO)₃, CpMn(CO)₃, Cp*Ir(CO)₂L ...
* CN9  ``HALF-SANDWICH-9`` η⁶-arene + 3 σ-donors                   (C3v)
                            — (η⁶-cymene)Ru(Cl)₂(PR₃),
                              (η⁶-mesitylene)Cr(CO)₃ ...

Each ``*_vertices()`` returns a deterministic ``(N, 3)`` unit-vector array.
The hapto vertices are NOT placed on the unit sphere at the M-centroid
position (which would re-introduce the very collapse this module fixes);
instead they are placed on the parallel ring at the geometrically
correct M-centroid distance / radius ratio for typical Cp / arene, then
RE-NORMALISED to unit norm at the call site (``ref_vectors_sandwich``)
so the assembler's standard ``ref[v] / np.linalg.norm(ref[v]) * md_distance``
contract is honoured.  The TRUE M-centroid distance is reasserted by the
assembler via :func:`hapto_modes.m_centroid_distance` when the new
polyhedra are dispatched.

Vertex ORDER contract (must match Pólya group generators):

  SANDWICH-10:    idx 0-4  = top Cp ring   (+z, 0°,72°,144°,216°,288°)
                  idx 5-9  = bottom Cp ring (-z, 36°,108°,180°,252°,324°)
                            (D5d staggered — typical CCDC ferrocene)
  PIANO-STOOL-8:  idx 0-4  = Cp ring        (+z, 0°,72°,144°,216°,288°)
                  idx 5-7  = σ-donor tripod (-z, 0°,120°,240°)
  HALF-SANDWICH-9: idx 0-5 = arene ring     (+z, 0°,60°,120°,180°,240°,300°)
                   idx 6-8 = σ-donor tripod (-z, 0°,120°,240°)

Env-gate
--------

``DELFIN_FFFREE_SANDWICH_POLYHEDRA=1`` (or ``DELFIN_FFFREE_PURE_TRACK3=1``)
activates dispatch from ``polyhedra.ref_vectors`` /
``polyhedra.geometries_for_cn`` / ``decompose._default_geometry``.  Default
OFF and **byte-identical** to the legacy path otherwise — none of the
SANDWICH-10 / PIANO-STOOL-8 / HALF-SANDWICH-9 keys are reachable from the
default catalogue without the flag.

Determinism contract
--------------------

* No randomness, no platform-dependent branches in the vertex builders.
* All trig calls use ``math.cos / math.sin / math.radians`` for byte-stability.
* Vertex index order is fixed and documented; downstream isomer-dedup groups
  in :mod:`delfin.fffree.polya_isomer_count` operate on the documented order.
"""
from __future__ import annotations

import math
import os
from typing import Optional

import numpy as np


# Empirical hapto-ring tilt / opening (radians).  The hapto vertices on
# the ideal-vertex polyhedron lie on a CONE around the M-centroid axis at
# this half-angle; this fixes the M-C distance / M-centroid distance ratio
# at the CCDC-mean ratio for η⁵-Cp (~ 1.27) and η⁶-arene (~ 1.21).  The
# absolute distance is set by the assembler via
# :func:`hapto_modes.m_centroid_distance` after the unit-vector polyhedron
# has been picked up by ``ref_vectors_sandwich``.
ETA5_CP_RING_HALF_ANGLE = math.radians(38.0)   # M-C : M-centroid ≈ tan(38°) → 1/0.788
ETA6_ARENE_RING_HALF_ANGLE = math.radians(35.0)
SIGMA_TRIPOD_OPEN_ANGLE = math.radians(125.0)  # below-plane tripod opens ~125° vs +z


def _norm_rows(V: np.ndarray) -> np.ndarray:
    n = np.linalg.norm(V, axis=1, keepdims=True)
    # Defence in depth (zero-norm impossible by construction).
    return V / np.where(n > 1e-12, n, 1.0)


# ---------------------------------------------------------------------------
# SANDWICH-10  ( bis-η⁵-Cp ; D5d staggered )
# ---------------------------------------------------------------------------


def sandwich_10_vertices() -> np.ndarray:
    """SANDWICH-10 (D5d staggered ferrocene-type) — 10 unit vectors.

    Top Cp ring on the +z hemisphere (azimuths 0°, 72°, 144°, 216°, 288°)
    at half-angle :data:`ETA5_CP_RING_HALF_ANGLE` from +z.
    Bottom Cp ring on the -z hemisphere (azimuths 36°, 108°, 180°, 252°,
    324° — staggered 36°) at the mirror-image half-angle.

    Index ORDER (must match Pólya / Burnside dedup in
    :mod:`polya_isomer_count`): top 0-4 then bottom 5-9.
    """
    s = math.sin(ETA5_CP_RING_HALF_ANGLE)
    c = math.cos(ETA5_CP_RING_HALF_ANGLE)
    top = [
        [s * math.cos(math.radians(72.0 * k)),
         s * math.sin(math.radians(72.0 * k)),
         +c]
        for k in range(5)
    ]
    bot = [
        [s * math.cos(math.radians(36.0 + 72.0 * k)),
         s * math.sin(math.radians(36.0 + 72.0 * k)),
         -c]
        for k in range(5)
    ]
    return _norm_rows(np.asarray(top + bot, dtype=float))


# ---------------------------------------------------------------------------
# PIANO-STOOL-8  ( η⁵-Cp + 3 σ-donors ; C3v )
# ---------------------------------------------------------------------------


def piano_stool_8_vertices() -> np.ndarray:
    """PIANO-STOOL-8 (C3v) — η⁵-Cp (5 vertices, +z) + σ-tripod (3, -z).

    Top Cp ring at half-angle :data:`ETA5_CP_RING_HALF_ANGLE` from +z,
    azimuths 0°, 72°, 144°, 216°, 288°.
    Bottom σ-donor tripod at half-angle :data:`SIGMA_TRIPOD_OPEN_ANGLE`
    from +z (i.e. opens downward, away from Cp), azimuths 0°, 120°, 240°.

    Index ORDER: Cp ring 0-4, σ-donor tripod 5-7.
    """
    s_cp = math.sin(ETA5_CP_RING_HALF_ANGLE)
    c_cp = math.cos(ETA5_CP_RING_HALF_ANGLE)
    cp = [
        [s_cp * math.cos(math.radians(72.0 * k)),
         s_cp * math.sin(math.radians(72.0 * k)),
         +c_cp]
        for k in range(5)
    ]
    s_t = math.sin(SIGMA_TRIPOD_OPEN_ANGLE)
    c_t = math.cos(SIGMA_TRIPOD_OPEN_ANGLE)
    tripod = [
        [s_t * math.cos(math.radians(120.0 * k)),
         s_t * math.sin(math.radians(120.0 * k)),
         c_t]
        for k in range(3)
    ]
    return _norm_rows(np.asarray(cp + tripod, dtype=float))


# ---------------------------------------------------------------------------
# HALF-SANDWICH-9  ( η⁶-arene + 3 σ-donors ; C3v )
# ---------------------------------------------------------------------------


def half_sandwich_9_vertices() -> np.ndarray:
    """HALF-SANDWICH-9 (C3v) — η⁶-arene (6 vertices, +z) + σ-tripod (3, -z).

    Top arene ring at half-angle :data:`ETA6_ARENE_RING_HALF_ANGLE`
    from +z, azimuths 0°, 60°, 120°, 180°, 240°, 300°.
    Bottom σ-donor tripod at half-angle :data:`SIGMA_TRIPOD_OPEN_ANGLE`
    from +z, azimuths 30°, 150°, 270° (offset 30° so the tripod legs
    sit between projected arene atoms, mirroring the typical CCDC
    staggered piano-stool).

    Index ORDER: arene ring 0-5, σ-donor tripod 6-8.
    """
    s_ar = math.sin(ETA6_ARENE_RING_HALF_ANGLE)
    c_ar = math.cos(ETA6_ARENE_RING_HALF_ANGLE)
    arene = [
        [s_ar * math.cos(math.radians(60.0 * k)),
         s_ar * math.sin(math.radians(60.0 * k)),
         +c_ar]
        for k in range(6)
    ]
    s_t = math.sin(SIGMA_TRIPOD_OPEN_ANGLE)
    c_t = math.cos(SIGMA_TRIPOD_OPEN_ANGLE)
    tripod = [
        [s_t * math.cos(math.radians(30.0 + 120.0 * k)),
         s_t * math.sin(math.radians(30.0 + 120.0 * k)),
         c_t]
        for k in range(3)
    ]
    return _norm_rows(np.asarray(arene + tripod, dtype=float))


# ---------------------------------------------------------------------------
# Registry + dispatch helpers
# ---------------------------------------------------------------------------


SANDWICH_GEOM_VERTICES = {
    "SANDWICH-10 bis-eta5-Cp":  sandwich_10_vertices,
    "PIANO-STOOL-8 eta5-Cp+L3": piano_stool_8_vertices,
    "HALF-SANDWICH-9 eta6+L3":  half_sandwich_9_vertices,
}

SANDWICH_GEOM_BY_CN = {
    8:  ["PIANO-STOOL-8 eta5-Cp+L3"],
    9:  ["HALF-SANDWICH-9 eta6+L3"],
    10: ["SANDWICH-10 bis-eta5-Cp"],
}

# Per-geometry hapto layout descriptor — used by the assembler to know
# which vertices belong to a hapto ring (and how many).  ``eta`` is the
# hapto count, ``ring_indices`` are the vertex indices belonging to the
# ring, ``cn_effective`` is the coordination-number count (ring = 1
# coordination site).
SANDWICH_HAPTO_LAYOUT = {
    "SANDWICH-10 bis-eta5-Cp": {
        "rings": [
            {"eta": 5, "ring_indices": [0, 1, 2, 3, 4], "axis_sign": +1},
            {"eta": 5, "ring_indices": [5, 6, 7, 8, 9], "axis_sign": -1},
        ],
        "tripod_indices": [],
        "cn_effective": 2,
    },
    "PIANO-STOOL-8 eta5-Cp+L3": {
        "rings": [
            {"eta": 5, "ring_indices": [0, 1, 2, 3, 4], "axis_sign": +1},
        ],
        "tripod_indices": [5, 6, 7],
        "cn_effective": 4,
    },
    "HALF-SANDWICH-9 eta6+L3": {
        "rings": [
            {"eta": 6, "ring_indices": [0, 1, 2, 3, 4, 5], "axis_sign": +1},
        ],
        "tripod_indices": [6, 7, 8],
        "cn_effective": 4,
    },
}


def is_sandwich_geometry(geometry: str) -> bool:
    """Return True iff ``geometry`` is one of the SANDWICH-10 /
    PIANO-STOOL-8 / HALF-SANDWICH-9 names registered here."""
    return geometry in SANDWICH_GEOM_VERTICES


def ref_vectors_sandwich(geometry: str) -> np.ndarray:
    """Return the unit-vector array for ``geometry``.  KeyError if the
    name is not one of the SANDWICH-10 / PIANO-STOOL-8 / HALF-SANDWICH-9
    keys — caller is expected to dispatch back to the legacy catalogue.

    The vectors are deterministic (no env / no platform dependency); two
    successive calls return byte-identical arrays.
    """
    try:
        return SANDWICH_GEOM_VERTICES[geometry]()
    except KeyError:
        raise KeyError(geometry)


def layout_for(geometry: str) -> Optional[dict]:
    """Hapto-ring layout descriptor for ``geometry`` (or ``None``)."""
    return SANDWICH_HAPTO_LAYOUT.get(geometry)


# ---------------------------------------------------------------------------
# Env-gate helper (single source of truth)
# ---------------------------------------------------------------------------


def env_active() -> bool:
    """``True`` iff the dispatch into SANDWICH-10 / PIANO-STOOL-8 /
    HALF-SANDWICH-9 geometries is enabled by the environment.  Honours
    ``DELFIN_FFFREE_SANDWICH_POLYHEDRA`` (primary) and
    ``DELFIN_FFFREE_PURE_TRACK3`` (umbrella)."""
    return (os.environ.get("DELFIN_FFFREE_SANDWICH_POLYHEDRA", "0") == "1"
            or os.environ.get("DELFIN_FFFREE_PURE_TRACK3", "0") == "1")
