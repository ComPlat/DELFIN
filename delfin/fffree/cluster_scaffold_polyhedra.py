"""delfin.fffree.cluster_scaffold_polyhedra — facade over the legacy
``delfin._cluster_scaffold`` PROTOTYPE so the multi-metal-assemble
orchestrator can import skeleton primitives from a single namespace.

Purpose
-------
The 844-LOC prototype in ``delfin/_cluster_scaffold.py`` defines a rich
set of metal-skeleton primitives (triangle, linear-M3, tetrahedron,
butterfly-M4, octahedron, star-M5, multihapto skeleton) plus donor
partition logic (terminal / edge-bridge / face-cap) that the multi-metal
assemble path needs.  Importing them from a stable namespace inside
``delfin.fffree`` keeps the orchestrator's dependency surface clean and
makes the prototype's promotion to production a one-line search.

This file does **not** wrap, reimplement or change behaviour.  It is a
re-export only — every callable below is a *direct* reference to the
underlying prototype symbol.  Default OFF byte-identical: importing
this module changes no runtime behaviour.

Env-gate (consumed by the orchestrator):
``DELFIN_FFFREE_MULTI_METAL=1`` (or ``DELFIN_FFFREE_PURE_TRACK3=1``).
This module itself is always-importable.
"""
from __future__ import annotations

from delfin._cluster_scaffold import (  # noqa: F401 — re-export
    MetalSkeleton,
    ClusterBuildPlan,
    is_cluster_complex,
    enumerate_skeletons_for_mol,
    classify_donors_for_skeleton,
    place_skeleton_with_ligands,
    proof_of_concept_run,
    is_multihapto_skeleton_candidate,
    build_multihapto_scaffold,
    _skeleton_triangle,
    _skeleton_linear_m3,
    _skeleton_tetrahedron,
    _skeleton_butterfly_m4,
    _skeleton_octahedron_m6,
    _skeleton_star_m5,
)

__all__ = [
    "MetalSkeleton",
    "ClusterBuildPlan",
    "is_cluster_complex",
    "enumerate_skeletons_for_mol",
    "classify_donors_for_skeleton",
    "place_skeleton_with_ligands",
    "proof_of_concept_run",
    "is_multihapto_skeleton_candidate",
    "build_multihapto_scaffold",
    "_skeleton_triangle",
    "_skeleton_linear_m3",
    "_skeleton_tetrahedron",
    "_skeleton_butterfly_m4",
    "_skeleton_octahedron_m6",
    "_skeleton_star_m5",
]
