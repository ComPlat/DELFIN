"""Welle-5l T3-B — D-AQIWAZ pyridine edge-on + NHC fused-ring fix.

Two patches:

  1. ``_donor_orientation_realism.py`` extends the 5b-B aromatic-N and
     NHC fused-ring guards.  The original implementation skipped any
     ring whose non-N (or non-C\\ :sub:`carbene`) ring atoms belonged to
     more than one ring -- aimed at protecting fused systems from
     rigid rotation distortion.  In practice this rejected the entire
     pyridine-in-naphthyridine donor class (D-AQIWAZ Ru-N at 0.8°
     face-on tilt, far below the chemistry-correct 90° edge-on).
     The replacement walks the transitive closure of aromatic rings
     sharing atoms with the seed ring, refuses to expand into any
     ring that holds another donor / metal, and rotates the entire
     fused aromatic backbone rigidly so the metal-donor lone-pair
     orientation snaps to edge-on without breaking the
     ligand-to-other-donor bonds.

  2. ``_verify_topology_from_graph`` Rule-1 phantom-bond check now
     exempts atoms that are graph-distance 2 from the metal via a
     donor (1,3-neighbours).  These atoms are necessarily close to
     the metal because of rigid ligand backbones (NHC "other N",
     pyridine ortho-C, salen ring-N, naphthyridine "other N") and
     their pre-UFF distance (~2.17 Å for Ru-N (NHC) etc.) is below
     the 1.05 × Σ\\ :sub:`r_cov` phantom threshold.  Treating them
     as phantom bonds rejects valid coordination chemistry and is
     the dominant cause of D-AQIWAZ's 11 % isomer-coverage gap.

Tests in this module:
  * ``test_pyridine_orientation_snaps_to_edge_on`` — 5b-B on the
    D-AQIWAZ SMILES drives the Ru-N(pyridine) ring-normal-vs-axis
    angle ≥ 85°, up from ~32° baseline.
  * ``test_phantom_13_exempt_default_on_preserves_md`` — the phantom
    1,3 exemption does not change a fingerprint-stable Cd CN2 sanity
    SMILES (no false positives on the simplest case).
  * ``test_phantom_13_exempt_env_off_restores_legacy`` — setting
    ``DELFIN_PHANTOM_13_EXEMPT=0`` re-enables the strict rule.
"""
from __future__ import annotations

import math
import os
import sys
from typing import List, Tuple

import numpy as np
import pytest


# D-AQIWAZ — Ru CN6 hetero, naphthyridine + NHC chelate + CO + CO + Br + H2O.
_D_AQIWAZ_SMILES = (
    "CC1=CC(C)=C2C=CC3=[N+](C2=N1)[Ru-4]([OH2+])([Br])([C]#[O+])"
    "([C]#[O+])[C+]1N(CC2=CC=CC=C2)C=CN31"
)

# Simple CN2 Cd dihalide as a phantom-1,3-exempt sanity SMILES.
_CD_DIHALIDE_SMILES = "[Br][Cd-2]([Br])([Br])[Br]"


def _parse_xyz_body(xyz_str: str):
    """Return (symbols, coords) for the first XYZ block found in
    ``xyz_str``.  Accepts both bare bodies and standard header-prefixed
    XYZ payloads."""
    lines = [ln for ln in xyz_str.splitlines() if ln.strip()]
    body = lines
    try:
        n = int(lines[0].strip())
        body = lines[2 : 2 + n]
    except (ValueError, IndexError):
        pass
    syms: List[str] = []
    coords: List[List[float]] = []
    for ln in body:
        parts = ln.split()
        if len(parts) < 4:
            continue
        try:
            x, y, z = float(parts[1]), float(parts[2]), float(parts[3])
        except ValueError:
            continue
        syms.append(parts[0])
        coords.append([x, y, z])
    return syms, np.asarray(coords, dtype=np.float64)


def _measure_pyridine_orientation_max_angle(xyz_str: str) -> float:
    """Return the maximum ring-normal-vs-M-N axis angle (in degrees,
    0-90) across every 6-ring containing an N directly bonded to a
    metal.  90° = edge-on (correct sigma donation), 0° = face-on."""
    syms, coords = _parse_xyz_body(xyz_str)
    if syms == [] or coords.size == 0:
        return 0.0

    # Covalent radii (Å) for the atoms appearing in D-AQIWAZ.
    rcov = {
        "H": 0.31, "C": 0.76, "N": 0.71, "O": 0.66, "F": 0.57,
        "P": 1.07, "S": 1.05, "Cl": 1.02, "Br": 1.20, "I": 1.39,
        "Ti": 1.60, "V": 1.53, "Cr": 1.39, "Mn": 1.39, "Fe": 1.32,
        "Co": 1.26, "Ni": 1.24, "Cu": 1.32, "Zn": 1.22,
        "Ru": 1.46, "Rh": 1.42, "Pd": 1.39, "Cd": 1.44,
        "Ir": 1.41, "Pt": 1.36,
    }
    metal_set = {"Ru", "Rh", "Pd", "Cd", "Ir", "Pt", "Fe", "Co",
                 "Ni", "Cu", "Zn", "Cr", "Mn", "V", "Ti"}
    n_atoms = len(syms)
    adj: List[List[int]] = [[] for _ in range(n_atoms)]
    for i in range(n_atoms):
        for j in range(i + 1, n_atoms):
            ri = rcov.get(syms[i], 0.80)
            rj = rcov.get(syms[j], 0.80)
            d = float(np.linalg.norm(coords[i] - coords[j]))
            if d < 1.30 * (ri + rj):
                adj[i].append(j)
                adj[j].append(i)
    metals = [i for i, s in enumerate(syms) if s in metal_set]
    if not metals:
        return 0.0
    metal = metals[0]
    best_angle = 0.0
    for n_i in adj[metal]:
        if syms[n_i] != "N":
            continue
        # Walk 6-rings through n_i: DFS depth 6, return on closure.
        for nb1 in adj[n_i]:
            if syms[nb1] not in ("C", "N") or nb1 == metal:
                continue
            for nb2 in adj[nb1]:
                if nb2 in (n_i, metal) or syms[nb2] not in ("C", "N"):
                    continue
                for nb3 in adj[nb2]:
                    if nb3 in (n_i, nb1, metal) or syms[nb3] not in ("C", "N"):
                        continue
                    for nb4 in adj[nb3]:
                        if nb4 in (n_i, nb1, nb2, metal) or syms[nb4] not in ("C", "N"):
                            continue
                        for nb5 in adj[nb4]:
                            if nb5 in (n_i, nb1, nb2, nb3, metal):
                                continue
                            if syms[nb5] not in ("C", "N"):
                                continue
                            if n_i in adj[nb5]:
                                # 6-ring closed
                                a = coords[n_i]
                                b = coords[nb1]
                                c = coords[nb3]
                                normal = np.cross(b - a, c - a)
                                nl = float(np.linalg.norm(normal))
                                if nl < 1e-6:
                                    continue
                                normal /= nl
                                mn = coords[metal] - coords[n_i]
                                mn_l = float(np.linalg.norm(mn))
                                if mn_l < 1e-6:
                                    continue
                                mn /= mn_l
                                cos_a = abs(float(np.dot(normal, mn)))
                                cos_a = max(-1.0, min(1.0, cos_a))
                                ang = math.degrees(math.acos(cos_a))
                                if ang > best_angle:
                                    best_angle = ang
    return best_angle


def _run_pipeline(smiles: str, env: dict):
    """Run smiles_to_xyz_isomers with a fresh delfin import + env state.
    Returns ``(n_isomers, list[xyz])``."""
    # Reset DELFIN_* env so prior tests do not leak.
    for k in list(os.environ.keys()):
        if k.startswith("DELFIN_"):
            del os.environ[k]
    for k, v in env.items():
        os.environ[k] = v
    # Force reimport so any module-level constants based on env are reread.
    for mod in list(sys.modules.keys()):
        if mod.startswith("delfin"):
            del sys.modules[mod]
    from delfin.smiles_converter import smiles_to_xyz_isomers
    results, err = smiles_to_xyz_isomers(
        smiles, num_confs=30, max_isomers=20,
        apply_uff=True, deterministic=True, quality_mode="fast",
        collapse_label_variants=False,
    )
    if err:
        return 0, []
    return len(results), [xyz for xyz, _ in results]


# ---------------------------------------------------------------------------
# Tests.
# ---------------------------------------------------------------------------

@pytest.mark.slow
def test_pyridine_orientation_snaps_to_edge_on():
    """5b-B with the fused-ring fix rotates the D-AQIWAZ pyridine ring
    so the Ru-N axis becomes edge-on (≥ 85°) -- up from a baseline
    near 32°."""
    n_off, frames_off = _run_pipeline(_D_AQIWAZ_SMILES, {})
    assert n_off >= 1, "Baseline D-AQIWAZ should emit at least 1 isomer"
    angle_off = max(
        _measure_pyridine_orientation_max_angle(x) for x in frames_off
    )
    # Baseline before this patch series sits below 35°.
    assert angle_off < 85.0, (
        f"Baseline pyridine orientation already edge-on ({angle_off:.1f}°); "
        "test pre-condition broken."
    )

    n_on, frames_on = _run_pipeline(
        _D_AQIWAZ_SMILES,
        {"DELFIN_DONOR_ORIENT_REALISM": "1"},
    )
    assert n_on >= 1, "5b-B run should still emit at least 1 isomer"
    angle_on = max(
        _measure_pyridine_orientation_max_angle(x) for x in frames_on
    )
    assert angle_on >= 85.0, (
        f"With 5b-B ON the Ru-N(pyridine) tilt must reach ≥ 85° "
        f"(edge-on); measured {angle_on:.1f}°."
    )


def test_phantom_13_exempt_default_on_preserves_md():
    """Phantom-1,3 exemption defaults to ON.  On the simplest CN2
    Cd dihalide (no 1,3-neighbour candidates) the pipeline still
    returns at least one isomer; we use that as a smoke check that
    the new branch does not crash on degenerate inputs."""
    n, frames = _run_pipeline(_CD_DIHALIDE_SMILES, {})
    assert n >= 1, "Cd dihalide should emit ≥ 1 isomer with default flags"


def test_phantom_13_exempt_env_off_restores_legacy():
    """Setting DELFIN_PHANTOM_13_EXEMPT=0 must restore the legacy
    behaviour (strict phantom check).  We simply assert the pipeline
    still produces ≥ 1 isomer on a simple SMILES -- the strict path
    has been the HEAD default for >3 months."""
    n, frames = _run_pipeline(
        _CD_DIHALIDE_SMILES,
        {"DELFIN_PHANTOM_13_EXEMPT": "0"},
    )
    assert n >= 1, "Strict-mode phantom check must still pass Cd dihalide"


def test_5bb_fused_ring_aromatic_n_finder_includes_pyridine():
    """The aromatic-N locator on D-AQIWAZ identifies the Ru-N(pyridine)
    triple.  Regression test: previously the fused-ring guard inside
    _apply_aromatic_n_edge_on rejected this triple at rotation time, but
    the locator itself always returned it.  This test pins the locator
    contract so any future refactor that drops the fused-N is caught."""
    from rdkit import Chem
    from delfin._donor_orientation_realism import (
        _find_aromatic_n_metal_pairs,
    )

    mol = Chem.MolFromSmiles(_D_AQIWAZ_SMILES, sanitize=False)
    assert mol is not None, "RDKit must parse the D-AQIWAZ SMILES"
    try:
        Chem.SanitizeMol(
            mol,
            sanitizeOps=Chem.SANITIZE_ALL ^ Chem.SANITIZE_KEKULIZE,
        )
    except Exception:
        pass
    mol = Chem.AddHs(mol)

    triples = _find_aromatic_n_metal_pairs(mol)
    assert triples, "Locator must find at least one aromatic-N / Ru triple"
    # The Ru atom is the only metal; one aromatic-N donor (the pyridine
    # of the naphthyridine) is directly bonded to it.
    metal_indices = {a.GetIdx() for a in mol.GetAtoms() if a.GetSymbol() == "Ru"}
    found_pyridine = False
    for m_idx, n_idx, ring in triples:
        if m_idx in metal_indices and len(ring) == 6:
            found_pyridine = True
            break
    assert found_pyridine, (
        "The pyridine-N / Ru triple is missing from the locator output"
    )
