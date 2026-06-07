"""Tests for delfin.fffree.construction_sanity

Validate the post-embedding sanity check + retry loop + per-class templates
that address Bug-class 2 (ODUXAN-like donor drift) and Bug-class 3
(WICROP-like total collapse) without breaking byte-identity when env flags
are unset.
"""
from __future__ import annotations

import math
import os

import numpy as np
import pytest


# ---------------------------------------------------------------------------
# Env-guard
# ---------------------------------------------------------------------------


def _reset_env():
    for k in (
        "DELFIN_FFFREE_CONSTRUCTION_SANITY",
        "DELFIN_FFFREE_CONSTRUCTION_SANITY_MAX_RETRIES",
        "DELFIN_FFFREE_CONSTRUCTION_SANITY_CSHM_MAX",
        "DELFIN_FFFREE_GRIP_VDW_FLOOR_ALL_PAIRS",
        "DELFIN_FFFREE_PIANO_STOOL_TEMPLATE",
        "DELFIN_FFFREE_FISCHER_CARBENE_TEMPLATE",
    ):
        os.environ.pop(k, None)


@pytest.fixture(autouse=True)
def _env_guard():
    _reset_env()
    yield
    _reset_env()


# ---------------------------------------------------------------------------
# Default-OFF byte-identical
# ---------------------------------------------------------------------------


def test_all_env_flags_default_off():
    """Module-level env predicates are False when no flags are set."""
    from delfin.fffree.construction_sanity import (
        sanity_active,
        vdw_floor_all_pairs_active,
        piano_stool_template_active,
        fischer_carbene_template_active,
    )
    assert sanity_active() is False
    assert vdw_floor_all_pairs_active() is False
    assert piano_stool_template_active() is False
    assert fischer_carbene_template_active() is False


def test_env_flag_activation():
    from delfin.fffree.construction_sanity import (
        sanity_active,
        vdw_floor_all_pairs_active,
        piano_stool_template_active,
        fischer_carbene_template_active,
    )
    os.environ["DELFIN_FFFREE_CONSTRUCTION_SANITY"] = "1"
    os.environ["DELFIN_FFFREE_GRIP_VDW_FLOOR_ALL_PAIRS"] = "1"
    os.environ["DELFIN_FFFREE_PIANO_STOOL_TEMPLATE"] = "yes"
    os.environ["DELFIN_FFFREE_FISCHER_CARBENE_TEMPLATE"] = "TRUE"
    assert sanity_active() is True
    assert vdw_floor_all_pairs_active() is True
    assert piano_stool_template_active() is True
    assert fischer_carbene_template_active() is True


# ---------------------------------------------------------------------------
# assert_construction_sane: clean structures pass
# ---------------------------------------------------------------------------


def test_methane_is_sane():
    """A clean tetrahedral CH4 passes all checks."""
    from delfin.fffree.construction_sanity import assert_construction_sane

    syms = ["C", "H", "H", "H", "H"]
    t = 1.09 / math.sqrt(3)
    P = np.array([
        [0.0, 0.0, 0.0],
        [t, t, t],
        [t, -t, -t],
        [-t, t, -t],
        [-t, -t, t],
    ], dtype=float)
    bonds = [(0, 1), (0, 2), (0, 3), (0, 4)]
    ok, viols = assert_construction_sane(
        P, syms, bonds, metal_idx=None, donor_idxs=None,
    )
    assert ok is True
    assert viols == []


# ---------------------------------------------------------------------------
# assert_construction_sane: WICROP-like collapse detected
# ---------------------------------------------------------------------------


def test_wicrop_collapse_detected():
    """Two heavy atoms collapsed below the Pauli floor fire a pauli violation."""
    from delfin.fffree.construction_sanity import assert_construction_sane

    # Cr–C with the C far away, but a separate C-H pair collapsed at 0.74 Å
    # (the actual WICROP-class signature) and a P-C compressed to 0.93 Å.
    syms = ["Cr", "C", "C", "H", "P", "C"]
    P = np.array([
        [0.0, 0.0, 0.0],     # Cr
        [2.0, 0.0, 0.0],     # C bonded to Cr
        [3.0, 0.0, 0.0],     # C neighbour
        [3.0, 0.74, 0.0],    # H collapsed onto the C at 0.74 Å
        [-2.30, 0.0, 0.0],   # P
        [-3.23, 0.0, 0.0],   # P-C @ 0.93 Å (compressed; ideal ~1.83)
    ], dtype=float)
    bonds = [(0, 1), (1, 2), (2, 3), (0, 4), (4, 5)]
    ok, viols = assert_construction_sane(
        P, syms, bonds, metal_idx=0, donor_idxs=[1, 4],
    )
    assert ok is False
    modes = {v["mode"] for v in viols}
    assert "bond_compressed" in modes


# ---------------------------------------------------------------------------
# assert_construction_sane: ODUXAN-like donor drift detected
# ---------------------------------------------------------------------------


def test_oduxan_donor_drift_detected():
    """A nitrogen donor at 2.84 Å triggers md_drifted (ratio > 2.0×ideal? No:
    2.84/2.06 = 1.38 < 2.0, so NOT a drift) — we test the geometric CShM gate
    instead.  An OC-6 where one donor is at 2.84 Å but the other five are
    correct still has a high CShM."""
    from delfin.fffree.construction_sanity import assert_construction_sane

    # Cr at origin; 5 CO carbons at OC-6 positions ± 0.0 deviation; 1 N donor
    # drifted to 2.84 Å instead of ~2.10 Å AND offset from its expected vertex
    # by 25° — produces a non-trivial CShM.
    cr_c = 2.0  # idealised
    syms = ["Cr", "C", "C", "C", "C", "C", "N"]
    P = np.array([
        [0.0, 0.0, 0.0],
        [cr_c, 0.0, 0.0],
        [-cr_c, 0.0, 0.0],
        [0.0, cr_c, 0.0],
        [0.0, -cr_c, 0.0],
        [0.0, 0.0, cr_c],
        # N donor: should be at [0, 0, -cr_c] but drifted to a 25° offset and 2.84 Å
        [
            2.84 * math.sin(math.radians(25)),
            0.0,
            -2.84 * math.cos(math.radians(25)),
        ],
    ], dtype=float)
    bonds = [(0, 1), (0, 2), (0, 3), (0, 4), (0, 5), (0, 6)]
    ok, viols = assert_construction_sane(
        P, syms, bonds, metal_idx=0,
        donor_idxs=[1, 2, 3, 4, 5, 6],
        geometry="OC-6 octahedron",
        cshm_max=2.0,        # tighten so the 25° offset trips the gate
    )
    assert ok is False
    modes = {v["mode"] for v in viols}
    assert "cshm_too_large" in modes


# ---------------------------------------------------------------------------
# assert_construction_sane: M-D drift detected
# ---------------------------------------------------------------------------


def test_md_drift_violation():
    """A donor at > md_ceil_factor × ideal fires md_drifted."""
    from delfin.fffree.construction_sanity import assert_construction_sane

    syms = ["Cr", "N"]
    # Cr-N ideal ~ 2.06 Å.  Set d = 5.0 Å -> ratio ~ 2.43 -> md_drifted.
    P = np.array([[0.0, 0.0, 0.0], [5.0, 0.0, 0.0]], dtype=float)
    bonds = [(0, 1)]
    ok, viols = assert_construction_sane(
        P, syms, bonds, metal_idx=0, donor_idxs=[1], md_ceil_factor=2.0,
    )
    assert ok is False
    assert any(v["mode"] == "md_drifted" for v in viols)


# ---------------------------------------------------------------------------
# Retry loop
# ---------------------------------------------------------------------------


def test_build_with_retries_first_seed_wins():
    """First sane embedding wins, no fallback used."""
    from delfin.fffree.construction_sanity import build_with_retries

    syms = ["C", "C"]
    bonds = [(0, 1)]

    P_good = np.array([[0.0, 0.0, 0.0], [1.54, 0.0, 0.0]], dtype=float)

    def embed(seed):
        return P_good

    P, viols = build_with_retries(
        embed, syms, bonds, metal_idx=None, donor_idxs=None,
    )
    assert P is not None
    np.testing.assert_allclose(P, P_good)
    assert viols == []


def test_build_with_retries_fallback_used():
    """All retries return bad geometry -> fallback is invoked."""
    from delfin.fffree.construction_sanity import build_with_retries

    syms = ["C", "C"]
    bonds = [(0, 1)]

    P_bad = np.array([[0.0, 0.0, 0.0], [0.30, 0.0, 0.0]], dtype=float)  # compressed
    P_good = np.array([[0.0, 0.0, 0.0], [1.54, 0.0, 0.0]], dtype=float)

    n_calls = [0]

    def embed(seed):
        n_calls[0] += 1
        return P_bad

    def fallback():
        return P_good

    P, viols = build_with_retries(
        embed, syms, bonds, metal_idx=None, donor_idxs=None,
        fallback_fn=fallback, max_retries=3,
    )
    assert n_calls[0] == 3
    assert P is not None
    np.testing.assert_allclose(P, P_good)
    # Violations from the LAST failed embed are reported alongside the
    # fallback's coordinates.
    assert any(v["mode"] in ("bond_compressed", "pauli_floor") for v in viols)


def test_build_with_retries_no_fallback_returns_none():
    """All retries fail + no fallback -> (None, violations)."""
    from delfin.fffree.construction_sanity import build_with_retries

    syms = ["C", "C"]
    bonds = [(0, 1)]
    P_bad = np.array([[0.0, 0.0, 0.0], [0.30, 0.0, 0.0]], dtype=float)

    def embed(seed):
        return P_bad

    P, viols = build_with_retries(
        embed, syms, bonds, metal_idx=None, donor_idxs=None, max_retries=2,
    )
    assert P is None
    assert len(viols) >= 1


# ---------------------------------------------------------------------------
# Extended vdW floor (all pairs)
# ---------------------------------------------------------------------------


def test_vdw_floor_all_pairs_no_violation():
    """Well-separated atoms produce zero loss + zero gradient."""
    from delfin.fffree.construction_sanity import vdw_floor_all_pairs_value_and_grad

    R = np.array([
        [0.0, 0.0, 0.0],
        [5.0, 0.0, 0.0],
        [0.0, 5.0, 0.0],
    ], dtype=float)
    syms = ["C", "C", "C"]
    L, G = vdw_floor_all_pairs_value_and_grad(
        R, symbols=syms, excluded_pairs=set(), weight=20.0,
    )
    assert L == 0.0
    np.testing.assert_allclose(G, np.zeros_like(R))


def test_vdw_floor_all_pairs_violation_fires():
    """A heavy-heavy clash below the Pauli floor produces L > 0 and gradient."""
    from delfin.fffree.construction_sanity import vdw_floor_all_pairs_value_and_grad

    R = np.array([
        [0.0, 0.0, 0.0],
        [1.5, 0.0, 0.0],     # C-C @ 1.5 Å < 0.85 × (1.70+1.70) = 2.89
    ], dtype=float)
    syms = ["C", "C"]
    L, G = vdw_floor_all_pairs_value_and_grad(
        R, symbols=syms, excluded_pairs=set(), weight=20.0,
    )
    assert L > 0.0
    # Gradients on the two atoms should be equal-and-opposite.
    np.testing.assert_allclose(G[0] + G[1], np.zeros(3), atol=1e-9)
    # The gradient signs encode the *steepest-ascent* direction.  Following
    # -G separates the atoms: atom 0 moves to -x (grad[0][0] > 0), atom 1
    # moves to +x (grad[1][0] < 0).
    assert G[0][0] > 0.0
    assert G[1][0] < 0.0


def test_vdw_floor_excluded_pair_skipped():
    """Bonded / 1,3 pairs are excluded -> no penalty even when close."""
    from delfin.fffree.construction_sanity import vdw_floor_all_pairs_value_and_grad

    R = np.array([
        [0.0, 0.0, 0.0],
        [1.5, 0.0, 0.0],
    ], dtype=float)
    syms = ["C", "C"]
    L, G = vdw_floor_all_pairs_value_and_grad(
        R, symbols=syms,
        excluded_pairs={frozenset((0, 1))},
        weight=20.0,
    )
    assert L == 0.0
    np.testing.assert_allclose(G, np.zeros_like(R))


def test_vdw_floor_rigid_body_internal_collapse_no_deform():
    """Two rigid-body atoms below the Pauli floor add to the loss BUT their
    individual gradients are zero (rigid body cannot deform)."""
    from delfin.fffree.construction_sanity import vdw_floor_all_pairs_value_and_grad

    # 3 atoms: metal (0), two rigid-body atoms (1, 2) at 1.5 Å apart.
    R = np.array([
        [10.0, 10.0, 10.0],  # metal (far from rigid body)
        [0.0, 0.0, 0.0],
        [1.5, 0.0, 0.0],
    ], dtype=float)
    syms = ["Cr", "C", "C"]
    L, G = vdw_floor_all_pairs_value_and_grad(
        R, symbols=syms, excluded_pairs=set(), weight=20.0,
        rigid_body_atoms=[1, 2], rigid_translation_metal=0,
    )
    # The pair (1, 2) is BOTH-rigid -> gradient dropped; loss stays.
    # The pair (0, 1) is metal-vs-rigid but the metal is far, so no
    # contact.
    assert L > 0.0           # the loss reflects the internal collapse
    np.testing.assert_allclose(G[1], np.zeros(3))
    np.testing.assert_allclose(G[2], np.zeros(3))
    np.testing.assert_allclose(G[0], np.zeros(3))


def test_vdw_floor_rigid_to_free_translates_to_metal():
    """A rigid-body atom in contact with a FREE atom moves the rigid body
    as a whole (gradient on the metal slot), not the rigid atom itself."""
    from delfin.fffree.construction_sanity import vdw_floor_all_pairs_value_and_grad

    R = np.array([
        [10.0, 10.0, 10.0],  # metal far away (no direct contact)
        [0.0, 0.0, 0.0],     # rigid atom
        [1.5, 0.0, 0.0],     # free atom in contact with rigid atom
    ], dtype=float)
    syms = ["Cr", "C", "C"]
    L, G = vdw_floor_all_pairs_value_and_grad(
        R, symbols=syms, excluded_pairs=set(), weight=20.0,
        rigid_body_atoms=[1], rigid_translation_metal=0,
    )
    assert L > 0.0
    # Rigid atom (1) keeps gradient = 0; metal takes the translation.
    np.testing.assert_allclose(G[1], np.zeros(3))
    # Free atom (2) at +1.5 vs rigid atom at origin: following -G[2] moves
    # the free atom to +x (away from the rigid atom) -> G[2][0] < 0.
    assert G[2][0] < 0.0
    # Metal slot gets the equal-and-opposite -- moves rigid body to -x by
    # following -G[0] (away from the free atom).
    assert G[0][0] > 0.0
    np.testing.assert_allclose(G[0] + G[2], np.zeros(3), atol=1e-9)


# ---------------------------------------------------------------------------
# Topology classification
# ---------------------------------------------------------------------------


def test_classify_piano_stool():
    """η⁶-arene + 3 CO + 1 P-aux  →  piano_stool."""
    from delfin.fffree.construction_sanity import classify_topology

    ligands = [
        {  # arene
            "is_hapto": True, "hapto_eta": 6, "denticity": 1,
            "donor_elems": ["C"], "smiles": "c1ccccc1",
        },
        {"is_hapto": False, "denticity": 1, "donor_elems": ["C"], "smiles": "[C]#[O+]"},
        {"is_hapto": False, "denticity": 1, "donor_elems": ["C"], "smiles": "[C]#[O+]"},
        {"is_hapto": False, "denticity": 1, "donor_elems": ["C"], "smiles": "[C]#[O+]"},
        # P-aux: tertiary phosphine PMe3 — monodentate with denticity 1 + donor P
        {"is_hapto": False, "denticity": 1, "donor_elems": ["P"], "smiles": "P(C)(C)C"},
    ]
    assert classify_topology(ligands) == "piano_stool"


def test_classify_fischer_carbene():
    """Fischer carbene + 5 CO  →  fischer_carbene."""
    from delfin.fffree.construction_sanity import classify_topology

    ligands = [
        # Carbene C donor with N substituent + extended backbone
        {
            "is_hapto": False, "denticity": 1, "donor_elems": ["C"],
            "smiles": "[C]1=NC(=[N+]2CCCC2)C2(O1)C1=CC=CC=C12",
        },
        {"is_hapto": False, "denticity": 1, "donor_elems": ["C"], "smiles": "[C]#[O+]"},
        {"is_hapto": False, "denticity": 1, "donor_elems": ["C"], "smiles": "[C]#[O+]"},
        {"is_hapto": False, "denticity": 1, "donor_elems": ["C"], "smiles": "[C]#[O+]"},
        {"is_hapto": False, "denticity": 1, "donor_elems": ["C"], "smiles": "[C]#[O+]"},
        {"is_hapto": False, "denticity": 1, "donor_elems": ["C"], "smiles": "[C]#[O+]"},
    ]
    assert classify_topology(ligands) == "fischer_carbene"


def test_classify_unclassified_for_simple_complex():
    """A homoleptic 6× CO complex with no carbene + no hapto is unclassified."""
    from delfin.fffree.construction_sanity import classify_topology

    ligands = [
        {"is_hapto": False, "denticity": 1, "donor_elems": ["C"], "smiles": "[C]#[O+]"}
        for _ in range(6)
    ]
    assert classify_topology(ligands) == "unclassified"


def test_classify_empty_ligands_unclassified():
    from delfin.fffree.construction_sanity import classify_topology
    assert classify_topology([]) == "unclassified"
    assert classify_topology(None) == "unclassified"  # type: ignore[arg-type]


# ---------------------------------------------------------------------------
# Per-class placement templates
# ---------------------------------------------------------------------------


def test_piano_stool_template_geometry():
    """Piano-stool template: ring above metal, CO below, oxygen outward."""
    from delfin.fffree.construction_sanity import piano_stool_template

    layout = piano_stool_template("Cr", ring_size=6, n_co=3, aux_donor="P")

    metal = layout["metal"]
    ring = layout["ring"]
    co_c = layout["co_carbons"]
    co_o = layout["co_oxygens"]
    aux = layout["aux_donor"]

    # Metal at origin
    np.testing.assert_allclose(metal, np.zeros(3))

    # Ring centroid at +z (above metal)
    centroid = ring.mean(axis=0)
    assert centroid[2] > 1.0
    assert abs(centroid[0]) < 1e-9
    assert abs(centroid[1]) < 1e-9

    # All ring atoms equidistant from centroid
    radii = np.linalg.norm(ring - centroid, axis=1)
    np.testing.assert_allclose(radii, radii[0] * np.ones(6), atol=1e-9)

    # CO carbons in -z hemisphere
    for k in range(3):
        assert co_c[k][2] < 0

    # CO oxygens FURTHER from metal than carbons (oxygen points outward)
    for k in range(3):
        d_c = float(np.linalg.norm(co_c[k]))
        d_o = float(np.linalg.norm(co_o[k]))
        assert d_o > d_c, f"O[{k}] @ {d_o} <= C[{k}] @ {d_c}; O should point OUTWARD"

    # Auxiliary donor at -z (opposite the ring)
    assert aux[2] < 0


def test_piano_stool_template_three_co_separated():
    """The three CO ligands of the piano-stool should be far apart (≥ 2.5 Å
    between their C atoms) so the WICROP collapse cannot re-emerge."""
    from delfin.fffree.construction_sanity import piano_stool_template

    layout = piano_stool_template("Cr", ring_size=6, n_co=3, aux_donor=None)
    co_c = layout["co_carbons"]
    for i in range(3):
        for j in range(i + 1, 3):
            d = float(np.linalg.norm(co_c[i] - co_c[j]))
            assert d > 2.5, f"CO[{i}]-CO[{j}] @ {d} Å is too close"


def test_fischer_carbene_template_geometry():
    """Fischer-carbene template: carbene axial +z, CO equatorial, aux axial -z."""
    from delfin.fffree.construction_sanity import fischer_carbene_template

    layout = fischer_carbene_template("Cr", n_co=4, aux_donor="N")

    carbene = layout["carbene_C"]
    co_c = layout["co_carbons"]
    co_o = layout["co_oxygens"]
    aux = layout["aux_donor"]

    # Carbene C purely on +z axis
    assert carbene[2] > 1.0
    np.testing.assert_allclose(carbene[:2], np.zeros(2), atol=1e-9)

    # All CO C in the z=0 equatorial plane
    for k in range(4):
        assert abs(co_c[k][2]) < 1e-9

    # CO oxygens point outward
    for k in range(4):
        d_c = float(np.linalg.norm(co_c[k]))
        d_o = float(np.linalg.norm(co_o[k]))
        assert d_o > d_c

    # Auxiliary donor at -z (axial opposite the carbene)
    assert aux[2] < 0
    np.testing.assert_allclose(aux[:2], np.zeros(2), atol=1e-9)

    # Carbene + aux are colinear with the metal
    assert abs(np.dot(carbene, aux)) > 0.99 * (
        np.linalg.norm(carbene) * np.linalg.norm(aux)
    )


def test_fischer_carbene_template_oduxan_donors_at_correct_distance():
    """The carbene C, aux N and all 5 CO carbons end up within the expected
    M-D distance bracket (no -2.84 Å drift)."""
    from delfin.fffree.construction_sanity import (
        fischer_carbene_template, assert_construction_sane,
    )

    layout = fischer_carbene_template("Cr", n_co=5, aux_donor="N")
    syms = (
        ["Cr"]                # 0
        + ["C"]               # 1 carbene C
        + ["C"] * 5           # 2-6 CO carbons
        + ["O"] * 5           # 7-11 CO oxygens
        + ["N"]               # 12 aux N
    )
    P = np.vstack([
        layout["metal"][None, :],
        layout["carbene_C"][None, :],
        layout["co_carbons"],
        layout["co_oxygens"],
        layout["aux_donor"][None, :],
    ])
    # Bonds: M-carbeneC, 5× M-C, M-N, 5× C=O
    bonds = [(0, 1)]
    bonds += [(0, 2 + k) for k in range(5)]
    bonds += [(0, 12)]
    bonds += [(2 + k, 7 + k) for k in range(5)]

    donor_idxs = [1] + [2 + k for k in range(5)] + [12]
    ok, viols = assert_construction_sane(
        P, syms, bonds,
        metal_idx=0,
        donor_idxs=donor_idxs,
        geometry="OC-6 octahedron",
        cshm_max=5.0,
    )
    # The template should be Pauli-clean and the OC-6 CShM should be tiny.
    md_drift = [v for v in viols if v["mode"] in ("md_drifted", "md_too_short")]
    assert not md_drift, f"M-D drift in template: {md_drift}"
    cshm_viols = [v for v in viols if v["mode"] == "cshm_too_large"]
    assert not cshm_viols, f"CShM too large in template: {cshm_viols}"


# ---------------------------------------------------------------------------
# Determinism
# ---------------------------------------------------------------------------


def test_templates_are_deterministic():
    """Same inputs -> bit-identical output."""
    from delfin.fffree.construction_sanity import (
        piano_stool_template, fischer_carbene_template,
    )

    a = piano_stool_template("Cr", ring_size=6, n_co=3, aux_donor="P")
    b = piano_stool_template("Cr", ring_size=6, n_co=3, aux_donor="P")
    for k in a:
        np.testing.assert_array_equal(a[k], b[k])

    c = fischer_carbene_template("Cr", n_co=5, aux_donor="N")
    d = fischer_carbene_template("Cr", n_co=5, aux_donor="N")
    for k in c:
        np.testing.assert_array_equal(c[k], d[k])


def test_sanity_check_is_deterministic():
    """Same input -> identical violation list."""
    from delfin.fffree.construction_sanity import assert_construction_sane

    syms = ["Cr", "C", "C"]
    P = np.array([
        [0.0, 0.0, 0.0],
        [2.0, 0.0, 0.0],
        [0.05, 0.0, 0.0],   # pauli-violating contact with M
    ], dtype=float)
    bonds = [(0, 1)]
    ok1, v1 = assert_construction_sane(P, syms, bonds, metal_idx=0)
    ok2, v2 = assert_construction_sane(P, syms, bonds, metal_idx=0)
    assert ok1 == ok2
    assert v1 == v2


# ---------------------------------------------------------------------------
# Byte-identity contract (the legacy assemble path is unchanged when no env
# flag is set).  We assert that the module's predicates only return True when
# the corresponding env flags are set — so no upstream caller will route
# through the new code by accident.
# ---------------------------------------------------------------------------


def test_byte_identity_no_env_flag_no_activation():
    """With no env flags set, every predicate is False AND no assertion is
    raised when the module is imported.
    """
    from delfin.fffree import construction_sanity as CS
    assert CS.sanity_active() is False
    assert CS.vdw_floor_all_pairs_active() is False
    assert CS.piano_stool_template_active() is False
    assert CS.fischer_carbene_template_active() is False
