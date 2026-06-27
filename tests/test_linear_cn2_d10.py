"""Tests for the linear-CN2 d10 special-case helpers (Wave 2026-05-13).

Covers ``delfin.manta._polyhedron_targets._is_d10_linear_cn2``,
``classify_geometry_from_cn_donors_with_metal`` and the new ``bent_2`` /
``linear_2`` ideal-vector tables.

Acceptance contract
-------------------

1. Au(I)/Ag(I)/Hg(II) CN=2 → ``linear_2`` (180°), regardless of formal
   charge value and regardless of the ``DELFIN_LINEAR_CN2`` env flag.
2. Cu(I)/Tl(I) CN=2 → ``linear_2`` only when formal charge matches the
   expected oxidation state (+1).  Other charges (Cu(II)/Tl(III)) fall
   back to ``bent_2`` when the env flag is enabled.
3. Non-d10 metals (Ti, V, Cr, Mn, Fe, …) at CN=2 → ``bent_2`` when env=1,
   otherwise ``linear_2`` (bit-exact legacy behaviour).
4. ``DELFIN_LINEAR_CN2`` unset / 0 → every CN=2 returns ``linear_2`` —
   identical to the pre-patch ``classify_geometry_from_cn_donors``.
5. CN!=2 ignores the env flag entirely.

The ideal-vector tables themselves are deterministic numpy arrays; we
check norm = 1 and the inter-vector angle (180° / ~104.5°).
"""
from __future__ import annotations

import math
import os

import numpy as np
import pytest

from delfin.manta._polyhedron_targets import (
    _D10_LINEAR_CN2_METALS_ALWAYS,
    _D10_LINEAR_CN2_METALS_CHARGED,
    _is_d10_linear_cn2,
    classify_geometry_from_cn_donors,
    classify_geometry_from_cn_donors_with_metal,
    get_ideal_donor_vectors,
    get_target_point_group,
)


# ---------------------------------------------------------------------------
# Helper – environment isolation
# ---------------------------------------------------------------------------


@pytest.fixture
def env_off(monkeypatch):
    """Default environment — DELFIN_LINEAR_CN2 unset → bit-exact OFF."""
    monkeypatch.delenv("DELFIN_LINEAR_CN2", raising=False)
    return None


@pytest.fixture
def env_on(monkeypatch):
    """DELFIN_LINEAR_CN2=1 → bent CN=2 enabled for non-d¹⁰ metals."""
    monkeypatch.setenv("DELFIN_LINEAR_CN2", "1")
    return None


# ---------------------------------------------------------------------------
# Section A — _is_d10_linear_cn2 predicate
# ---------------------------------------------------------------------------


def test_au_cn2_is_d10_linear():
    """Au CN=2 → d10 linear regardless of formal charge."""
    assert _is_d10_linear_cn2("Au", cn=2) is True
    assert _is_d10_linear_cn2("Au", cn=2, metal_formal_charge=+1) is True
    assert _is_d10_linear_cn2("Au", cn=2, metal_formal_charge=-1) is True
    assert _is_d10_linear_cn2("Au", cn=2, metal_formal_charge=0) is True


def test_ag_hg_cn2_is_d10_linear():
    """Ag and Hg at CN=2 are always d10 linear."""
    assert _is_d10_linear_cn2("Ag", cn=2) is True
    assert _is_d10_linear_cn2("Hg", cn=2) is True


def test_cu_i_cn2_is_d10_linear():
    """Cu(I) (+1 formal charge) is d10 → linear."""
    assert _is_d10_linear_cn2("Cu", cn=2, metal_formal_charge=+1) is True


def test_cu_ii_cn2_not_d10_linear():
    """Cu(II) is d9 → not in the linear set."""
    assert _is_d10_linear_cn2("Cu", cn=2, metal_formal_charge=+2) is False


def test_cu_unknown_charge_cn2_not_d10_linear():
    """Cu without explicit charge → conservative False (cannot prove d10)."""
    assert _is_d10_linear_cn2("Cu", cn=2) is False
    assert _is_d10_linear_cn2("Cu", cn=2, metal_formal_charge=None) is False


def test_tl_i_cn2_is_d10_linear():
    """Tl(I) (s2 post-transition) routes to linear too."""
    assert _is_d10_linear_cn2("Tl", cn=2, metal_formal_charge=+1) is True


def test_ti_cn2_not_d10_linear():
    """Ti is d0/d1/d2 — never in the linear CN=2 set."""
    assert _is_d10_linear_cn2("Ti", cn=2) is False
    assert _is_d10_linear_cn2("Ti", cn=2, metal_formal_charge=+4) is False


def test_cn_other_than_2_returns_false():
    """Predicate is CN=2 specific — CN=3/4/6 etc never matches."""
    for cn in (1, 3, 4, 5, 6, 7, 8, 9, 12):
        assert _is_d10_linear_cn2("Au", cn=cn) is False
        assert _is_d10_linear_cn2("Cu", cn=cn, metal_formal_charge=+1) is False


def test_empty_or_none_metal_symbol():
    """Missing metal symbol → False, no crash."""
    assert _is_d10_linear_cn2(None, cn=2) is False
    assert _is_d10_linear_cn2("", cn=2) is False


def test_metal_set_membership_documented():
    """Constants reflect the d10 metal set documented in the audit."""
    assert "Au" in _D10_LINEAR_CN2_METALS_ALWAYS
    assert "Ag" in _D10_LINEAR_CN2_METALS_ALWAYS
    assert "Hg" in _D10_LINEAR_CN2_METALS_ALWAYS
    assert _D10_LINEAR_CN2_METALS_CHARGED["Cu"] == +1
    assert _D10_LINEAR_CN2_METALS_CHARGED["Tl"] == +1


# ---------------------------------------------------------------------------
# Section B — classify_geometry_from_cn_donors_with_metal
# ---------------------------------------------------------------------------


def test_classifier_default_off_bit_exact_linear(env_off):
    """Default OFF: every CN=2 metal returns linear_2 — matches legacy."""
    for metal, q in [
        ("Au", +1), ("Ag", +1), ("Cu", +1), ("Cu", +2),
        ("Hg", +2), ("Tl", +1), ("Ti", +4), ("V", +5),
        ("Cr", +3), ("Pd", +2), ("Fe", +2), ("Ni", +2),
    ]:
        geom, iso = classify_geometry_from_cn_donors_with_metal(
            2, ["N", "N"], metal, q,
        )
        assert geom == "linear_2", (metal, q, geom)
        assert iso is None
        # bit-exact identical to legacy when env=0
        legacy_geom, legacy_iso = classify_geometry_from_cn_donors(
            2, ["N", "N"],
        )
        assert geom == legacy_geom
        assert iso == legacy_iso


def test_classifier_env_on_d10_metals_still_linear(env_on):
    """env=1: d10 metals still linear (Au, Ag, Hg unconditional; Cu(I)/Tl(I))."""
    for metal, q in [("Au", +1), ("Ag", +1), ("Hg", +2)]:
        geom, _ = classify_geometry_from_cn_donors_with_metal(
            2, ["N", "N"], metal, q,
        )
        assert geom == "linear_2", (metal, q, geom)
    # Cu(I) explicit
    geom, _ = classify_geometry_from_cn_donors_with_metal(
        2, ["N", "C"], "Cu", +1,
    )
    assert geom == "linear_2"
    # Tl(I)
    geom, _ = classify_geometry_from_cn_donors_with_metal(
        2, ["O", "O"], "Tl", +1,
    )
    assert geom == "linear_2"


def test_classifier_env_on_non_d10_bent(env_on):
    """env=1: non-d10 metals at CN=2 → bent_2 (~104°)."""
    for metal, q in [("Ti", +4), ("V", +3), ("Cr", +3), ("Mn", +2),
                     ("Fe", +3), ("Co", +2), ("Ni", +2), ("Pd", +2)]:
        geom, iso = classify_geometry_from_cn_donors_with_metal(
            2, ["N", "N"], metal, q,
        )
        assert geom == "bent_2", (metal, q, geom)
        assert iso is None


def test_classifier_env_on_cu_ii_bent(env_on):
    """env=1: Cu(II) is d9, NOT in d10 set → bent_2."""
    geom, _ = classify_geometry_from_cn_donors_with_metal(
        2, ["N", "N"], "Cu", +2,
    )
    assert geom == "bent_2"


def test_classifier_env_on_cu_unknown_charge_bent(env_on):
    """env=1 + Cu without charge info → conservative bent (cannot prove d10)."""
    geom, _ = classify_geometry_from_cn_donors_with_metal(
        2, ["N", "N"], "Cu", None,
    )
    assert geom == "bent_2"


def test_classifier_env_on_does_not_affect_higher_cn(env_on):
    """env=1 only changes CN=2 — higher CN unchanged."""
    geom, _ = classify_geometry_from_cn_donors_with_metal(
        4, ["N"] * 4, "Au", +3,
    )
    assert geom == "Td"
    geom, _ = classify_geometry_from_cn_donors_with_metal(
        6, ["N"] * 6, "Cu", +2,
    )
    assert geom == "Oh"


# ---------------------------------------------------------------------------
# Section C — ideal-vector tables (linear_2 vs bent_2)
# ---------------------------------------------------------------------------


def test_linear_2_vectors_are_180_degrees():
    """Linear CN=2 vectors must be exactly antiparallel (cos = -1, 180°)."""
    v = get_ideal_donor_vectors(2, "linear_2")
    assert v.shape == (2, 3)
    np.testing.assert_allclose(np.linalg.norm(v, axis=1), 1.0, atol=1e-9)
    cos = float(np.dot(v[0], v[1]))
    assert math.isclose(cos, -1.0, abs_tol=1e-9)


def test_bent_2_vectors_are_104p5_degrees():
    """Bent CN=2 (C2v) vectors at H2O-like ~104.5° aperture."""
    v = get_ideal_donor_vectors(2, "bent_2")
    assert v.shape == (2, 3)
    np.testing.assert_allclose(np.linalg.norm(v, axis=1), 1.0, atol=1e-9)
    cos = float(np.dot(v[0], v[1]))
    angle = math.degrees(math.acos(max(-1.0, min(1.0, cos))))
    assert math.isclose(angle, 104.5, abs_tol=1e-3)


def test_bent_2_alias_resolution():
    """Aliases ``bent`` / ``BENT`` resolve to the same table as ``bent_2``."""
    v1 = get_ideal_donor_vectors(2, "bent_2")
    v2 = get_ideal_donor_vectors(2, "bent")
    v3 = get_ideal_donor_vectors(2, "BENT")
    np.testing.assert_allclose(v1, v2, atol=1e-12)
    np.testing.assert_allclose(v1, v3, atol=1e-12)


def test_linear_2_alias_resolution():
    """Aliases ``linear`` / ``LIN`` resolve to ``linear_2`` table."""
    v1 = get_ideal_donor_vectors(2, "linear_2")
    v2 = get_ideal_donor_vectors(2, "linear")
    v3 = get_ideal_donor_vectors(2, "LIN")
    np.testing.assert_allclose(v1, v2, atol=1e-12)
    np.testing.assert_allclose(v1, v3, atol=1e-12)


def test_linear_vs_bent_point_groups():
    """Point-group classifications: D∞h linear, C2v bent."""
    assert get_target_point_group(2, "linear_2") == "D_inf_h"
    assert get_target_point_group(2, "bent_2") == "C2v"


# ---------------------------------------------------------------------------
# Section D — backward compatibility / default-OFF bit exactness
# ---------------------------------------------------------------------------


def test_legacy_classifier_unchanged(env_off):
    """Legacy classifier must always return linear_2 for CN=2."""
    for donors in [["N", "N"], ["C", "C"], ["O", "O"], ["P", "P"], ["N", "Cl"]]:
        geom, iso = classify_geometry_from_cn_donors(2, donors)
        assert geom == "linear_2"
        assert iso is None


def test_metal_aware_falls_through_to_legacy_for_cn3_to_12(env_on):
    """Metal-aware classifier delegates for CN != 2 — even with env=1."""
    for cn, donors in [
        (3, ["O"] * 3),
        (4, ["N"] * 4),
        (5, ["N"] * 5),
        (6, ["N"] * 6),
        (7, ["O"] * 7),
        (8, ["O"] * 8),
        (9, ["O"] * 9),
    ]:
        legacy = classify_geometry_from_cn_donors(cn, donors)
        metal_aware = classify_geometry_from_cn_donors_with_metal(
            cn, donors, "Au", +1,
        )
        assert legacy == metal_aware, (cn, legacy, metal_aware)


def test_default_off_for_real_pool_metals(env_off):
    """Bit-exact OFF on every CN=2 metal observed in master_v3 pool."""
    # Metals observed at CN=2 in
    # /home/qmchem_max/agent_workspace/quality_framework/results/
    # smiles_master_v3__DELFIN__20260513_122130.jsonl
    observed_metals = [
        "Au", "Ag", "Cu", "Hg", "Zn", "Cd", "Fe", "Ru", "V", "Pt",
        "Ti", "Mn", "Cr", "Pd", "Zr", "Sn", "Re", "Rh", "Ir", "Co", "Ni",
    ]
    for m in observed_metals:
        for q in (None, 0, +1, +2, +3, +4, +5):
            geom, _ = classify_geometry_from_cn_donors_with_metal(
                2, ["N", "N"], m, q,
            )
            assert geom == "linear_2", (m, q, geom)


# ---------------------------------------------------------------------------
# Section E — Integration with Tier-A donor target builder
# ---------------------------------------------------------------------------


def test_tier_a_integration_au_cn2_returns_linear():
    """Tier-A builder produces 180° targets for an Au CN=2 complex."""
    rdkit = pytest.importorskip("rdkit")
    from rdkit import Chem
    from rdkit.Chem import AllChem

    # Au(NH3)2+ analog written as bare donors to keep RDKit happy.
    smiles = "[NH3][Au+][NH3]"
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        pytest.skip("RDKit failed to parse Au(NH3)2 surrogate")
    mol_h = Chem.AddHs(mol)
    if AllChem.EmbedMolecule(mol_h, randomSeed=42) != 0:
        pytest.skip("Embedding failed in this RDKit build")
    # Force a deliberately bent geometry → expect Tier-A to push back to 180°.
    conf = mol_h.GetConformer()
    coords = conf.GetPositions().copy()
    metal_idx = next(
        a.GetIdx() for a in mol_h.GetAtoms() if a.GetSymbol() == "Au"
    )
    metal_pos = coords[metal_idx]
    donor_idxs = [
        n.GetIdx() for n in mol_h.GetAtomWithIdx(metal_idx).GetNeighbors()
    ]
    # Set donors at +x and +y (90° apart) to confirm we are NOT
    # accidentally reading the input geometry as 180°.
    coords[donor_idxs[0]] = metal_pos + np.array([1.5, 0.0, 0.0])
    coords[donor_idxs[1]] = metal_pos + np.array([0.0, 1.5, 0.0])

    from delfin.manta._symmetry_detection import hungarian_assign_donors_to_slots

    targets = hungarian_assign_donors_to_slots(coords, mol_h, metal_idx)
    assert len(targets) == 2
    # Compute vectors metal → target and verify 180° aperture.
    tgt_vecs = np.array([targets[d] - metal_pos for d in donor_idxs])
    norms = np.linalg.norm(tgt_vecs, axis=1)
    assert (norms > 0.5).all()
    unit = tgt_vecs / norms[:, None]
    cos = float(np.dot(unit[0], unit[1]))
    # Target unit vectors should be antiparallel within numerical noise.
    assert cos < -0.99, f"Au CN=2 Tier-A target not linear (cos={cos:.4f})"
