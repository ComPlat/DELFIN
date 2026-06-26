"""Unit tests for ``delfin/_e6761e4_soft_donor.py`` (Welle-5l Track-4).

Covers the four contract points:

1. ``detect_sp3_c_donor_feature`` flags M-CH₂-X (WUXQAK pattern).
2. ``detect_terminal_ch3_on_donor_feature`` flags N-Me (29-Ni pattern).
3. ``is_feature_gate_enabled`` is default-OFF (bit-exact preservation).
4. ``donor_weight_for_atom`` clamps to [0, 1] and honours env override.
"""

from __future__ import annotations

import os
from unittest import mock

import pytest


pytest.importorskip("rdkit")


def _mol(smiles: str):
    from rdkit import Chem
    m = Chem.MolFromSmiles(smiles)
    assert m is not None, f"RDKit failed to parse {smiles!r}"
    return Chem.AddHs(m)


def test_detect_sp3_c_donor_feature_positive_wuxqak_like():
    """An Hf bearing three CH₂-CH₃ σ-donors must be flagged.

    Mirrors the WUXQAK pattern (Hf with three sp3-C donors, M-C-X linear
    collapse) documented in feedback_sp3_c_donor_linear.
    """
    from delfin._e6761e4_soft_donor import detect_sp3_c_donor_feature
    # Tetra-ethyl hafnium (Hf-CH₂-CH₃ ×4): all four σ-donors are sp3-C.
    m = _mol("[Hf](CC)(CC)(CC)CC")
    assert detect_sp3_c_donor_feature(m) is True


def test_detect_sp3_c_donor_feature_negative_pure_halide():
    """A purely halide-coordinated complex must NOT be flagged.

    No σ-donor here is a carbon, so the feature gate stays OFF — this
    pins the negative branch of ``detect_sp3_c_donor_feature``.
    """
    from delfin._e6761e4_soft_donor import detect_sp3_c_donor_feature
    # cisplatin-like Pt(NH3)2Cl2: donors are N (sp³) and Cl — no C donor.
    m = _mol("[Pt](Cl)(Cl)(N)N")
    assert detect_sp3_c_donor_feature(m) is False


def test_detect_terminal_ch3_on_donor_feature_positive_29ni_like():
    """An N-methyl-imidazolyl Ni σ-donor must be flagged.

    Mirrors the 29-Ni pincer-tBu-imid gold-test (CH₃ umbrella inversion
    in 12/12 frames per User test-case archive).
    """
    from delfin._e6761e4_soft_donor import (
        detect_terminal_ch3_on_donor_feature,
    )
    # Ni-NMe₃ (trimethylamine σ-donor): N donor bonded to three terminal CH₃.
    m = _mol("[Ni](Cl)(Cl)N(C)(C)C")
    assert detect_terminal_ch3_on_donor_feature(m) is True


def test_is_feature_gate_enabled_default_on_sigma_class():
    """Iter-19 (2026-05-19): default flipped 0 → 1 for {sigma, multi_sigma}.

    Feature-positive sigma-class molecule, no env-overrides → True
    (per cross-archive analysis: e6761e4 wins F19/F24/pi_planar/funcgrp_bond).
    """
    from delfin._e6761e4_soft_donor import is_feature_gate_enabled
    with mock.patch.dict(
        os.environ,
        {k: v for k, v in os.environ.items()
         if not k.startswith("DELFIN_5L_T4_")},
        clear=True,
    ):
        # Feature-positive sigma class → True (default ON for sigma).
        m = _mol("[Hf](CC)(CC)(CC)CC")
        assert is_feature_gate_enabled(m) is True

    # Explicit GATE=0 still disables entirely.
    with mock.patch.dict(
        os.environ,
        {"DELFIN_5L_T4_SOFT_DONOR_FEATURE_GATE": "0"},
        clear=False,
    ):
        m = _mol("[Hf](CC)(CC)(CC)CC")
        assert is_feature_gate_enabled(m) is False


def test_is_feature_gate_enabled_on_with_feature():
    """With flag ON and feature present, the gate fires."""
    from delfin._e6761e4_soft_donor import is_feature_gate_enabled
    with mock.patch.dict(
        os.environ,
        {"DELFIN_5L_T4_SOFT_DONOR_FEATURE_GATE": "1"},
        clear=False,
    ):
        m = _mol("[Hf](CC)(CC)(CC)CC")
        assert is_feature_gate_enabled(m) is True


def test_donor_weight_for_atom_default_and_clamp():
    """``donor_weight_for_atom`` defaults to 0.3 and clamps to [0, 1]."""
    from delfin._e6761e4_soft_donor import donor_weight_for_atom
    m = _mol("[Hf](CC)(CC)(CC)CC")
    # Default base weight = 0.3.
    with mock.patch.dict(os.environ, {}, clear=False):
        # Strip any pre-set override
        os.environ.pop("DELFIN_5L_T4_SOFT_DONOR_WEIGHT", None)
        w = donor_weight_for_atom(m, 0)
        assert abs(w - 0.3) < 1e-9, f"default weight 0.3, got {w}"
    # Env override valid.
    with mock.patch.dict(
        os.environ,
        {"DELFIN_5L_T4_SOFT_DONOR_WEIGHT": "0.7"},
        clear=False,
    ):
        w = donor_weight_for_atom(m, 0)
        assert abs(w - 0.7) < 1e-9, f"env override 0.7, got {w}"
    # Above 1.0 clamps to 1.0.
    with mock.patch.dict(
        os.environ,
        {"DELFIN_5L_T4_SOFT_DONOR_WEIGHT": "1.5"},
        clear=False,
    ):
        w = donor_weight_for_atom(m, 0)
        assert abs(w - 1.0) < 1e-9, f"clamp upper 1.0, got {w}"
    # Below 0.0 clamps to 0.0.
    with mock.patch.dict(
        os.environ,
        {"DELFIN_5L_T4_SOFT_DONOR_WEIGHT": "-0.2"},
        clear=False,
    ):
        w = donor_weight_for_atom(m, 0)
        assert abs(w - 0.0) < 1e-9, f"clamp lower 0.0, got {w}"
