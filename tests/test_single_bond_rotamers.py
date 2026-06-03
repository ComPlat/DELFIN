"""Tests for delfin.fffree.single_bond_rotamers — deterministic single-bond
rotamer enumeration (gauche+/anti/gauche-).
"""
from __future__ import annotations

import os

import numpy as np
import pytest


def _build_mol(smiles: str):
    """Helper: SMILES -> AddHs + ETKDG embed, returns (mol, coords)."""
    from rdkit import Chem
    from rdkit.Chem import AllChem
    m = Chem.AddHs(Chem.MolFromSmiles(smiles))
    res = AllChem.EmbedMolecule(m, randomSeed=42)
    assert res == 0, f"Embed failed for {smiles}"
    P = np.array(m.GetConformer().GetPositions(), dtype=float)
    return m, P


# ----------------------------------------------------------------------
# find_rotatable_bonds
# ----------------------------------------------------------------------


def test_find_rotatable_bonds_butane_one_rotor():
    """Butane (CCCC) has a single rotatable bond (the central C-C)."""
    from delfin.fffree.single_bond_rotamers import find_rotatable_bonds
    m, _ = _build_mol("CCCC")
    rotors = find_rotatable_bonds(m)
    assert len(rotors) == 1, f"butane expected 1 rotor, got {len(rotors)}"


def test_find_rotatable_bonds_pentane_two_rotors():
    """Pentane (CCCCC) -> two rotatable central C-C bonds."""
    from delfin.fffree.single_bond_rotamers import find_rotatable_bonds
    m, _ = _build_mol("CCCCC")
    rotors = find_rotatable_bonds(m)
    assert len(rotors) == 2, f"pentane expected 2 rotors, got {len(rotors)}"


def test_find_rotatable_bonds_methane_none():
    """Methane has no heavy-heavy bonds -> zero rotors."""
    from delfin.fffree.single_bond_rotamers import find_rotatable_bonds
    m, _ = _build_mol("C")
    rotors = find_rotatable_bonds(m)
    assert rotors == [], "methane should have no rotatable bonds"


def test_find_rotatable_bonds_benzene_none_ring():
    """Benzene: aromatic + in-ring -> all bonds excluded."""
    from delfin.fffree.single_bond_rotamers import find_rotatable_bonds
    m, _ = _build_mol("c1ccccc1")
    rotors = find_rotatable_bonds(m)
    assert rotors == [], "benzene should have no rotatable bonds"


def test_find_rotatable_bonds_ethane_excluded():
    """Ethane: both endpoints terminal (no heavy neighbours other than
    each other) -> excluded."""
    from delfin.fffree.single_bond_rotamers import find_rotatable_bonds
    m, _ = _build_mol("CC")
    rotors = find_rotatable_bonds(m)
    assert rotors == [], "ethane should be excluded (no distinct rotamer)"


def test_find_rotatable_bonds_deterministic_ordering():
    """find_rotatable_bonds must return a sorted list (deterministic)."""
    from delfin.fffree.single_bond_rotamers import find_rotatable_bonds
    m, _ = _build_mol("CCCCCC")  # hexane
    rotors = find_rotatable_bonds(m)
    assert rotors == sorted(rotors), "rotor list not sorted"


def test_find_rotatable_bonds_metal_excluded():
    """Synthetic Co-C-C-C: the Co-C bond must NOT count as rotatable.
    Constructed via RDKit RWMol because SMILES for metals is brittle."""
    from rdkit import Chem
    from delfin.fffree.single_bond_rotamers import find_rotatable_bonds
    rw = Chem.RWMol()
    rw.AddAtom(Chem.Atom(27))   # Co  (idx 0)
    rw.AddAtom(Chem.Atom(6))    # C   (idx 1)
    rw.AddAtom(Chem.Atom(6))    # C   (idx 2)
    rw.AddAtom(Chem.Atom(6))    # C   (idx 3)
    rw.AddBond(0, 1, Chem.BondType.SINGLE)
    rw.AddBond(1, 2, Chem.BondType.SINGLE)
    rw.AddBond(2, 3, Chem.BondType.SINGLE)
    mol = rw.GetMol()
    rotors = find_rotatable_bonds(mol)
    # Allowed: (1,2), (2,3) is terminal -> excluded;
    # (0,1) metal -> excluded.  Only (1,2) survives.
    assert (0, 1) not in rotors, "metal-incident bond must be excluded"


# ----------------------------------------------------------------------
# enumerate_rotamer_configs
# ----------------------------------------------------------------------


def test_enumerate_rotamer_configs_zero_rotors():
    from delfin.fffree.single_bond_rotamers import enumerate_rotamer_configs
    configs = list(enumerate_rotamer_configs(0))
    assert configs == [()], "zero rotors should yield exactly identity tuple"


def test_enumerate_rotamer_configs_count_3_rotors():
    """3 rotors x 3 states = 27 configs."""
    from delfin.fffree.single_bond_rotamers import enumerate_rotamer_configs
    configs = list(enumerate_rotamer_configs(3, n_states=3, max_configs=243))
    assert len(configs) == 27, f"expected 27 got {len(configs)}"


def test_enumerate_rotamer_configs_first_is_identity():
    """The first config MUST be all-zero offsets (= identity = as-built)."""
    from delfin.fffree.single_bond_rotamers import enumerate_rotamer_configs
    cfg0 = next(enumerate_rotamer_configs(4, n_states=3, max_configs=243))
    assert cfg0 == (0.0, 0.0, 0.0, 0.0), f"first config must be identity, got {cfg0}"


def test_enumerate_rotamer_configs_cap_honored():
    from delfin.fffree.single_bond_rotamers import enumerate_rotamer_configs
    configs = list(enumerate_rotamer_configs(5, n_states=3, max_configs=10))
    assert len(configs) == 10


def test_enumerate_rotamer_configs_deterministic():
    """Same call twice -> identical sequence."""
    from delfin.fffree.single_bond_rotamers import enumerate_rotamer_configs
    a = list(enumerate_rotamer_configs(3, n_states=3))
    b = list(enumerate_rotamer_configs(3, n_states=3))
    assert a == b


# ----------------------------------------------------------------------
# enumerate_single_bond_rotamers (full pipeline)
# ----------------------------------------------------------------------


def test_enumerate_full_butane_first_is_input():
    """First yielded variant equals the input coordinates."""
    from delfin.fffree.single_bond_rotamers import enumerate_single_bond_rotamers
    m, P = _build_mol("CCCC")
    out = list(enumerate_single_bond_rotamers(m, P, max_configs=27))
    assert len(out) >= 1
    P0_variant, lab0 = out[0]
    assert np.allclose(P0_variant, P, atol=1e-8), \
        "first variant must equal input (identity offset)"
    assert "0" in lab0


def test_enumerate_full_pentane_count():
    """Pentane: 2 rotors x 3 states = 9 variants (no cap).

    Pinned to n_states=3 (legacy ±120° step) so this regression test
    survives the subagent #129 finer-grid default change to ±60° (6 states).
    """
    from delfin.fffree.single_bond_rotamers import enumerate_single_bond_rotamers
    m, P = _build_mol("CCCCC")
    out = list(enumerate_single_bond_rotamers(m, P, max_configs=243, n_states=3))
    assert len(out) == 9, f"pentane expected 9 variants, got {len(out)}"


def test_enumerate_full_pentane_count_default_step():
    """Pentane: 2 rotors x 6 states = 36 variants under the new default
    ±60° step (subagent #129 finer-grid default).
    """
    from delfin.fffree.single_bond_rotamers import enumerate_single_bond_rotamers
    m, P = _build_mol("CCCCC")
    out = list(enumerate_single_bond_rotamers(m, P, max_configs=243))
    assert len(out) == 36, f"pentane default step expected 36 variants, got {len(out)}"


def test_enumerate_full_methane_yields_identity():
    """No rotors -> yield just the input as identity."""
    from delfin.fffree.single_bond_rotamers import enumerate_single_bond_rotamers
    m, P = _build_mol("C")
    out = list(enumerate_single_bond_rotamers(m, P))
    assert len(out) == 1
    assert np.allclose(out[0][0], P, atol=1e-8)


def test_enumerate_full_changes_geometry():
    """A non-identity offset must actually MOVE atoms."""
    from delfin.fffree.single_bond_rotamers import enumerate_single_bond_rotamers
    m, P = _build_mol("CCCCC")
    out = list(enumerate_single_bond_rotamers(m, P, max_configs=243))
    # The last variant should differ from the first (unless degenerate).
    Plast, _ = out[-1]
    assert not np.allclose(Plast, P, atol=1e-3), \
        "non-trivial rotamer should not equal input"


def test_enumerate_full_max_rotors_cap():
    """max_rotors=1 on hexane (2 rotors) should yield only 3 variants under
    the legacy n_states=3 (±120°) step.  Subagent #129 default ±60° step
    -> 6 variants; the n_states=3 pin keeps this test focused on the
    max_rotors-cap behaviour.
    """
    from delfin.fffree.single_bond_rotamers import enumerate_single_bond_rotamers
    m, P = _build_mol("CCCCCC")
    out = list(enumerate_single_bond_rotamers(
        m, P, max_rotors=1, max_configs=243, n_states=3,
    ))
    assert len(out) == 3, f"expected 3 with 1 rotor, got {len(out)}"


def test_enumerate_full_max_rotors_cap_default_step():
    """max_rotors=1 on hexane under the new default ±60° step -> 6 variants."""
    from delfin.fffree.single_bond_rotamers import enumerate_single_bond_rotamers
    m, P = _build_mol("CCCCCC")
    out = list(enumerate_single_bond_rotamers(m, P, max_rotors=1, max_configs=243))
    assert len(out) == 6, f"expected 6 with 1 rotor (default step), got {len(out)}"


def test_enumerate_full_deterministic():
    """Two calls with same input -> identical variant sequence."""
    from delfin.fffree.single_bond_rotamers import enumerate_single_bond_rotamers
    m, P = _build_mol("CCCC")
    a = list(enumerate_single_bond_rotamers(m, P, max_configs=27))
    b = list(enumerate_single_bond_rotamers(m, P, max_configs=27))
    assert len(a) == len(b)
    for (Pa, la), (Pb, lb) in zip(a, b):
        assert la == lb
        assert np.allclose(Pa, Pb, atol=1e-9)


def test_enumerate_full_labels_distinct():
    """Each yielded variant has a distinct rotamer label."""
    from delfin.fffree.single_bond_rotamers import enumerate_single_bond_rotamers
    m, P = _build_mol("CCCCC")
    out = list(enumerate_single_bond_rotamers(m, P, max_configs=243))
    labels = [lab for _, lab in out]
    assert len(labels) == len(set(labels)), "labels must be unique"


# ----------------------------------------------------------------------
# Subagent #129 follow-up: step-degrees env-flag + new defaults
# ----------------------------------------------------------------------


def test_step_degrees_env_default_60():
    """Env-flag DELFIN_FFFREE_ROTAMER_STEP_DEGREES default = 60.0."""
    from delfin.fffree.single_bond_rotamers import _step_degrees, _n_states
    saved = os.environ.pop("DELFIN_FFFREE_ROTAMER_STEP_DEGREES", None)
    os.environ.pop("DELFIN_FFFREE_ENUMERATE_ROTAMERS_STATES", None)
    try:
        assert _step_degrees() == 60.0
        assert _n_states() == 6
    finally:
        if saved is not None:
            os.environ["DELFIN_FFFREE_ROTAMER_STEP_DEGREES"] = saved


def test_step_degrees_env_override():
    """Setting the env-flag changes the derived n_states."""
    from delfin.fffree.single_bond_rotamers import _step_degrees, _n_states
    saved = os.environ.get("DELFIN_FFFREE_ROTAMER_STEP_DEGREES")
    saved_states = os.environ.pop("DELFIN_FFFREE_ENUMERATE_ROTAMERS_STATES", None)
    try:
        os.environ["DELFIN_FFFREE_ROTAMER_STEP_DEGREES"] = "30"
        assert abs(_step_degrees() - 30.0) < 1e-9
        assert _n_states() == 12
        os.environ["DELFIN_FFFREE_ROTAMER_STEP_DEGREES"] = "120"
        assert _n_states() == 3
    finally:
        if saved is not None:
            os.environ["DELFIN_FFFREE_ROTAMER_STEP_DEGREES"] = saved
        else:
            os.environ.pop("DELFIN_FFFREE_ROTAMER_STEP_DEGREES", None)
        if saved_states is not None:
            os.environ["DELFIN_FFFREE_ENUMERATE_ROTAMERS_STATES"] = saved_states


def test_max_configs_env_default_144():
    """Env-flag DELFIN_FFFREE_ENUMERATE_ROTAMERS_MAX_CONFIGS default = 144."""
    from delfin.fffree.single_bond_rotamers import _max_configs
    saved = os.environ.pop("DELFIN_FFFREE_ENUMERATE_ROTAMERS_MAX_CONFIGS", None)
    try:
        assert _max_configs() == 144
    finally:
        if saved is not None:
            os.environ["DELFIN_FFFREE_ENUMERATE_ROTAMERS_MAX_CONFIGS"] = saved


def test_enumerate_step_degrees_kwarg_overrides_env():
    """Explicit step_degrees kwarg trumps the env-default."""
    from delfin.fffree.single_bond_rotamers import enumerate_single_bond_rotamers
    m, P = _build_mol("CCCC")  # 1 rotor
    out120 = list(enumerate_single_bond_rotamers(
        m, P, step_degrees=120.0, max_configs=243,
    ))
    out60 = list(enumerate_single_bond_rotamers(
        m, P, step_degrees=60.0, max_configs=243,
    ))
    assert len(out120) == 3
    assert len(out60) == 6


def test_enumerate_six_states_pentane_36():
    """Pentane 2-rotor under 60° step -> 6**2 = 36 variants."""
    from delfin.fffree.single_bond_rotamers import enumerate_single_bond_rotamers
    m, P = _build_mol("CCCCC")
    out = list(enumerate_single_bond_rotamers(
        m, P, step_degrees=60.0, max_configs=144,
    ))
    assert len(out) == 36


def test_enumerate_step_degrees_invalid_falls_back():
    """Malformed env -> fall back to default step (no crash)."""
    from delfin.fffree.single_bond_rotamers import _step_degrees
    saved = os.environ.get("DELFIN_FFFREE_ROTAMER_STEP_DEGREES")
    try:
        os.environ["DELFIN_FFFREE_ROTAMER_STEP_DEGREES"] = "garbage"
        # Should not raise.
        v = _step_degrees()
        assert v == 60.0
    finally:
        if saved is not None:
            os.environ["DELFIN_FFFREE_ROTAMER_STEP_DEGREES"] = saved
        else:
            os.environ.pop("DELFIN_FFFREE_ROTAMER_STEP_DEGREES", None)


def test_enumerate_max_configs_cap_caps_product():
    """With ridiculously low max_configs, output is truncated."""
    from delfin.fffree.single_bond_rotamers import enumerate_single_bond_rotamers
    m, P = _build_mol("CCCCC")
    out = list(enumerate_single_bond_rotamers(
        m, P, step_degrees=60.0, max_configs=5,
    ))
    assert len(out) == 5
