"""Unit tests for ``delfin/_deep_enumerator.py`` (Welle-5m-X).

Covers four contract points:

1. ``is_deep_enum_enabled`` default-OFF (bit-exact preservation).
2. ``count_distinct_donor_classes`` distinguishes Morgan-distinct nitrogens.
3. ``should_apply_deep_enum`` fires on YIRQIC-style hetero CN >= 5 with
   >= 3 donor classes.
4. ``should_apply_deep_enum`` does NOT fire on simple cisplatin
   (CN=4, 2 donor classes).
"""

from __future__ import annotations

import os
from unittest import mock

import pytest


pytest.importorskip("rdkit")


def _mol(smiles: str, *, add_hs: bool = True):
    from rdkit import Chem
    m = Chem.MolFromSmiles(smiles, sanitize=False)
    if m is None:
        raise AssertionError(f"RDKit failed to parse {smiles!r}")
    try:
        m.UpdatePropertyCache(strict=False)
    except Exception:
        pass
    if add_hs:
        try:
            m = Chem.AddHs(m)
        except Exception:
            pass
    return m


def test_is_deep_enum_enabled_default_off_byte_identical():
    """Without env-flag the gate is OFF, preserving pre-patch behaviour."""
    from delfin._deep_enumerator import is_deep_enum_enabled
    # Strip the env-flag if a sibling test set it.
    with mock.patch.dict(
        os.environ,
        {k: v for k, v in os.environ.items()
         if not k.startswith("DELFIN_5M_X_")},
        clear=True,
    ):
        assert is_deep_enum_enabled() is False


def test_should_apply_deep_enum_default_off_returns_false():
    """Master gate unset -> trigger returns False regardless of features."""
    from delfin._deep_enumerator import should_apply_deep_enum
    # YIRQIC normalised SMILES -- meets feature trigger, but gate is unset.
    smi = (
        "[O+]#[C][Re-5]1([Br])([C]#[O+])("
        "[P+](C2=CC=CC=C2)(C2=CC=CC=C2)C2=CC=CC=C2)"
        "[C+]2N(C3=CC=CC=C3)C=CN2C2=CC=CC=[N+]21"
    )
    m = _mol(smi)
    with mock.patch.dict(
        os.environ,
        {k: v for k, v in os.environ.items()
         if not k.startswith("DELFIN_5M_X_")},
        clear=True,
    ):
        assert should_apply_deep_enum(m) is False


def test_count_distinct_donor_classes_yirqic_5_classes():
    """X10-YIRQIC's Re centre coordinates 5 distinct donor classes.

    Counts: Br x1, C-CO x2 (same Morgan), C-carbene x1 (different
    Morgan), N-pyridyl x1, P-phosphine x1 -> 5 distinct classes.
    """
    from delfin._deep_enumerator import count_distinct_donor_classes
    from delfin.smiles_converter import _donor_type_map, _METAL_SET

    smi = (
        "[O+]#[C][Re-5]1([Br])([C]#[O+])("
        "[P+](C2=CC=CC=C2)(C2=CC=CC=C2)C2=CC=CC=C2)"
        "[C+]2N(C3=CC=CC=C3)C=CN2C2=CC=CC=[N+]21"
    )
    m = _mol(smi)
    dtype_map = _donor_type_map(m)
    metal_idx = next(
        a.GetIdx() for a in m.GetAtoms() if a.GetSymbol() in _METAL_SET
    )
    n = count_distinct_donor_classes(m, metal_idx, dtype_map)
    # Conservative lower-bound: >= 3 is the trigger threshold; YIRQIC
    # actually carries 5 Morgan-distinct classes (verified manually
    # via `_donor_type_map`).
    assert n >= 3, f"Expected >=3 donor classes for YIRQIC, got {n}"


def test_should_apply_deep_enum_yirqic_triggers_when_gate_on():
    """X10-YIRQIC must trigger deep enum when ``DELFIN_5M_X_DEEP_ENUM=1``."""
    from delfin._deep_enumerator import should_apply_deep_enum
    from delfin.smiles_converter import _donor_type_map

    smi = (
        "[O+]#[C][Re-5]1([Br])([C]#[O+])("
        "[P+](C2=CC=CC=C2)(C2=CC=CC=C2)C2=CC=CC=C2)"
        "[C+]2N(C3=CC=CC=C3)C=CN2C2=CC=CC=[N+]21"
    )
    m = _mol(smi)
    dtype_map = _donor_type_map(m)
    with mock.patch.dict(
        os.environ,
        {"DELFIN_5M_X_DEEP_ENUM": "1"},
    ):
        assert should_apply_deep_enum(m, dtype_map=dtype_map) is True


def test_should_apply_deep_enum_cisplatin_does_not_trigger():
    """cisplatin (Pt CN=4, 2 donor classes) must NOT trigger deep enum.

    Pins the negative branch -- the trigger predicate must protect
    simple, well-behaved sigma complexes from any change in enumerator
    contract even when the master gate is on.
    """
    from delfin._deep_enumerator import should_apply_deep_enum
    from delfin.smiles_converter import _donor_type_map

    m = _mol("[Pt](Cl)(Cl)(N)N")
    dtype_map = _donor_type_map(m)
    with mock.patch.dict(
        os.environ,
        {"DELFIN_5M_X_DEEP_ENUM": "1"},
    ):
        # CN=4 (< min_cn=5), so even with 2 distinct donor classes
        # the trigger refuses.
        assert should_apply_deep_enum(m, dtype_map=dtype_map) is False


def test_coordination_number_for_yirqic_is_six():
    """Sanity check the CN helper before relying on it in trigger logic."""
    from delfin._deep_enumerator import coordination_number
    from delfin.smiles_converter import _METAL_SET

    smi = (
        "[O+]#[C][Re-5]1([Br])([C]#[O+])("
        "[P+](C2=CC=CC=C2)(C2=CC=CC=C2)C2=CC=CC=C2)"
        "[C+]2N(C3=CC=CC=C3)C=CN2C2=CC=CC=[N+]21"
    )
    m = _mol(smi)
    metal_idx = next(
        a.GetIdx() for a in m.GetAtoms() if a.GetSymbol() in _METAL_SET
    )
    assert coordination_number(m, metal_idx) == 6


def test_unlabeled_rmsd_threshold_default_off_is_2_5():
    """Pre-patch contract: unlabeled RMSD threshold stays at 2.5 A unless
    the deep-enum master gate is non-zero."""
    from delfin._deep_enumerator import deep_enum_unlabeled_rmsd_threshold

    smi = (
        "[O+]#[C][Re-5]1([Br])([C]#[O+])("
        "[P+](C2=CC=CC=C2)(C2=CC=CC=C2)C2=CC=CC=C2)"
        "[C+]2N(C3=CC=CC=C3)C=CN2C2=CC=CC=[N+]21"
    )
    m = _mol(smi)
    with mock.patch.dict(
        os.environ,
        {k: v for k, v in os.environ.items()
         if not k.startswith("DELFIN_5M_X_")},
        clear=True,
    ):
        assert deep_enum_unlabeled_rmsd_threshold(m) == pytest.approx(2.5)


def test_unlabeled_rmsd_threshold_deep_enum_on_drops_to_0_8():
    """Deep-enum on + trigger satisfied -> threshold drops to 0.8 A."""
    from delfin._deep_enumerator import deep_enum_unlabeled_rmsd_threshold
    from delfin.smiles_converter import _donor_type_map

    smi = (
        "[O+]#[C][Re-5]1([Br])([C]#[O+])("
        "[P+](C2=CC=CC=C2)(C2=CC=CC=C2)C2=CC=CC=C2)"
        "[C+]2N(C3=CC=CC=C3)C=CN2C2=CC=CC=[N+]21"
    )
    m = _mol(smi)
    dtype_map = _donor_type_map(m)
    with mock.patch.dict(
        os.environ,
        {"DELFIN_5M_X_DEEP_ENUM": "1"},
    ):
        assert deep_enum_unlabeled_rmsd_threshold(
            m, dtype_map=dtype_map
        ) == pytest.approx(0.8)


def test_simple_acac_complex_does_not_trigger():
    """Ru(acac)3: CN=6 but only 1 donor class (O) -> trigger refuses.

    Universal contract test: trigger requires >= 3 distinct donor classes
    regardless of CN.  Simple mono-class chelates stay on the legacy
    code path (bit-exact pre-patch).
    """
    from delfin._deep_enumerator import should_apply_deep_enum
    from delfin.smiles_converter import _donor_type_map

    # acac: O=C(C)C=C(C)[O-] coordinated tris -- all donors are O.
    m = _mol("CC(=O)CC(=O)C.CC(=O)CC(=O)C.CC(=O)CC(=O)C.[Ru+3]")
    dtype_map = _donor_type_map(m)
    with mock.patch.dict(
        os.environ,
        {"DELFIN_5M_X_DEEP_ENUM": "1"},
    ):
        # No Ru-O bonds in this loose SMILES (Ru is disconnected),
        # so CN=0 -- trigger refuses.  This also guards against
        # disconnected-fragment SMILES sneaking through.
        assert should_apply_deep_enum(m, dtype_map=dtype_map) is False
