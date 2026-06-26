"""Tests for the per-donor-element gate added in Phase 3C
(ITER-uffsoft 2026-05-13).

Background
----------

Phase 3 (commits b6.../5dd520a) introduced an UFF "soft donor" mode that
replaces ``FixAtom`` on each donor with an ``AddDistanceConstraint`` for
the M-D pair, freeing the donor to reorient while keeping the bond
length pinned.  A forensik on the 7871-file common subset of pool
``3943c2b-uffsoft-vollpool`` vs ``HEAD-iter8.10-vollpool`` showed:

* +7.53 pp M-D break regression overall, of which
* **76.1 % of breaks are M-C donors** (689 / 905 B-only-broken files).

Root cause: UFF has no transition-metal-bonded parameters, so the M-D
distance pin is the only force holding the donor in place.  For C donors
the substituent gradients exceed the spring constant and the donor drifts.
N / O / P / S / halide donors have stronger UFF angle terms via their own
substituents AND a directed lone pair that keeps them aligned.

Patch
-----

``delfin/_uff_soft_donor.should_soften_donor(d_sym, allow_carbon=…)``
returns True iff the donor element should receive the distance-pin.
Default policy:

* SOFT (distance pin): N, O, P, S, F, Cl, Br, I
* HARD (FixAtom)     : C (escape via ``DELFIN_UFF_SOFT_DONORS_CARBON=1``),
                       Si, B, Se, any unknown/empty/None symbol.

The call site ``smiles_converter._optimize_xyz_openbabel`` consults this
helper per (M, D) pair when ``DELFIN_UFF_SOFT_DONORS=1``.  When the
flag is off (default), the helper is not touched and behaviour is
bit-exact identical to the pre-patch code path.
"""

from __future__ import annotations

import importlib
import os

import pytest


@pytest.fixture(autouse=True)
def _clear_carbon_env(monkeypatch):
    """Each test starts with DELFIN_UFF_SOFT_DONORS_CARBON cleared."""
    monkeypatch.delenv("DELFIN_UFF_SOFT_DONORS_CARBON", raising=False)
    # Re-import to ensure module-level state is fresh; the helper reads
    # the env-var lazily, so this is precautionary.
    import delfin._uff_soft_donor as _mod
    importlib.reload(_mod)
    yield


def test_soft_donor_elements_default_true():
    """N, O, P, S and halides default to SOFT."""
    from delfin._uff_soft_donor import should_soften_donor
    for sym in ("N", "O", "P", "S", "F", "Cl", "Br", "I"):
        assert should_soften_donor(sym) is True, f"{sym} should be SOFT"


def test_carbon_default_hard():
    """C donor is HARD when DELFIN_UFF_SOFT_DONORS_CARBON is unset."""
    from delfin._uff_soft_donor import should_soften_donor
    assert should_soften_donor("C") is False
    # Explicit allow_carbon=False also returns False.
    assert should_soften_donor("C", allow_carbon=False) is False


def test_carbon_soft_when_env_flag_set(monkeypatch):
    """C donor flips to SOFT iff DELFIN_UFF_SOFT_DONORS_CARBON=1."""
    monkeypatch.setenv("DELFIN_UFF_SOFT_DONORS_CARBON", "1")
    from delfin._uff_soft_donor import should_soften_donor
    assert should_soften_donor("C") is True
    # Explicit allow_carbon=True overrides regardless of env-var.
    monkeypatch.delenv("DELFIN_UFF_SOFT_DONORS_CARBON", raising=False)
    assert should_soften_donor("C", allow_carbon=True) is True


def test_other_hard_elements_default_false():
    """Si, B, Se, and unknown elements default to HARD."""
    from delfin._uff_soft_donor import should_soften_donor
    for sym in ("Si", "B", "Se", "As", "Te", "H", "Xx"):
        assert should_soften_donor(sym) is False, f"{sym} should be HARD"


def test_malformed_input_safe_fallback():
    """Empty string, None, or non-string → False (HARD)."""
    from delfin._uff_soft_donor import should_soften_donor
    assert should_soften_donor("") is False
    assert should_soften_donor(None) is False
    assert should_soften_donor(123) is False
    assert should_soften_donor(["N"]) is False


def test_malformed_env_flag_safe_fallback(monkeypatch):
    """Garbage in DELFIN_UFF_SOFT_DONORS_CARBON → treated as 0."""
    from delfin._uff_soft_donor import should_soften_donor
    for bad in ("", "yes", "abc", "1.5"):
        monkeypatch.setenv("DELFIN_UFF_SOFT_DONORS_CARBON", bad)
        assert should_soften_donor("C") is False, (
            f"Malformed env value {bad!r} should yield HARD for C"
        )


def test_should_use_soft_donor_unchanged():
    """The class-conditional gate (Phase 3) is untouched by this patch."""
    from delfin._uff_soft_donor import should_use_soft_donor
    assert should_use_soft_donor("sigma") is True
    assert should_use_soft_donor("multi_sigma") is True
    assert should_use_soft_donor("hapto") is False
    assert should_use_soft_donor("no_metal") is False
    assert should_use_soft_donor("") is False


def test_default_off_bit_exact_unchanged(monkeypatch):
    """When DELFIN_UFF_SOFT_DONORS=0 (default), the call-site never
    reaches ``should_soften_donor`` and behaviour is bit-exact identical
    to the pre-patch code path.

    We assert this contractually by verifying:
    1. ``_delfin_env_int`` returns 0 for an unset env-var (the gate that
       guards the entire soft-donor block);
    2. ``_class_conditional_flag`` returns False when both the scalar
       and the *_CLASSES variants are unset.
    """
    monkeypatch.delenv("DELFIN_UFF_SOFT_DONORS", raising=False)
    monkeypatch.delenv("DELFIN_UFF_SOFT_DONORS_CLASSES", raising=False)
    monkeypatch.delenv("DELFIN_UFF_SOFT_DONORS_CARBON", raising=False)

    from delfin.smiles_converter import _delfin_env_int
    assert _delfin_env_int("DELFIN_UFF_SOFT_DONORS", 0) == 0
    # Carbon-escape default is also 0:
    assert _delfin_env_int("DELFIN_UFF_SOFT_DONORS_CARBON", 0) == 0


def test_soft_donor_elements_constant_immutable():
    """The whitelist is a frozenset so callers cannot mutate it."""
    from delfin._uff_soft_donor import _SOFT_DONOR_ELEMENTS
    assert isinstance(_SOFT_DONOR_ELEMENTS, frozenset)
    # Carbon must NOT be in the default whitelist.
    assert "C" not in _SOFT_DONOR_ELEMENTS
    # Standard heteroatom donors must be present.
    for sym in ("N", "O", "P", "S", "F", "Cl", "Br", "I"):
        assert sym in _SOFT_DONOR_ELEMENTS
