"""Reviewer s15 adversarial tests for package ML1 (ML potentials / ASE).

These test the PUBLIC path only: the registered adapters'
validate_params + the calculator factory in delfin.mlp_tools.calculators.
No torchani/aimnet/mace install is required: every case here must be
decided by validation, before any model is created.
"""

from __future__ import annotations

import pytest

from delfin.tools.adapters.mlp import (
    MlpOptimizeAdapter,
    MlpRankAdapter,
    MlpSinglePointAdapter,
)

SP = MlpSinglePointAdapter()
OPT = MlpOptimizeAdapter()
RANK = MlpRankAdapter()


# ---------------------------------------------------------------- A1: element coverage must fail closed
class TestA1ElementCoverage:
    """ANI-2x covers C,H,N,O,S,F,Cl only. A metal or Br must be refused."""

    def test_unknown_backend_is_refused(self):
        with pytest.raises(ValueError):
            SP.validate_params(backend="magic_potential")

    def test_backend_must_be_required(self):
        with pytest.raises(ValueError):
            SP.validate_params()
        with pytest.raises(ValueError):
            OPT.validate_params()
        with pytest.raises(ValueError):
            RANK.validate_params(backend="ani2x")

    def test_metal_on_ani2x_declared_unsupported(self):
        # A supported-element declaration that does not exist yet must at
        # least make the limitation findable/declared somewhere.
        from delfin.tools.adapters import mlp as mlp_mod

        tables = [
            getattr(mlp_mod, name, None)
            for name in dir(mlp_mod)
            if "element" in name.lower() or "supported" in name.lower()
        ]
        tables = [t for t in tables if t is not None]
        assert tables, (
            "MLP adapter declares no supported-elements table: "
            "the physical limit (ANI-2x = C,H,N,O,S,F,Cl) is undeclared"
        )

    def test_ani2x_table_excludes_metals_and_bromine(self):
        from delfin.tools.adapters import mlp as mlp_mod

        names = [
            n
            for n in dir(mlp_mod)
            if "element" in n.lower() or "supported" in n.lower()
        ]
        for n in names:
            t = getattr(mlp_mod, n)
            if isinstance(t, (set, frozenset, tuple, list)) and all(
                isinstance(x, str) and len(x) <= 3 for x in t
            ):
                assert "Fe" not in t and "Zn" not in t and "Br" not in t
                assert "C" in t and "Cl" in t


# ---------------------------------------------------------------- A2: charge/mult honored or refused
class TestA2ChargeMultiplicity:
    """7 of 8 backends ignore charge/mult. The adapter must refuse them."""

    def test_nonzero_charge_on_non_aimnet2_backend_is_refused(self):
        # Currently this SILENTLY passes: the defect under test.
        with pytest.raises(ValueError):
            SP.validate_params(backend="ani2x", charge=1)

    def test_nonzero_mult_on_non_aimnet2_backend_is_refused(self):
        with pytest.raises(ValueError):
            SP.validate_params(backend="ani2x", mult=3)

    def test_aimnet2_may_carry_charge_and_mult(self):
        SP.validate_params(backend="aimnet2", charge=1, mult=2)


# ---------------------------------------------------------------- A4: aliases must not bypass the enum
class TestA4AliasBypass:
    def test_alias_mace_is_accepted(self):
        # 'mace' is a documented alias for mace_off; enum-level acceptance
        # (or a normalized failure) is fine, a crash is not.
        try:
            SP.validate_params(backend="mace")
        except ValueError as exc:
            assert "mace" in str(exc).lower()


# ---------------------------------------------------------------- A5: edge cases
class TestA5EdgeCases:
    def test_empty_backend_string_refused(self):
        with pytest.raises(ValueError):
            SP.validate_params(backend="")

    def test_none_backend_refused(self):
        with pytest.raises(ValueError):
            SP.validate_params(backend=None)

    def test_case_insensitive_backend(self):
        # calculators.py lowercases; the adapter enum check must not be
        # stricter than the factory in one direction or the other.
        SP.validate_params(backend="ANI2X")
