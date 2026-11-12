"""ML1 builder control (nacht-s13): physical-realism enforcement in the MLP
adapters must fail closed through the PUBLIC path.

These cases drive the real call site (adapter validate_params + execute),
not a shadow copy.  The element-refusal cases go through ``execute`` because
that is the only place a structure is available; they must fail BEFORE any
ML model is created (no torchani/aimnet install needed).

Sources for every declared limit (verified 2026-10-24):
- ANI-2x elements {H,C,N,O,S,F,Cl}: Devereux et al., J. Chem. Theory Comput.
  16 (2020) 4192, DOI 10.1021/acs.jctc.0c00121.
- MACE-OFF elements {H,C,N,O,F,P,S,Cl,Br,I}: Kovacs et al.,
  J. Am. Chem. Soc. 147 (2025) 9254, DOI 10.1021/jacs.4c07099.
- AIMNet2 elements {H,B,C,N,O,F,Si,P,S,Cl,As,Se,Br,I} and charge support:
  Anstine, Zubatyuk, Isayev, Chem. Sci. 16 (2025) 10228,
  DOI 10.1039/D4SC08572H.
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


class TestBackendValueChecked:
    """The backend must be checked against the enum VALUE, not presence."""

    def test_unknown_backend_refused(self):
        with pytest.raises(ValueError):
            SP.validate_params(backend="magic_potential")

    def test_empty_backend_refused(self):
        with pytest.raises(ValueError):
            SP.validate_params(backend="")

    def test_none_backend_refused(self):
        with pytest.raises(ValueError):
            SP.validate_params(backend=None)

    def test_backend_required_all_adapters(self):
        with pytest.raises(ValueError):
            SP.validate_params()
        with pytest.raises(ValueError):
            OPT.validate_params()
        with pytest.raises(ValueError):
            RANK.validate_params(ensemble="x.xyz")

    def test_case_insensitive_accepted(self):
        SP.validate_params(backend="ANI2X")

    def test_alias_mace_accepted(self):
        SP.validate_params(backend="mace")


class TestElementsDeclaredAndEnforced:
    """A declared per-backend element table must exist and be enforced."""

    def test_element_table_declared(self):
        from delfin.tools.adapters import mlp as mlp_mod

        assert hasattr(mlp_mod, "SUPPORTED_ELEMENTS")
        table = mlp_mod.SUPPORTED_ELEMENTS
        assert isinstance(table, dict) and table, "empty element table"
        assert "ani2x" in table
        assert set(table["ani2x"]) == {"H", "C", "N", "O", "S", "F", "Cl"}
        # only element-scoped backends carry a verified DOI source
        for name, scope in table.items():
            if scope is not None:
                assert name in mlp_mod.ELEMENT_TABLE_SOURCES

    def test_metal_refused_via_real_execute_path(self, tmp_path):
        geometry = tmp_path / "fe_molecule.xyz"
        geometry.write_text(
            "2\n\nFe 0.0 0.0 0.0\nCl 1.5 0.0 0.0\n", encoding="utf-8"
        )
        result = SP.execute(
            tmp_path, geometry=geometry, backend="ani2x"
        )
        assert result.status.name == "FAILED"
        assert "Fe" in (result.error or "")

    def test_supported_organic_passes_element_check(self, tmp_path):
        geometry = tmp_path / "chloroethane.xyz"
        geometry.write_text(
            "2\n\nC 0.0 0.0 0.0\nCl 1.8 0.0 0.0\n", encoding="utf-8"
        )
        # element check passes; execution may then fail on a missing model
        # install, but NEVER with an element-refusal error
        result = SP.execute(
            tmp_path, geometry=geometry, backend="ani2x"
        )
        error = result.error or ""
        if result.status.name == "FAILED":
            assert "not in the supported-elements set" not in error


class TestChargeMultiplicityFailClosed:
    """Non-AIMNet2 backends are charge-neutral closed-shell models."""

    def test_charge_on_ani2x_refused(self):
        with pytest.raises(ValueError):
            SP.validate_params(backend="ani2x", charge=1)

    def test_mult_on_ani2x_refused(self):
        with pytest.raises(ValueError):
            SP.validate_params(backend="ani2x", mult=3)

    def test_aimnet2_may_carry_charge_and_mult(self):
        SP.validate_params(backend="aimnet2", charge=1, mult=2)

    def test_charge_enforced_in_execute_path(self, tmp_path):
        geometry = tmp_path / "ethane.xyz"
        geometry.write_text("2\n\nC 0.0 0.0 0.0\nH 1.1 0.0 0.0\n", encoding="utf-8")
        result = SP.execute(
            tmp_path, geometry=geometry, backend="ani2x", charge=1
        )
        assert result.status.name == "FAILED"
        assert "charge" in (result.error or "").lower()


class TestValidityDeclared:
    """Every backend states what it is valid for and what it is not."""

    def test_validity_statement_per_backend(self):
        from delfin.tools.adapters import mlp as mlp_mod

        assert hasattr(mlp_mod, "VALIDITY")
        for backend in mlp_mod._MLP_BACKENDS:
            assert backend in mlp_mod.VALIDITY, backend
            assert mlp_mod.VALIDITY[backend].strip(), backend

    def test_ani2x_not_valid_for_periodic_systems(self):
        from delfin.tools.adapters import mlp as mlp_mod

        assert "periodic" in mlp_mod.VALIDITY["ani2x"].lower()
