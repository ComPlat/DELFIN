"""ML1 P2 builder control (nacht-s13): physical realism in the ASE adapters.

ASE's calculators accept any Atoms object; a wrong element or a charge/mult
the calculator cannot model fails at best with an obscure error, at worst
silently wrong (EMT raises; LennardJones and the MACE/ANI wrappers compute a
number for any element).  The adapter must declare each calculator's limits
and enforce them fail-closed on the ACTUAL structure BEFORE the calculator
runs, through the public path.

Behavior sources (verified 2026-10-24):
- xtb GFN: parametrized for all spd-block elements Z=1-86; Grimme,
  Bannwarth, Shushkov, J. Chem. Theory Comput. 13 (2017) 1989-2009,
  DOI 10.1021/acs.jctc.7b00118 (xtb repo assets/references.bib).
  Charge/mult pass through to xtb (-c / -u) — supported.
- orca: full periodic-table scope, charge/mult passed to ORCA — supported.
- emt: ASE calculator docstring (installed ase.calculators.emt):
  standard metals Al, Cu, Ag, Au, Ni, Pd, Pt; additionally H, C, N, O
  with parameters "not for any serious use".
- lj / mace / ani: compute energies for ANY element (generic wrappers) —
  element scope None (model-defined), so no element refusal.
- mace / ani / lj: no charge/mult handling in the adapter's calculator
  constructor (ase_adapter.py:72-79) — silently ignored, so refused.
"""

from __future__ import annotations

import pytest

from delfin.tools.adapters.ase_adapter import (
    AseMdAdapter,
    AseOptimizeAdapter,
    AseSinglePointAdapter,
)

OPT = AseOptimizeAdapter()
SP = AseSinglePointAdapter()
MD = AseMdAdapter()


class TestDeclaredTables:
    def test_element_and_validity_tables_declared(self):
        from delfin.tools.adapters import ase_adapter as mod

        assert hasattr(mod, "ASE_SUPPORTED_ELEMENTS")
        assert hasattr(mod, "ASE_VALIDITY")
        assert hasattr(mod, "ASE_SOURCES")
        for calc in mod._ASE_CALCS:
            assert calc in mod.ASE_SUPPORTED_ELEMENTS, calc
            assert calc in mod.ASE_VALIDITY, calc
            assert mod.ASE_VALIDITY[calc].strip(), calc

    def test_emt_table_matches_ase_docstring(self):
        from delfin.tools.adapters import ase_adapter as mod

        assert set(mod.ASE_SUPPORTED_ELEMENTS["emt"]) == (
            frozenset({"Al", "Cu", "Ag", "Au", "Ni", "Pd", "Pt"})
            | frozenset({"H", "C", "N", "O"})
        )

    def test_xtb_table_is_full_spd_block(self):
        from delfin.tools.adapters import ase_adapter as mod

        xtb = mod.ASE_SUPPORTED_ELEMENTS["xtb"]
        # spd-block = main groups + all d-block up to Z=86: spot-checks
        for el in ("H", "C", "Fe", "Zn", "I", "Xe"):
            assert el in xtb, el
        # beyond Z=86 must be excluded
        for el in ("Rn", "Ra", "Ac"):
            assert el not in xtb, el


class TestChargeMultFailClosed:
    def test_charge_on_mace_refused(self):
        with pytest.raises(ValueError):
            SP.validate_params(calculator="mace", charge=1)

    def test_mult_on_mace_refused(self):
        with pytest.raises(ValueError):
            SP.validate_params(calculator="mace", mult=3)

    def test_charge_on_xtb_supported(self):
        SP.validate_params(calculator="xtb", charge=1, mult=2)

    def test_charge_on_orca_supported(self):
        SP.validate_params(calculator="orca", charge=1, mult=2)

    def test_charge_on_lj_refused(self):
        with pytest.raises(ValueError):
            SP.validate_params(calculator="lj", charge=1)


class TestElementEnforcementOnRealPath:
    def test_gold_on_emt_passes_element_check(self, tmp_path):
        geom = tmp_path / "au.xyz"
        geom.write_text("1\n\nAu 0.0 0.0 0.0\n", encoding="utf-8")
        result = SP.execute(tmp_path, geometry=geom, calculator="emt")
        error = result.error or ""
        if result.status.name == "FAILED":
            assert "unsupported element" not in error

    def test_tin_on_emt_refused_before_calculation(self, tmp_path):
        geom = tmp_path / "sn.xyz"
        geom.write_text("1\n\nSn 0.0 0.0 0.0\n", encoding="utf-8")
        result = SP.execute(tmp_path, geometry=geom, calculator="emt")
        assert result.status.name == "FAILED"
        assert "Sn" in (result.error or "")

    def test_uranium_on_mace_element_and_charge_refused(self, tmp_path):
        # mace has model-defined element scope (no refusal) but does NOT
        # model charge — a mace run with charge must fail closed.
        geom = tmp_path / "u.xyz"
        geom.write_text("1\n\nU 0.0 0.0 0.0\n", encoding="utf-8")
        result = SP.execute(tmp_path, geometry=geom, calculator="mace",
                            charge=1)
        assert result.status.name == "FAILED"
        assert "charge" in (result.error or "").lower()

    def test_unknown_calculator_value_refused(self):
        with pytest.raises(ValueError):
            SP.validate_params(calculator="magic_force_field")
