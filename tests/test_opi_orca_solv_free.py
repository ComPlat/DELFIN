"""Tests for delfin.tools.adapters.orca_solv_free_opi.

An OPI building block for solvation free energies: gas-phase and solvated
single points at the same geometry, dG_solv, and a documented thermodynamic
cycle for pKa and the redox potential E0 against a stated reference
electrode.  These are red-first tests that pin the public surface of the new
module; they require the OPI venv (DELFIN_OPI_REQUIRED=1).
"""

import math

import pytest

pytest.importorskip("opi")

from delfin.tools.adapters import orca_solv_free_opi as solvfree
from delfin.tools._registry import get as registry_get

# A tiny closed-shell molecule: water, 10 electrons, singlet.
H2O_ELEMENTS = ["O", "H", "H"]
H2O_COORDINATES = [
    [0.000, 0.000, 0.117],
    [0.000, 0.755, -0.467],
    [0.000, -0.755, -0.467],
]

# O2, 16 electrons: triplet (mult 3) is valid, singlet (mult 1) is valid,
# mult 2 (an even multiplicity on an even electron count) is physical garbage.
O2_ELEMENTS = ["O", "O"]
O2_COORDINATES = [
    [0.0, 0.0, 0.0],
    [0.0, 0.0, 1.2],
]

_EH_TO_KJ = 2625.499639479947  # Hartree -> kJ/mol, DELFIN's conv constant


# ---------------------------------------------------------------------------
# electron count / spin validation
# ---------------------------------------------------------------------------
def test_electron_count_water():
    assert solvfree._electron_count(H2O_ELEMENTS, 0) == 10


def test_electron_count_anion():
    # hydroxide-ish test charge: neutral water minus a proton
    assert solvfree._electron_count(H2O_ELEMENTS, -1) == 11


def test_spin_valid_singlet_water():
    n = solvfree._electron_count(H2O_ELEMENTS, 0)
    assert solvfree._validate_spin(0, 1, n) is None


def test_spin_valid_triplet_o2():
    n = solvfree._electron_count(O2_ELEMENTS, 0)
    assert solvfree._validate_spin(0, 3, n) is None


def test_spin_refused_even_mult_on_even_electron_count():
    n = solvfree._electron_count(O2_ELEMENTS, 0)
    err = solvfree._validate_spin(0, 2, n)
    assert err is not None
    assert "electron" in err.lower()


def test_spin_refused_mult_below_one():
    n = solvfree._electron_count(H2O_ELEMENTS, 0)
    assert solvfree._validate_spin(0, 0, n) is not None


# ---------------------------------------------------------------------------
# input build: gas and solvated legs on the SAME geometry level
# ---------------------------------------------------------------------------
def test_build_gas_leg_has_no_solvation_keyword(tmp_path):
    gas_inp, _solv_inp = solvfree.build_solv_free_inputs(
        tmp_path,
        elements=H2O_ELEMENTS,
        coordinates=H2O_COORDINATES,
        charge=0,
        mult=1,
        functional="B3LYP",
        basis="def2-SVP",
        ri=False,
        solvent_model="CPCM",
        solvent="Water",
        maxcore=None,
        ncores=1,
    )
    text = gas_inp.read_text()
    # OPI writes !-keywords lowercased; the gas leg carries no solvation.
    low = text.lower()
    assert "cpcm(" not in low and "smd(" not in low
    assert "b3lyp" in low and "def2-svp" in low


def test_build_solv_leg_cpcm_water(tmp_path):
    _gas_inp, solv_inp = solvfree.build_solv_free_inputs(
        tmp_path,
        elements=H2O_ELEMENTS,
        coordinates=H2O_COORDINATES,
        charge=0,
        mult=1,
        functional="B3LYP",
        basis="def2-SVP",
        ri=False,
        solvent_model="CPCM",
        solvent="Water",
        maxcore=None,
        ncores=1,
    )
    text = solv_inp.read_text()
    low = text.lower()
    assert "cpcm(" in low
    assert "water" in low


def test_build_solv_leg_smd_water(tmp_path):
    _gas_inp, solv_inp = solvfree.build_solv_free_inputs(
        tmp_path,
        elements=H2O_ELEMENTS,
        coordinates=H2O_COORDINATES,
        charge=0,
        mult=1,
        functional="B3LYP",
        basis="def2-SVP",
        ri=False,
        solvent_model="SMD",
        solvent="Water",
        maxcore=None,
        ncores=1,
    )
    text = solv_inp.read_text()
    assert "smd(" in text.lower()


def test_legs_share_geometry_and_method(tmp_path):
    gas_inp, solv_inp = solvfree.build_solv_free_inputs(
        tmp_path,
        elements=H2O_ELEMENTS,
        coordinates=H2O_COORDINATES,
        charge=0,
        mult=1,
        functional="B3LYP",
        basis="def2-SVP",
        ri=False,
        solvent_model="CPCM",
        solvent="Water",
        maxcore=None,
        ncores=1,
    )
    g = gas_inp.read_text()
    s = solv_inp.read_text()
    # same functional, basis, charge, multiplicity, coordinates (OPI lowercases)
    for token in ("b3lyp", "def2-svp", "* xyz 0 1"):
        assert token in g and token in s
    # the solvated leg alone carries the solvation keyword
    assert "cpcm(" in s and "cpcm(" not in g


# ---------------------------------------------------------------------------
# refusals
# ---------------------------------------------------------------------------
def test_validate_unknown_solvent_refused():
    errors, _warnings = solvfree.validate_solv_free(
        elements=H2O_ELEMENTS, charge=0, mult=1,
        functional="B3LYP", basis="def2-SVP", ri=False,
        solvent_model="CPCM", solvent="NotASolvent",
    )
    assert any("solvent" in e.lower() for e in errors)


def test_validate_unknown_model_refused():
    errors, _warnings = solvfree.validate_solv_free(
        elements=H2O_ELEMENTS, charge=0, mult=1,
        functional="B3LYP", basis="def2-SVP", ri=False,
        solvent_model="PCM", solvent="Water",
    )
    assert any("model" in e.lower() or "PCM" in e for e in errors)


def test_validate_alpb_refused_outside_xtb():
    errors, _warnings = solvfree.validate_solv_free(
        elements=H2O_ELEMENTS, charge=0, mult=1,
        functional="B3LYP", basis="def2-SVP", ri=False,
        solvent_model="ALPB", solvent="Water",
    )
    assert errors  # ALPB is an xTB-only model; this block is DFT single points


def test_validate_wrong_spin_refused():
    errors, _warnings = solvfree.validate_solv_free(
        elements=O2_ELEMENTS, charge=0, mult=2,
        functional="B3LYP", basis="def2-SVP", ri=False,
        solvent_model="CPCM", solvent="Water",
    )
    assert any("multiplicit" in e.lower() for e in errors)


def test_validate_missing_coordinates_refused():
    errors, _warnings = solvfree.validate_solv_free(
        elements=H2O_ELEMENTS, charge=0, mult=1,
        functional="B3LYP", basis="def2-SVP", ri=False,
        solvent_model="CPCM", solvent="Water",
    )
    # no coordinates supplied: the block must refuse rather than guess
    assert any("coordinate" in e.lower() for e in errors)


def test_validate_ok_returns_no_errors():
    errors, warnings = solvfree.validate_solv_free(
        elements=H2O_ELEMENTS, charge=0, mult=1,
        functional="B3LYP", basis="def2-SVP", ri=False,
        solvent_model="CPCM", solvent="Water",
        coordinates=H2O_COORDINATES,
    )
    assert not errors


def test_validate_mixing_levels_between_legs_refused():
    # the two legs must share one method/geometry set; asking for different
    # functionals on the gas and solvated sides must be refused out loud
    errors, _warnings = solvfree.validate_solv_free(
        elements=H2O_ELEMENTS, charge=0, mult=1,
        functional="B3LYP", basis="def2-SVP", ri=False,
        solvent_model="CPCM", solvent="Water",
        coordinates=H2O_COORDINATES,
        solv_functional="PBE0",  # a different method on the solvated leg
    )
    assert any("solv_functional" in e or "different" in e.lower() or "level" in e.lower()
               for e in errors)


# ---------------------------------------------------------------------------
# reading: dG_solv with units, from synthetic energies (no real ORCA output)
# ---------------------------------------------------------------------------
def test_read_dg_solv_value_and_units():
    e_gas = -76.4080  # Eh
    e_solv = -76.4241  # Eh
    res = solvfree.read_solv_free(e_gas, e_solv, solvent_model="CPCM",
                                  solvent="Water")
    assert res.ok and res.meaningful
    assert res.dg_solv_Eh == pytest.approx(e_solv - e_gas)
    assert res.dg_solv_kJ_per_mol == pytest.approx((e_solv - e_gas) * _EH_TO_KJ)


def test_read_failed_leg_marks_not_meaningful():
    res = solvfree.read_solv_free(None, -76.4, solvent_model="CPCM",
                                  solvent="Water")
    assert not res.meaningful
    assert not res.ok


def test_read_named_reference_electrode():
    res = solvfree.read_solv_free(-76.4, -76.42, solvent_model="SMD",
                                  solvent="Water")
    assert res.reference_electrode is not None
    value, name = res.reference_electrode
    assert name  # the reference is always stated
    assert isinstance(value, float)


# ---------------------------------------------------------------------------
# thermodynamic cycle: redox potential E0 vs a stated reference electrode
# ---------------------------------------------------------------------------
def test_redox_e0_oxidation_numeric():
    # E_ox = -((G0_valence - G_+1) * conv) / F - E_ref  (DELFIN/cli_calculations)
    g0 = -100.0
    gox = -99.5
    e_ref = 4.345  # V vs SHE
    conv = _EH_TO_KJ
    f = 96.4853321233100184  # kJ/(V mol)
    expected = -((g0 - gox) * conv) / f - e_ref
    got = solvfree.e0_ox(g0, gox, n=1, e_ref_v=e_ref)
    assert got == pytest.approx(expected)


def test_redox_e0_reduction_numeric():
    # E_red = ((G0_valence - G_-1) * conv) / F - E_ref
    g0 = -100.0
    gred = -99.4
    e_ref = 4.345
    conv = _EH_TO_KJ
    f = 96.4853321233100184
    expected = ((g0 - gred) * conv) / f - e_ref
    got = solvfree.e0_red(g0, gred, n=1, e_ref_v=e_ref)
    assert got == pytest.approx(expected)


# ---------------------------------------------------------------------------
# thermodynamic cycle: pKa via isodesmic pKa_from_cycle anchor
# ---------------------------------------------------------------------------
def test_pka_cycle_imported_and_numeric(tmp_path):
    from delfin.pka import pka_from_cycle
    # isodesmic cycle, four solution-phase free energies; acetic acid anchor
    gibbs = {
        "target_HA": -100.0,
        "target_A": -99.2,
        "reference_HA": -90.0,
        "reference_A": -89.5,
    }
    pka = pka_from_cycle(gibbs, reference_pka=4.756, temperature_k=298.15)
    assert math.isfinite(pka)


def test_pka_cycle_incomplete_raises(tmp_path):
    from delfin.pka import pka_from_cycle
    # missing reference_A must raise rather than guess
    with pytest.raises(ValueError):
        pka_from_cycle({"target_HA": -100.0, "target_A": -99.2,
                        "reference_HA": -90.0}, reference_pka=4.756)


# ---------------------------------------------------------------------------
# registration in the tool registry
# ---------------------------------------------------------------------------
def test_adapter_registered():
    adapter = registry_get("orca_solv_free_opi")
    assert adapter is not None
    assert adapter.name == "orca_solv_free_opi"


# ---------------------------------------------------------------------------
# probe: ORCA's own parser must accept the built solvated input
# ---------------------------------------------------------------------------
def test_orca_parser_accepts_solvated_input(tmp_path):
    """Hand the OPI-built solvated input to ORCA's own parser via
    ``delfin.orca_knowledge.probe.probe_input`` (<= 2 cores, small molecule,
    parser check only — the probe guard keeps it a probe, not a calculation).
    Skipped with a reason when the ORCA binary is absent."""
    from delfin.orca_knowledge.probe import probe_input

    _gas_inp, solv_inp = solvfree.build_solv_free_inputs(
        tmp_path,
        elements=H2O_ELEMENTS,
        coordinates=H2O_COORDINATES,
        charge=0,
        mult=1,
        functional="B3LYP",
        basis="def2-SVP",
        ri=False,
        solvent_model="CPCM",
        solvent="Water",
        maxcore=1000,
        ncores=1,
    )
    result = probe_input(solv_inp.read_text(), timeout_s=300)
    if result["verdict"] == "unavailable":
        pytest.skip("ORCA binary not available for probe")
    assert result["verdict"] in ("parsed", "finished"), (
        "ORCA rejected the OPI-built solvated input; messages: "
        f"{result['messages']}"
    )


def test_redox_e0_divides_by_the_electron_count():
    """Real-run finding: n was accepted but ignored, so a two-electron
    potential came out equal to the one-electron one."""
    g0, gox, e_ref = -1.0, -0.9, 4.345
    one = solvfree.e0_ox(g0, gox, n=1, e_ref_v=e_ref) + e_ref
    two = solvfree.e0_ox(g0, gox, n=2, e_ref_v=e_ref) + e_ref
    assert two == pytest.approx(one / 2)
    one = solvfree.e0_red(g0, -1.1, n=1, e_ref_v=e_ref) + e_ref
    two = solvfree.e0_red(g0, -1.1, n=2, e_ref_v=e_ref) + e_ref
    assert two == pytest.approx(one / 2)
    with pytest.raises(ValueError):
        solvfree.e0_ox(g0, gox, n=0)
