"""Equivalence tests: DELFIN's thermochemistry readers vs OPI (orca-pi).

For the four small ORCA 6.1.1 fixtures in tests/fixtures/orca611/ (H2O and O2
only, committed), assert that the existing DELFIN functions and OPI
(get_free_energy, get_enthalpy, get_entropy, get_zpe, get_inner_energy) return
the same value within a stated tolerance, and that both behave sensibly on the
not-converged run.

Where they differ, the test states which is right and why (measured, not
guessed) -- see the per-function commentary below.

OPI background (measured on ORCA 6.1.1):
  * get_free_energy  = ORCA "Final Gibbs free energy ... X Eh"         (G)
  * get_enthalpy     = ORCA "Total Enthalpy     ... X Eh"              (H = U + kB*T)
  * get_entropy      = ORCA "Final entropy term ... X Eh"   (== -T*S, the
                       entropy multiplied by temperature, in Hartree --
                       NOT molar entropy in cal/mol/K)
  * get_zpe          = ORCA "Zero point energy ... X Eh"               (ZPE)
  * get_inner_energy = ORCA "Total thermal energy ... X Eh"            (U)
  All in atomic units (Hartree).

Units: DELFIN's readers are all in Hartree too (find_gibbs_energy,
find_ZPE, find_electronic_energy in delfin/energies.py; extract_thermochem
in delfin/api.py). The DELFIN values are read at ORCA's printed precision
(8 decimals after the point), OPI at full double precision; the measured
differences are < 1e-8 Hartree. A tolerance of 1e-6 Hartree is therefore
generous and stable, while far smaller than any thermochemistry difference
that would change a conclusion.

Fixture semantics (measured):
  * h2o_optfreq       -- converged Opt+Freq: the ONLY fixture with a
                         thermochemistry block. All five thermo values exist.
  * o2_triplet, h2o_tddft -- converged single-point / TDDFT, NO thermochemistry:
                         ORCA prints no Thermodynamics block, so both DELFIN
                         and OPI return None for every thermo quantity.
  * scf_not_converged -- SCF failed after 2 cycles: no thermochemistry and no
                         final energy. OPI scf_converged()/terminated_normally()
                         and DELFIN agree the run failed, and both return None.
                         (This is safer by construction than, e.g., the
                         HOMO-LUMO gap, which OPI would return from
                         unconverged orbitals; the Thermodynamics block only
                         exists after a completed Freq print.)

OPI needs the orca_2json executable (part of an ORCA >= 6.1.1 install) on
PATH for create_property_json(). The tests skip when OPI is not installed or
orca_2json is not on PATH -- these are machine requirements, not DELFIN
behaviour. create_property_json() writes files, so every run's .out and
.property.txt are copied into a per-test tmp_path first; the committed
fixture folder is never written to.
"""
import os
import shutil

import pytest

from delfin.api import extract_thermochem
from delfin.energies import (
    find_electronic_energy,
    find_gibbs_energy,
    find_ZPE,
)

FIXTURE_DIR = os.path.join(os.path.dirname(__file__), "fixtures", "orca611")
FIXTURES = ["h2o_optfreq", "o2_triplet", "h2o_tddft", "scf_not_converged"]

# In Hartree (atomic units). ORCA prints 8 decimals after the point; OPI
# resolves the same number to full double precision. Measured deviations on
# h2o_optfreq are 1.7e-10 .. 3.7e-9 Hartree, so this tolerance is robust.
TOL_EH = 1e-6


def _opi_available() -> bool:
    """OPI is importable and ORCA's orca_2json is on PATH (machine req.)."""
    try:
        import opi  # noqa: F401
    except Exception:
        return False
    return shutil.which("orca_2json") is not None


# `pytest` may run this module on a machine without OPI: guard the import.
from opi.output.core import Output  # noqa: E402  (only after the check below)

if not _opi_available():
    pytest.skip(
        "OPI (orca-pi) or ORCA orca_2json not available; "
        "OPI-equivalence tests skipped",
        allow_module_level=True,
    )


def _opi_result(tmp_path, name):
    """An OPI Output over a fresh copy of fixture `name`, in a tmp_path so
    create_property_json() never writes into the committed fixture folder."""
    work = tmp_path / name
    work.mkdir()
    for ext in (".out", ".property.txt"):
        shutil.copy(os.path.join(FIXTURE_DIR, name + ext), work)
    out = Output(basename=name, working_dir=work)
    out.create_property_json()
    out.parse()
    return out


def _delfin_thermo_full(out_path):
    """extract_thermochem on the .out FILE (not a folder the caller may parse
    as 'first thermochemistry in dir')."""
    return extract_thermochem(str(out_path))


def _collector_inner_energy(out_path):
    """U (inner energy) as the reporting collector reads it: the finished run's
    'Total thermal energy ... X Eh' line, via parse_thermochemistry."""
    from delfin.reporting.delfin_collector import parse_thermochemistry

    thermo = parse_thermochemistry(out_path)
    if thermo is None:
        return None
    return thermo.get("total_thermal_energy_hartree")


def _fixture_out(name):
    return os.path.join(FIXTURE_DIR, name + ".out")


def test_converged_freq_matches(tmp_path):
    """On the converged Opt+Freq fixture every thermo reader agrees with OPI."""
    name = "h2o_optfreq"
    # Enforce that the fixture really carries a thermochemistry block, so a
    # regression that silently drops it is caught as an infra error, not as a
    # masked comparison of two None values.
    assert any("THERMOCHEMISTRY" in line
               for line in open(_fixture_out(name), errors="replace"))
    out = _opi_result(tmp_path, name)
    copy = tmp_path / name / (name + ".out")
    d = _delfin_thermo_full(copy)
    assert out.scf_converged() is True
    assert out.terminated_normally() is True

    # G: get_free_energy == find_gibbs_energy == extract_thermochem.final_gibbs
    assert out.get_free_energy() == pytest.approx(
        find_gibbs_energy(_fixture_out(name)), abs=TOL_EH)
    assert out.get_free_energy() == pytest.approx(d.final_gibbs, abs=TOL_EH)
    # ZPE: get_zpe == find_ZPE == extract_thermochem.zpe
    assert out.get_zpe() == pytest.approx(
        find_ZPE(_fixture_out(name)), abs=TOL_EH)
    assert out.get_zpe() == pytest.approx(d.zpe, abs=TOL_EH)
    # H: get_enthalpy == extract_thermochem.enthalpy_corr (both are the TOTAL
    # enthalpy; 'enthalpy_corr' is misnamed in DELFIN -- it reads the total).
    assert out.get_enthalpy() == pytest.approx(d.enthalpy_corr, abs=TOL_EH)
    # S: get_entropy == extract_thermochem.entropy_total; BOTH are the entropy
    # term T*S in Hartree, NOT molar entropy in cal/mol/K.
    assert out.get_entropy() == pytest.approx(d.entropy_total, abs=TOL_EH)
    # U (inner): get_inner_energy == the collector's total_thermal_energy.
    assert out.get_inner_energy() == pytest.approx(
        _collector_inner_energy(copy), abs=TOL_EH)


@pytest.mark.parametrize("name", ["o2_triplet", "h2o_tddft"])
def test_no_thermochemistry_returns_none(tmp_path, name):
    """Converged but no Frequency block: no Thermodynamics in the output, so
    both DELFIN and OPI return None for every thermo quantity -- while the
    electronic energy still exists and agrees."""
    assert any("ORCA TERMINATED NORMALLY" in line
               for line in open(_fixture_out(name), errors="replace"))
    out = _opi_result(tmp_path, name)
    copy = tmp_path / name / (name + ".out")
    d = _delfin_thermo_full(copy)
    assert out.scf_converged() is True
    assert out.terminated_normally() is True
    assert out.get_free_energy() is None
    assert out.get_enthalpy() is None
    assert out.get_entropy() is None
    assert out.get_zpe() is None
    assert out.get_inner_energy() is None
    # DELFIN agrees: no Gibbs, no ZPE, no thermochemistry block at all.
    assert find_gibbs_energy(_fixture_out(name)) is None
    assert find_ZPE(_fixture_out(name)) is None
    assert d.final_gibbs is None and d.zpe is None and d.enthalpy_corr is None
    assert _collector_inner_energy(copy) is None
    # Electronic energy is still read and agrees.
    assert out.get_final_energy() == pytest.approx(
        find_electronic_energy(_fixture_out(name)), abs=TOL_EH)


def test_not_converged_returns_none_not_a_value(tmp_path):
    """SCF failed: both readers agree the run failed and hand back None --
    OPI must not report an unconverged thermochemistry value (contrast the
    HOMO-LUMO gap, which OPI would still return from unconverged orbitals)."""
    name = "scf_not_converged"
    assert any("SCF NOT CONVERGED" in line
               for line in open(_fixture_out(name), errors="replace"))
    out = _opi_result(tmp_path, name)
    copy = tmp_path / name / (name + ".out")
    d = _delfin_thermo_full(copy)
    # OPI flags the failure.
    assert out.scf_converged() is False
    assert out.terminated_normally() is False
    # Every thermo value is None, not a plausible-looking number.
    for fn in (out.get_free_energy, out.get_enthalpy, out.get_entropy,
               out.get_zpe, out.get_inner_energy, out.get_final_energy):
        assert fn() is None
    # DELFIN agrees on None and on the failure.
    assert find_gibbs_energy(_fixture_out(name)) is None
    assert find_ZPE(_fixture_out(name)) is None
    assert d.final_gibbs is None and d.zpe is None and d.enthalpy_corr is None
    assert _collector_inner_energy(copy) is None
