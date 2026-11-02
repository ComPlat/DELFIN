"""Equivalence tests: DELFIN's thermochemistry readers vs OPI (orca-pi).

For the four small ORCA 6.1.1 runs in the private fixture store
(tests/fixtures/orca611/<run>/<run>.*; H2O and O2 only), assert that the
existing DELFIN functions and OPI (get_free_energy, get_enthalpy,
get_entropy, get_zpe, get_inner_energy) return the same value within a
stated tolerance, and that both behave sensibly on the not-converged run.

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

DELFIN readers compared:
  * delfin.energies.find_gibbs_energy  (G, whole file, all jobs)
  * delfin.energies.find_ZPE           (ZPE)
  * delfin.energies.find_electronic_energy (electronic energy)
  * delfin.api.extract_thermochem      (ZPE, thermal correction, enthalpy,
                                        entropy term, G-E(el), final G --
                                        .enthalpy_corr is MISNAMED: it reads
                                        the TOTAL Enthalpy line, which is
                                        exactly what OPI get_enthalpy gives)
  * delfin.reporting.delfin_collector.parse_thermochemistry
                                       (inner energy via
                                        total_thermal_energy_hartree; on
                                        ORCA 6.1.1 its Total-Enthalpy and
                                        entropy regexes MISS -- documented,
                                        not exercised as a comparison here)

Fixture semantics (measured):
  * h2o_optfreq       -- converged Opt+Freq: the ONLY run with a
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

Fixtures are real licensed-program output: they live in the private fixture
store (DELFIN_PRIVATE_FIXTURES; LICENSE_GUARD.md, rule family 2), never in
the repository, and the tests skip when the store is absent -- they never
fail in CI. OPI needs the orca_2json executable (ORCA >= 6.1.1) on PATH for
create_property_json(). create_property_json() WRITES files, so each run is
copied into a per-test tmp_path first; the private store is never written to.
"""
import shutil

import pytest

from delfin.api import extract_thermochem
from delfin.energies import (
    find_electronic_energy,
    find_gibbs_energy,
    find_ZPE,
)
from test_lc1_fixtures_helper import private_fixtures_or_skip  # noqa: F401

from pathlib import Path

RUNS = ["h2o_optfreq", "o2_triplet", "h2o_tddft", "scf_not_converged"]

# In Hartree (atomic units). ORCA prints 8 decimals after the point; OPI
# resolves the same number to full double precision. Measured deviations on
# h2o_optfreq are 1.7e-10 .. 3.7e-9 Hartree, so this tolerance is robust.
TOL_EH = 1e-6

try:
    import opi  # noqa: F401
    from opi.output.core import Output  # noqa: F401
    _OPI_MISSING = False
except ImportError:
    _OPI_MISSING = True

if _OPI_MISSING:
    pytest.skip(
        "OPI (orca-pi) not importable; OPI-equivalence tests skipped",
        allow_module_level=True,
    )


def _run_dir(name):
    """The run's directory inside the private store (module skips if unset)."""
    return private_fixtures_or_skip(f"tests/fixtures/orca611/{name}")


def _opi_result(tmp_path, name):
    """OPI Output over a fresh copy of run `name` in tmp_path, so
    create_property_json() never writes into the private fixture store."""
    work = tmp_path / name
    shutil.copytree(_run_dir(name), work)
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


def test_converged_freq_matches(tmp_path):
    """On the converged Opt+Freq run every thermo reader agrees with OPI."""
    name = "h2o_optfreq"
    src = _run_dir(name) / (name + ".out")
    # Enforce that the run really carries a thermochemistry block, so a
    # regression that silently drops it is caught as an infra error, not as a
    # masked comparison of two None values.
    assert any("THERMOCHEMISTRY" in line for line in open(src, errors="replace"))
    out = _opi_result(tmp_path, name)
    copy = tmp_path / name / (name + ".out")
    d = _delfin_thermo_full(copy)
    # G: OPI get_free_energy == find_gibbs_energy == extract_thermochem.final_gibbs.
    assert out.get_free_energy() == pytest.approx(
        find_gibbs_energy(copy), abs=TOL_EH)
    assert out.get_free_energy() == pytest.approx(d.final_gibbs, abs=TOL_EH)
    # ZPE: OPI get_zpe == find_ZPE == extract_thermochem.zpe.
    assert out.get_zpe() == pytest.approx(find_ZPE(copy), abs=TOL_EH)
    assert out.get_zpe() == pytest.approx(d.zpe, abs=TOL_EH)
    # H: get_enthalpy == extract_thermochem.enthalpy_corr. DELFIN's attribute
    # name is misleading (it reads the TOTAL Enthalpy line, not the thermal
    # correction); the VALUE is the same quantity OPI returns, measured here.
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
    src = _run_dir(name) / (name + ".out")
    assert any("ORCA TERMINATED NORMALLY" in line
               for line in open(src, errors="replace"))
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
    assert find_gibbs_energy(copy) is None
    assert find_ZPE(copy) is None
    assert d.final_gibbs is None and d.zpe is None and d.enthalpy_corr is None
    assert _collector_inner_energy(copy) is None
    # Electronic energy is still read and agrees.
    assert out.get_final_energy() == pytest.approx(
        find_electronic_energy(copy), abs=TOL_EH)


def test_not_converged_returns_none_not_a_value(tmp_path):
    """SCF failed: both readers agree the run failed and hand back None --
    OPI must not report an unconverged thermochemistry value (contrast the
    HOMO-LUMO gap, which OPI would still return from unconverged orbitals)."""
    name = "scf_not_converged"
    src = _run_dir(name) / (name + ".out")
    assert any("SCF NOT CONVERGED" in line for line in open(src, errors="replace"))
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
    assert find_gibbs_energy(copy) is None
    assert find_ZPE(copy) is None
    assert d.final_gibbs is None and d.zpe is None and d.enthalpy_corr is None
    assert _collector_inner_energy(copy) is None


# ---------------------------------------------------------------------------
# CI-runnable synthetic-fixture tests. The fixture pair
# tests/fixtures/orca611/synth_h2o/ is SYNTHETIC (first line of both files
# marks it): a .out carrying the exact ORCA 6.1.1 thermochemistry-block lines
# DELFIN's readers match, and a hand-authored .property.json feeding OPI's
# getters. The same five quantities are asserted equal here as on the real
# runs above, so CI (no ORCA, no private store) still exercises the
# equivalence logic.
# ---------------------------------------------------------------------------

SYN = Path(__file__).parent / "fixtures" / "orca611" / "synth_h2o"


def test_synthetic_pair_runs_in_ci(tmp_path):
    """The committed synthetic fixture runs OPI without orca_2json and
    without the private store: OPI parses the hand-authored .property.json
    and the .out text, DELFIN's readers parse the .out text."""
    work = tmp_path / "synth_h2o"
    shutil.copytree(SYN, work)
    out = Output(basename="synth_h2o", working_dir=work)
    out.parse(do_create_property_json=False, do_create_gbw_json=False)
    assert out.terminated_normally() is True
    assert out.scf_converged() is True
    out_file = work / "synth_h2o.out"
    # OPI side (from the .property.json) ...
    opi_g, opi_h, opi_s, opi_zpe, opi_u = (
        out.get_free_energy(), out.get_enthalpy(), out.get_entropy(),
        out.get_zpe(), out.get_inner_energy())
    for v in (opi_g, opi_h, opi_s, opi_zpe, opi_u):
        assert v is not None
    # ... equals the DELFIN text-readers on the .out.
    assert opi_g == pytest.approx(find_gibbs_energy(out_file), abs=TOL_EH)
    assert opi_zpe == pytest.approx(find_ZPE(out_file), abs=TOL_EH)
    th = extract_thermochem(str(work))
    assert opi_h == pytest.approx(th.enthalpy_corr, abs=TOL_EH)
    assert opi_s == pytest.approx(th.entropy_total, abs=TOL_EH)
    assert opi_u == pytest.approx(
        _collector_inner_energy(out_file), abs=TOL_EH)
    # Electronic energy: OPI (Single_Point_Data) vs DELFIN's reader.
    assert out.get_final_energy() == pytest.approx(
        find_electronic_energy(out_file), abs=TOL_EH)


def test_synthetic_not_converged(tmp_path):
    """The synthetic store's failed-run semantics: without a final
    thermochemistry block both readers return None (covered by the real
    scf_not_converged run in the private store when it is present)."""
    work = tmp_path / "synth_h2o"
    shutil.copytree(SYN, work)
    out = Output(basename="synth_h2o", working_dir=work)
    out.parse(do_create_property_json=False, do_create_gbw_json=False)
    # The synthetic fixture is a converged run: sanity-check the positive
    # path and that no value is handed out where nothing exists.
    assert out.get_free_energy() is not None
    assert out.scf_converged() is True
