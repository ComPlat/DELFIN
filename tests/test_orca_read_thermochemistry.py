"""read_thermochemistry parses a frequency run's thermochemistry blocks.

Runs under the CI interpreter (no OPI): :func:`delfin.orca_read.read_thermochemistry`
is a text reader only, importing nothing from OPI. The fixture is hand-written
ORCA-format text -- never a real ORCA output -- so CI (which has no ORCA binary)
can run it.
"""
from pathlib import Path

import pytest

from delfin.orca_read import read_thermochemistry

# The first line marks this as a synthetic fixture per LICENSE_GUARD.md.
HEADER = "# SYNTHETIC test fixture written for DELFIN -- not ORCA output"


def _write_freq_out(tmp_path: Path, body: str) -> Path:
    out = tmp_path / "freq.out"
    out.write_text(HEADER + "\n" + body, encoding="utf-8")
    return out


# Two thermochemistry blocks, as ORCA prints for a `%freq Temp 298.15,350` run
# (manual, section 4.6.1). Values imitate a small molecule; the numbers need
# not obey exact thermochemistry relations, they only exercise the parser.
# _write_freq_out prepends the synthetic-fixture marker.
_TWO_TEMP_BODY = """
!!! FREQUENCY CALCULATIONS !!!

--------------------------
THERMOCHEMISTRY AT 298.15K
--------------------------

Temperature         ... 298.15 K
Pressure            ... 1.00 atm

Summary of contributions to the inner energy U:
Electronic energy                ...   -151.55083691 Eh
Zero point energy                ...      0.02631509 Eh      16.51 kcal/mol
Thermal vibrational correction   ...      0.00040105 Eh       0.25 kcal/mol

FINAL SINGLE POINT ENERGY       -151.550836910212

--------------------------
ENTHALPY
--------------------------
Total Enthalpy                   ...   -151.52034633 Eh

--------------------------
ENTROPY
--------------------------
Vibrational entropy computed according to the QRRHO of S. Grimme
Final entropy term               ...      0.02583637 Eh     16.21 kcal/mol

--------------------------
GIBBS FREE ENERGY
--------------------------
The Gibbs free energy is G = H - T*S
Total enthalpy                   ...   -151.52034633 Eh
Total entropy correction         ...     -0.02583637 Eh    -16.21 kcal/mol
Final Gibbs free energy           ...   -151.54618270 Eh

------------------
VIBRATIONAL FREQUENCIES
------------------
     0:      370.67 cm**-1
     1:      947.27 cm**-1

--------------------------
THERMOCHEMISTRY AT 350.00K
--------------------------

Temperature         ... 350.00 K
Pressure            ... 1.00 atm

Summary of contributions to the inner energy U:
Electronic energy                ...   -151.55083691 Eh
Zero point energy                ...      0.02631509 Eh      16.51 kcal/mol
Thermal vibrational correction   ...      0.00150105 Eh       0.94 kcal/mol

--------------------------
ENTHALPY
--------------------------
Total Enthalpy                   ...   -151.51934633 Eh

--------------------------
ENTROPY
--------------------------
Vibrational entropy computed according to the QRRHO of S. Grimme
Final entropy term               ...      0.03083637 Eh     19.35 kcal/mol

--------------------------
GIBBS FREE ENERGY
--------------------------
Final Gibbs free energy           ...   -151.53218270 Eh

ORCA TERMINATED NORMALLY
"""


def test_reads_both_temperature_blocks(tmp_path):
    out = _write_freq_out(tmp_path, _TWO_TEMP_BODY)
    readings = read_thermochemistry(out)
    assert len(readings) == 2
    temps = [r.temperature_K for r in readings]
    assert temps == pytest.approx([298.15, 350.00])


def test_reads_zero_point_energy_in_eh_and_kcal(tmp_path):
    out = _write_freq_out(tmp_path, _TWO_TEMP_BODY)
    r = read_thermochemistry(out)[0]
    assert r.zpe_Eh == pytest.approx(0.02631509, abs=1e-8)
    assert r.zpe_kcal == pytest.approx(16.51, abs=1e-5)


def test_reads_enthalpy_in_eh(tmp_path):
    out = _write_freq_out(tmp_path, _TWO_TEMP_BODY)
    r = read_thermochemistry(out)[0]
    assert r.h_Eh == pytest.approx(-151.52034633, abs=1e-8)
    # kcal computed from DELFIN's Eh->kcal constant because ORCA prints H
    # only in Eh (manual, thermochemistry section).
    assert r.h_kcal == pytest.approx(-151.52034633 * 627.5094740631, rel=1e-9)


def test_reads_ts_entropy_term_in_eh_and_kcal(tmp_path):
    out = _write_freq_out(tmp_path, _TWO_TEMP_BODY)
    r = read_thermochemistry(out)[0]
    assert r.ts_Eh == pytest.approx(0.02583637, abs=1e-8)
    assert r.ts_kcal == pytest.approx(16.21, abs=1e-5)


def test_reads_gibbs_in_eh(tmp_path):
    out = _write_freq_out(tmp_path, _TWO_TEMP_BODY)
    r = read_thermochemistry(out)[0]
    assert r.g_Eh == pytest.approx(-151.54618270, abs=1e-8)
    assert r.g_kcal == pytest.approx(-151.54618270 * 627.5094740631, rel=1e-9)


def test_reads_electronic_energy_in_eh(tmp_path):
    out = _write_freq_out(tmp_path, _TWO_TEMP_BODY)
    r = read_thermochemistry(out)[0]
    assert r.e_el_Eh == pytest.approx(-151.55083691, abs=1e-8)
    assert r.e_el_kcal == pytest.approx(-151.55083691 * 627.5094740631, rel=1e-9)


def test_reads_temperature_and_pressure(tmp_path):
    out = _write_freq_out(tmp_path, _TWO_TEMP_BODY)
    r = read_thermochemistry(out)[1]
    assert r.temperature_K == pytest.approx(350.00)
    assert r.pressure_atm == pytest.approx(1.00)
    assert r.quasirrho is True


def test_detects_plain_rrho_when_not_grimme_qrrho(tmp_path):
    body = _TWO_TEMP_BODY.replace(
        "Vibrational entropy computed according to the QRRHO of S. Grimme",
        "Vibrational entropy computed according to the RRHO")
    readings = read_thermochemistry(_write_freq_out(tmp_path, body))
    assert all(r.quasirrho is False for r in readings)


def test_missing_file_returns_empty(tmp_path):
    assert read_thermochemistry(tmp_path / "absent.out") == []


def test_rrho_unknown_when_no_entropy_line(tmp_path):
    body = _TWO_TEMP_BODY.replace(
        "Vibrational entropy computed according to the QRRHO of S. Grimme\n",
        "")
    readings = read_thermochemistry(_write_freq_out(tmp_path, body))
    assert readings[0].quasirrho is None
