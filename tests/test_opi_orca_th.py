"""Tests for the ``orca_th_opi`` adapter: reading the thermochemistry of an
existing frequency run, refusing imaginary modes for a minimum, and building a
re-evaluation input at another temperature/pressure through OPI.

Run with the OPI interpreter (``DELFIN_OPI_REQUIRED=1``); skipped without OPI.
All fixtures below are hand-written, synthetic ORCA-format text — never real
ORCA output — and carry the mandated first line.
"""

import pytest

pytest.importorskip("opi")

from pathlib import Path

from delfin.imag import collect_imaginary_modes
from delfin.orca_read import read_thermochemistry
from delfin.tools import _registry

_HEADER = "# SYNTHETIC test fixture written for DELFIN -- not ORCA output"


def _out_block(temp: float, *, quasirrho: bool = True) -> str:
    """One synthetic ``THERMOCHEMISTRY AT <temp>K`` block the reader targets."""
    qr_line = (
        "Vibrational entropy computed according to the QRRHO of S. Grimme "
        "(The Journal of Chemical Physics, 2012, 136, 9955)"
        if quasirrho
        else "Vibrational entropy computed according to the standard RRHO of an ideal gas"
    )
    return f"""------------------
THERMOCHEMISTRY AT {temp}K
------------------
Temperature                         ...     {temp} K
Pressure                            ...      1.00 atm
Bonding analysis:
Lowest frequency mode                18.2 cm-1
--------------------------------------------------
INNER ENERGY
--------------------------------------------------
Electronic energy                    ...    -151.55083691 Eh
Zero point energy                    ...       0.02631509 Eh      16.51 kcal/mol

--------------------------------------------------
ENTHALPY
--------------------------------------------------
Thermal Enthalpy correction          ...       0.00094421 Eh
Total Enthalpy                       ...    -151.52034633 Eh

--------------------------------------------------
ENTROPY
--------------------------------------------------
Thermal entropy correction           ...       0.00072596 Eh
{qr_line}
Total entropy correction             ...       0.02583637 Eh      16.21 kcal/mol

--------------------------------------------------
GIBBS FREE ENERGY
--------------------------------------------------
Thermal Gibbs free energy correction ...       0.00465890 Eh
Total Gibbs free energy correction   ...       0.00156606 Eh
Final Gibbs free energy              ...    -151.54618270 Eh
"""

_TWO_TEMP_OUT = (
    _HEADER
    + "\n"
    + _out_block(298.15, quasirrho=True)
    + "\n"
    + _out_block(350.0, quasirrho=True)
)


def _out_with_imaginary() -> str:
    """A freq output whose geometry sits on a saddle point: one imaginary mode."""
    imag = """------------------
THERMOCHEMISTRY AT 298.15K
------------------
--------------------------------------------------
IR SPECTRUM
--------------------------------------------------
  0:       -134.20 cm**-1   ***imaginary mode***
"""
    return _HEADER + "\n" + imag


def _write_out(tmp_path: Path, body: str) -> Path:
    p = tmp_path / "freq.out"
    p.write_text(body)
    return p


def _adapter():
    return _registry.get("orca_th_opi")


# --------------------------------------------------------------------------- #
# registration in the shared registry
# --------------------------------------------------------------------------- #


def test_registered_in_registry():
    """The adapter is discoverable by other blocks via the shared registry."""
    assert _adapter() is not None


# --------------------------------------------------------------------------- #
# reading an existing frequency output
# --------------------------------------------------------------------------- #


def test_reads_thermochemistry_from_existing_output(tmp_path):
    out = _write_out(tmp_path, _TWO_TEMP_OUT)
    result = _adapter().execute(tmp_path, output_file=str(out))
    assert result.ok
    d = result.data
    assert d["temperature_K"] == pytest.approx(298.15)
    assert d["pressure_atm"] == pytest.approx(1.0)
    assert d["quasirrho"] is True
    assert d["E_el_Eh"] == pytest.approx(-151.55083691)
    assert d["ZPE_Eh"] == pytest.approx(0.02631509)
    assert d["ZPE_kcal"] == pytest.approx(16.51)
    assert d["H_Eh"] == pytest.approx(-151.52034633)
    assert d["TS_Eh"] == pytest.approx(0.02583637)
    assert d["TS_kcal"] == pytest.approx(16.21)
    assert d["G_Eh"] == pytest.approx(-151.54618270)
    # G and H are printed by ORCA only in Eh; the adapter prints kcal too.
    assert "G_kcal" in d
    assert "H_kcal" in d


def test_selects_requested_temperature_block(tmp_path):
    """A %freq Temp list run carries several blocks; pick the chosen T."""
    out = _write_out(tmp_path, _TWO_TEMP_OUT)
    result = _adapter().execute(tmp_path, output_file=str(out), thermo_temp=350.0)
    assert result.ok
    assert result.data["temperature_K"] == pytest.approx(350.0)


def test_reads_rrho_when_not_grimme_qrrho(tmp_path):
    body = _HEADER + "\n" + _out_block(298.15, quasirrho=False)
    out = _write_out(tmp_path, body)
    result = _adapter().execute(tmp_path, output_file=str(out))
    assert result.ok
    assert result.data["quasirrho"] is False


# --------------------------------------------------------------------------- #
# refusals
# --------------------------------------------------------------------------- #


def test_refuses_imaginary_modes_for_a_minimum(tmp_path):
    out = _write_out(tmp_path, _out_with_imaginary())
    result = _adapter().execute(tmp_path, output_file=str(out))
    assert not result.ok
    assert "imaginary" in result.error.lower()


def test_allows_imaginary_modes_for_a_stated_ts(tmp_path):
    out = _write_out(tmp_path, _out_with_imaginary())
    result = _adapter().execute(tmp_path, output_file=str(out), is_ts=True)
    assert result.ok
    assert result.data["Nimag"] == 1
    assert len(result.data["imaginary_freqs_cm"]) == 1


def test_refuses_absent_temperature_not_in_output(tmp_path):
    """Requesting a T ORCA did not compute cannot be answered from the .out."""
    out = _write_out(tmp_path, _TWO_TEMP_OUT)
    result = _adapter().execute(tmp_path, output_file=str(out), thermo_temp=500.0)
    assert not result.ok
    assert "50" in result.error  # the absent temperature appears in the message


def test_refuses_non_positive_temperature(tmp_path):
    out = _write_out(tmp_path, _TWO_TEMP_OUT)
    result = _adapter().execute(tmp_path, output_file=str(out), re_eval_temp=0.0)
    assert not result.ok


def test_refuses_missing_output_file(tmp_path):
    result = _adapter().execute(tmp_path, output_file=str(tmp_path / "nope.out"))
    assert not result.ok


# --------------------------------------------------------------------------- #
# re-evaluation input (OPI typed objects)
# --------------------------------------------------------------------------- #


def test_builds_re_eval_input_with_orca_own_rules(tmp_path):
    """A new T/p that ORCA did not print is requested: build a %freq rerun that
    reuses the existing Hessian, so ORCA only repeats the thermochemistry
    (manual 4.6.1 / 4.7)."""
    out = _write_out(tmp_path, _TWO_TEMP_OUT)
    inp = tmp_path / "reval.inp"
    result = _adapter().execute(
        tmp_path,
        output_file=str(out),
        re_eval_temp=500.0,
        re_eval_pressure=2.0,
        hessian=str(out.parent / "freq.hess"),
        input_path=inp,
    )
    # The reading still succeeds for the existing block; the re-eval input is
    # written out for ORCA to run.
    assert result.ok
    text = inp.read_text()
    # OPI writes the keyword and block option in lower case (`!freq`,
    # `temp`); assert case-insensitively so the check targets the
    # directive and value, not OPI's spelling.
    assert "freq" in text.lower()
    assert "temp" in text.lower()
    assert "500" in text
    assert "2" in text


# --------------------------------------------------------------------------- #
# imaginary-mode helper is the same one the adapter uses
# --------------------------------------------------------------------------- #


def test_imaginary_collector_reuse(tmp_path):
    out = _write_out(tmp_path, _out_with_imaginary())
    assert collect_imaginary_modes(str(out)) == [(0, -134.20)]
