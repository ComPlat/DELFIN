"""The ESD report must carry every rate the parser fills (reviewer F4).

The parser fills summary.fluor / summary.phosp from the FLUOR/PHOSP
outputs (esd_results.py:417-425) and counts them in has_data
(esd_results.py:165-172), but generate_esd_report wrote only the
states/ISC/IC sections -- a fluorescence or phosphorescence rate,
DELFIN's product, vanished from the human-readable report, whose body
then said "No ESD data available." while has_data was True.
"""

from __future__ import annotations

from pathlib import Path

from delfin.esd_results import (
    ESDSummary,
    FluorResult,
    ICResult,
    ISCResult,
    PhospResult,
    StateResult,
)
from delfin.reporting.esd_report import generate_esd_report


def _summary_fluor_only() -> ESDSummary:
    summary = ESDSummary()
    summary.fluor["S1>S0"] = FluorResult(
        rate=1.23e8,
        temperature=298.15,
        delta_cm1=35301.0,
        source=Path("S1_S0_FLUOR.out"),
    )
    return summary


def _summary_phosp_only() -> ESDSummary:
    summary = ESDSummary()
    summary.phosp["T1>S0"] = PhospResult(
        sublevel_rates=(1.0e3, 1.0e3, 1.0e3),
        rate_mean=1.0e3,
        temperature=298.15,
        delta_cm1=28185.0,
        source=Path("T1_S0_PHOSP.out"),
    )
    return summary


def test_report_carries_the_fluorescence_rate(tmp_path):
    text = Path(tmp_path / "ESD.txt")
    generate_esd_report(_summary_fluor_only(), text)
    body = text.read_text()
    assert "No ESD data available." not in body
    assert "FLUOR rate constant" in body
    assert "1.230000e+08" in body
    assert "[NOT A RESULT" not in body


def test_report_carries_the_phosphorescence_rate(tmp_path):
    text = Path(tmp_path / "ESD.txt")
    generate_esd_report(_summary_phosp_only(), text)
    body = text.read_text()
    assert "No ESD data available." not in body
    assert "PHOSP rate constant" in body
    assert "1.000000e+03" in body


def test_report_marks_a_flawed_fluor_rate_not_a_result(tmp_path):
    summary = _summary_fluor_only()
    summary.fluor["S1>S0"].problem = "negative rate"
    text = Path(tmp_path / "ESD.txt")
    generate_esd_report(summary, text)
    body = text.read_text()
    assert "[NOT A RESULT: negative rate]" in body


def test_report_still_says_no_data_when_the_parser_found_nothing(tmp_path):
    text = Path(tmp_path / "ESD.txt")
    generate_esd_report(ESDSummary(), text)
    assert "No ESD data available." in text.read_text()


def test_isc_and_ic_sections_are_unchanged(tmp_path):
    summary = ESDSummary()
    summary.isc["S1>T1(Ms=0)"] = ISCResult(
        rate=2.46e7, temperature=298.15, delta_cm1=1000.0,
        soc=(50.0, 0.0), fc_percent=60.0, ht_percent=40.0,
        source=Path("S1_T1_ISC_ms0.out"),
    )
    summary.ic["S1>S0"] = ICResult(
        rate=3.79e-3, temperature=298.15, delta_cm1=35301.0,
        source=Path("S1_S0_IC.out"),
    )
    text = Path(tmp_path / "ESD.txt")
    generate_esd_report(summary, text)
    body = text.read_text()
    assert "ISC rate constants (s^-1):" in body
    assert "2.460000e+07" in body
    assert "IC rate constants:" in body
    assert "3.790000e-03" in body
