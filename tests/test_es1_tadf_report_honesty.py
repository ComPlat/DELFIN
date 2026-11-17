"""The TADF-xTB report says what an xTB screening may claim.

sTDA-xTB / sTD-DFT-xTB excitation energies and Delta E_ST values are
screening estimates.  Grimme & Bannwarth's sTDA-xTB paper (J. Chem.
Phys. 145, 054103 (2016), DOI 10.1063/1.4959605 — accuracy typically
~0.5 eV for excitation energies; XTB-based screening, not a
production method) and the sTD-DFT-xTB extension target ranking, not
reference-grade numbers.  A referee accepts a TADF screening when the
report says: rank/order (S1 vs T1, sign of Delta E_ST), not absolute
energies; confirm at TD-DFT (or STEOM-CCSD / SCS double-hybrid) level
before quoting rates or absolute values.

The report carries this where a consumer reads it:

* the printed one-line result names the screening estimate as an
  estimate and points to the confirming level;
* the JSON payload carries the same note as data, so any consumer
  (dashboard, downstream workflow) can rely on it, not just parse
  prose.
"""

from __future__ import annotations

import json

import pytest

from delfin.tadf_xtb import TadfXtbResult, _print_result, screening_note


def _result(**overrides) -> TadfXtbResult:
    fields = {
        "label": "mol1",
        "smiles": "c1ccccc1",
        "charge": 0,
        "multiplicity": 1,
        "xtb_method": "XTB2",
        "excited_method": "stda",
        "workdir": "tadf_xtb_runs/mol1",
        "initial_xyz": "mol1_start.xyz",
        "s0_xyz": "mol1_S0.xyz",
        "crest_xyz": None,
        "goat_xyz": None,
        "crest_output": None,
        "goat_output": None,
        "s0_opt_output": "S0_opt.out",
        "s0_xtb_energy_eh": -1.0,
        "s0_sp_output": "S0_sp.out",
        "t1_opt_xyz": None,
        "t1_vertical_xtb_output": None,
        "t1_opt_output": None,
        "t1_xtb_energy_eh": None,
        "t1_vertical_xtb_ev": None,
        "t1_adiabatic_ev": None,
        "t1_relaxation_ev": None,
        "stokes_shift_est_ev": None,
        "s1_relaxed_est_ev": None,
        "t1_relaxed_ev": None,
        "singlet_output": "singlet.out",
        "triplet_output": "triplet.out",
        "tda_singlet": "stda_s.out",
        "tda_triplet": "stda_t.out",
        "s1_state": 1,
        "s1_ev": 4.8,
        "s1_nm": 258.3,
        "lambda_abs_nm": 258.3,
        "lambda_em_s1_nm": 280.0,
        "lambda_pl_est_nm": 280.0,
        "s1_f_osc": 0.01,
        "first_allowed_singlet_state": 1,
        "first_allowed_singlet_ev": 4.8,
        "first_allowed_singlet_nm": 258.3,
        "first_allowed_singlet_f_osc": 0.01,
        "brightest_singlet_state": 1,
        "brightest_singlet_ev": 4.8,
        "brightest_singlet_nm": 258.3,
        "brightest_singlet_f_osc": 0.01,
        "t1_ev": 3.9,
        "t1_nm": 317.9,
        "delta_est_ev": 0.9,
        "delta_est_relaxed_est_ev": None,
    }
    fields.update(overrides)
    return TadfXtbResult(**fields)


def test_screening_note_names_limit_and_confirmation_level():
    note = screening_note()
    assert "estimate" in note.lower()
    assert "TD-DFT" in note or "STEOM-CCSD" in note


def test_printed_result_carries_the_note(capsys):
    _print_result(_result())
    out = capsys.readouterr().out.lower()
    assert "screening estimate" in out
    assert "td-dft" in out or "steom-ccsd" in out


def test_payload_carries_the_note_as_data():
    """The JSON payload includes the note as a field, not only in prose."""
    from delfin.tadf_xtb import _payload_with_notes

    payload = _payload_with_notes([_result()], {})
    note = payload["screening_note"]
    assert isinstance(note, str)
    assert "estimate" in note.lower()
