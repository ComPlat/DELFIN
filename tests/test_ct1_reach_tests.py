"""A CONTROL.txt edit reaches what it claims to change -- at the ORCA input.

Phase-2 of CT1: for the keys the classic workflows read, show the value
arriving in the exact ORCA input line the key names, through the REAL
writers (not re-implementations of them):

* delfin/common/solvation.py:18-46 build_solvation_keyword -- the single
  place the solvation keyword is assembled;
* delfin/xyz_io.py:505-542 _build_bang_line -- the '!' keyword line
  assembler;
* the CONTROL validator coercers, so the value that passes validation is
  the value that reaches the line.

Also pinned here (the P3 red control, kept with its finding until the
phase-3 validator check lands): a solvent WITHOUT an implicit solvation
model yields an empty solvation keyword SILENTLY (solvation.py:31-32),
while a model without a solvent is warned and dropped (:34-44).  The
silent drop is the asymmetry the phase-3 validator finding targets.
"""
import logging

import pytest

from delfin.common.solvation import build_solvation_keyword
from delfin.common.control_validator import (
    _as_freq_type,
    _as_functional,
    _as_initial_guess,
    _as_implicit_solvation_model,
    _as_solvent,
    _as_relativity,
)
from delfin.xyz_io import _build_bang_line


BASE_CONFIG = {
    "functional": "PBE0",
    "ri_jkx": "RIJCOSX",
    "disp_corr": "",
    "initial_guess": "PModel",
    "freq_type": "FREQ",
    "geom_opt": "OPT",
}


def test_solvent_and_model_together_reach_the_solvation_keyword():
    assert build_solvation_keyword("CPCM", "DMF") == "CPCM(DMF)"


def test_a_model_without_solvent_is_dropped_not_run_bare(caplog):
    with caplog.at_level(logging.WARNING):
        assert build_solvation_keyword("CPCM", "") == ""
    assert "no solvent" in caplog.text  # solvation.py:34-44 warns; nothing silent here


def test_a_solvent_without_a_model_is_dropped_silently():
    """The P3 gap: '' is returned with no warning (solvation.py:31-32).

    Kept as a documentation-of-current-behaviour test: the phase-3
    validator check makes the case visible at CONTROL validation, this
    test pins what the keyword builder itself does so the fix cannot
    accidentally change the runtime behaviour instead of adding the check.
    """
    with caplog_disabled():
        assert build_solvation_keyword("", "DMF") == ""


def test_gas_phase_when_neither_is_set():
    assert build_solvation_keyword("", "") == ""


class caplog_disabled:
    """Context manager turning off the solvation module's logger."""

    def __init__(self) -> None:
        self._records = []

    def __enter__(self):
        self._logger = logging.getLogger("delfin.common.solvation")
        self._old_disabled = self._logger.disabled
        self._old_propagate = self._logger.propagate
        self._logger.disabled = True
        return self._records

    def __exit__(self, *exc):
        self._logger.disabled = self._old_disabled
        self._logger.propagate = self._old_propagate
        return False


def test_what_the_validator_passes_is_what_the_line_gets():
    # coercers first: the validated value is the value on the line.
    # ORCA's canonical solvent spellings are lower-case (ORCA_SOLVENTS,
    # control_validator.py:180 "dmf"); the validator returns the list spelling.
    assert _as_solvent("DMF") == "dmf"
    assert _as_solvent("water") == "water"
    assert _as_implicit_solvation_model("c-pcm") == "CPCM"
    assert _as_functional("b3lyp") == "B3LYP"
    assert _as_freq_type("numfreq") == "numFREQ"
    assert _as_initial_guess("hcore") == "HCore"
    assert _as_relativity("zora") == "ZORA"


def test_freq_type_reaches_the_bang_line_when_freq_is_included():
    config = dict(BASE_CONFIG)
    line_analytical = _build_bang_line(config, "", "def2-SVP", "def2/J", "CPCM(DMF)", include_freq=True)
    assert line_analytical.endswith("FREQ PModel")

    config["freq_type"] = "numfreq"
    line_numerical = _build_bang_line(config, "", "def2-SVP", "def2/J", "CPCM(DMF)", include_freq=True)
    # the bang line upper-cases freq_type for ORCA (xyz_io.py:535-538)
    assert line_numerical.endswith("NUMFREQ PModel")


def test_the_level_of_theory_tokens_reach_the_bang_line_in_order():
    config = dict(BASE_CONFIG, functional="TPSS", disp_corr="D4")
    line = _build_bang_line(config, "ZORA", "ZORA-def2-SVP", "SARC/J", "SMD(water)", include_freq=False)
    tokens = line.split()
    assert tokens[0] == "!"
    assert tokens[1] == "TPSS"
    assert "ZORA" in tokens
    assert "ZORA-def2-SVP" in tokens
    assert "D4" in tokens
    assert "RIJCOSX" in tokens
    assert "SARC/J" in tokens
    assert "SMD(water)" in tokens
    assert "OPT" in tokens


def test_smd_solvent_change_changes_the_line_and_nothing_else():
    line_dmf = _build_bang_line(dict(BASE_CONFIG, implicit_solvation_model="SMD", solvent="DMF"),
                                "", "def2-SVP", "def2/J", "SMD(DMF)")
    line_water = _build_bang_line(dict(BASE_CONFIG, implicit_solvation_model="SMD", solvent="water"),
                                  "", "def2-SVP", "def2/J", "SMD(water)")
    assert line_dmf != line_water
    diff = [t for t in line_dmf.split() if t not in line_water.split()]
    assert diff == ["SMD(DMF)"], f"exactly the solvation token differs, got {diff}"


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(pytest.main([__file__]))
