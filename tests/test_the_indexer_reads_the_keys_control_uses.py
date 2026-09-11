"""An operator answering "wie funktioniert aktuell der co2 coordinator"
from the code noticed that the calculation indexer reads a key named
``CO2_coordinator`` while CONTROL calls it ``co2_coordination``
(2026-09-11). The same held for ``tadf_xTB`` and ``hyperpol_xTB``: the
indexer compared exact, case-sensitive names that no CONTROL file
writes, so no indexed run ever carried those modules. It now compares
keys the way config.py does, and takes "on" for what define.py
documents.
"""

from __future__ import annotations

from delfin.doc_server import calc_indexer as ci


_CONTROL = """\
NAME=cu_complex
method=OCCUPIER
functional=PBE0
main_basisset=def2-SVP
solvent=water
co2_coordination=on
co2_species_delta=0
tadf_xTB=yes
hyperpol_xTB=no
"""


def test_the_module_flags_are_read_by_the_keys_control_uses():
    rec = ci._extract_from_control_txt(_CONTROL)
    mods = rec["modules"]
    assert "CO2" in mods, mods
    assert "TADF_xTB" in mods, mods
    assert "OCCUPIER" in mods, "method=OCCUPIER is how a run says it"
    assert "hyperpol_xtb" not in mods


def test_the_legacy_spellings_still_count():
    rec = ci._extract_from_control_txt("CO2_coordinator=yes\nTADF_xTB_module=yes\nhyperpol_xtb_module=yes\n")
    assert {"CO2", "TADF_xTB", "hyperpol_xtb"} <= set(rec["modules"])


def test_key_comparison_ignores_case_and_separators():
    rec = ci._extract_from_control_txt("Co2-Coordination = ON\n")
    assert "CO2" in rec["modules"]
    assert ci._norm_key("CO2_coordination") == ci._norm_key("co2 coordination") == "co2coordination"


def test_on_is_a_yes_and_off_is_not():
    assert ci._is_yes("on") and ci._is_yes("ON") and ci._is_yes("yes")
    assert not ci._is_yes("off") and not ci._is_yes("no") and not ci._is_yes("")


def test_the_flag_helper_returns_the_first_present_value():
    kv = {"tadf_xTB": "yes", "TADF_xTB_module": "no"}
    assert ci._flag(kv, "tadf_xTB", "TADF_xTB_module") == "yes"
    assert ci._flag({}, "tadf_xTB") == ""
