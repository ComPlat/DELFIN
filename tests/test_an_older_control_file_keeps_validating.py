"""CONTROL files written from older templates validate as they did.

Checked against every archived CONTROL file, old validator against new: 264
files that validate now failed before, 254 of them Jerome's, on two lines an
older template wrote itself --

* ``states=[S0,S1,T1,...]``: the template's help line listed S0 among the
  states, and S0 is computed whenever ESD is on anyway;
* ``OCCUPIER_method=auto|manually``: an older template shipped the choice,
  and the template reader takes its first option.

And the state list, like the transition lists, is only read by the ESD
module: with ESD off it is inert.
"""

from __future__ import annotations

from delfin.config import set_control_value, validate_control_text
from delfin.common.control_validator import validate_control_config
from delfin.define import TEMPLATE


def _control(**values: str) -> str:
    text = set_control_value(set_control_value(TEMPLATE, "charge", "0"), "solvent", "water")
    for key, value in values.items():
        text = set_control_value(text, key, value)
    return text


def _errors(text: str) -> list[str]:
    return [e for e in validate_control_text(text, converts_smiles=False)]


_REQUIRED = {"charge": 0, "method": "classic"}


def test_s0_in_the_state_list_is_the_state_esd_computes_anyway():
    listed = _control(ESD_modul="yes", states="[S0,S1,T1,T2]")

    assert _errors(listed) == _errors(_control(ESD_modul="yes", states="[S1,T1,T2]"))
    assert validate_control_config({**_REQUIRED, "ESD_modul": "yes", "states": "[S0,S1,T1]"})["states"] == ["S1", "T1"]


def test_the_state_list_is_inert_without_esd():
    assert _errors(_control(ESD_modul="no", states="[S0,S1,Q7]")) == _errors(_control(ESD_modul="no"))


def test_a_template_that_shipped_the_choice_reads_its_first_option():
    assert _errors(_control(OCCUPIER_method="auto|manually")) == _errors(_control(OCCUPIER_method="auto"))
    assert validate_control_config({**_REQUIRED, "OCCUPIER_method": "auto|manually"})["OCCUPIER_method"] == "auto"
