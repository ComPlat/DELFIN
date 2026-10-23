"""The four ESD lists are read the same way everywhere, brackets or not.

The template shows them in brackets (states=[S1,T1,S2,T2], ISCs=[S1>T1,T1>S1],
ICs=[S1>S0], emission_rates=[f,p]).  states and ICs were read with the
brackets stripped; ISCs and emission_rates were split at the commas with the
brackets still on ("[S1>T1", "T1>S1]", "[f", "p]"), so a file that switched
ESD on computed the states and the IC and left out every ISC and emission
rate without a word -- while the Submit tab said the lists were "not set".
The 1 103 archived ESD runs write these lists without brackets or leave them
empty; those read exactly as before.
"""

from __future__ import annotations

import logging

import pytest

from delfin.common.control_validator import esd_list_items, validate_control_config
from delfin.config import get_esd_hints, set_control_value, validate_control_text
from delfin.define import TEMPLATE
from delfin.esd_module import parse_emission_rates, parse_esd_config


@pytest.mark.parametrize("iscs, rates", [
    ("[S1>T1,T1>S1]", "[f,p]"),
    ("S1>T1,T1>S1", "f,p"),
    (["S1>T1", "T1>S1"], ["f", "p"]),
    ("['S1>T1', 'T1 > S1']", "f p"),
])
def test_every_spelling_gives_the_same_jobs(iscs, rates):
    config = {"ESD_modul": "yes", "states": "[S1,T1]", "ISCs": iscs, "ICs": "[S1>S0]", "emission_rates": rates}

    enabled, states, isc_list, ics = parse_esd_config(config)

    assert enabled
    assert states == ["S0", "S1", "T1"]
    assert isc_list == ["S1>T1", "T1>S1"]
    assert ics == ["S1>S0"]
    assert parse_emission_rates(config) == {"f", "p"}


def test_a_transition_that_is_not_an_isc_is_skipped_with_its_reason(caplog):
    config = {"ESD_modul": "yes", "ISCs": "S1>T1,S1>S2,S1", "emission_rates": "f,x"}
    with caplog.at_level(logging.WARNING, logger="delfin.esd_module"):
        _, _, iscs, _ = parse_esd_config(config)
        rates = parse_emission_rates(config)

    assert iscs == ["S1>T1"]
    assert rates == {"f"}
    assert "keeps its spin" in caplog.text and "is not a transition" in caplog.text
    assert "'x'" in caplog.text


def test_the_validator_reads_brackets_and_refuses_only_spelling():
    template = set_control_value(TEMPLATE, "ESD_modul", "yes")
    assert not [e for e in validate_control_text(template) if "ISC" in e or "emission" in e]

    bad = set_control_value(template, "ISCs", "[S1-T1]")
    assert any("ISCs" in e and "not a transition" in e for e in validate_control_text(bad))
    bad = set_control_value(template, "emission_rates", "[f,q]")
    assert any("emission_rates" in e for e in validate_control_text(bad))


def test_with_esd_off_the_lists_are_inert():
    validated = validate_control_config({"charge": 0, "method": "classic", "ESD_modul": "no",
                                         "ISCs": "[junk]", "emission_rates": "zz"})

    assert validated["ISCs"] == "" and validated["emission_rates"] == ""


def test_the_submit_tab_says_what_will_be_computed():
    hints = get_esd_hints(set_control_value(TEMPLATE, "ESD_modul", "yes"))

    assert "ESD_modul=yes computes: states S0,S1,T1,S2,T2; ISC S1>T1,T1>S1; IC S1>S0; emission f,p" in hints
    assert not [h for h in hints if "not set" in h and ("ISCs" in h or "emission_rates" in h)]


def test_an_archived_plain_list_reads_as_before():
    config = {"ESD_modul": "yes", "states": "S1,T1", "ISCs": "S1>T1,T1>S1", "ICs": "", "emission_rates": "f,p"}

    assert parse_esd_config(config) == (True, ["S0", "S1", "T1"], ["S1>T1", "T1>S1"], [])
    assert parse_emission_rates(config) == {"f", "p"}
    assert esd_list_items("") == [] and esd_list_items("[]") == [] and esd_list_items(None) == []
