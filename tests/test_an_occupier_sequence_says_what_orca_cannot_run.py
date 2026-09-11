"""An OCCUPIER sequence entry ORCA cannot run is named before the run, not after it.

``BrokenSym NA,NB`` converges the high-spin state with NA+NB unpaired
electrons and flips it to the broken-symmetry state of multiplicity
|NA-NB|+1 (ORCA 6.1.1 manual, 5.27), so NA+NB must have the parity of the
electron count.  An older sequence profile had BrokenSym 3,2 and 5,2 in the
even sequence and 4,2 and 6,2 in the odd one; 146 archived OCCUPIER runs in
58 jobs stopped on "multiplicity (4) is even and number of electrons (242) is
even -> impossible".  Of the 3 724 archived OCCUPIER CONTROL files, 260 carry
such entries, and every one of the 58 jobs is among them.  The current
template has none.  A hint, not an error: those files ran, and still run.
"""

from __future__ import annotations

import logging

from delfin.config import get_occupier_hints, read_control_file, set_control_value, validate_control_text
from delfin.define import TEMPLATE

OLD_PROFILE = """
even electron number:
even_seq = [
  {"index": 1, "m": 1, "BS": ""},
  {"index": 2, "m": 1, "BS": "1,1", "from": 1},
  {"index": 3, "m": 1, "BS": "2,2", "from": 2},
  {"index": 4, "m": 3, "BS": "", "from": 1},
  {"index": 5, "m": 3, "BS": "3,1", "from": 4},
  {"index": 6, "m": 3, "BS": "3,2", "from": 5},
  {"index": 7, "m": 5, "BS": "", "from": 1},
  {"index": 8, "m": 5, "BS": "5,1", "from": 7},
  {"index": 9, "m": 5, "BS": "5,2", "from": 8}
]
odd electron number:
odd_seq = [
  {"index": 1, "m": 2, "BS": ""},
  {"index": 2, "m": 2, "BS": "2,1", "from": 1},
  {"index": 3, "m": 2, "BS": "3,2", "from": 2},
  {"index": 4, "m": 4, "BS": "", "from": 1},
  {"index": 5, "m": 4, "BS": "4,1", "from": 4},
  {"index": 6, "m": 4, "BS": "4,2", "from": 5},
  {"index": 7, "m": 6, "BS": "", "from": 1},
  {"index": 8, "m": 6, "BS": "6,1", "from": 7},
  {"index": 9, "m": 6, "BS": "6,2", "from": 8}
]
"""


def _occupier(text: str) -> str:
    return set_control_value(text, "method", "OCCUPIER")


def test_the_template_sequences_are_all_runnable():
    assert get_occupier_hints(_occupier(TEMPLATE)) == []


def test_the_old_profile_is_named_entry_by_entry():
    hints = get_occupier_hints(_occupier("method=OCCUPIER\ncharge=0\n" + OLD_PROFILE))

    assert len(hints) == 4
    assert any("even_seq index 6: BrokenSym 3,2 has 5 unpaired electrons" in h for h in hints)
    assert any("even_seq index 9: BrokenSym 5,2" in h for h in hints)
    assert any("odd_seq index 6: BrokenSym 4,2 has 6 unpaired electrons" in h for h in hints)
    assert any("odd_seq index 9: BrokenSym 6,2" in h for h in hints)


def test_a_broken_symmetry_state_listed_as_another_multiplicity_is_named():
    text = "method=OCCUPIER\neven_seq = [\n  {\"index\": 1, \"m\": 3, \"BS\": \"2,2\"}\n]\n"

    assert get_occupier_hints(text) == [
        "sequence, even_seq index 1: BrokenSym 2,2 is a multiplicity-1 state, listed as m=3"]


def test_it_is_a_hint_and_the_file_still_validates():
    old = _occupier(TEMPLATE) + OLD_PROFILE

    assert validate_control_text(old) == validate_control_text(_occupier(TEMPLATE))


def test_without_occupier_the_sequences_are_not_used_and_not_named():
    assert get_occupier_hints("method=classic\n" + OLD_PROFILE) == []


def test_a_run_says_it_in_its_log(tmp_path, caplog):
    control = tmp_path / "CONTROL.txt"
    control.write_text(set_control_value(set_control_value(_occupier(TEMPLATE), "charge", "0"),
                                         "solvent", "water") + OLD_PROFILE)
    with caplog.at_level(logging.WARNING, logger="delfin.config"):
        read_control_file(str(control))

    assert "BrokenSym 3,2 has 5 unpaired electrons" in caplog.text
