"""Grounded and wrong — the case none of the other guards could express.

Field report 20260902-080635. DELFIN had already computed beta_HRS for
404 folders and stored it. The agent computed its own with a guessed
formula, compared the two, found a mean deviation of 17.98% (min -2.33%,
max +47.04%), wrote "deutet auf eine unterschiedliche Formel-Konvention
hin", and shipped ITS values in the CSV under the header `beta_HRS_au`.

Every existing guard passed it, correctly. The numbers came out of the
agent's own script, so they were in the observed-number pool and the
quantity scanner grounded them. The citation checker had real paths. The
answer even carried the caveat — as prose, at the bottom, framed as a
question of convention.

What no guard could see is that two values existed for one field of one
record. That needs no model judgement: it is two floats and a tolerance.

    stored in DELFIN_Data.json     171232.01 au
    delivered in the CSV           180721.43 au

The deviation is also the refutation the agent missed. A convention is a
constant factor; -2.33% to +47.04% is not constant. The formula omitted
the cross terms of the real orientational average (parser.py:575), whose
sign and size depend on the tensor, which is exactly why it varies.
"""

from __future__ import annotations

import pytest

from delfin.agent import verify_guard as vg


@pytest.fixture(autouse=True)
def _fresh_turn():
    vg.reset_keyed_values()
    yield
    vg.reset_keyed_values()


def _saw(output: str, source: str) -> None:
    vg.record_keyed_values(output, source=source)


# The two witnesses, as they actually reached the harness.
STORED = ('{"ground_state_S0": {"hyperpolarizability": {'
          '"beta_HRS_au": 171232.0148, "beta_zzz_au": -294746.6558,'
          ' "beta_tot_au": 410891.8621}}}')
STORED_FROM = ('{"path": "/w/archive/ADMA-2021-BAF1_C_119/DELFIN_Data.json"}')

DELIVERED = (
    "Folder,beta_HRS_au,beta_HRS_esu_30,beta_zzz_au,beta_tot_au\n"
    "ADMA-2021-BAF1_C_119,180721.4333,1561.2886,-294746.6558,410891.8621\n")
DELIVERED_FROM = '{"command": "head -10 /w/archive/beta_hrs_summary.csv"}'


# ---------------------------------------------------------------------------
# The case it exists for
# ---------------------------------------------------------------------------

def test_a_delivered_figure_that_fights_the_stored_one_is_caught():
    _saw(STORED, STORED_FROM)
    _saw(DELIVERED, DELIVERED_FROM)
    flags = vg.scan_for_conflicting_figures()
    assert [f.field for f in flags] == ["beta_HRS_au"]
    assert flags[0].values == (171232.0148, 180721.4333)


def test_the_fields_that_agree_stay_silent():
    """beta_zzz_au and beta_tot_au are identical in both witnesses. A
    guard that flagged them would bury the one that matters."""
    _saw(STORED, STORED_FROM)
    _saw(DELIVERED, DELIVERED_FROM)
    fields = {f.field for f in vg.scan_for_conflicting_figures()}
    assert "beta_zzz_au" not in fields and "beta_tot_au" not in fields


def test_the_feedback_names_both_numbers_and_who_wins():
    _saw(STORED, STORED_FROM)
    _saw(DELIVERED, DELIVERED_FROM)
    msg = vg.conflicting_figure_feedback(vg.scan_for_conflicting_figures())
    assert "171232" in msg and "180721" in msg
    assert "the source wins" in msg
    # The sentence that turns the field case from judgement into arithmetic.
    assert "constant factors" in msg


# ---------------------------------------------------------------------------
# Resolution cannot be talked into existence
# ---------------------------------------------------------------------------

def _one_conflict():
    _saw(STORED, STORED_FROM)
    _saw(DELIVERED, DELIVERED_FROM)
    return vg.scan_for_conflicting_figures()[0]


def test_the_answer_that_shipped_does_not_count_as_resolved():
    """Verbatim from the field report."""
    assert not vg.conflict_is_addressed(_one_conflict(),
        "Die berechneten Werte weichen um durchschnittlich ~18% von den "
        "vor-berechneten ab. Dies deutet auf eine unterschiedliche "
        "Konventions- oder Gewichtungsfaktor-Verwendung hin.")


@pytest.mark.parametrize("answer", [
    "Ich liefere den gespeicherten Wert 171232.01 au.",
    "Der DELFIN-Wert ist 171232.0148 und den nehme ich.",
    "Meine Rechnung ergab 180721.43; ich verwerfe sie zugunsten des Datei-Werts.",
])
def test_naming_one_of_the_two_numbers_resolves_it(answer):
    assert vg.conflict_is_addressed(_one_conflict(), answer)


def test_an_empty_correction_resolves_nothing():
    assert not vg.conflict_is_addressed(_one_conflict(), "")


# ---------------------------------------------------------------------------
# The COST side: what must NOT fire
# ---------------------------------------------------------------------------

CONFIG = '{"path": "/w/project/config.json"}'


def test_a_setting_changed_and_read_back_is_not_a_conflict():
    """The most ordinary thing an agent does. Before the witness rule
    this fired: `timeout` 30 then 60, one file, reported as a
    disagreement. Two readings of one witness are a history."""
    _saw('{"timeout": 30, "retries": 3}', CONFIG)
    _saw('{"timeout": 60, "retries": 3}', CONFIG)
    assert vg.scan_for_conflicting_figures() == []


def test_a_float_setting_changed_in_place_is_not_one_either():
    _saw('{"ratio": 0.50}', CONFIG)
    _saw('{"ratio": 0.75}', CONFIG)
    assert vg.scan_for_conflicting_figures() == []


def test_the_same_value_seen_twice_is_not_a_conflict():
    _saw('{"energy_au": -1705.197734}', '{"path": "/w/calc/run_7/d.json"}')
    _saw("Folder,energy_au\nrun_7,-1705.197734\n",
         '{"command": "cat /w/calc/table.csv"}')
    assert vg.scan_for_conflicting_figures() == []


def test_many_records_with_their_own_values_are_not_conflicts():
    """404 folders each holding a different beta_HRS is the normal shape
    of this work, not 404 contradictions."""
    for i in range(60):
        _saw(f'{{"beta_HRS_au": {1000 + i}.5}}',
             f'{{"path": "/w/calc/sys_{i}/DELFIN_Data.json"}}')
    assert vg.scan_for_conflicting_figures() == []


def test_rounding_is_not_disagreement():
    """0.5% tolerance: the field case was caught at 5.3% and would
    survive an order of magnitude tighter."""
    _saw('{"gap_ev": 1.554000}', '{"path": "/w/calc/sys_2/d.json"}')
    _saw("Folder,gap_ev\nsys_2,1.554100\n", '{"command": "cat t.csv"}')
    assert vg.scan_for_conflicting_figures() == []


def test_nothing_recorded_means_nothing_claimed():
    assert vg.scan_for_conflicting_figures() == []
    assert vg.conflicting_figure_feedback([]) == ""
    assert vg.conflicting_figure_caveat([]) == ""


def test_the_ledger_is_per_turn():
    _saw(STORED, STORED_FROM)
    _saw(DELIVERED, DELIVERED_FROM)
    assert vg.scan_for_conflicting_figures()
    vg.reset_keyed_values()
    assert vg.scan_for_conflicting_figures() == []


def test_a_garbled_result_records_nothing_rather_than_raising():
    assert vg.record_keyed_values(b"\x00\x01not text", source="") == 0
    assert vg.record_keyed_values(None, source=None) == 0
    assert vg.scan_for_conflicting_figures() == []
