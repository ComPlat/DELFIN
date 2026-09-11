"""Three runs under one folder were one record to the verifier.

The record a figure belongs to is read from the tool call that produced
it, and that call arrives rendered -- "{'folder': '/tmp/x/spectra/dye_a'}".
Split on "/" alone, the deepest segment is "dye_a'}", which is no path
segment, so the parent "spectra" was taken instead: every run under it
shared one record, and the verifier told the model its three honest
answers were one record contradicting itself. It re-read all three
files to resolve a conflict that did not exist (2026-09-11).
"""

from __future__ import annotations

import pytest

from delfin.agent import verify_guard as V


@pytest.mark.parametrize("source,record", [
    ("{'folder': '/tmp/x/spectra/dye_a'}", "dye_a"),
    ("{'path': '/tmp/x/spectra/dye_a/run.out'}", "dye_a"),
    ("{'folder': '/tmp/x/spectra/dye_a', 'n': 3}", "dye_a"),
    ('{"folders": "/tmp/x/spectra/dye_a,/tmp/x/spectra/dye_b"}', "dye_b"),
    ("/tmp/x/spectra/dye_a", "dye_a"),
    ("/tmp/x/spectra/dye_a/run.out", "dye_a"),
    ("/tmp/proj-x/calc/", "proj-x"),
])
def test_the_record_is_the_deepest_folder_of_the_calls_path(source, record):
    assert V._record_id(source) == record


def test_two_runs_under_one_folder_do_not_contradict_each_other():
    V.reset_keyed_values()
    V.record_keyed_values('{"gap_ev": 3.54, "homo_ev": -5.85}', source="{'folder': '/tmp/x/spectra/dye_a'}")
    V.record_keyed_values('{"gap_ev": 4.90, "homo_ev": -6.53}', source="{'folder': '/tmp/x/spectra/dye_b'}")
    assert V.scan_for_conflicting_figures() == []


def test_one_run_read_twice_with_different_values_still_is_a_conflict():
    V.reset_keyed_values()
    V.record_keyed_values('{"gap_ev": 3.54}', source="{'folder': '/tmp/x/spectra/dye_a'}")
    V.record_keyed_values('{"gap_ev": 3.91}', source="{'path': '/tmp/x/spectra/dye_a/DELFIN_Data.json'}")
    flags = V.scan_for_conflicting_figures()
    assert flags and flags[0].record == "dye_a" and flags[0].field == "gap_ev"
