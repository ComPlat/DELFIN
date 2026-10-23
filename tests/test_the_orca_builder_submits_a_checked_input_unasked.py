"""`/orca submit` sends a checked ORCA Builder input without the dialog.

Input: the ORCA Builder's input preview, PAL, MaxCore and Slurm time, as
the agent's `/orca submit` action reads them (``_orca_builder_verdict``).

`/orca submit` clicks the builder's submit button, and the builder holds
the whole input. The same rule as the other routes applies (user decision
2026-10-22): an input whose charge and multiplicity fit its electrons,
with no validator error, whose PAL and time fit the host or the
partition, goes without the dialog; anything else is asked about with the
reason. `/submit` (a DELFIN pipeline folder) is not read and keeps its
dialog.
"""
from __future__ import annotations

from types import SimpleNamespace as W

import pytest

from delfin.agent import calc_policy as C
from delfin.dashboard.tab_agent import _orca_builder_verdict

WATER = "O 0 0 0\nH 0 0.76 0.59\nH 0 -0.76 0.59\n"


def _refs(mult=1, pal=4, time="02:00:00", inp_pal=""):
    text = (f"! r2SCAN-3c Opt {inp_pal}\n%maxcore 2000\n"
            f"* xyz 0 {mult}\n{WATER}*\n")
    return {"orca_preview": W(value=text), "orca_pal": W(value=pal),
            "orca_maxcore": W(value=2000), "orca_slurm_time": W(value=time)}


@pytest.fixture(autouse=True)
def cluster(monkeypatch):
    monkeypatch.setattr(C, "host_kind", lambda env=None: "login_node")
    monkeypatch.setattr(C, "partition_limits", lambda name: {
        "PartitionName": "cpu", "MaxTime": "3-00:00:00",
        "MaxCPUsPerNode": "64"})
    monkeypatch.setattr(C, "queued_jobs", lambda: 0)


def test_a_checked_input_goes_unasked():
    v = _orca_builder_verdict(_refs())
    assert v.calculation and v.admit, v.reason


def test_a_multiplicity_the_electrons_cannot_have_is_asked_about():
    v = _orca_builder_verdict(_refs(mult=2))
    assert not v.admit and "multiplicity 2" in v.reason


def test_a_time_beyond_the_partition_is_asked_about():
    v = _orca_builder_verdict(_refs(time="5-00:00:00"))
    assert not v.admit and "partition limit" in v.reason


def test_the_input_s_own_pal_counts(monkeypatch):
    v = _orca_builder_verdict(_refs(pal=4, inp_pal="PAL128"))
    assert not v.admit and "128" in v.reason


def test_a_builder_without_a_preview_is_not_admitted():
    assert _orca_builder_verdict({}) is None
