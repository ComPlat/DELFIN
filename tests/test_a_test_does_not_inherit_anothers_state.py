"""A test reads its own evidence, not the process's.

Second class of machine dependence, and one the universality ratchet in
test_a_test_is_bound_to_delfin_not_to_a_machine.py does NOT see. That
file reads skip CONDITIONS for host lookups; this failure has no
`skipif` at all. It reads shared process state and therefore depends on
what ran before it on this machine.

Measured: `test_the_calc_ledger_reaches_the_completion_check.py::
test_the_provider_serves_the_ledger_the_executor_fills` failed the full
suite twice in a row, on two different branches, with

    {'folder': '/home/<user>/calc/opt_freq_energy_.../01_xtb_preopt',
     'outcome': 'unknown (no run log and no exit code)', 'worst': 'warn'}

in a list it expected to hold only its own entry. Alone it passed at the
baseline and on both branches, twice each. `api_client._doc_executor` is
a module singleton, so that list is whatever the process has
accumulated.

The fix is isolation in conftest (`_isolate_shared_ledgers`), beside the
fixtures that already do this for trust caches and subagent state --
not a change to the test that read it. A test that assumes it is the
only writer is reasonable; a singleton that silently spans the session
is the defect.
"""

from __future__ import annotations

import importlib
import pathlib

import pytest


def _conftest_source() -> str:
    return (pathlib.Path(__file__).resolve().parent
            / "conftest.py").read_text(encoding="utf-8")


def test_every_shared_ledger_is_named_in_one_place():
    """Adding a shared singleton has to be a deliberate act, so the list
    is explicit rather than discovered by reflection."""
    src = _conftest_source()
    assert "_SHARED_LEDGERS = (" in src
    assert '("delfin.agent.api_client", "_doc_executor", "_calc_evidence")' in src


def test_the_ledger_this_test_sees_is_empty():
    """The fixture runs for every test, so this one starts clean even
    though earlier tests in this very file and session wrote to it."""
    from delfin.agent import api_client as A

    assert A._doc_executor._calc_evidence == []


def test_what_this_test_writes_does_not_reach_the_next_one():
    from delfin.agent import api_client as A

    A._doc_executor._calc_evidence.append(
        {"ts": 1.0, "folder": "calc/left_behind",
         "outcome": "succeeded", "worst": "ok"})
    # ... and test_the_ledger_is_still_empty_after below proves it is gone.


def test_the_ledger_is_still_empty_after():
    from delfin.agent import api_client as A

    assert A._doc_executor._calc_evidence == [], (
        "the previous test's entry survived; the isolation does not hold "
        "and every assertion about this list is about the whole session")


def test_entries_are_restored_and_not_discarded():
    """A test's own view is clean, but the process's record must come
    back afterwards: the entries belong to whoever recorded them, and
    dropping them would move the surprise to the next reader.

    Driven through the fixture itself rather than asserted about its
    source, because "clears" and "restores" look identical in a body
    that only clears.
    """
    from delfin.agent import api_client as A

    ledger = A._doc_executor._calc_evidence
    outside = {"ts": 2.0, "folder": "calc/recorded_earlier",
               "outcome": "succeeded", "worst": "ok"}
    ledger.append(outside)
    try:
        gen = _isolate()
        next(gen)                      # enters: the list is cleared
        assert ledger == [], "the fixture did not isolate"
        ledger.append({"ts": 3.0, "folder": "calc/inner",
                       "outcome": "succeeded", "worst": "ok"})
        try:
            next(gen)                  # exits: the list is restored
        except StopIteration:
            pass
        assert ledger == [outside], (
            "the fixture cleared instead of restoring; the entry recorded "
            "before it is gone")
    finally:
        ledger.clear()


def _isolate():
    """The conftest fixture's underlying generator."""
    conftest = importlib.import_module("conftest")
    fixture = conftest._isolate_shared_ledgers
    func = getattr(fixture, "__wrapped__", fixture)
    return func()


def test_the_universality_ratchet_names_this_second_class():
    """The ratchet file claims to be about tests bound to a machine. It
    reads skip conditions only, so this class walked straight past it.
    It has to say so, or the next reader trusts it for more than it
    checks."""
    text = (pathlib.Path(__file__).resolve().parent
            / "test_a_test_is_bound_to_delfin_not_to_a_machine.py"
            ).read_text(encoding="utf-8")
    assert "shared process state" in text, (
        "the ratchet must name the class it does NOT cover")
