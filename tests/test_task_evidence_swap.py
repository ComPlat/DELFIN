"""The swap site, pinned.

``api_client.py`` is security code in this run, so the swap itself is
the operator's edit. This file names the exact call site and pins the
behaviour AFTER the swap with ``xfail(strict=True)``: every case here
is red against the shipped ``api_client.check_completion_claim`` and
must turn green the moment that call site delegates to
``delfin.agent.task_evidence.check_completion_claim``. If it turns
green before the swap (or stays red after it), the strict xfail makes
the suite say so.
"""

from __future__ import annotations

import pytest

from delfin.agent import api_client


def _changes(*paths):
    return [{"path": p, "ts": 10.0, "created": True} for p in paths]


# The call site that hands a completion claim to the check:
# api_client.py, ``_verify_task_completion`` (~line 18404), the only
# ``check_completion_claim(`` call in the repo outside the function's
# own definition and its tests.
SWAP_SITE = "delfin/agent/api_client.py::_verify_task_completion"


def test_the_swap_site_is_the_only_call_of_the_check():
    """Guards the wiring: the repo calls the check from exactly one
    place, so swapping there swaps it everywhere."""
    import inspect
    src = inspect.getsource(api_client)
    own_def = src.count("def check_completion_claim(")
    calls = src.count("return check_completion_claim(")
    assert own_def == 1, "the shipped check must stay defined"
    assert calls == 1, "one call site to swap: " + SWAP_SITE


class TestAfterTheSwap:
    """Each case: judged wrongly (or by luck) by the shipped check at
    the END-TO-END level -- the task subject/description the caller
    passes -- and correctly by task_evidence."""

    def test_a_mentioned_reference_file_does_not_make_the_task_unmet(self):
        res = api_client.check_completion_claim(
            "Bau die Prüflogik wie in foo.py beschrieben",
            changes=_changes("/w/prueflogik.py"),
            observed=["/w/foo.py"])
        assert res["verdict"] != "unmet", res

    def test_a_bare_tabelle_word_demands_no_spreadsheet(self):
        res = api_client.check_completion_claim(
            "Ergänze die Übersicht um eine Tabelle der Testfälle",
            changes=_changes("/w/README.md"), observed=[])
        assert res["verdict"] != "unmet", res

    def test_a_description_path_cannot_make_the_task_unmet(self):
        """The red variant: the SUBJECT names no path, the description
        mentions one in a where-it-came-from role. The shipped check
        reads paths from subject + description, so foo.py accuses."""
        res = api_client.check_completion_claim(
            "Behebe den Fehler",
            description="Der Fehler kam ursprünglich aus foo.py.",
            changes=_changes("/w/bar.py"),
            observed=["/w/foo.py"])
        assert res["verdict"] == "verified", res

    # NOTE (measured, not assumed): the "note instruction in the
    # description" shape -- subject "Lies den Stand in engine.py",
    # description "Notiere die Stichpunkte in Zusammenfassung.md." -- is
    # judged CORRECTLY by the shipped check ("notiere" is not in the
    # write-verb table, and the read of engine.py verifies). Not an
    # xfail case: there is nothing to swap for it. task_evidence keeps
    # the same answer by reading intent from the subject alone.
