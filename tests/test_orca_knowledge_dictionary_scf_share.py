"""The approximation/SCF/task share: honest status, dash-tolerant lookup.

A grounded record states what its keyword does; a placeholder (an OPI
catalog docstring, "see manual page") is a draft. The input spelling
without dashes finds the record the manual writes with them.
"""
from __future__ import annotations

from delfin.orca_knowledge.dictionary import explain, records

_PLACEHOLDERS = ("OPI catalog:", "see manual page", "see page",
                 "Not documented under this exact spelling")


def test_the_share_is_loaded():
    assert explain("stabmaxdim") is not None


def test_a_grounded_record_states_its_purpose():
    for rec in records():
        if rec.get("status") != "grounded":
            continue
        purpose = rec.get("purpose") or ""
        assert purpose and not any(p in purpose for p in _PLACEHOLDERS), \
            rec.get("keyword")


def test_a_dashed_manual_spelling_is_found_without_dashes():
    dashed = next((r["keyword"] for r in records()
                   if "-" in r["keyword"]), None)
    assert dashed is not None
    found = explain(dashed.replace("-", ""))
    assert found is not None
