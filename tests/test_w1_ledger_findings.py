"""W1 red controls for reviewer findings F1 + F2 (nacht-s33, 647e2279).

F1: two URL-less primaries must NOT count as two independent voices.
F2: a negated claim vs its unnegated form must be a conflict even when
    only one side carries a polarity tag.
"""
from __future__ import annotations

import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))

from delfin.agent.source_ledger import (  # noqa: E402
    SourceLedger,
    review_conflicts,
)


def test_f1_urlless_primaries_cannot_be_two_voices():
    led = SourceLedger()
    led.record("overlap is -0.42", sources=[
        {"url": "", "kind": "primary", "retrieved": "2026-10-18", "quote": ""},
        {"url": "", "kind": "primary", "retrieved": "2026-10-18", "quote": ""},
    ])
    issues = led.review()
    types = {i["type"] for i in issues}
    assert "few_sources" in types


def test_f1_url_required_for_a_primary_voice():
    led = SourceLedger()
    led.record("overlap is -0.42", sources=[
        {"url": "https://a.example/x", "kind": "primary",
         "retrieved": "2026-10-18", "quote": ""},
        {"url": "", "kind": "primary", "retrieved": "2026-10-18", "quote": ""},
    ])
    issues = led.review()
    types = {i["type"] for i in issues}
    assert "few_sources" in types


def test_f2_negated_vs_unnegated_claim_is_a_conflict():
    records = [
        ("the method does not proceed",
         [{"url": "https://a.example/x", "kind": "primary",
           "retrieved": "2026-10-18", "quote": ""}]),
        ("the method does proceed",
         [{"url": "https://b.example/y", "kind": "primary",
           "retrieved": "2026-10-18", "quote": ""}]),
    ]
    conflicts = review_conflicts(records)
    assert conflicts, "negation-vs-None polarity must fire"


def test_f2_agreed_behaviour_still_green():
    # Two claims that are merely worded differently but agree in polarity
    # are NOT a conflict (conservative rule stays).
    records = [
        ("the overlap is -0.42",
         [{"url": "https://a.example/x", "kind": "primary",
           "retrieved": "2026-10-18", "quote": ""}]),
        ("the overlap is -0.42",  # duplicate wording, different source
         [{"url": "https://b.example/y", "kind": "primary",
           "retrieved": "2026-10-18", "quote": ""}]),
    ]
    assert review_conflicts(records) == []
