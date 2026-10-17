"""Adversarial tests for source_ledger.py (reviewer nacht-s33, package W1).

Attacks the two phrase-list logics in the builder's implementation:
  1. independence (>=2 distinct primaries; same URL counts once; a
     secondary never adds a voice),
  2. conflict detection (polarity flips on a shared stem).

Red = the code answers wrongly. Findings go to the builder with test id.
Green cases pin behaviour I agree with, so a later change that breaks
them is visible.
"""
from __future__ import annotations

import pytest

from delfin.agent import source_ledger as sl


# ---------- independence ----------

def test_same_url_repeated_does_not_make_two_voices():
    # SAME URL twice with kind primary: one voice -> few_sources must fire
    src = [{"url": "https://x.example/a", "kind": "primary",
            "retrieved": "2026-10-17"}] * 2
    issues = sl.review_sources("claim", src)
    types = [i["type"] for i in issues]
    assert "few_sources" in types


def test_same_url_different_case_counts_once():
    # URL comparison lower-cases both sides, so case variants collide.
    src = [
        {"url": "HTTPS://X.EXAMPLE/a", "kind": "primary",
         "retrieved": "2026-10-17"},
        {"url": "https://x.example/A", "kind": "primary",
         "retrieved": "2026-10-17"},
    ]
    issues = sl.review_sources("claim", src)
    types = [i["type"] for i in issues]
    assert "few_sources" in types


# FINDING to nacht-s32 (test id below): source_ledger.py:46-53 counts
# two empty-URL primaries as two voices (empty url never enters
# seen_urls). xfail until fixed.
@pytest.mark.xfail(reason="finding: url-less primaries double-count as two voices, source_ledger.py:46-53", strict=True)
def test_urlless_primaries_cannot_be_two_independent_voices():
    # A URL-less primary cannot be distinguished from another URL-less
    # primary: attacker model = the model citing one paper twice, no URL
    # either time. Two empty-url primaries must NOT count as two voices.
    src = [{"url": "", "kind": "primary", "retrieved": "2026-10-17"}] * 2
    issues = sl.review_sources("claim", src)
    types = [i["type"] for i in issues]
    assert "few_sources" in types, (
        "two URL-less 'primary' entries are indistinguishable from one "
        "source cited twice; they must NOT count as two voices")


def test_secondary_only_claim_never_passes_the_two_voice_gate():
    src = [
        {"url": "https://x.example/a", "kind": "secondary",
         "retrieved": "2026-10-17"},
        {"url": "https://y.example/b", "kind": "secondary",
         "retrieved": "2026-10-17"},
    ]
    issues = sl.review_sources("claim", src)
    types = [i["type"] for i in issues]
    assert "few_sources" in types


def test_unicode_lookalike_urls_are_distinct_byte_strings():
    # Cyrillic 'а' (U+0430) vs Latin 'a': different byte strings, so the
    # deduplication does NOT merge them and two voices is defensible.
    src = [
        {"url": "https://x.example/\u0430", "kind": "primary",
         "retrieved": "2026-10-17"},
        {"url": "https://x.example/a", "kind": "primary",
         "retrieved": "2026-10-17"},
    ]
    issues = sl.review_sources("claim", src)
    types = [i["type"] for i in issues]
    assert "few_sources" not in types


def test_no_sources_at_all_is_flagged():
    issues = sl.review_sources("claim", [])
    assert [i["type"] for i in issues] == ["few_sources"]
    assert issues[0]["n_sources"] == 0


def test_missing_retrieved_date_is_flagged_per_source():
    src = [
        {"url": "https://x.example/a", "kind": "primary"},
        {"url": "https://y.example/b", "kind": "primary",
         "retrieved": "2026-10-17"},
    ]
    issues = sl.review_sources("claim", src)
    types = [i["type"] for i in issues]
    assert "no_retrieved_date" in types
    # The dateless source is named in the issue.
    dateless = [i for i in issues if i["type"] == "no_retrieved_date"]
    assert dateless[0]["url"] == "https://x.example/a"


# ---------- conflicts ----------

def test_conflict_positive_vs_negative_is_found():
    a = ("overlap is negative",
         [{"url": "u1", "kind": "primary", "retrieved": "d"}])
    b = ("overlap is positive",
         [{"url": "u2", "kind": "primary", "retrieved": "d"}])
    assert sl.review_conflicts([a, b]) != []


# FINDING to nacht-s32: source_ledger.py:108-110 needs TWO truthy
# polarity tags; negated-vs-affirmed (neg vs None) never conflicts,
# although the docstring at :79-82 claims this shape is covered.
@pytest.mark.xfail(reason="finding: negation-vs-affirmed claims never conflict, source_ledger.py:108-110", strict=True)
def test_conflict_negation_forms_are_found():
    a = ("the reaction does not proceed",
         [{"url": "u1", "kind": "primary", "retrieved": "d"}])
    b = ("the reaction does proceed",
         [{"url": "u2", "kind": "primary", "retrieved": "d"}])
    assert sl.review_conflicts([a, b]) != []


def test_no_conflict_when_wording_differs_but_polarity_same():
    a = ("the overlap is large",
         [{"url": "u1", "kind": "primary", "retrieved": "d"}])
    b = ("the overlap is enormous",
         [{"url": "u2", "kind": "primary", "retrieved": "d"}])
    assert sl.review_conflicts([a, b]) == []


def test_identical_claims_are_not_a_conflict():
    t = "the band gap is 2.1 eV"
    a = (t, [{"url": "u1", "kind": "primary", "retrieved": "d"}])
    b = (t, [{"url": "u2", "kind": "primary", "retrieved": "d"}])
    assert sl.review_conflicts([a, b]) == []


def test_conflict_scan_survives_empty_list():
    assert sl.review_conflicts([]) == []


def test_review_aggregates_both_issue_families():
    # One claim with too few sources, one pair in conflict: review()
    # must return both families in one list.
    weak = ("weak claim", [{"url": "u1", "kind": "secondary",
                            "retrieved": "d"}])
    a = ("spin density is negative",
         [{"url": "u1", "kind": "primary", "retrieved": "d"}])
    b = ("spin density is positive",
         [{"url": "u2", "kind": "primary", "retrieved": "d"}])
    issues = sl.review([weak, a, b])
    types = {i["type"] for i in issues}
    assert "few_sources" in types
    assert "conflict" in types


def test_source_ledger_class_round_trip():
    ledger = sl.SourceLedger()
    ledger.record("gap is 2.1 eV",
                  [{"url": "u1", "kind": "primary", "retrieved": "d"}])
    assert ledger.claims() == 1
    issues = ledger.review()
    assert [i["type"] for i in issues] == ["few_sources"]
