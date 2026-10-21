"""W1 red control: the source ledger flags weakly-sourced claims.

Runs against delfin/agent/source_ledger.py (owned file). Every case here
must FAIL on the commit before source_ledger.py exists and PASS after.
"""
from __future__ import annotations

import sys
from datetime import date
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))

from delfin.agent.source_ledger import (  # noqa: E402
    SourceLedger,
    check_ledger,
)


def _src(url: str, kind: str = "primary", retrieved: str = "2026-10-18",
         quote: str = ""):
    return {"url": url, "kind": kind, "retrieved": retrieved, "quote": quote}


def test_two_independent_primary_sources_is_clean():
    led = SourceLedger()
    led.record("B3LYP overestimates frequencies by ~1.5%",
               sources=[_src("https://a.example/x"), _src("https://b.example/y")])
    issues = led.review()
    assert issues == [], issues


def test_single_source_is_flagged():
    led = SourceLedger()
    led.record("CAM-B3LYP doubles the CT excitation energy",
               sources=[_src("https://a.example/x")])
    issues = led.review()
    assert len(issues) == 1
    issue = issues[0]
    assert issue["type"] == "few_sources"
    assert issue["n_sources"] == 1


def test_same_url_twice_is_not_two_sources():
    led = SourceLedger()
    led.record("Delta-SCF captures S1 correctly",
               sources=[_src("https://a.example/x"), _src("https://a.example/x")])
    issues = led.review()
    assert issues and issues[0]["type"] == "few_sources"


def test_secondary_backing_primary_is_not_independence():
    led = SourceLedger()
    led.record("EOM-CCSD vertical errors are ~0.3 eV",
               sources=[_src("https://a.example/paper"),
                        _src("https://b.example/blog", kind="secondary")])
    issues = led.review()
    # One primary only: a secondary restating the primary is the same voice.
    assert issues and issues[0]["type"] == "few_sources"


def test_missing_retrieved_date_is_flagged():
    led = SourceLedger()
    led.record("solvent parameters from paper X",
               sources=[_src("https://a.example/x", retrieved=""),
                        _src("https://b.example/y", retrieved="")])
    issues = led.review()
    types = {i["type"] for i in issues}
    assert "no_retrieved_date" in types


def test_conflicting_claims_are_flagged():
    led = SourceLedger()
    led.record("SCS overlap is negative",
               sources=[_src("https://a.example/x")])
    led.record("SCS overlap is positive",
               sources=[_src("https://b.example/y")])
    issues = led.review()
    types = {i["type"] for i in issues}
    assert "conflict" in types


def test_check_ledger_reports_everything_without_recording():
    raw = [
        ("claim A", [_src("https://a.example/x")]),
        ("claim B", [_src("https://b.example/y", retrieved="")]),
    ]
    report = check_ledger(raw)
    assert report["claims"] == 2
    types = {i["type"] for i in report["issues"]}
    assert "few_sources" in types and "no_retrieved_date" in types
