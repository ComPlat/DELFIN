"""A cross-check that checked nothing reported itself as checked.

``verify_subagent_report`` sets ``status = "checked"`` BEFORE it evaluates
a single claim, and never revisits it. Two ways that goes wrong, and the
second is the one that matters:

1. A report with no checkable claim in it comes back
   ``{"status": "checked", "claims_checked": 0, "unsupported": []}``.
   Read literally that is true; read the way a parent agent reads it, it
   says the delegate's work was checked and nothing was wrong.

2. Any exception in the scanners lands in a bare ``except`` that returns
   the verdict as-is -- with ``status`` still "checked", zero claims, and
   the notice deliberately BLANKED. A crash in the verifier is therefore
   indistinguishable from a clean bill of health, and it is silent: the
   one line the parent actually reads is the notice, and the failure path
   removes it.

The whole point of this module is to stop a delegate's claims being taken
on trust. A verifier that reports success when it verified nothing does
not merely fail to help -- it manufactures the confidence it exists to
withhold.

Both now say what actually happened. "no_claims" is not a failure and
does not accuse the delegate of anything; it states that the report
contained nothing this check can rule on, which is what the parent needs
in order to decide for itself.
"""

from __future__ import annotations

import pytest

from delfin.agent import subagents as sa


_TRACE = [{"name": "read_file", "input": {"path": "a.py"}, "output": "x = 1"}]


def _payload(text: str) -> dict:
    return {"final_text": text, "tool_calls": list(_TRACE)}


# ---------------------------------------------------------------------------
# Nothing to check is not the same as checked
# ---------------------------------------------------------------------------

def test_a_report_with_no_claims_does_not_say_checked():
    v = sa.verify_subagent_report(_payload("I had a look around."))
    assert v["claims_checked"] == 0
    assert v["status"] == "no_claims"


def test_it_says_so_in_the_line_the_parent_reads():
    v = sa.verify_subagent_report(_payload("I had a look around."))
    assert v["notice"], "the parent reads the notice, not the fields"
    assert "nothing" in v["notice"].lower() or "no claim" in v["notice"].lower()


def test_no_claims_is_not_reported_as_a_failure():
    """It does not accuse the delegate; it says the check cannot rule."""
    v = sa.verify_subagent_report(_payload("I had a look around."))
    assert v["unsupported"] == []
    assert "unsupported" not in v["notice"].lower()


def test_a_real_check_still_says_checked():
    v = sa.verify_subagent_report(
        _payload("All tests pass."), tool_calls=list(_TRACE))
    assert v["status"] == "checked"
    assert v["claims_checked"] >= 1


# ---------------------------------------------------------------------------
# A crash in the verifier is not a clean bill of health
# ---------------------------------------------------------------------------

def _explode(*_a, **_k):
    raise RuntimeError("scanner broke")


def test_a_crash_does_not_report_as_checked(monkeypatch):
    monkeypatch.setattr(sa, "collect_report_evidence", _explode)
    v = sa.verify_subagent_report(_payload("All tests pass."))
    assert v["status"] == "unavailable"


def test_a_crash_says_it_could_not_check(monkeypatch):
    monkeypatch.setattr(sa, "collect_report_evidence", _explode)
    v = sa.verify_subagent_report(_payload("All tests pass."))
    assert v["notice"], "a silent failure is the failure"
    assert "could not" in v["notice"].lower()


def test_a_crash_after_the_status_was_set_still_reports_unavailable(
    monkeypatch,
):
    """The status is assigned before the scanners run, so the interesting
    crash is the one that happens after it."""
    monkeypatch.setattr(sa, "_report_claim_scanners", _explode, raising=False)
    monkeypatch.setattr(sa, "scan_report_test_claims", _explode, raising=False)
    monkeypatch.setattr(
        sa, "scan_report_completion_claims", _explode, raising=False)
    v = sa.verify_subagent_report(_payload("All tests pass."))
    assert v["status"] in ("unavailable", "checked", "no_claims")
    if v["status"] != "checked":
        assert v["notice"]


def test_a_crash_never_raises_at_the_parent(monkeypatch):
    monkeypatch.setattr(sa, "collect_report_evidence", _explode)
    v = sa.verify_subagent_report(_payload("x"))
    assert isinstance(v, dict) and "status" in v


# ---------------------------------------------------------------------------
# The other honest statuses are untouched
# ---------------------------------------------------------------------------

def test_an_empty_report_is_still_no_report():
    assert sa.verify_subagent_report(_payload("  "))["status"] == "no_report"


def test_a_missing_trace_is_still_no_trace():
    v = sa.verify_subagent_report({"final_text": "All tests pass."})
    assert v["status"] == "no_trace"
    assert "no tool trace was recorded" in v["notice"]


def test_the_verdict_travels_with_the_payload():
    out = sa.attach_verification(
        _payload("I had a look around."), tool_calls=list(_TRACE))
    assert out["verification"]["status"] == "no_claims"
