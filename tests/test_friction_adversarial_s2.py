"""Adversarial tests for package E friction metrics (reviewer s2).

These attack the false-positive and silent-failure surfaces of
``compute_friction`` / ``friction_summary`` (delfin/agent/session_report.py,
builder commit 819bd811). Each test states the DESIRED behavior; a red run is
a finding for the builder, not a spec change.

Entry format is the real trace record (delfin/agent/tool_trace.py:97-107).
"""

from __future__ import annotations

import json

from delfin.agent import session_report as sr


def _bash(ts, command, ok=True, error=""):
    return {
        "ts": ts, "tool": "bash", "input": json.dumps({"command": command}),
        "output": "", "duration_ms": 0, "ok": ok, "error": error,
    }


def _tool(ts, tool, payload, ok=True, error=""):
    return {
        "ts": ts, "tool": tool, "input": payload,
        "output": "", "duration_ms": 0, "ok": ok, "error": error,
    }


_DENIED = "not on the auto-allow list"


# ---------------------------------------------------------------------------
# A1: non-bash denials must not match on empty commands.
# _command() returns "" for every non-bash entry (session_report.py:622-623),
# so the identical-retry branch compares "" == "" and reports an approval
# wait for two UNRELATED calls of the same tool.
# ---------------------------------------------------------------------------

def test_a1_unrelated_non_bash_retry_is_not_an_approval_wait():
    entries = [
        _tool(100.0, "read_file", '{"path": "a.py"}', ok=False, error=_DENIED),
        _tool(160.0, "read_file", '{"path": "b.py"}', ok=True),
    ]
    items = sr.friction_summary(entries)
    waits = [i for i in items if i["kind"] == "approval_wait"]
    assert waits == [], (
        "two unrelated read_file calls matched on empty commands: "
        f"{waits!r}")


def test_a1b_same_file_non_bash_retry_still_counts():
    # The legitimate case must keep working: same tool, same file, later ok.
    entries = [
        _tool(100.0, "read_file", '{"path": "a.py"}', ok=False, error=_DENIED),
        _tool(160.0, "read_file", '{"path": "a.py"}', ok=True),
    ]
    items = sr.friction_summary(entries)
    assert any(i["kind"] == "approval_wait" for i in items) or items == []


# ---------------------------------------------------------------------------
# A2: a later command sharing only a basename is not automatically a rework.
# _object_tokens adds every basename (session_report.py:644) and the rework
# branch fires on ANY shared token (session_report.py:791), so inspecting the
# denied file afterwards (or touching an unrelated file with the same name)
# is reported as "the denial was reworked".
# ---------------------------------------------------------------------------

def test_a2_inspecting_the_denied_file_is_not_a_rework():
    entries = [
        _bash(200.0, "rm scratch.tmp", ok=False, error=_DENIED),
        # agent looks at the file it was refused to delete: not a rework
        _bash(240.0, "sed -n '1,10p' /work/scratch.tmp", ok=True),
    ]
    items = sr.friction_summary(entries)
    reworks = [i for i in items if i["kind"] == "denial_reworked"]
    assert reworks == [], f"inspection misread as rework: {reworks!r}"


def test_a2b_different_file_with_same_basename_is_not_a_rework():
    entries = [
        _bash(200.0, "rm a.py", ok=False, error=_DENIED),
        _bash(300.0, "python other/a.py", ok=True),
    ]
    items = sr.friction_summary(entries)
    reworks = [i for i in items if i["kind"] == "denial_reworked"]
    assert reworks == [], f"basename collision misread as rework: {reworks!r}"


# ---------------------------------------------------------------------------
# A3: one hostile timestamp must not blank the whole report.
# A NaN/Inf gap passes the `gap < idle_min_s` test (NaN comparisons are
# False), then round(nan) raises (session_report.py:703-713) and
# friction_summary's blanket except (session_report.py:833-834) silently
# returns [] — losing every real friction in the session.
# ---------------------------------------------------------------------------

def test_a3_nan_timestamp_does_not_blank_the_report():
    entries = [
        {"ts": float("nan"), "tool": "bash", "input": "{}",
         "output": "", "duration_ms": 0, "ok": True},
        _bash(100.0, "sed -n '1,50p' big.out"),
        _bash(101.0, "sed -n '50,100p' big.out"),
        _bash(102.0, "sed -n '100,150p' big.out"),
    ]
    items = sr.friction_summary(entries)
    assert items, "NaN entry blanked the whole friction report"
    assert any(i["kind"] == "repeated_reads" for i in items)


def test_a3b_inf_timestamp_does_not_blank_the_report():
    entries = [
        {"ts": float("inf"), "tool": "bash", "input": "{}",
         "output": "", "duration_ms": 0, "ok": True},
        _bash(100.0, "sed -n '1,50p' big.out"),
        _bash(101.0, "sed -n '50,100p' big.out"),
        _bash(102.0, "sed -n '100,150p' big.out"),
    ]
    items = sr.friction_summary(entries)
    assert items, "Inf entry blanked the whole friction report"


# ---------------------------------------------------------------------------
# A4: a redirected slice still targets the file it reads.
# _extract_target takes the LAST path-like token (session_report.py:659-667),
# so `sed -n 1,50p big.out > slice.txt` is grouped under slice.txt and real
# repeated slicing of big.out is never reported.
# ---------------------------------------------------------------------------

def test_a4_redirected_slices_are_grouped_under_the_read_file():
    entries = [
        _bash(10.0, "sed -n '1,50p' big.out > s1.txt"),
        _bash(20.0, "sed -n '50,100p' big.out > s2.txt"),
        _bash(30.0, "sed -n '100,150p' big.out > s3.txt"),
    ]
    items = sr.friction_summary(entries, read_min_count=3)
    reads = [i for i in items if i["kind"] == "repeated_reads"]
    assert len(reads) == 1, "redirect hid three slices of big.out"
    assert reads[0]["detail"] == "big.out"


# ---------------------------------------------------------------------------
# A5: robustness contract — malformed entries never raise and never poison
# the patterns around them (compute_friction must be total).
# ---------------------------------------------------------------------------

def test_a5_malformed_entries_never_raise():
    entries = [
        None,
        {},
        {"ts": "not-a-number", "tool": "bash", "input": "{}"},
        {"ts": 5.0, "tool": 42, "input": None, "ok": "yes"},
        _bash(10.0, "sed -n '1,50p' big.out"),
        _bash(11.0, "sed -n '50,100p' big.out"),
        _bash(12.0, "sed -n '100,150p' big.out"),
    ]
    items = sr.friction_summary(entries)
    assert any(i["kind"] == "repeated_reads" for i in items)
