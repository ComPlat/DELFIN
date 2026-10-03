"""Friction metrics for package E — computed from the recorded tool trace only.

Each of the four wave-10 friction patterns has one test built on a fabricated
sample trace (no live session, no wall-clock): lost time on an approval wait,
a denied action reworked in another form, repeated partial reads of one file,
and idle time after a job-starting call. A clean trace must show no friction.

Entry format is the real per-session trace record
(``delfin.agent.tool_trace:97-107``): ``{ts, tool, input, output,
duration_ms, ok, error}``; the bash command is read from the JSON ``input``
like ``tool_trace.command_of`` does.
"""

from __future__ import annotations

import json

import pytest

from delfin.agent import session_report as sr


def _bash(ts, command, ok=True, error=""):
    """One bash trace entry with a JSON ``{"command": ...}`` input."""
    return {
        "ts": ts, "tool": "bash", "input": json.dumps({"command": command}),
        "output": "", "duration_ms": 0, "ok": ok, "error": error,
    }


def _other(ts, tool, payload, ok=True, error=""):
    """A non-bash entry; payload is the literal JSON the input carried."""
    return {
        "ts": ts, "tool": tool, "input": payload,
        "output": "", "duration_ms": 0, "ok": ok, "error": error,
    }


# A refusal denial uses the same error text the permission gate records.
_DENIED = 'not on the auto-allow list'


# ---------------------------------------------------------------------------
# approval_wait: a refusal followed later by the identical successful action.
# ---------------------------------------------------------------------------

def _approval_wait_trace():
    return [
        _bash(100.0, "sed -n '1,50p' job.log", ok=False, error=_DENIED),
        _bash(115.0, "grep Energy job.log", ok=True),
        _bash(4700.0, "sed -n '1,50p' job.log", ok=True),  # retried after a wait
    ]


def test_approval_wait_counts_the_wait_and_wait_time():
    items = sr.friction_summary(_approval_wait_trace())
    waits = [i for i in items if i["kind"] == "approval_wait"]
    # the gap 4700 - 100 = 4600s is the wait between the denied call and its
    # successful identical retry
    assert len(waits) == 1
    w = waits[0]
    assert w["detail"] == "sed -n '1,50p' job.log"
    assert w["unit"] == "s"
    assert w["magnitude"] == 4600


# ---------------------------------------------------------------------------
# denial_reworked: a refusal, then the SAME end retried in another form.
# ---------------------------------------------------------------------------

def _reworked_trace():
    return [
        # denied: cd-wrapped rm is refused
        _bash(200.0, "cd /work && rm scratch.tmp", ok=False, error=_DENIED),
        # reworked: rm via a different wrapper, same target
        _bash(240.0, "rm -f /work/scratch.tmp", ok=True),
    ]


def test_denial_reworked_detects_another_form_of_the_same_action():
    items = sr.friction_summary(_reworked_trace())
    reworks = [i for i in items if i["kind"] == "denial_reworked"]
    assert len(reworks) == 1
    r = reworks[0]
    assert "/work/scratch.tmp" in r["detail"]
    assert r["unit"] == "n"


# ---------------------------------------------------------------------------
# repeated_reads: the same file sliced many times instead of one read.
# ---------------------------------------------------------------------------

def _repeated_reads_trace():
    cmds = [
        "sed -n '1,50p' big.out",
        "sed -n '50,100p' big.out",
        "head -n 20 big.out",
        "tail -n 30 big.out",
    ]
    return [_bash(10.0 + i, c) for i, c in enumerate(cmds)]


def test_repeated_partial_reads_of_one_file_are_reported():
    items = sr.friction_summary(_repeated_reads_trace(), read_min_count=3)
    reads = [i for i in items if i["kind"] == "repeated_reads"]
    assert len(reads) == 1
    r = reads[0]
    assert r["detail"] == "big.out"
    assert r["magnitude"] >= 3
    assert r["unit"] == "reads"


# ---------------------------------------------------------------------------
# idle_gap: a long dead stretch, here after a job-starting call.
# ---------------------------------------------------------------------------

def _idle_after_job_trace():
    return [
        _other(400.0, "bash_background",
               '{"command": "run_long_calc"}'),
        _other(3700.0, "bash_status", '{"job_id": 1}'),  # >idle_min_s later
    ]


def test_idle_after_a_job_start_is_reported():
    items = sr.friction_summary(_idle_after_job_trace(), idle_min_s=300.0)
    idle = [i for i in items if i["kind"] == "idle_gap"]
    assert len(idle) == 1
    g = idle[0]
    assert g["unit"] == "s"
    assert g["magnitude"] == 3300


# ---------------------------------------------------------------------------
# A clean trace shows no friction at all.
# ---------------------------------------------------------------------------

def _clean_trace():
    return [
        _bash(0.0, "git log --oneline -3"),
        _bash(5.0, "pytest -q tests/test_x.py"),
        _bash(12.0, "git diff --stat"),
    ]


def test_clean_trace_shows_no_friction():
    items = sr.friction_summary(_clean_trace(), idle_min_s=300.0)
    assert items == []


# ---------------------------------------------------------------------------
# The summary is bounded to the top N and carries the causing command.
# ---------------------------------------------------------------------------

def _busy_trace():
    # approval wait (4600s) + repeated reads (4 of big.out) + an idle gap after
    # a job start at ts 6000/10000 — three DISTINCT frictions. The idle gap sits
    # outside the approval-wait window [100,4700] so it is not double-counted.
    idle = [_other(6000.0, "bash_background", '{"command": "submit_job"}'),
            _other(10000.0, "bash_status", '{"job_id": 7}')]
    return _approval_wait_trace() + _repeated_reads_trace() + idle


def test_summary_is_bounded_to_top_three():
    items = sr.friction_summary(_busy_trace(), top=3)
    assert len(items) == 3
    for i in items:
        assert "kind" in i and "magnitude" in i and "command" in i
