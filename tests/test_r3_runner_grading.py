"""The graded long-horizon run carries the derailment counts (package R3).

The runner is the real path a live run grades through (runner.run_task),
not a helper called directly: a threshold the yaml names must reject a
transcript that derailed and accept one that did not, through the same
grading the CLI performs.
"""
from __future__ import annotations

import pytest

from delfin.agent.longhorizon import runner


def _transcript(n_failing_repeats: int):
    """A tool record in tool_trace.record's accepted field set: three
    identical failing calls for n_failing_repeats=2, plus a clean tail."""
    events = [
        {"kind": "tool", "tool": "bash", "input": "xtb s.xyz",
         "output": "command not found", "ok": False}
        for _ in range(n_failing_repeats + 1)
    ]
    events.append({"kind": "tool", "tool": "bash",
                   "input": "echo done", "output": "done", "ok": True})
    events.append({"kind": "turn", "stopped": True,
                   "error": "", "tool_calls": 3,
                   "input_tokens": 100, "output_tokens": 20})
    return events


def test_repeated_failing_calls_bound_rejects_a_derailed_run():
    # Fails if the threshold did not bind: 3 identical failing calls must
    # score 2 and exceed a max of 1.
    payload = runner.run_task(
        "lh-code-change-fact", _transcript(2), thresholds={
            "repeated_failing_calls_max": 1})
    assert payload["accept"] is False
    assert "repeated_failing_calls=2" in payload["verdict"]


def test_zero_bound_accepts_a_clean_run():
    # The control: the same grading rule accepts a transcript with no
    # repeats — fails if the check always rejected.
    payload = runner.run_task(
        "lh-code-change-fact", _transcript(0), thresholds={
            "repeated_failing_calls_max": 0})
    assert payload["accept"] is True
    assert payload["metrics"]["repeated_failing_calls"] == 0


def test_no_threshold_entry_keeps_the_old_verdict():
    # A task yaml that names no derailment bound grades byte-for-byte as
    # before — fails if the new check made a bound implicit.
    clean = _transcript(0)
    derailed = _transcript(2)
    for t in (clean, derailed):
        payload = runner.run_task("lh-code-change-fact", t, thresholds={})
        assert payload["accept"] is True, payload["verdict"]


def test_report_carries_the_new_lines():
    import tempfile
    from pathlib import Path
    with tempfile.TemporaryDirectory() as td:
        payload = runner.run_task(
            "lh-code-change-fact", _transcript(2), thresholds={},
            workspace=Path(td))
        text = (Path(td) / "report.md").read_text("utf-8")
        assert "repeated_failing_calls: 2" in text
        assert "unproductive_loops: 0" in text
