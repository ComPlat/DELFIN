"""Tests for the tool-result truncation that protects the context window
from MCP outputs that exceed a few kB.

The truncation uses ``api_client._smart_truncate`` which keeps a HEAD
slice + a TAIL slice + a "chars omitted" marker so Python tracebacks
(usually at the end of the output) always survive.
"""

from __future__ import annotations

from delfin.agent.api_client import _smart_truncate


CAP = 5000
LABEL = "tool_result"


def test_short_output_passes_through_unchanged():
    s = "a" * 100
    assert _smart_truncate(s, cap=CAP, label=LABEL) == s


def test_exact_cap_passes_through_unchanged():
    s = "a" * CAP
    assert _smart_truncate(s, cap=CAP, label=LABEL) == s


def test_oversized_output_gets_truncated():
    s = "a" * (CAP * 3)
    out = _smart_truncate(s, cap=CAP, label=LABEL)
    assert len(out) < len(s)
    assert "truncated" in out
    assert "head and tail preserved" in out


def test_truncation_preserves_traceback_tail():
    """The most useful part of an error is the last few lines.  Make
    sure they survive the truncation."""
    body = "noise " * 5000
    traceback_tail = (
        "Traceback (most recent call last):\n"
        '  File "x.py", line 1, in <module>\n'
        "    raise RuntimeError(\"the actual error\")\n"
        "RuntimeError: the actual error"
    )
    s = body + traceback_tail
    out = _smart_truncate(s, cap=CAP, label=LABEL)
    assert "the actual error" in out
    assert "RuntimeError" in out


def test_empty_input_returns_empty():
    assert _smart_truncate("", cap=CAP, label=LABEL) == ""


def test_marker_records_actual_omitted_count():
    s = "x" * 10_000
    out = _smart_truncate(s, cap=CAP, label=LABEL)
    # Exact omitted = original_len - (head_len + tail_len)
    omitted_in_marker = next(
        int(p.split()[0])
        for p in out.split("(tool_result truncated, ")
        if p[:1].isdigit()
    )
    assert omitted_in_marker > 0
    assert omitted_in_marker < len(s)


def test_label_appears_in_marker():
    out = _smart_truncate("z" * 9000, cap=CAP, label="cargo-test")
    assert "(cargo-test truncated" in out
