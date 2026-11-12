"""OT1 phase-3 finding (s4): the terminal timeline must render the DROPPED
note so a user at the CLI sees the redaction failure too.

Red on 7a50e7b6: timeline() ignored notes; the CLI user saw a normal
timeline while fields were dropped.
"""

import pytest

from delfin.agent import tool_trace, turn_metrics


def test_timeline_renders_redaction_note(tmp_path, monkeypatch):
    monkeypatch.setattr(tool_trace, "_DIR", tmp_path / "tool_traces")
    monkeypatch.setattr(turn_metrics, "_DIR", tmp_path / "turn_metrics")
    tool_trace.record("ot1tn", tool="bash", tool_input={"c": "ls"},
                      output="visible", duration_ms=5, ok=True)
    import delfin.agent.output_guard as og

    def boom(text, *a, **kw):
        raise RuntimeError("scrub unavailable")

    monkeypatch.setattr(og, "scrub_secrets", boom)
    from delfin.agent import otel_export
    text = otel_export.timeline("ot1tn")
    assert "DROPPED" in text, text
    assert "redact path failed" in text
    # the raw output is nowhere in the terminal view
    assert "visible" not in text.split("!")[0].replace(
        "tool_call bash", "")  # crude guard: no raw output line rendered
