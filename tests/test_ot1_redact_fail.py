"""OT1 fix control (F1, QS-ruled shape): a failing redact path must never
export raw text; the field content is DROPPED and the failure is visible.

Red against 82363ca3: monkeypatches output_guard.scrub_secrets to raise.
On the previous commit spans() silently returned the raw text. Required
here: the field reads "[redaction unavailable]" and the notes (spans(...,
notes=True)) name the dropped fields. s4's adversarial test
(test_ot1_adversarial_otel.py::test_redaction_failure_is_visible_not_silent)
asserts the same rule from its side.

Also controlled: normal export is unchanged — a clean session still
exports its fields (nothing new refused).
"""

import json

import pytest

from delfin.agent import tool_trace, turn_metrics


def test_redact_failure_drops_content_and_notes(tmp_path, monkeypatch):
    monkeypatch.setattr(tool_trace, "_DIR", tmp_path / "tool_traces")
    monkeypatch.setattr(turn_metrics, "_DIR", tmp_path / "turn_metrics")
    tool_trace.record("ot1f1", tool="bash", tool_input={"command": "ls"},
                      output="hello world", duration_ms=5, ok=True)
    import delfin.agent.output_guard as og

    def boom(text, *a, **kw):
        raise RuntimeError("scrub unavailable")

    monkeypatch.setattr(og, "scrub_secrets", boom)
    from delfin.agent import otel_export
    spans, notes = otel_export.spans("ot1f1", notes=True)
    blob = json.dumps(spans)
    assert "hello world" not in blob, "raw text exported on redact failure"
    tools = [s for s in spans if s["name"] == "tool_call"]
    assert tools[0]["attributes"]["tool.output"] == "[redaction unavailable]"
    assert any("DROPPED" in n for n in notes), notes


def test_normal_export_unchanged(tmp_path, monkeypatch):
    monkeypatch.setattr(tool_trace, "_DIR", tmp_path / "tool_traces")
    monkeypatch.setattr(turn_metrics, "_DIR", tmp_path / "turn_metrics")
    tool_trace.record("ot1f2", tool="read_file", tool_input={"path": "a"},
                      output="plain content", duration_ms=5, ok=True)
    from delfin.agent import otel_export
    spans, notes = otel_export.spans("ot1f2", notes=True)
    tools = [s for s in spans if s["name"] == "tool_call"]
    assert tools[0]["attributes"]["tool.output"] == "plain content"
    assert not any("DROPPED" in n for n in notes)
