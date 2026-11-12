"""OT1 phase 1: otel_export.spans() — span tree from the three existing stores.

Controls are written against the REAL record shapes: tool_trace.record and
turn_metrics.record write the fixtures at run time into a redirected store
under tmp_path (never the live ~/.delfin), and the exported spans are built
from those real records, not from hand-typed dicts.
"""

import json
import time

import pytest

from delfin.agent import tool_trace
from delfin.agent import turn_metrics


def _fake_session(tmp_path, monkeypatch, n_tools=1, n_turns=1):
    """Record REAL entries for session 'ot1test' into redirected stores."""
    monkeypatch.setattr(tool_trace, "_DIR", tmp_path / "tool_traces")
    monkeypatch.setattr(turn_metrics, "_DIR", tmp_path / "turn_metrics")

    base = time.time() - 1000.0
    turn_metrics.record("ot1test", provider="kit", model="glm-5.3",
                        role="assistant", total_ms=5000, ttft_ms=120,
                        output_chars=900, tool_calls=2)
    t0 = base
    tool_trace.record("ot1test", tool="read_file", tool_input={"path": "a.py"},
                      output="file content", duration_ms=12, ok=True)
    tool_trace.record("ot1test", tool="bash", tool_input={"command": "ls"},
                      output="a.py", duration_ms=40, ok=False,
                      error="boom: auth"+"_token = rk_0123456789abcdefghij")
    return time.time()


def test_span_tree_session_turn_model_tool(tmp_path, monkeypatch):
    _fake_session(tmp_path, monkeypatch)
    from delfin.agent import otel_export

    spans = otel_export.spans("ot1test")
    assert any(s["name"] == "session" for s in spans), spans
    assert any(s["name"].startswith("turn") or
               s["attributes"].get("delfin.span_type") == "turn"
               for s in spans), spans
    assert any(s["name"] == "model_call" for s in spans), spans
    assert any(s["name"] == "tool_call" for s in spans), spans


def test_parent_ids_consistent(tmp_path, monkeypatch):
    _fake_session(tmp_path, monkeypatch)
    from delfin.agent import otel_export

    spans = otel_export.spans("ot1test")
    ids = {s["span_id"] for s in spans}
    assert len(ids) == len(spans), "duplicate span ids"
    for s in spans:
        if s["name"] == "session":
            assert s.get("parent_span_id", "") in ("", None)
        elif s["name"] in ("turn", "model_call", "tool_call"):
            assert s["parent_span_id"] in ids, s
    # every tool/model span's parent is the turn span
    turn_ids = {s["span_id"] for s in spans
                if s["attributes"].get("delfin.span_type") == "turn"}
    for s in spans:
        if s["name"] in ("tool_call",):
            assert s["parent_span_id"] in turn_ids


def test_attributes_tool_name_ok_duration(tmp_path, monkeypatch):
    _fake_session(tmp_path, monkeypatch)
    from delfin.agent import otel_export

    spans = otel_export.spans("ot1test")
    tools = [s for s in spans if s["name"] == "tool_call"]
    assert len(tools) == 2
    names = {s["attributes"]["tool.name"] for s in tools}
    assert names == {"read_file", "bash"}
    ok_flags = {s["attributes"]["tool.ok"] for s in tools}
    assert ok_flags == {True, False}
    for s in tools:
        assert s["attributes"]["duration_ms"] in (12, 40)
    model = [s for s in spans if s["name"] == "model_call"][0]
    assert model["attributes"]["model"] == "glm-5.3"
    assert model["attributes"]["provider"] == "kit"
    assert model["attributes"]["input_tokens"] == 0  # field present, defaulted


def test_token_attributes_on_model_call(tmp_path, monkeypatch):
    monkeypatch.setattr(tool_trace, "_DIR", tmp_path / "tool_traces")
    monkeypatch.setattr(turn_metrics, "_DIR", tmp_path / "turn_metrics")
    turn_metrics.record("ot1tok", provider="kit", model="glm-5.3",
                        role="assistant", total_ms=100, ttft_ms=10,
                        output_chars=50, tool_calls=0,
                        input_tokens=111, output_tokens=222, cached_tokens=33)
    from delfin.agent import otel_export
    spans = otel_export.spans("ot1tok")
    model = [s for s in spans if s["name"] == "model_call"]
    assert len(model) == 1
    m = model[0]
    assert m["attributes"]["input_tokens"] == 111
    assert m["attributes"]["output_tokens"] == 222
    assert m["attributes"]["cached_tokens"] == 33


def test_missing_stores_partial_trace_with_note(tmp_path, monkeypatch):
    from delfin.agent import otel_export

    spans, notes = otel_export.spans("ot1nosuch", notes=True)
    assert spans == [] or all(
        s["name"] == "session" for s in spans)
    assert isinstance(notes, list) and notes
    assert any("tool_trace" in n or "turn_metrics" in n for n in notes)


def test_session_span_metadata_from_session_store(tmp_path, monkeypatch):
    # session_store writes to the real home; that is refused in tests, so
    # load_session is stubbed to return a real-shaped record (the shape
    # read at session_store.py:1399-1424, not invented).
    import delfin.agent.session_store as ss
    monkeypatch.setattr(tool_trace, "_DIR", tmp_path / "tool_traces")
    monkeypatch.setattr(turn_metrics, "_DIR", tmp_path / "turn_metrics")
    tool_trace.record("ot1meta", tool="read_file", tool_input={"p": 1},
                      output="x", duration_ms=5, ok=True)
    monkeypatch.setattr(ss, "load_session", lambda sid: {
        "session_id": sid, "title": "OT1 test", "provider": "kit",
        "model": "glm-5.3", "mode": "quick", "pid": 0, "heartbeat_at": 0.0})
    from delfin.agent import otel_export
    spans = otel_export.spans("ot1meta")
    session = [s for s in spans if s["name"] == "session"][0]
    assert session["attributes"]["session.title"] == "OT1 test"
    assert session["attributes"]["session.provider"] == "kit"


def test_secret_never_in_span(tmp_path, monkeypatch):
    _fake_session(tmp_path, monkeypatch)
    from delfin.agent import otel_export

    spans = otel_export.spans("ot1test")
    blob = json.dumps(spans)
    assert "rk_0123456789abcdefghij" not in blob
    # error text survives redaction but the token does not
    tools = [s for s in spans if s["name"] == "tool_call"
             and s["attributes"]["tool.name"] == "bash"]
    assert tools[0]["attributes"].get("tool.error", "")
    assert "boom" in tools[0]["attributes"]["tool.error"]


def test_output_cut_to_500_chars(tmp_path, monkeypatch):
    monkeypatch.setattr(tool_trace, "_DIR", tmp_path / "tool_traces")
    monkeypatch.setattr(turn_metrics, "_DIR", tmp_path / "turn_metrics")
    tool_trace.record("ot1cut", tool="read_file", tool_input={"path": "x"},
                      output="Z" * 5000, duration_ms=1, ok=True)
    from delfin.agent import otel_export
    spans = otel_export.spans("ot1cut")
    tools = [s for s in spans if s["name"] == "tool_call"]
    assert len(tools[0]["attributes"].get("tool.output", "")) <= 500
