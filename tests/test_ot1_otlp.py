"""OT1 phase 2: OTLP/JSON file export, text timeline, and the perf gate.

Real records again: fixtures write through tool_trace.record /
turn_metrics.record into redirected stores. The OTLP file is validated
structurally (resourceSpans -> scopeSpans -> spans), not against a schema
library: the point is that Jaeger/Tempo can import it, so field names and
nesting must match the OTLP/JSON shape.
"""

import json
import time

from delfin.agent import tool_trace
from delfin.agent import turn_metrics


def _session(tmp_path, monkeypatch, name="ot1otlp", n_tools=2):
    monkeypatch.setattr(tool_trace, "_DIR", tmp_path / "tool_traces")
    monkeypatch.setattr(turn_metrics, "_DIR", tmp_path / "turn_metrics")
    turn_metrics.record(name, provider="kit", model="glm-5.3",
                        role="assistant", total_ms=3000, ttft_ms=100,
                        output_chars=400, tool_calls=n_tools)
    for i in range(n_tools):
        tool_trace.record(name, tool=f"tool{i}", tool_input={"i": i},
                          output=f"out{i}", duration_ms=10 + i, ok=True)


def test_write_otlp_json_structure(tmp_path, monkeypatch):
    _session(tmp_path, monkeypatch)
    from delfin.agent import otel_export
    out = tmp_path / "trace.json"
    p = otel_export.write_otlp_json("ot1otlp", out)
    assert str(p) == str(out)
    doc = json.loads(out.read_text())
    rs = doc["resourceSpans"]
    assert isinstance(rs, list) and rs
    spans = []
    for r in rs:
        for ss in r.get("scopeSpans", []):
            spans.extend(ss.get("spans", []))
    names = {s["name"] for s in spans}
    assert {"session", "model_call", "tool_call"} <= names
    for s in spans:
        assert s["traceId"]
        assert s["spanId"]
        assert "startTimeUnixNano" in s and "endTimeUnixNano" in s
        assert "attributes" in s


def test_write_otlp_no_network(tmp_path, monkeypatch):
    _session(tmp_path, monkeypatch)
    from delfin.agent import otel_export
    out = tmp_path / "trace.json"
    otel_export.write_otlp_json("ot1otlp", out)
    assert out.exists() and out.stat().st_size > 0


def test_text_timeline_marks_longest(tmp_path, monkeypatch):
    _session(tmp_path, monkeypatch, n_tools=3)
    from delfin.agent import otel_export
    text = otel_export.timeline("ot1otlp")
    lines = text.splitlines()
    assert any("session" in ln for ln in lines)
    assert any("model_call" in ln for ln in lines)
    assert any("tool_call" in ln for ln in lines)
    # exactly three longest marked
    marks = [ln for ln in lines if ln.strip().startswith("*")]
    assert len(marks) == 3


def test_perf_10k_tool_calls(tmp_path, monkeypatch):
    import delfin.agent.otel_export as oe
    monkeypatch.setattr(tool_trace, "_DIR", tmp_path / "tool_traces")
    monkeypatch.setattr(turn_metrics, "_DIR", tmp_path / "turn_metrics")
    # write the store directly in the recorded format: the record() call
    # loop itself is the thing that would dominate, not the exporter.
    store = tmp_path / "tool_traces"
    store.mkdir()
    base = time.time() - 100.0
    with open(store / "ot1perf.jsonl", "w") as f:
        for i in range(10_000):
            f.write(json.dumps({
                "ts": base + i * 0.001, "tool": f"t{i % 7}",
                "input": "{}", "output": "x" * 40,
                "duration_ms": 5, "ok": True, "error": ""}) + "\n")
    turn_metrics.record("ot1perf", provider="kit", model="m",
                        role="assistant", total_ms=600_000, ttft_ms=10,
                        output_chars=1, tool_calls=10_000)
    from delfin.agent import otel_export
    t0 = time.monotonic()
    spans = otel_export.spans("ot1perf")
    elapsed = time.monotonic() - t0
    assert len(spans) >= 10_000
    assert elapsed < 5.0, elapsed


def test_perf_10k_write_otlp(tmp_path, monkeypatch):
    import delfin.agent.otel_export as oe
    monkeypatch.setattr(tool_trace, "_DIR", tmp_path / "tool_traces")
    monkeypatch.setattr(turn_metrics, "_DIR", tmp_path / "turn_metrics")
    store = tmp_path / "tool_traces"
    store.mkdir()
    base = time.time() - 100.0
    with open(store / "ot1perf2.jsonl", "w") as f:
        for i in range(10_000):
            f.write(json.dumps({
                "ts": base + i * 0.001, "tool": f"t{i % 7}",
                "input": "{}", "output": "x" * 40,
                "duration_ms": 5, "ok": True, "error": ""}) + "\n")
    turn_metrics.record("ot1perf2", provider="kit", model="m",
                        role="assistant", total_ms=600_000, ttft_ms=10,
                        output_chars=1, tool_calls=10_000)
    from delfin.agent import otel_export
    t0 = time.monotonic()
    otel_export.write_otlp_json("ot1perf2", tmp_path / "big.json")
    elapsed = time.monotonic() - t0
    assert elapsed < 5.0, elapsed
