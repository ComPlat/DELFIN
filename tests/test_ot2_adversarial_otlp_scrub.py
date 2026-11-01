"""OT2 adversarial tests (reviewer nacht-s2): how often does the export
resolve the held-credential values?

The defect (PR #507 CI): ``otel_export`` scrubs every span field through
``output_guard.scrub_secrets``, which calls ``web_tools._known_secret_values()``
on EVERY call (file + keyring + env lookups). 10 000 spans -> ~10 000-30 000
lookups per export, so the 5.0 s perf test is hostage to how slow a
predecessor test made that lookup.

Contract under test (the fix must satisfy): one export calls the held-secret
lookup at most once, and the redaction result is byte-identical to today.

Records are written through the REAL writers (tool_trace.record) into
redirected stores; the export runs through the real public path
``otel_export.write_otlp_json``. ``_known_secret_values`` is monkeypatched as
a black box (never touched the credentials module), with a counting wrapper.
"""

from __future__ import annotations

import json

from delfin.agent import otel_export
from delfin.agent import tool_trace
from delfin.agent import turn_metrics
from delfin.agent import web_tools


def _record(tmp_path, monkeypatch, name, n_tools):
    monkeypatch.setattr(tool_trace, "_DIR", tmp_path / "tool_traces")
    monkeypatch.setattr(turn_metrics, "_DIR", tmp_path / "turn_metrics")
    turn_metrics.record(name, provider="kit", model="glm",
                        input_tokens=1, output_tokens=2, total_ms=10.0)
    for i in range(n_tools):
        tool_trace.record(name, tool=f"tool{i}",
                          tool_input={"i": i},
                          output=f"result number {i} with text",
                          ok=True)
    return name


def _count_lookups(monkeypatch):
    """Replace web_tools._known_secret_values with a counting wrapper."""
    calls = {"n": 0}

    def fake():
        calls["n"] += 1
        return ["super-secret-held-value-0001"]

    monkeypatch.setattr(web_tools, "_known_secret_values", fake,
                        raising=False)
    return calls


def test_export_resolves_held_secrets_at_most_once(tmp_path, monkeypatch):
    """10 tool spans -> 30 scrub calls on main -> 30 lookups. The export
    must resolve the held values once, not once per field."""
    _record(tmp_path, monkeypatch, "ot2adv", n_tools=10)
    calls = _count_lookups(monkeypatch)
    otel_export.write_otlp_json("ot2adv", tmp_path / "trace.json")
    assert calls["n"] <= 1, (
        f"held-secret lookup ran {calls['n']} times for one export of "
        "10 tool spans; it must be resolved once per export")


def test_export_redaction_bytes_unchanged(tmp_path, monkeypatch):
    """Byte-for-byte: a held secret, a shape-only secret and plain text
    come out the same as today through the real export path."""
    held = "super-secret-held-value-0001"
    shaped = "ghp_" + "A" * 36          # shape scan only, not held
    plain = "no secrets here at all, just a result string"
    monkeypatch.setattr(tool_trace, "_DIR", tmp_path / "tool_traces")
    monkeypatch.setattr(turn_metrics, "_DIR", tmp_path / "turn_metrics")
    turn_metrics.record("ot2bytes", provider="kit", model="glm",
                        input_tokens=1, output_tokens=2, total_ms=10.0)
    tool_trace.record("ot2bytes", tool="t0", tool_input={"text": held},
                      output=shaped, ok=True)
    tool_trace.record("ot2bytes", tool="t1", tool_input={"text": plain},
                      output=plain, ok=True)
    calls = _count_lookups(monkeypatch)
    p = otel_export.write_otlp_json("ot2bytes", tmp_path / "trace.json")
    doc = json.loads((tmp_path / "trace.json").read_text(encoding="utf-8"))
    spans = doc["resourceSpans"][0]["scopeSpans"][0]["spans"]
    tool_spans = [s for s in spans if s["name"] == "tool_call"]
    texts = []
    for s in tool_spans:
        for a in s["attributes"]:
            v = a["value"]
            texts.append(v.get("stringValue", ""))
    blob = "\n".join(texts)
    assert held not in blob, "held credential leaked into the export"
    assert "[redacted:credential]" in blob, "held secret not replaced"
    assert shaped not in blob, "shaped secret leaked into the export"
    assert "[redacted:github-token]" in blob, "shape scan not applied"
    assert plain in blob, "plain text was wrongly redacted"
    # the fix must also keep the lookup count small for this export
    assert calls["n"] <= 1, (
        f"lookup ran {calls['n']} times on a 2-span export")
    assert p == str(tmp_path / "trace.json")
