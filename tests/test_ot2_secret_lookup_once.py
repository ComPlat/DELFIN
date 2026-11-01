"""OT2: the OTLP export scrubs secrets ONCE per export, not once per span.

The defect (operator, 01.11., PR #507 blocked): ``otel_export._redact`` calls
``output_guard.scrub_secrets`` for every span field, and ``scrub_secrets``
calls ``web_tools._known_secret_values()`` on every call. If a predecessor
test makes that lookup slow, the export is 10 000 x slow lookups.

Red test: wrap ``_known_secret_values`` in a counting shim and show that one
``write_otlp_json`` over 100 tool spans triggers many lookups (3 per span on
main). After the fix: exactly one lookup per export call.

Note: the secret-shaped literals below are assembled by concatenation so
they survive DELFIN's own redaction of tool-call text; they are the classic
public AWS-docs example and a dummy string, never a real credential.
"""

import json

from delfin.agent import otel_export
from delfin.agent import output_guard
from delfin.agent import tool_trace
from delfin.agent import turn_metrics

# shape-test values (public dummy examples, not real credentials)
_AWS = "AKI" + "AIOSFODNN7EXAMPLE"          # aws-access-key shape
_SK = "sk-" + "proj-abcdefghijklmnop"       # provider-api-key shape
_HELD = "kit-secret-value-abc123456789"     # fake "held" exact secret


def _counting_secrets(monkeypatch, delay_s=0.0002):
    """Replace the real lookup with a counted, slightly slow stand-in."""
    calls = {"n": 0}

    def fake_values():
        calls["n"] += 1
        import time as _t
        _t.sleep(delay_s)
        return [_HELD]

    monkeypatch.setattr("delfin.agent.web_tools._known_secret_values",
                        fake_values)
    return calls


def _write_store(tmp_path, monkeypatch, name, n_spans):
    monkeypatch.setattr(tool_trace, "_DIR", tmp_path / "tool_traces")
    monkeypatch.setattr(turn_metrics, "_DIR", tmp_path / "turn_metrics")
    store = tmp_path / "tool_traces"
    store.mkdir()
    with open(store / f"{name}.jsonl", "w") as f:
        for i in range(n_spans):
            f.write(json.dumps({
                "ts": 100.0 + i * 0.001, "tool": f"t{i % 5}",
                "input": json.dumps({"path": f"/x/{i}"}),
                "output": f"result {i}",
                "duration_ms": 5, "ok": True, "error": ""}) + "\n")
    turn_metrics.record(name, provider="kit", model="m", role="assistant",
                        total_ms=50_000, ttft_ms=10, output_chars=1,
                        tool_calls=n_spans)


def test_export_looks_up_held_secrets_once(tmp_path, monkeypatch):
    calls = _counting_secrets(monkeypatch)
    _write_store(tmp_path, monkeypatch, "ot2once", 100)
    out = tmp_path / "trace.json"
    otel_export.write_otlp_json("ot2once", out)
    # One lookup per export, not one per span field (100 spans = up to 300).
    assert calls["n"] <= 1, (
        f"scrub_secrets fetched the held secrets {calls['n']} times for one "
        f"export over 100 spans — must be computed once per export")


def test_export_redaction_result_is_identical(tmp_path, monkeypatch):
    """Same redaction result with the once-per-export path as with today's
    per-call lookup: held secret, shaped secrets and plain text."""
    monkeypatch.setattr("delfin.agent.web_tools._known_secret_values",
                        lambda: [_HELD])
    _write_store(tmp_path, monkeypatch, "ot2pin", 3)
    # overwrite the store with content hitting all three classes
    store = tmp_path / "tool_traces"
    rows = [
        {"ts": 100.0, "tool": "t0", "input": json.dumps({"k": _HELD}),
         "output": f"leak {_HELD}", "duration_ms": 5, "ok": True,
         "error": ""},
        {"ts": 100.1, "tool": "t1", "input": _SK,
         "output": _AWS, "duration_ms": 5, "ok": True, "error": ""},
        {"ts": 100.2, "tool": "t2", "input": "plain nothing here",
         "output": "benzene CC(=O)Oc1ccccc1C", "duration_ms": 5,
         "ok": True, "error": ""},
    ]
    with open(store / "ot2pin.jsonl", "w") as f:
        for r in rows:
            f.write(json.dumps(r) + "\n")
    out = tmp_path / "trace.json"
    otel_export.write_otlp_json("ot2pin", out)
    doc = json.loads(out.read_text())
    # OTLP attributes are [{key, value}, ...] lists (_otlp_attr); flatten
    # them back to a dict per span, keyed by spanId.
    tool_spans = {}
    for r in doc["resourceSpans"]:
        for sc in r["scopeSpans"]:
            for s in sc["spans"]:
                attrs = {a["key"]: a["value"].get("stringValue", "")
                         for a in s["attributes"]}
                if attrs.get("delfin.span_type") == "tool_call":
                    tool_spans[s["spanId"]] = attrs
    assert len(tool_spans) == 3
    texts = []
    for i in range(3):
        a = tool_spans[otel_export._id("ot2pin", f"tool{i}")]
        texts += [a.get("tool.input", ""), a.get("tool.output", "")]
    blob = "\n".join(texts)
    # exact held secret redacted wherever it appears
    assert _HELD not in blob and "[redacted:credential]" in blob
    # shape checks still applied
    assert _AWS not in blob
    assert _SK not in blob
    # plain chemistry / ordinary text untouched
    assert "benzene CC(=O)Oc1ccccc1C" in blob
    assert "plain nothing here" in blob


def test_lookup_failure_degrades_like_today(tmp_path, monkeypatch):
    """A raising lookup must not change the redaction outcome: today the
    exception is swallowed inside scrub_secrets and the shape scan still
    runs. The once-per-export path must keep exactly that behaviour."""
    def boom():
        raise RuntimeError("keyring unavailable")

    monkeypatch.setattr("delfin.agent.web_tools._known_secret_values", boom)
    _write_store(tmp_path, monkeypatch, "ot2fail", 2)
    # shape-only redaction must still happen (exact pass degraded, as today)
    assert output_guard.scrub_secrets(_AWS) != _AWS
    out = tmp_path / "trace.json"
    otel_export.write_otlp_json("ot2fail", out)
    blob = out.read_text()
    assert len(json.loads(blob)["resourceSpans"]) == 1
    assert '"tool_call"' in blob
    assert _AWS not in blob
    assert "[redaction unavailable]" not in blob
