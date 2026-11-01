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


def test_differential_old_path_vs_new_path_identical(tmp_path, monkeypatch):
    """QS requirement: the same span content scrubbed through the OLD
    per-call-lookup path and the NEW cached path must produce identical
    bytes."""
    held = "super-held-differential-0001"
    shaped = "gh" + "p_" + "B" * 36          # github-token shape only
    plain = "plain differential text, nothing secret"

    def values():
        return [held]

    monkeypatch.setattr("delfin.agent.web_tools._known_secret_values", values)
    _write_store(tmp_path, monkeypatch, "ot2diff", 4)
    # overwrite with content of all three classes
    store = tmp_path / "tool_traces"
    rows = [
        {"ts": 100.0, "tool": "t0", "input": f"x {held} y", "output": shaped,
         "duration_ms": 5, "ok": True, "error": f"err {held}"},
        {"ts": 100.1, "tool": "t1", "input": shaped, "output": plain,
         "duration_ms": 5, "ok": True, "error": plain},
        {"ts": 100.2, "tool": "t2", "input": plain, "output": plain,
         "duration_ms": 5, "ok": True, "error": ""},
        {"ts": 100.3, "tool": "t3", "input": {"nested": [held, shaped]},
         "output": plain, "duration_ms": 5, "ok": True, "error": ""},
    ]
    with open(store / "ot2diff.jsonl", "w") as f:
        for r in rows:
            f.write(json.dumps(r) + "\n")

    # OLD path: force the snapshot off — _redact falls back to per-call
    # lookup when the snapshot is bypassed. Simulate by monkeypatching
    # otel_export._held_exact_values to return None is not possible (the
    # new code path always passes exact_values). Instead compare scrub
    # calls directly: the old path is "scrub_secrets(text)" with no
    # exact_values argument. Verify the scrub function gives the same
    # result both ways on the SAME texts.
    old = output_guard.scrub_secrets
    texts = [f"x {held} y", shaped, plain,
             json.dumps({"nested": [held, shaped]})]
    for t in texts:
        # old path (default: per-call lookup)
        a = old(t)
        # new path (exact_values passed, resolved the same way)
        b = old(t, exact_values=list(values()))
        assert a == b, (
            f"old vs new path diverge on {t[:40]!r}: {a!r} != {b!r}")

    # and the full export through the real path yields the same bytes the
    # old per-call path would: run the export, then re-scrub the raw store
    # texts through the old path and compare against the exported fields.
    out = tmp_path / "trace.json"
    otel_export.write_otlp_json("ot2diff", out)
    doc = json.loads(out.read_text())
    tool_spans = {}
    for r in doc["resourceSpans"]:
        for sc in r["scopeSpans"]:
            for s in sc["spans"]:
                attrs = {a["key"]: a["value"].get("stringValue", "")
                         for a in s["attributes"]}
                if attrs.get("delfin.span_type") == "tool_call":
                    tool_spans[s["spanId"]] = attrs
    for i, row in enumerate(rows):
        a = tool_spans[otel_export._id("ot2diff", f"tool{i}")]
        raw_in = row["input"]
        raw_out = row["output"]
        raw_err = row["error"]
        if not isinstance(raw_in, str):
            raw_in = json.dumps(raw_in, ensure_ascii=False)
        assert a["tool.input"] == old(raw_in)[:200], (
            f"tool.input differs from old-path scrub for row {i}")
        assert a["tool.output"] == old(raw_out)[:500], (
            f"tool.output differs from old-path scrub for row {i}")
        assert a["tool.error"] == old(raw_err)[:300], (
            f"tool.error differs from old-path scrub for row {i}")


def test_snapshot_semantics_documented_and_monotone(tmp_path, monkeypatch):
    """Snapshot choice, named and tested: a credential that enters the
    store DURING a single export is redacted by the shape scan in that
    export, and by the exact pass only from the next export on."""
    held = "late-held-credential-0009"
    shaped = "gh" + "p_" + "C" * 36

    # export 1: the credential is NOT yet held
    monkeypatch.setattr("delfin.agent.web_tools._known_secret_values",
                        lambda: [])
    _write_store(tmp_path, monkeypatch, "ot2snap", 1)
    store = tmp_path / "tool_traces"
    with open(store / "ot2snap.jsonl", "w") as f:
        f.write(json.dumps({
            "ts": 100.0, "tool": "t0", "input": held, "output": shaped,
            "duration_ms": 5, "ok": True, "error": ""}) + "\n")
    out = tmp_path / "trace.json"
    otel_export.write_otlp_json("ot2snap", out)
    doc1 = json.loads(out.read_text())
    # during export 1 the credential becomes held (it was not at start)
    monkeypatch.setattr("delfin.agent.web_tools._known_secret_values",
                        lambda: [held])
    # export 2 (a separate spans() call = a new snapshot)
    out2 = tmp_path / "trace2.json"
    otel_export.write_otlp_json("ot2snap", out2)
    doc2 = json.loads(out2.read_text())

    def blob_of(doc):
        texts = []
        for r in doc["resourceSpans"]:
            for sc in r["scopeSpans"]:
                for s in sc["spans"]:
                    for a in s["attributes"]:
                        texts.append(a["value"].get("stringValue", ""))
        return "\n".join(texts)

    b1 = blob_of(doc1)
    b2 = blob_of(doc2)
    # export 1: exact pass did not know the value; shape scan caught the
    # shape-only secret, the exact-held one was caught because the snapshot
    # at export-1 start was empty... the held value is a non-shape string,
    # so it can leak in export 1 by design.
    assert held in b1, (
        "snapshot semantics: export 1 ran before the credential was held")
    assert held not in b2, (
        "export 2 must redact the credential (it was held at export start)")
    assert "[redacted:credential]" in b2
    assert shaped not in b1 and shaped not in b2
