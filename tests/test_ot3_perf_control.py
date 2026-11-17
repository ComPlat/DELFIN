"""OT3 Part A: the perf tests must measure the EXPORTER, not the runner.

Observed CI condition (#507 shard 1, operator): after ~110 predecessor test
files the wall clock for the 10k-span export is 6.4-6.5 s while the exporter
alone takes 1.9 s. On a contended runner wall clock and process CPU time
diverge; wall time inflates although the exporter's work did not change.

That condition is injected here deterministically: ``time.monotonic`` is
scaled (factor 100) exactly like wall-vs-CPU divergence on a busy runner.
Under that injection:

* ``test_ot1_otlp.test_perf_10k_write_otlp`` must still PASS — it measures
  process CPU time, which the injected condition does not touch. RED on the
  previous commit, where it measured bare wall clock (AssertionError: 140 s).
* The exporter's own CPU cost stays far below the bound while the wall
  measurement breaches it — the measurement was the defect, not the exporter.

Control provenance (red on unchanged code, seen through the gate): probe
tests/test_zz_ot3_perf_red.py on ea7d62774 gave ".F" — the real perf test
failed under the injection with elapsed 140.4 s while the exporter's process
CPU baseline stayed < 1.5 s.
"""

import gc
import importlib.util
import json
import pathlib
import time


def _load_test_ot1_otlp():
    here = pathlib.Path(__file__).resolve().parent
    spec = importlib.util.spec_from_file_location(
        "test_ot1_otlp_ot3_control", here / "test_ot1_otlp.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def _scale_monotonic(monkeypatch, factor=100.0):
    real = time.monotonic

    def scaled():
        return real() * factor

    monkeypatch.setattr(time, "monotonic", scaled)


def _cpu_of_export(session, tmp_path):
    from delfin.agent import otel_export
    gc.collect()
    t0 = time.process_time()
    otel_export.write_otlp_json(session, tmp_path / "out.json")
    return time.process_time() - t0


def test_perf_test_survives_injected_wall_divergence(tmp_path, monkeypatch):
    """The real perf test passes under an injected contended-runner clock."""
    mod = _load_test_ot1_otlp()
    _scale_monotonic(monkeypatch)
    mod.test_perf_10k_write_otlp(tmp_path, monkeypatch)


def test_perf_sibling_survives_injected_wall_divergence(tmp_path, monkeypatch):
    """The sibling spans() perf test also survives the injected condition."""
    mod = _load_test_ot1_otlp()
    _scale_monotonic(monkeypatch)
    mod.test_perf_10k_tool_calls(tmp_path, monkeypatch)


def test_exporter_cpu_far_below_bound_but_wall_breaches_under_injection(
        tmp_path, monkeypatch):
    """The injection reproduces the CI condition on THIS measurement:
    wall clock (old style) breaches 5 s while the exporter's process CPU
    time stays a small fraction of that — the exporter is not the defect."""
    from delfin.agent import otel_export
    from delfin.agent import tool_trace
    from delfin.agent import turn_metrics
    monkeypatch.setattr(tool_trace, "_DIR", tmp_path / "tool_traces")
    monkeypatch.setattr(turn_metrics, "_DIR", tmp_path / "turn_metrics")
    store = tmp_path / "tool_traces"
    store.mkdir()
    base = time.time() - 100.0
    with open(store / "ot3ctl.jsonl", "w") as f:
        for i in range(10_000):
            f.write(json.dumps({
                "ts": base + i * 0.001, "tool": f"t{i % 7}",
                "input": "{}", "output": "x" * 40,
                "duration_ms": 5, "ok": True, "error": ""}) + "\n")
    turn_metrics.record("ot3ctl", provider="kit", model="m",
                        role="assistant", total_ms=600_000, ttft_ms=10,
                        output_chars=1, tool_calls=10_000)

    cpu = _cpu_of_export("ot3ctl", tmp_path)

    _scale_monotonic(monkeypatch)
    wall_start = time.monotonic()
    otel_export.write_otlp_json("ot3ctl", tmp_path / "wall.json")
    wall_old_style = time.monotonic() - wall_start
    monkeypatch.undo()

    # Old-style wall measurement breaches the 5.0 s bound under the
    # injected condition ...
    assert wall_old_style >= 5.0, wall_old_style
    # ... while the exporter's process CPU time stays far below it.
    assert cpu < 5.0, cpu
    assert wall_old_style > 2.0 * max(cpu, 1.0), (wall_old_style, cpu)
