"""Tests for delfin.agent.report_costs — fabricated metrics only."""

from __future__ import annotations

import json
import time

from delfin.agent import report_costs


def _write_log(tmp_path, records):
    p = tmp_path / "agent_metrics.jsonl"
    p.write_text(
        "".join(json.dumps(r) + "\n" for r in records), encoding="utf-8")
    return p


def _rec(model, ts, cost=0.0, in_tok=0, out_tok=0, deleg=0.0, **extra):
    r = {
        "model": model, "ts": ts, "cost_usd": cost,
        "input_tokens": in_tok, "output_tokens": out_tok,
        "delegated_cost_usd": deleg,
    }
    r.update(extra)
    return r


def test_collect_empty_log_never_raises(tmp_path):
    p = tmp_path / "missing.jsonl"
    data = report_costs.collect(path=p)
    assert data["models"] == {}
    assert "No agent turn records" in report_costs.format_text(data)


def test_collect_per_model_sums(tmp_path):
    now = time.time()
    p = _write_log(tmp_path, [
        _rec("glm-5.3", now - 100, cost=0.10, in_tok=1000, out_tok=200),
        _rec("glm-5.3", now - 50, cost=0.30, in_tok=2000, out_tok=400,
             deleg=0.05, delegate_count=1),
        _rec("gpt-5.4", now - 10, cost=1.00, in_tok=500, out_tok=900),
    ])
    data = report_costs.collect(path=p)
    m = data["models"]
    assert m["glm-5.3"]["n_turns"] == 2
    assert m["glm-5.3"]["input_tokens"] == 3000
    assert m["glm-5.3"]["output_tokens"] == 600
    assert abs(m["glm-5.3"]["total_cost_usd"] - 0.40) < 1e-9
    assert abs(m["glm-5.3"]["total_delegated_cost_usd"] - 0.05) < 1e-9
    assert abs(m["glm-5.3"]["avg_cost_usd"] - 0.20) < 1e-9
    assert m["gpt-5.4"]["n_turns"] == 1
    assert abs(m["gpt-5.4"]["total_cost_usd"] - 1.00) < 1e-9


def test_collect_respects_window(tmp_path):
    now = time.time()
    old = now - 30 * 86_400
    p = _write_log(tmp_path, [
        _rec("glm-5.3", old, cost=5.0),
        _rec("glm-5.3", now - 60, cost=0.20),
    ])
    data = report_costs.collect(path=p)
    assert data["models"]["glm-5.3"]["n_turns"] == 1
    assert abs(data["models"]["glm-5.3"]["total_cost_usd"] - 0.20) < 1e-9


def test_collect_handles_corrupt_lines_and_missing_fields(tmp_path):
    now = time.time()
    p = tmp_path / "agent_metrics.jsonl"
    p.write_text(
        "not json\n"
        + json.dumps({"model": "m1", "ts": now}) + "\n"
        + json.dumps({"model": "m1", "ts": now, "cost_usd": "x",
                      "input_tokens": None}) + "\n",
        encoding="utf-8",
    )
    data = report_costs.collect(path=p)
    # The corrupt line is skipped by read_turns; the record with
    # cost_usd="x" is DROPPED by report_costs' sanitizer (it would
    # raise ValueError inside aggregate_by_model). Nothing raises.
    assert data["models"]["m1"]["n_turns"] == 1
    assert data["models"]["m1"]["total_cost_usd"] == 0.0


def test_window_comparison_included_when_enough_turns(tmp_path, monkeypatch):
    now = time.time()
    # 6 old turns (>5) in the week before yesterday, 6 new turns since.
    recs = []
    for i in range(6):
        recs.append(_rec("glm-5.3", now - 2 * 86_400 - 100 + i, cost=0.50))
    for i in range(6):
        recs.append(_rec("glm-5.3", now - 100 + i, cost=0.10))
    p = _write_log(tmp_path, recs)
    data = report_costs.collect(path=p)
    cmp_data = data["window_comparison"].get("glm-5.3")
    assert cmp_data is not None
    assert cmp_data["avg_cost_usd"]["delta"] < 0
    assert cmp_data["avg_cost_usd"]["improved"] is True

    text = report_costs.format_text(data)
    assert "avg_cost_usd" in text
    assert "(improved)" in text
    assert "glm-5.3" in text
    assert "$3.60" in text  # all 12 turns are inside the 7-day window


def test_format_text_notes_missing_comparison(tmp_path):
    now = time.time()
    p = _write_log(tmp_path, [_rec("m1", now, cost=0.1)])
    data = report_costs.collect(path=p)
    text = report_costs.format_text(data)
    assert "no recent-vs-earlier comparison" in text
    assert "Total direct cost: $0.10" in text


def test_main_smoke(tmp_path, capsys, monkeypatch):
    # Point the module's default log path at an empty temp file so the
    # CLI entrypoint runs without touching the user's real metrics.
    empty = tmp_path / "empty.jsonl"
    empty.write_text("", encoding="utf-8")

    class FakePath:
        def exists(self):
            return True

        def read_text(self, encoding="utf-8"):
            return ""

    import delfin.agent.agent_metrics as am
    monkeypatch.setattr(am, "_LOG_PATH", empty)
    rc = report_costs.main()
    out = capsys.readouterr().out
    assert rc == 0
    assert "No agent turn records" in out
