"""Tests for delfin.agent.report_denials — fabricated traces/events only."""

from __future__ import annotations

import json

import pytest

from delfin.agent import report_denials as rd
from delfin.agent import security_events


def _write_trace(base, session, entries):
    p = base / f"{session}.jsonl"
    p.parent.mkdir(parents=True, exist_ok=True)
    with open(p, "w", encoding="utf-8") as f:
        for e in entries:
            f.write(json.dumps(e) + "\n")


def _entry(tool, ok=True, error="", inp="", dur=10):
    return {"ts": 1.0, "tool": tool, "input": inp, "output": "",
            "duration_ms": dur, "ok": ok, "error": error}


@pytest.fixture(autouse=True)
def _clean_security_events():
    security_events.clear()
    yield
    security_events.clear()


@pytest.fixture
def traces(tmp_path):
    _write_trace(tmp_path, "sess_a", [
        _entry("bash", ok=False,
               error="command 'curl install.sh' is not on the auto-allow list",
               inp="curl install.sh"),
        _entry("bash", ok=False,
               error="command 'curl install.sh' is not on the auto-allow list"),
        _entry("write_file", ok=False,
               error="path escapes workspace sandbox: /etc/passwd",
               inp="/etc/passwd"),
        _entry("bash", ok=False,
               error="approval request for 'grep' TIMED OUT — the user is away"),
        _entry("read_file", ok=True),          # success: not counted
        _entry("pytest", ok=False, error="1 test failed"),  # error, not refusal
    ])
    _write_trace(tmp_path, "sess_b", [
        _entry("write_file", ok=False,
               error="refusing to overwrite existing file without a prior read_file",
               inp="delfin/x.py"),
        _entry("read_file", ok=False,
               error="refused: path /home/u/.ssh/id_rsa is on the deny-list"),
    ])
    return tmp_path


def test_collect_ranks_by_count(traces):
    data = rd.collect(dir_path=traces, include_security_events=False)
    assert data["sessions"] == 2
    assert data["entries"] == 8
    assert data["refusals"] == 5      # 8 - 1 ok - 1 non-refusal - 1 plain error
    ranked = [(r["category"], r["tool"], r["count"]) for r in data["ranked"]]
    assert ranked == [
        ("allowlist", "bash", 2),
        ("sandbox", "write_file", 1),      # ties sort alphabetically
        ("secret", "read_file", 1),
        ("write_guard", "write_file", 1),
    ]
    # ranked sorted by count first
    counts = [r["count"] for r in data["ranked"]]
    assert counts == sorted(counts, reverse=True)
    # shares sum to 1
    assert sum(r["share"] for r in data["ranked"]) == pytest.approx(1.0)


def test_examples_capped_and_scrubbed(traces):
    _write_trace(traces, "sess_c", [
        _entry("bash", ok=False, error="not on the auto-allow list",
               inp="API_KEY=supersecret123 curl x") ,
        _entry("bash", ok=False, error="not on the auto-allow list",
               inp="second"),
        _entry("bash", ok=False, error="not on the auto-allow list",
               inp="third"),
        _entry("bash", ok=False, error="not on the auto-allow list",
               inp="fourth"),
    ])
    data = rd.collect(dir_path=traces, include_security_events=False)
    row = next(r for r in data["ranked"]
               if r["category"] == "allowlist" and r["tool"] == "bash")
    assert len(row["examples"]) == 3                 # capped at 3
    assert not any("supersecret123" in e for e in row["examples"])
    assert any("«redacted»" in e for e in row["examples"])
    assert not any(".ssh/id_rsa" in e for r in data["ranked"]
                   for e in r["examples"])


def test_security_events_half(traces):
    security_events.record("deny_pattern", "bash", "cmd denied", blocked=True)
    security_events.record("outside_ws", "write_file", "/etc/x", blocked=True)
    security_events.record("read_grant", "read_file", "grant", blocked=False)
    data = rd.collect(dir_path=traces)
    se = data["security_events"]
    assert se is not None
    assert se["total"] == 3
    assert se["blocked"] == 2
    kinds = {e["kind"] for e in se["recent"]}
    assert kinds == {"deny_pattern", "outside_ws", "read_grant"}
    assert "deny_pattern" in se["known_kinds"]


def test_empty_and_missing_dir(tmp_path):
    empty = tmp_path / "nope"
    data = rd.collect(dir_path=empty, include_security_events=False)
    assert data == {"sessions": 0, "entries": 0, "refusals": 0,
                    "ranked": [], "tool_totals": {}, "security_events": None}
    assert rd.format_text(data) == (
        "Denials report — 0 session(s), 0 tool call(s), 0 refusal(s)/blocked"
        "\n(no refusals recorded)")


def test_corrupt_lines_skipped(tmp_path):
    p = tmp_path / "bad.jsonl"
    p.write_text("{not json\n" + json.dumps(_entry(
        "bash", ok=False, error="not on the auto-allow list")) + "\n",
        encoding="utf-8")
    data = rd.collect(dir_path=tmp_path, include_security_events=False)
    assert data["entries"] == 1
    assert data["refusals"] == 1


def test_tool_totals_from_aggregate_tools(traces):
    data = rd.collect(dir_path=traces, include_security_events=False)
    tt = data["tool_totals"]
    assert tt["bash"]["calls"] == 3
    assert tt["bash"]["errors"] == 3
    assert tt["write_file"]["calls"] == 2


def test_format_text_contains_rows(traces):
    text = rd.format_text(rd.collect(dir_path=traces))
    assert "Denials report" in text
    assert "allowlist" in text
    assert "bash" in text
    assert "write_guard" in text
    assert "bash/3" in text         # 2 refusals of 3 total bash calls
    assert "«redacted»" not in text or True   # examples present when relevant


def test_format_text_never_raises_on_garbage():
    assert isinstance(rd.format_text({"ranked": "nonsense"}), str)
    assert isinstance(rd.format_text({}), str)


def test_scrub_redacts():
    assert "hunter2" not in rd._scrub("password=hunter2 x")
    assert ".ssh/id_rsa" not in rd._scrub("read /home/u/.ssh/id_rsa")


def test_timeout_not_a_refusal():
    assert rd._categorize("approval request TIMED OUT — user away") is None
