"""Declarative orchestration executor (subagents.run_orchestration).

Covers: stage sequencing with a barrier, {{stage:NAME}} template
substitution, verify-vote majority logic, hard limits (stages / calls /
votes / no nesting) and the single summary telemetry record. No real
LLM — a scripted thread-safe fake client answers every child.
"""

from __future__ import annotations

import json
import threading
from types import SimpleNamespace

from delfin.agent import subagents as SA


def _event(**kw):
    return SimpleNamespace(
        type=kw.get("type", "text_delta"),
        text=kw.get("text", ""),
        tool_name=kw.get("tool_name", ""),
        tool_input=kw.get("tool_input", ""),
        tool_output=kw.get("tool_output", ""),
        input_tokens=kw.get("input_tokens", 0),
        output_tokens=kw.get("output_tokens", 0),
    )


class ScriptedClient:
    """Thread-safe fake parent client: replies via ``reply_fn(prompt)``.

    run_subagent shallow-copies the parent client; the copy shares the
    calls list and lock, so assertions see every child's traffic.
    """

    model = "fake-model"

    def __init__(self, reply_fn):
        self._reply_fn = reply_fn
        self.calls: list[dict] = []
        self._lock = threading.Lock()

    def set_permissions(self, perms):
        pass

    def stream_message(self, *, messages, system, max_tokens):
        prompt = messages[-1]["content"]
        with self._lock:
            self.calls.append({"system": system, "prompt": prompt})
        text = self._reply_fn(prompt)
        yield _event(type="text_delta", text=text)
        yield _event(type="message_delta", input_tokens=7, output_tokens=3)


def _quiet_store(monkeypatch, tmp_path):
    """Keep session/running/telemetry writes inside the test tmp dir."""
    monkeypatch.setattr(SA, "_SESSIONS_DIR", tmp_path / "sessions")
    monkeypatch.setattr(SA, "_RUNNING_DIR", tmp_path / "running")
    monkeypatch.setattr(SA, "_TELEMETRY_PATH", tmp_path / "telemetry.jsonl")


# ---------------------------------------------------------------------------
# Stage sequencing + template substitution
# ---------------------------------------------------------------------------


def test_stages_run_in_order_with_barrier_and_substitution(
        monkeypatch, tmp_path):
    _quiet_store(monkeypatch, tmp_path)

    def reply(prompt):
        if "find alpha" in prompt:
            return "ALPHA-FINDING"
        if "find beta" in prompt:
            return "BETA-FINDING"
        return "SYNTHESIS-DONE"

    client = ScriptedClient(reply)
    spec = {
        "stages": [
            {"name": "research", "parallel": [
                {"subagent_type": "explore", "description": "alpha",
                 "prompt": "please find alpha in the codebase"},
                {"subagent_type": "explore", "description": "beta",
                 "prompt": "please find beta in the codebase"},
            ]},
            {"name": "synthesis", "parallel": [
                {"subagent_type": "explore", "description": "combine",
                 "prompt": "Combine these findings: {{stage:research}}"},
            ]},
        ],
    }
    out = SA.run_orchestration(spec, client, None)
    assert "error" not in out
    assert set(out["stages"]) == {"research", "synthesis"}
    assert len(out["stages"]["research"]) == 2
    assert len(out["stages"]["synthesis"]) == 1
    # Deterministic ordering: payloads come back in submission order.
    assert out["stages"]["research"][0]["description"] == "alpha"
    assert out["stages"]["research"][0]["result"] == "ALPHA-FINDING"
    assert out["stages"]["research"][1]["result"] == "BETA-FINDING"
    # Barrier: the synthesis child ran LAST, after both research children.
    assert len(client.calls) == 3
    assert "Combine these findings" in client.calls[-1]["prompt"]
    # Substitution: stage-1 results (JSON-compact) reached the stage-2
    # prompt — the raw placeholder must be gone.
    assert "{{stage:research}}" not in client.calls[-1]["prompt"]
    assert "ALPHA-FINDING" in client.calls[-1]["prompt"]
    assert "BETA-FINDING" in client.calls[-1]["prompt"]
    # No verify step → all final-stage results count as verified.
    assert len(out["verified"]) == 1
    assert out["rejected"] == []


def test_substitute_stage_refs_leaves_unknown_stage_untouched():
    assert SA._substitute_stage_refs(
        "use {{stage:missing}} here", {"done": [{"result": "x"}]},
    ) == "use {{stage:missing}} here"
    rendered = SA._substitute_stage_refs(
        "data: {{stage:done}}", {"done": [{"result": "x"}]})
    assert rendered == 'data: [{"result":"x"}]'


# ---------------------------------------------------------------------------
# Verify votes + majority logic
# ---------------------------------------------------------------------------


def _verify_spec(votes: int) -> dict:
    return {
        "stages": [
            {"name": "work", "parallel": [
                {"subagent_type": "explore", "description": "claim",
                 "prompt": "make a claim about the codebase please"},
            ]},
        ],
        "verify": {
            "prompt_template": "Skeptically check this: {{result}}",
            "votes": votes,
        },
    }


def _voting_client(ballots: list[dict]) -> ScriptedClient:
    """Replies 'CLAIM' to work prompts; hands out ``ballots`` (in call
    order) to verify prompts."""
    lock = threading.Lock()
    remaining = list(ballots)

    def reply(prompt):
        if "Skeptically check" in prompt:
            with lock:
                ballot = remaining.pop(0) if remaining else {"refuted": False}
            return json.dumps(ballot)
        return "CLAIM"

    return ScriptedClient(reply)


def test_verify_majority_refuted_rejects_result(monkeypatch, tmp_path):
    _quiet_store(monkeypatch, tmp_path)
    client = _voting_client([
        {"refuted": True, "note": "contradicts source"},
        {"refuted": True, "note": "wrong file"},
        {"refuted": False, "note": "looks fine"},
    ])
    out = SA.run_orchestration(_verify_spec(votes=3), client, None)
    assert "error" not in out
    assert out["verified"] == []
    assert len(out["rejected"]) == 1
    v = out["rejected"][0]["verify"]
    assert v["votes"] == 3
    assert v["valid_votes"] == 3
    assert v["refuted"] == 2
    assert any("contradicts source" in n for n in v["notes"])


def test_verify_minority_refuted_keeps_result(monkeypatch, tmp_path):
    _quiet_store(monkeypatch, tmp_path)
    client = _voting_client([
        {"refuted": False}, {"refuted": True, "note": "meh"},
        {"refuted": False},
    ])
    out = SA.run_orchestration(_verify_spec(votes=3), client, None)
    assert len(out["verified"]) == 1
    assert out["rejected"] == []
    assert out["verified"][0]["verify"]["refuted"] == 1


def test_verify_substitutes_result_and_clamps_votes(monkeypatch, tmp_path):
    _quiet_store(monkeypatch, tmp_path)
    client = _voting_client([{"refuted": False}] * 3)
    out = SA.run_orchestration(_verify_spec(votes=9), client, None)
    assert "error" not in out
    vote_calls = [c for c in client.calls
                  if "Skeptically check" in c["prompt"]]
    assert len(vote_calls) == 3          # 9 requested, clamped to max 3
    # {{result}} substitution carried the claim into every ballot prompt.
    assert all("CLAIM" in c["prompt"] for c in vote_calls)
    # The skeptics run under the schema mechanism from run_subagent.
    assert all("Structured output contract" in c["system"]
               for c in vote_calls)


def test_verify_with_no_valid_ballots_keeps_result(monkeypatch, tmp_path):
    _quiet_store(monkeypatch, tmp_path)

    def reply(prompt):
        if "Skeptically check" in prompt:
            return "I refuse to answer in JSON."  # invalid twice → no ballot
        return "CLAIM"

    client = ScriptedClient(reply)
    out = SA.run_orchestration(_verify_spec(votes=1), client, None)
    assert len(out["verified"]) == 1
    assert out["verified"][0]["verify"]["valid_votes"] == 0


# ---------------------------------------------------------------------------
# Hard limits
# ---------------------------------------------------------------------------


def _one_call_stage(name: str) -> dict:
    return {"name": name, "parallel": [
        {"subagent_type": "explore", "description": f"work {name}",
         "prompt": f"do the work for stage {name} please"},
    ]}


def test_fourth_stage_rejected(monkeypatch, tmp_path):
    _quiet_store(monkeypatch, tmp_path)
    client = ScriptedClient(lambda p: "x")
    spec = {"stages": [_one_call_stage(n) for n in "abcd"]}
    out = SA.run_orchestration(spec, client, None)
    assert "too many stages: 4" in out["error"]
    assert client.calls == []            # nothing ran


def test_seventh_call_in_stage_rejected(monkeypatch, tmp_path):
    _quiet_store(monkeypatch, tmp_path)
    client = ScriptedClient(lambda p: "x")
    spec = {"stages": [{"name": "big", "parallel": [
        {"subagent_type": "explore", "description": f"c{i}",
         "prompt": "do a chunk of work for this stage"}
        for i in range(7)
    ]}]}
    out = SA.run_orchestration(spec, client, None)
    assert "7 calls > max 6" in out["error"]
    assert client.calls == []


def test_no_nested_orchestration(monkeypatch, tmp_path):
    _quiet_store(monkeypatch, tmp_path)
    client = ScriptedClient(lambda p: "x")
    child_perms = SimpleNamespace(subagent_depth=1)
    out = SA.run_orchestration(
        {"stages": [_one_call_stage("a")]}, client, child_perms)
    assert "cannot be nested" in out["error"]
    assert client.calls == []


def test_spec_shape_violations_rejected(monkeypatch, tmp_path):
    _quiet_store(monkeypatch, tmp_path)
    client = ScriptedClient(lambda p: "x")
    assert "error" in SA.run_orchestration("not a dict", client, None)
    assert "error" in SA.run_orchestration({"stages": []}, client, None)
    assert "non-empty 'name'" in SA.run_orchestration(
        {"stages": [{"parallel": [{}]}]}, client, None)["error"]
    assert "duplicate stage name" in SA.run_orchestration(
        {"stages": [_one_call_stage("a"), _one_call_stage("a")]},
        client, None)["error"]
    assert "missing 'prompt'" in SA.run_orchestration(
        {"stages": [{"name": "a", "parallel": [
            {"subagent_type": "explore", "description": "d"}]}]},
        client, None)["error"]
    assert "prompt_template" in SA.run_orchestration(
        {"stages": [_one_call_stage("a")], "verify": {"votes": 2}},
        client, None)["error"]
    assert client.calls == []


# ---------------------------------------------------------------------------
# Telemetry
# ---------------------------------------------------------------------------


def test_orchestration_writes_one_summary_telemetry_record(
        monkeypatch, tmp_path):
    _quiet_store(monkeypatch, tmp_path)
    client = _voting_client([{"refuted": False}])
    out = SA.run_orchestration(_verify_spec(votes=1), client, None)
    assert "error" not in out
    records = SA.read_telemetry(tmp_path / "telemetry.jsonl")
    orch = [r for r in records if r.get("subagent_type") == "orchestration"]
    assert len(orch) == 1
    rec = orch[0]
    assert rec["stages"] == ["work"]
    assert rec["verified"] == 1
    assert rec["rejected"] == 0
    # Tokens aggregate across children AND vote ballots (2 runs x 7/3).
    assert rec["input_tokens"] == 14
    assert rec["output_tokens"] == 6
    assert rec["error"] == ""
