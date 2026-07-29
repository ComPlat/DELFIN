"""Tests for the delegate-report verification pass in delfin.agent.subagents.

Two of the cases below reproduce real failure modes: a delegate that reports
a green test suite it never ran, and a delegate that reports an unchanged
contract without having looked at anything. Both were confident, both were
wrong, and both would have been passed on as fact by the parent.
"""

from __future__ import annotations

import json

from delfin.agent import subagents as SA


def _payload(text: str, calls=None, **extra) -> dict:
    """A minimal sub-agent result payload as the parent receives it."""
    p = {
        "subagent_type": "general-purpose",
        "description": "delegated task",
        "result": text,
    }
    if calls is not None:
        p["tool_calls"] = calls
    p.update(extra)
    return p


def _read(path: str) -> dict:
    return {"name": "read_file", "input": json.dumps({"path": path}),
            "output": "file body"}


def _bash(cmd: str, output: str = "") -> dict:
    return {"name": "bash", "input": json.dumps({"command": cmd}),
            "output": output}


def _write(path: str) -> dict:
    return {"name": "edit_file", "input": json.dumps({"path": path}),
            "output": "ok"}


# ---------------------------------------------------------------------------
# Field case 1 — "all tests pass" with no test-running tool call
# ---------------------------------------------------------------------------


def test_tests_pass_claim_without_test_tool_call_is_unsupported():
    verdict = SA.verify_subagent_report(_payload(
        "I regenerated the catalogue. All 60 tests pass and the suite is "
        "green.",
        [_read("delfin/agent/api_client.py")]))
    assert verdict["status"] == "checked"
    kinds = [u["kind"] for u in verdict["unsupported"]]
    assert "test_result" in kinds
    flagged = next(u for u in verdict["unsupported"]
                   if u["kind"] == "test_result")
    assert "no test-running tool call" in flagged["why"]
    assert "run the tests yourself" in flagged["recheck"]


def test_tests_pass_claim_with_a_real_test_run_is_supported():
    verdict = SA.verify_subagent_report(_payload(
        "All tests pass — the suite is green.",
        [_bash("python -m pytest tests/ -q", "60 passed in 4.2s")]))
    assert verdict["unsupported"] == []
    assert any(s["kind"] == "test_result" for s in verdict["supported"])
    assert verdict["evidence"]["test_runs"] == 1


def test_claimed_test_count_absent_from_the_recorded_output_is_unsupported():
    verdict = SA.verify_subagent_report(_payload(
        "The suite is green: 60 passed.",
        [_bash("python -m pytest tests/ -q", "12 passed in 1.1s")]))
    flagged = [u for u in verdict["unsupported"] if u["kind"] == "test_result"]
    assert flagged, verdict
    assert "60" in flagged[0]["why"]


def test_hedged_test_claim_is_not_flagged():
    verdict = SA.verify_subagent_report(_payload(
        "I could not run them; presumably all tests pass.",
        [_read("delfin/agent/engine.py")]))
    assert [u for u in verdict["unsupported"]
            if u["kind"] == "test_result"] == []


# ---------------------------------------------------------------------------
# Field case 2 — a completion/"contract unchanged" claim with zero evidence
# ---------------------------------------------------------------------------


def test_completion_claim_with_zero_tool_calls_is_unsupported():
    verdict = SA.verify_subagent_report(_payload(
        "Done. Parameter names, types, enums and required lists are "
        "bit-identical; the contract is unchanged.",
        []))
    assert verdict["status"] == "checked"
    assert verdict["evidence"]["tool_calls"] == 0
    completion = [u for u in verdict["unsupported"]
                  if u["kind"] == "completion"]
    assert len(completion) >= 2
    assert all("no tool calls at all" in u["why"] for u in completion)
    assert "bit-identical" in verdict["notice"]


def test_mutation_claim_without_any_write_is_unsupported():
    verdict = SA.verify_subagent_report(_payload(
        "I implemented the missing guard and fixed the enum.",
        [_read("delfin/agent/subagents.py")]))
    flagged = [u for u in verdict["unsupported"] if u["kind"] == "completion"]
    assert flagged, verdict
    assert "no file-writing tool call" in flagged[0]["why"]
    assert "diff the files" in flagged[0]["recheck"]


def test_finding_about_existing_code_is_not_read_as_a_mutation_claim():
    # A read-only delegate REPORTING that an enum was added to a parameter
    # is describing the code, not claiming it changed anything.
    verdict = SA.verify_subagent_report(_payload(
        "The two module versions differ: an enum was added to the "
        "'mode' parameter of delfin/agent/subagents.py.",
        [_read("delfin/agent/subagents.py")]))
    assert [u for u in verdict["unsupported"]
            if u["kind"] == "completion"] == []


def test_german_report_claims_are_classified():
    verdict = SA.verify_subagent_report(_payload(
        "Alle Tests sind grün. Fertig.", [_read("delfin/agent/engine.py")]))
    flagged = [u for u in verdict["unsupported"] if u["kind"] == "test_result"]
    assert flagged, verdict
    assert "Alle Tests sind grün" in verdict["notice"]


# ---------------------------------------------------------------------------
# File / location claims
# ---------------------------------------------------------------------------


def test_file_the_delegate_never_opened_is_unsupported():
    verdict = SA.verify_subagent_report(_payload(
        "The handler lives in delfin/agent/engine.py:281.",
        [_read("delfin/agent/subagents.py")]))
    refs = [u for u in verdict["unsupported"]
            if u["kind"] in ("file_reference", "location")]
    assert refs, verdict
    assert any("engine.py" in u["claim"] for u in refs)


def test_file_the_delegate_read_is_supported():
    verdict = SA.verify_subagent_report(_payload(
        "The runner is wired in delfin/agent/subagents.py.",
        [_read("delfin/agent/subagents.py")]))
    assert [u for u in verdict["unsupported"]
            if u["kind"] in ("file_reference", "location")] == []
    assert any(s["kind"] == "file_reference" for s in verdict["supported"])


# ---------------------------------------------------------------------------
# Clean reports, honesty about missing bookkeeping, run status
# ---------------------------------------------------------------------------


def test_well_evidenced_report_passes_clean():
    verdict = SA.verify_subagent_report(_payload(
        "I implemented the verification pass in delfin/agent/subagents.py "
        "and ran the suite: 12 passed.",
        [_read("delfin/agent/subagents.py"),
         _write("delfin/agent/subagents.py"),
         _bash("python -m pytest tests/test_subagents.py -q",
               "12 passed in 0.9s")]))
    assert verdict["unsupported"] == []
    assert verdict["claims_checked"] >= 3
    assert "all are backed by its own evidence" in verdict["notice"]
    # The notice never sells the check as a guarantee.
    assert "consistency check only" in verdict["notice"]


def test_missing_trace_is_not_treated_as_fabrication():
    # No tool_calls key at all (an externally attached runner): absence of
    # bookkeeping must not be reported as an unsupported claim.
    verdict = SA.verify_subagent_report(_payload("All tests pass. Done."))
    assert verdict["status"] == "no_trace"
    assert verdict["unsupported"] == []
    assert "none of its claims could be cross-checked" in verdict["notice"]


def test_empty_report_yields_no_verdict_noise():
    verdict = SA.verify_subagent_report(_payload("", []))
    assert verdict["status"] == "no_report"
    assert verdict["notice"] == ""


def test_failed_run_is_surfaced_in_the_notice():
    verdict = SA.verify_subagent_report(_payload(
        "Everything is done.", [], error="wall-clock budget exhausted (300s)"))
    assert verdict["run_incomplete"].startswith("wall-clock")
    assert "did not finish cleanly" in verdict["notice"]


# ---------------------------------------------------------------------------
# Never raise, never alter the delegate's own text
# ---------------------------------------------------------------------------


def test_garbage_payloads_never_raise():
    for bad in (None, "", [], 3, {}, {"result": None},
                {"result": 5, "tool_calls": "not a list"},
                {"result": "x", "tool_calls": [None, 7, {"name": None},
                                               {"input": {"path": 3}}]},
                {"result": "x", "tool_calls": [{"name": "bash",
                                                "input": "{broken"}]},
                {"result": "x", "worktree": "not a dict"}):
        verdict = SA.verify_subagent_report(bad)
        assert isinstance(verdict, dict)
        assert verdict["status"] in ("checked", "no_report", "no_trace",
                                     "unavailable")
        assert isinstance(SA.collect_report_evidence(bad), dict)
        # Annotation degrades to returning the input untouched.
        SA.attach_verification(bad)


def test_scanners_never_raise_on_odd_text():
    for text in ("", "   ", "```\nall tests pass\n```", "x" * 50_000,
                 "\n".join(["- Done."] * 500)):
        assert isinstance(SA.scan_report_test_claims(text), list)
        assert isinstance(SA.scan_report_completion_claims(text), list)


def test_delegate_text_is_never_altered_and_input_not_mutated():
    text = "All tests pass.\n\n- Implemented everything.\n"
    payload = _payload(text, [])
    before = json.dumps(payload, sort_keys=True)
    annotated = SA.attach_verification(payload)
    assert annotated["result"] == text
    assert json.dumps(payload, sort_keys=True) == before   # input untouched
    assert annotated is not payload
    for key, value in payload.items():
        assert annotated[key] == value


def test_annotation_is_idempotent():
    payload = SA.attach_verification(_payload("All tests pass.", []))
    again = SA.attach_verification(payload)
    assert again == payload
    assert list(again.keys()) == list(payload.keys())


# ---------------------------------------------------------------------------
# Wiring: the parent cannot miss the verdict
# ---------------------------------------------------------------------------


def test_to_payload_carries_the_verdict_inline():
    res = SA.SubagentResult(
        subagent_type="general-purpose",
        description="regenerate the catalogue",
        final_text="All 60 tests pass; the contract is unchanged.",
        tool_calls=[{"name": "read_file",
                     "input": '{"path": "delfin/agent/api_client.py"}',
                     "output": "..."}],
        sa_id="ab12cd34",
    )
    payload = res.to_payload()
    # The delegate's own text survives verbatim ...
    assert payload["result"] == res.final_text
    # ... and the notice heads the JSON the parent reads.
    serialized = json.dumps(payload, ensure_ascii=False)
    assert serialized.startswith('{"verification_notice":')
    assert "[subagent-verify]" in serialized
    assert "All 60 tests pass" in payload["verification_notice"]
    assert payload["verification"]["unsupported"]


def test_to_payload_of_a_clean_run_still_states_the_limit():
    res = SA.SubagentResult(
        subagent_type="explore",
        description="locate the handler",
        final_text="The handler is in delfin/agent/subagents.py.",
        tool_calls=[{"name": "read_file",
                     "input": '{"path": "delfin/agent/subagents.py"}',
                     "output": "..."}],
    )
    payload = res.to_payload()
    assert payload["verification"]["unsupported"] == []
    assert "consistency check only" in payload["verification"]["limits"]


def test_collected_background_result_is_verified_too(tmp_path, monkeypatch):
    monkeypatch.setattr(SA, "_SESSIONS_DIR", tmp_path / "sessions")
    monkeypatch.setattr(SA, "_RUNNING_DIR", tmp_path / "running")
    (tmp_path / "sessions").mkdir()
    (tmp_path / "sessions" / "bg01.json").write_text(json.dumps({
        "sa_id": "bg01",
        "subagent_type": "general-purpose",
        "description": "background job",
        "messages": [{"role": "user", "content": "do it"},
                     {"role": "assistant",
                      "content": "Done — all tests pass."}],
        "interactions": [{"name": "read_file",
                          "input": '{"path": "delfin/agent/engine.py"}',
                          "output": "..."}],
        "error": "",
    }), encoding="utf-8")
    out = SA.get_subagent_result("bg01")
    assert out["status"] == "finished"
    assert out["final_text"] == "Done — all tests pass."
    assert "[subagent-verify]" in out["verification_notice"]
    assert any(u["kind"] == "test_result"
               for u in out["verification"]["unsupported"])


def test_evidence_summary_stays_compact():
    calls = [_read(f"pkg/module_{i}.py") for i in range(50)]
    verdict = SA.verify_subagent_report(_payload("Done.", calls))
    ev = verdict["evidence"]
    assert ev["tool_calls"] == 50
    assert len(ev["files_touched"]) <= 8
    assert len(ev["tools"]) <= 12
    assert len(json.dumps(verdict)) < 4000
