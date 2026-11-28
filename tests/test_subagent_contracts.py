"""Schema-validated subagent returns (run_subagent output_schema).

Covers the dependency-free JSON-Schema subset validator, JSON extraction
from prose/fences, and the run_subagent contract: valid final message →
structured_output; invalid → ONE correction round; still invalid →
structured_output=None + schema_error with the text kept. No real LLM —
everything streams from scripted fake clients.
"""

from __future__ import annotations

import json
from types import SimpleNamespace
from unittest.mock import MagicMock

from delfin.agent import subagents as SA


def _fake_event(**kwargs):
    return SimpleNamespace(
        type=kwargs.get("type", "text_delta"),
        text=kwargs.get("text", ""),
        tool_name=kwargs.get("tool_name", ""),
        tool_input=kwargs.get("tool_input", ""),
        tool_output=kwargs.get("tool_output", ""),
        input_tokens=kwargs.get("input_tokens", 0),
        output_tokens=kwargs.get("output_tokens", 0),
    )


def _round(text: str):
    """One streamed model round: a text delta + a usage delta."""
    return iter([
        _fake_event(type="text_delta", text=text),
        _fake_event(type="message_delta", input_tokens=10, output_tokens=5),
    ])


def _client_rounds(*texts: str):
    """A fake client whose successive stream_message calls yield *texts*."""
    client = MagicMock()
    client.stream_message = MagicMock(
        side_effect=[_round(t) for t in texts])
    client._permissions = None
    client.set_permissions = MagicMock()
    return client


_SCHEMA = {
    "type": "object",
    "properties": {
        "answer": {"type": "integer"},
        "confidence": {"type": "string",
                       "enum": ["low", "medium", "high"]},
        "sources": {"type": "array", "items": {"type": "string"}},
    },
    "required": ["answer", "confidence"],
}


# ---------------------------------------------------------------------------
# Validator subset unit tests
# ---------------------------------------------------------------------------


def test_validator_accepts_conforming_object():
    val = {"answer": 3, "confidence": "high", "sources": ["a.py", "b.py"]}
    assert SA.validate_json_schema(val, _SCHEMA) == []


def test_validator_rejects_wrong_types():
    errs = SA.validate_json_schema(
        {"answer": "three", "confidence": "high"}, _SCHEMA)
    assert any("expected integer" in e for e in errs)


def test_validator_rejects_missing_required():
    errs = SA.validate_json_schema({"answer": 1}, _SCHEMA)
    assert any("missing required property 'confidence'" in e for e in errs)


def test_validator_rejects_bad_enum():
    errs = SA.validate_json_schema(
        {"answer": 1, "confidence": "certain"}, _SCHEMA)
    assert any("not in enum" in e for e in errs)


def test_validator_checks_array_items():
    errs = SA.validate_json_schema(
        {"answer": 1, "confidence": "low", "sources": ["ok", 7]}, _SCHEMA)
    assert any("sources[1]" in e and "expected string" in e for e in errs)


def test_validator_bool_is_not_integer_and_vice_versa():
    assert SA.validate_json_schema(True, {"type": "integer"})
    assert SA.validate_json_schema(True, {"type": "number"})
    assert SA.validate_json_schema(1, {"type": "boolean"})
    assert SA.validate_json_schema(True, {"type": "boolean"}) == []


def test_validator_top_level_type_mismatch():
    errs = SA.validate_json_schema(["not", "an", "object"], _SCHEMA)
    assert errs and "expected object" in errs[0]


def test_validator_ignores_unsupported_keywords():
    # Subset contract: unknown keywords are ignored, never enforced.
    schema = {"type": "string", "minLength": 99, "pattern": "^zzz"}
    assert SA.validate_json_schema("hi", schema) == []


def test_extract_json_object_tolerates_fences_and_prose():
    text = "Here you go:\n```json\n{\"a\": 1, \"b\": [true]}\n```\nDone."
    assert SA.extract_json_object(text) == {"a": 1, "b": [True]}


def test_extract_json_object_skips_non_object_json():
    assert SA.extract_json_object("[1, 2, 3] then {\"k\": 5}") == {"k": 5}
    assert SA.extract_json_object("no json here at all") is None
    assert SA.extract_json_object("{broken: json") is None


# ---------------------------------------------------------------------------
# run_subagent contract
# ---------------------------------------------------------------------------


def test_schema_success_populates_structured_output():
    good = json.dumps({"answer": 42, "confidence": "high"})
    client = _client_rounds(good)
    res = SA.run_subagent(
        subagent_type="explore",
        description="structured answer",
        prompt="Report the answer as structured JSON please.",
        parent_client=client,
        parent_perms=None,
        output_schema=_SCHEMA,
    )
    assert res.error == ""
    assert res.schema_error == ""
    assert res.structured_output == {"answer": 42, "confidence": "high"}
    payload = res.to_payload()
    assert payload["structured_output"] == {"answer": 42,
                                            "confidence": "high"}
    assert "schema_error" not in payload
    # Only one model round — no correction needed.
    assert client.stream_message.call_count == 1


def test_schema_instruction_lands_in_system_prompt():
    good = json.dumps({"answer": 1, "confidence": "low"})
    client = _client_rounds(good)
    SA.run_subagent(
        subagent_type="explore",
        description="structured answer",
        prompt="Report the answer as structured JSON please.",
        parent_client=client,
        parent_perms=None,
        output_schema=_SCHEMA,
    )
    sys_prompt = client.stream_message.call_args_list[0].kwargs["system"]
    assert "Structured output contract" in sys_prompt
    assert '"required":["answer","confidence"]' in sys_prompt


def test_schema_mismatch_gets_one_correction_round():
    bad = "The answer is 42 and I am very confident."
    good = "```json\n" + json.dumps(
        {"answer": 42, "confidence": "high"}) + "\n```"
    client = _client_rounds(bad, good)
    res = SA.run_subagent(
        subagent_type="explore",
        description="structured answer",
        prompt="Report the answer as structured JSON please.",
        parent_client=client,
        parent_perms=None,
        output_schema=_SCHEMA,
    )
    assert client.stream_message.call_count == 2
    assert res.structured_output == {"answer": 42, "confidence": "high"}
    assert res.schema_error == ""
    # The corrected message becomes the final report.
    assert "answer" in res.final_text
    # The correction turn carried the validation errors to the child.
    corr_messages = client.stream_message.call_args_list[1].kwargs["messages"]
    assert "did not satisfy the required output schema" in \
        corr_messages[-1]["content"]
    assert corr_messages[-2] == {"role": "assistant", "content": bad}


def test_persistently_invalid_yields_schema_error_and_keeps_text():
    bad1 = "No JSON at all, just prose."
    bad2 = json.dumps({"answer": "still a string", "confidence": "high"})
    client = _client_rounds(bad1, bad2)
    res = SA.run_subagent(
        subagent_type="explore",
        description="structured answer",
        prompt="Report the answer as structured JSON please.",
        parent_client=client,
        parent_perms=None,
        output_schema=_SCHEMA,
    )
    assert client.stream_message.call_count == 2  # exactly ONE correction
    assert res.structured_output is None
    assert "expected integer" in res.schema_error
    assert res.final_text  # text stays available
    payload = res.to_payload()
    assert "structured_output" not in payload
    assert payload["schema_error"] == res.schema_error


def test_empty_correction_round_reports_first_errors():
    bad = "Prose only, no object."
    client = _client_rounds(bad, "")
    res = SA.run_subagent(
        subagent_type="explore",
        description="structured answer",
        prompt="Report the answer as structured JSON please.",
        parent_client=client,
        parent_perms=None,
        output_schema=_SCHEMA,
    )
    assert res.structured_output is None
    assert "no JSON object found" in res.schema_error
    assert res.final_text == bad


def test_no_output_schema_keeps_legacy_behaviour():
    client = _client_rounds("plain prose report")
    res = SA.run_subagent(
        subagent_type="explore",
        description="plain",
        prompt="Just report back in plain prose please.",
        parent_client=client,
        parent_perms=None,
    )
    assert res.structured_output is None
    assert res.schema_error == ""
    assert "structured_output" not in res.to_payload()
    assert client.stream_message.call_count == 1


def test_errored_run_skips_validation():
    # Stream raises → error path; validation must not run or overwrite it.
    client = MagicMock()
    client.stream_message = MagicMock(side_effect=RuntimeError("boom"))
    client._permissions = None
    client.set_permissions = MagicMock()
    res = SA.run_subagent(
        subagent_type="explore",
        description="failing child",
        prompt="This child is going to crash mid-stream.",
        parent_client=client,
        parent_perms=None,
        output_schema=_SCHEMA,
    )
    assert "sub-agent stream raised" in res.error
    assert res.structured_output is None
    assert res.schema_error == ""


# ---------------------------------------------------------------------------
# The verdict handoff: the preset's own output_schema, no caller schema
# ---------------------------------------------------------------------------

# A conforming verdict for a read-only reviewer: it runs no tests, so an
# honest answer reports 0/0/0 and ran=false. The caller passes NO
# output_schema — the preset's own is the contract under test.
_VERDICT_OK = {
    "verdict": "comparable: same level of theory for both",
    "evidence": ["a.inp:1 !B3LYP def2-TZVP"],
    "files_written": [],
    "open_questions": [],
    "live_run": {
        "ran": False,
        "where_ran": "the two calc folders read",
        "head_sha": "",
        "test_counts": {"passed": 0, "failed": 0, "skipped": 0},
    },
}
# The failure this exists to catch: a verdict with no live_run block. It
# looks fine as prose, but the run it rests on is invisible — the caller
# cannot tell whether a check actually happened.
_VERDICT_MISSING_RUN = {"verdict": "comparable", "evidence": ["x"]}


def _verdict_schema():
    return SA.SUBAGENT_PRESETS["method-reviewer"].output_schema


def test_the_verdict_presets_declare_a_schema_the_builtin_do_not():
    """The shipped verdict-returners carry the contract; the free-form ones
    do not, so their behaviour is untouched (never-worse)."""
    for name in ("verifier", "method-reviewer", "chemistry-reviewer",
                 "security-reviewer"):
        schema = SA.SUBAGENT_PRESETS[name].output_schema
        assert isinstance(schema, dict), f"{name} declares no output_schema"
        assert "live_run" in schema.get("required", []), (
            f"{name}'s schema does not make the run visible")
        assert ("skipped"
                in schema["properties"]["live_run"]["properties"]
                   ["test_counts"]["required"]), (
            f"{name}'s schema lets a verdict hide its skipped tests")
    # ``code-reviewer`` returns a DIFFERENT contract — a finding list, not a
    # verdict. It belongs to neither group: it must not silently drop back to
    # free prose (that is what made a review a checklist nobody could check),
    # and it must not acquire the verdict contract, which is about a run and
    # its test counts and a read-only reviewer has neither.
    review_schema = SA.SUBAGENT_PRESETS["code-reviewer"].output_schema
    assert review_schema is SA.REVIEW_FINDINGS_SCHEMA
    assert "live_run" not in review_schema.get("required", [])
    for name in ("explore", "plan", "general-purpose"):
        assert SA.SUBAGENT_PRESETS[name].output_schema is None, (
            f"{name} kept a schema it never had")


def test_a_conforming_verdict_returns_structured_with_no_caller_schema():
    """The live path: the caller names no schema, the preset's own is read
    and the returned verdict is validated against it — no correction round
    needed for a conforming one."""
    client = _client_rounds(json.dumps(_VERDICT_OK))
    res = SA.run_subagent(
        subagent_type="method-reviewer",
        description="compare two levels",
        prompt="Are these two calculations at the same level of theory?",
        parent_client=client,
        parent_perms=None,
        # NO output_schema — the preset's own is the contract under test.
    )
    assert res.error == ""
    assert res.schema_error == ""
    assert res.structured_output == _VERDICT_OK
    assert client.stream_message.call_count == 1


def test_a_verdict_without_a_live_run_gets_one_correction_round():
    """The check that catches failure #1-4: a verdict that omits the run it
    rests on does not pass silently — it is corrected once, and if the
    delegate still omits it, the error is the one the caller reads."""
    # Round one omits live_run; round two (the correction) supplies it.
    client = _client_rounds(json.dumps(_VERDICT_MISSING_RUN),
                            json.dumps(_VERDICT_OK))
    res = SA.run_subagent(
        subagent_type="method-reviewer",
        description="compare two levels",
        prompt="Are these two calculations at the same level of theory?",
        parent_client=client,
        parent_perms=None,
    )
    assert client.stream_message.call_count == 2
    assert res.structured_output == _VERDICT_OK
    # The correction carried the missing field to the delegate.
    corr = client.stream_message.call_args_list[1].kwargs["messages"]
    assert "live_run" in corr[-1]["content"]


def test_a_persistently_missing_run_is_a_schema_error_not_a_pass():
    """If the delegate never supplies the run, the verdict is NOT returned
    as structured output: structured_output is None and the caller gets the
    reason. That is the difference between a check that fires and prose."""
    client = _client_rounds(json.dumps(_VERDICT_MISSING_RUN),
                            json.dumps(_VERDICT_MISSING_RUN))
    res = SA.run_subagent(
        subagent_type="method-reviewer",
        description="compare two levels",
        prompt="Are these two calculations at the same level of theory?",
        parent_client=client,
        parent_perms=None,
    )
    assert res.structured_output is None
    assert "live_run" in res.schema_error
    assert res.final_text  # the text stays available
