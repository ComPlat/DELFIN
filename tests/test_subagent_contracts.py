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
