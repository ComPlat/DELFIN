"""Structured output as a harness service (structured_output module).

Covers request_structured (happy path, retry-with-errors, exhausted
retries, transport failure), the re-export contract that keeps the old
subagents import paths valid, and the memory-distill adoption (structured
setting on → schema-validated facts; off → legacy free-text path
byte-identical). No real LLM — scripted fakes only.
"""

from __future__ import annotations

import json
from types import SimpleNamespace

from delfin.agent import memory_distill as md
from delfin.agent import structured_output as so


class FakeClient:
    """Scripted streaming client: successive calls yield queued replies."""

    def __init__(self, *replies: str):
        self.replies = list(replies)
        self.calls: list[dict] = []

    def stream_message(self, *, messages, system="", max_tokens=0):
        self.calls.append({
            "prompt": str((messages or [{}])[-1].get("content", "")),
            "system": system,
            "max_tokens": max_tokens,
        })
        text = self.replies.pop(0)
        yield SimpleNamespace(type="text_delta", text=text)
        yield SimpleNamespace(type="message_delta", input_tokens=5,
                              output_tokens=3)


class RaisingClient:
    def stream_message(self, **kwargs):
        raise RuntimeError("socket down")
        yield  # pragma: no cover


_SCHEMA = {
    "type": "object",
    "required": ["answer"],
    "properties": {"answer": {"type": "integer"}},
}


# ---------------------------------------------------------------------------
# request_structured
# ---------------------------------------------------------------------------


def test_happy_path_with_fenced_json():
    client = FakeClient(
        "Here you go:\n```json\n" + json.dumps({"answer": 42}) + "\n```\n")
    res = so.request_structured(
        client, prompt="What is the answer?", schema=_SCHEMA)
    assert res["data"] == {"answer": 42}
    assert res["error"] == ""
    assert res["attempts"] == 1
    assert "42" in res["raw"]
    assert len(client.calls) == 1
    # The schema contract is appended to the system prompt.
    assert "Structured output contract" in client.calls[0]["system"]
    assert '"required":["answer"]' in client.calls[0]["system"]


def test_system_prefix_is_preserved():
    client = FakeClient(json.dumps({"answer": 1}))
    so.request_structured(client, prompt="p", schema=_SCHEMA,
                          system="You are a distiller.")
    assert client.calls[0]["system"].startswith("You are a distiller.")
    assert "Structured output contract" in client.calls[0]["system"]


def test_retry_on_invalid_then_success():
    client = FakeClient(
        json.dumps({"answer": "forty-two"}),          # invalid: wrong type
        json.dumps({"answer": 42}),                   # corrected
    )
    res = so.request_structured(
        client, prompt="What is the answer?", schema=_SCHEMA, retries=1)
    assert res["attempts"] == 2
    assert res["data"] == {"answer": 42}
    assert res["error"] == ""
    # The correction prompt carried the validation error back to the model.
    second_prompt = client.calls[1]["prompt"]
    assert "expected integer" in second_prompt
    assert "What is the answer?" in second_prompt
    # The first prompt did not contain correction chatter.
    assert "did not satisfy" not in client.calls[0]["prompt"]


def test_exhausted_retries_keeps_error_and_raw():
    bad = "no JSON here, just prose"
    client = FakeClient(bad, bad)
    res = so.request_structured(
        client, prompt="p", schema=_SCHEMA, retries=1)
    assert res["data"] is None
    assert "no JSON object found" in res["error"]
    assert res["raw"] == bad
    assert res["attempts"] == 2
    assert len(client.calls) == 2


def test_zero_retries_means_single_attempt():
    client = FakeClient("prose only")
    res = so.request_structured(client, prompt="p", schema=_SCHEMA, retries=0)
    assert res["attempts"] == 1
    assert res["data"] is None
    assert len(client.calls) == 1


def test_transport_failure_reports_error_without_raising():
    res = so.request_structured(
        RaisingClient(), prompt="p", schema=_SCHEMA, retries=3)
    assert res["data"] is None
    assert "socket down" in res["error"]
    assert res["attempts"] == 1          # broken client: no blind retries


def test_invalid_schema_argument_is_an_error_not_a_raise():
    client = FakeClient("{}")
    res = so.request_structured(client, prompt="p", schema=None)  # type: ignore[arg-type]
    assert res["data"] is None
    assert "schema" in res["error"]
    assert client.calls == []            # no tokens burned


# ---------------------------------------------------------------------------
# Re-export contract: the old subagents import paths stay valid
# ---------------------------------------------------------------------------


def test_subagents_reexports_are_the_same_objects():
    from delfin.agent import subagents as SA
    assert SA.validate_json_schema is so.validate_json_schema
    assert SA.extract_json_object is so.extract_json_object
    assert SA._schema_instruction is so.schema_instruction
    assert "validate_json_schema" in SA.__all__
    assert "extract_json_object" in SA.__all__


def test_moved_helpers_still_behave():
    from delfin.agent.subagents import extract_json_object, validate_json_schema
    assert extract_json_object("x {\"a\": 1} y") == {"a": 1}
    assert validate_json_schema({"answer": 1}, _SCHEMA) == []
    assert validate_json_schema({}, _SCHEMA)


# ---------------------------------------------------------------------------
# memory_distill adoption (agent.auto_memory.structured)
# ---------------------------------------------------------------------------


_CHAT = [
    {"role": "user", "content": "please always use English strings in code"},
    {"role": "assistant", "content": "Understood, English strings only."},
    {"role": "user", "content": "and prefer the cheap-tier model for chores"},
    {"role": "assistant", "content": "OK."},
    {"role": "user", "content": "the benchmark suite lives under tests/bench"},
]


def _capture_store(monkeypatch):
    saved: list[str] = []
    monkeypatch.setattr("delfin.agent.memory_store.load_memories", lambda: [])
    monkeypatch.setattr("delfin.agent.memory_store.save_memory",
                        lambda text, source="": saved.append(text))
    return saved


def test_distill_structured_path_saves_validated_facts(monkeypatch):
    saved = _capture_store(monkeypatch)
    systems: list[str] = []

    def llm(prompt, system, settings):
        systems.append(system)
        return json.dumps(
            {"facts": ["feedback: keep code strings in English"]})

    n = md.distill_and_save(
        _CHAT,
        settings={"agent": {"auto_memory": {"enabled": True,
                                            "structured": True}}},
        llm_fn=llm,
    )
    assert n == 1
    assert saved == ["feedback: keep code strings in English"]
    assert len(systems) == 1             # one structured call, no fallback
    assert "Structured output contract" in systems[0]
    assert '"facts"' in systems[0]


def test_distill_structured_falls_back_to_legacy_parse(monkeypatch):
    saved = _capture_store(monkeypatch)
    systems: list[str] = []

    def llm(prompt, system, settings):
        systems.append(system)
        # Free text, never a JSON object: the structured call (attempt +
        # one retry) fails, then the legacy parse handles the same reply.
        return "feedback: always run the linter before pushing"

    n = md.distill_and_save(
        _CHAT,
        settings={"agent": {"auto_memory": {"enabled": True,
                                            "structured": True}}},
        llm_fn=llm,
    )
    assert n == 1
    assert saved == ["feedback: always run the linter before pushing"]
    # structured attempt + retry + legacy fallback call
    assert len(systems) == 3
    assert "Structured output contract" in systems[0]
    assert "Structured output contract" not in systems[-1]


def test_distill_default_path_unchanged_with_setting_off(monkeypatch):
    saved = _capture_store(monkeypatch)
    systems: list[str] = []

    def llm(prompt, system, settings):
        systems.append(system)
        return "feedback: keep commit messages plain"

    n = md.distill_and_save(
        _CHAT,
        settings={"agent": {"auto_memory": {"enabled": True}}},
        llm_fn=llm,
    )
    assert n == 1
    assert saved == ["feedback: keep commit messages plain"]
    assert len(systems) == 1             # exactly one legacy call
    assert "Structured output contract" not in systems[0]
    assert "One fact per line" in systems[0]   # the legacy briefing


def test_auto_memory_settings_structured_defaults_off():
    assert md.auto_memory_settings({})["structured"] is False
    assert md.auto_memory_settings(
        {"agent": {"auto_memory": {"structured": True}}})["structured"] is True
