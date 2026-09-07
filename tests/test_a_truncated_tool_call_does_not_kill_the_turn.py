"""A cut-off tool call must not be replayed to the endpoint.

An assistant tool_call is sent back on every later round of the turn, and
vLLM parses its arguments on the way IN. So a call the model truncated —
cut off mid-string while writing a file — does not fail once and get
retried. It makes every later request a 400 and the turn dies with a
backend error instead of a tool result the model could act on.

Observed 2026-09-07 on kit.deepseek-v4-flash, asked to build a module in a
user project. Its own words: "Ich baue das Skript und teste den reinen
Export." What came back was

    litellm.BadRequestError: Hosted_vllmException —
    "Unterminated string starting at: line 1 column 34 (char 33)"

and the task scored 0/3 across three samples with no write tool called at
all. The same defect was fixed once for the id field beside this one.
"""

import json

import pytest

from delfin.agent.api_client import _sendable_tool_args


def _is_object(text: str) -> bool:
    return isinstance(json.loads(text), dict)


def test_a_truncated_call_becomes_something_the_endpoint_accepts():
    """The exact shape from the field: a write_file call cut off inside
    the content string."""
    raw = '{"file_path": "export.py", "content": "import csv\\n\\ndef mai'
    with pytest.raises(json.JSONDecodeError):
        json.loads(raw)
    assert _is_object(_sendable_tool_args(raw))


def test_valid_arguments_are_passed_through_byte_for_byte():
    raw = '{"file_path": "export.py", "content": "x"}'
    assert _sendable_tool_args(raw) == raw


def test_near_json_is_repaired_rather_than_discarded():
    """A weak model's near-JSON carries real arguments; dropping them
    would turn a recoverable call into an empty one."""
    out = _sendable_tool_args("{'file_path': 'a.py', 'content': 'x',}")
    assert json.loads(out) == {"file_path": "a.py", "content": "x"}


@pytest.mark.parametrize("raw", ["", "   ", None, 42, "true", "[1,2]",
                                 "not json at all"])
def test_anything_unusable_becomes_an_empty_object(raw):
    out = _sendable_tool_args(raw)
    assert _is_object(out)


def test_a_dict_is_serialised():
    assert json.loads(_sendable_tool_args({"a": 1})) == {"a": 1}


def test_the_wire_copy_never_carries_unparseable_text():
    """The property that matters: whatever the model emitted, what goes
    back is a JSON object."""
    for raw in ('{"a": "b', "{", "}{", '{"a": }', '{"a": "b"} trailing',
                '```json\n{"a": 1}\n```'):
        assert _is_object(_sendable_tool_args(raw)), raw


def test_the_model_still_learns_the_call_was_broken():
    """Sanitising the wire copy must not silence the error. The executor
    reads the RAW emission, kept beside the wire copy for exactly this."""
    import inspect

    from delfin.agent import api_client

    src = inspect.getsource(api_client.OpenAIClient.stream_message)
    assert "_raw_args_by_id[_call_id] = _raw_emitted" in src
    assert '_raw_args = _raw_args_by_id.get(' in src
    # ... and the wire copy is the sanitised one, not the raw string.
    assert '"arguments": _sendable_tool_args(_raw_emitted)' in src
