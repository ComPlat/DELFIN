"""The single failure modes from the task brief, case by case.

Each case names the tool schema it is about, so the test reads like
the brief: missing required field, wrong name, wrong type, unknown
extra field, nested multi_edit shape, safe repairs, truncation.
"""

import json

import pytest

from delfin.agent import api_client as A
from delfin.agent import tool_args


def _schema(name):
    for tool in A._DOC_TOOLS_OPENAI:
        fn = tool.get("function") or {}
        if fn.get("name") == name:
            return fn.get("parameters") or {}
    raise KeyError(name)


def test_missing_required_field_names_it_and_the_types():
    res = tool_args.check(_schema("write_file"), {"path": "a.txt"})
    assert not res.ok
    assert "'content'" in res.error_text
    # the message teaches the correct call: required fields with type
    assert "string" in res.error_text


def test_wrong_field_name_gives_hint_not_silent_rename():
    # `file_path` instead of `path`: a hint, and the args are NOT
    # rewritten -- executing on a guess would be a silent rename.
    res = tool_args.check(_schema("read_file"), {"file_path": "a.txt"})
    assert not res.ok
    assert "file_path" in res.error_text
    assert "path" in res.error_text  # the suggestion
    assert "did you mean" in res.error_text
    assert res.args == {"file_path": "a.txt"}  # unchanged


def test_hint_comes_from_similarity_not_a_list():
    # No hand-maintained table maps "paths" -> "path"; similarity must.
    res = tool_args.check(_schema("read_file"), {"paths": "a.txt"})
    assert "did you mean 'path'?" in res.error_text
    # ... but a field nothing like any known field gets no invented hint
    res = tool_args.check(_schema("read_file"), {"zzzzqqqq": 1})
    assert "did you mean" not in res.error_text
    assert "known fields" in res.error_text


def test_number_as_string_is_repaired_and_recorded():
    res = tool_args.check(
        _schema("read_file"), {"path": "a.txt", "limit": "3"})
    assert res.ok
    assert res.args["limit"] == 3
    assert res.repairs == ["limit: '3' -> integer"]


def test_list_as_json_string_is_repaired_and_recorded():
    edits = json.dumps([
        {"old_string": "a", "new_string": "b"},
        {"old_string": "c", "new_string": "d"},
    ])
    res = tool_args.check(
        _schema("multi_edit"), {"path": "a.txt", "edits": edits})
    assert res.ok
    assert res.args["edits"][1]["new_string"] == "d"
    assert any("edits" in r for r in res.repairs)


def test_quoted_float_is_not_repaired():
    # only integers are safely unquoted; "1.5" must not become 1.5
    schema = {"type": "object",
              "properties": {"n": {"type": "number"}},
              "required": ["n"]}
    res = tool_args.check(schema, {"n": "1.5"})
    assert not res.ok
    assert "expected number" in res.error_text


def test_nested_multi_edit_shape_is_checked():
    # an edit object missing new_string: the error must reach INTO the
    # array element and name the field there.
    res = tool_args.check(_schema("multi_edit"),
                          {"path": "a.txt",
                           "edits": [{"old_string": "a"}]})
    assert not res.ok
    assert "'new_string'" in res.error_text
    assert "edits[0]" in res.error_text


def test_edits_as_string_is_rejected_with_shape():
    res = tool_args.check(
        _schema("multi_edit"), {"path": "a.txt", "edits": "not a list"})
    assert not res.ok
    assert "edits" in res.error_text
    assert "array" in res.error_text


def test_single_edit_as_json_string_is_repaired():
    # one array element arrived as a JSON string holding an object --
    # the same safe repair as at the property level, applied per item
    edit = json.dumps({"old_string": "a", "new_string": "b"})
    res = tool_args.check(
        _schema("multi_edit"), {"path": "a.txt", "edits": [edit]})
    assert res.ok
    assert res.args["edits"][0]["new_string"] == "b"
    assert any("edits[0]" in r for r in res.repairs)


def test_unknown_extra_field_is_rejected():
    res = tool_args.check(
        _schema("read_file"), {"path": "a.txt", "chunks": 4})
    assert not res.ok
    assert "chunks" in res.error_text
    assert "unknown field" in res.error_text


def test_boolean_where_integer_expected_is_rejected():
    # bool is an int subclass in Python; the checker must not be fooled
    schema = {"type": "object",
              "properties": {"n": {"type": "integer"}},
              "required": ["n"]}
    res = tool_args.check(schema, {"n": True})
    assert not res.ok
    assert "expected integer" in res.error_text


def test_enum_violation_names_the_allowed_values():
    schema = {"type": "object",
              "properties": {"mode": {"type": "string",
                                      "enum": ["a", "b"]}},
              "required": ["mode"]}
    res = tool_args.check(schema, {"mode": "c"})
    assert not res.ok
    assert '["a", "b"]' in res.error_text


def test_truncated_json_arguments_are_named_as_cut_off():
    blob = '{"path": "a.txt", "content": "one\\ntwo\\nthr'
    res = tool_args.check(_schema("write_file"), blob)
    assert not res.ok
    assert "cut off" in res.error_text
    assert "output limit" in res.error_text
    assert "smaller" in res.error_text


def test_malformed_but_not_truncated_json_is_not_called_cut_off():
    blob = '{"path": "a.txt", "content" "missing colon"}'
    res = tool_args.check(_schema("write_file"), blob)
    assert not res.ok
    assert "cut off" not in res.error_text


def test_arguments_as_non_dict_are_reported_not_crashed():
    for bad in ([], 3, None):
        res = tool_args.check(_schema("read_file"), bad)
        assert not res.ok
        assert "JSON object" in res.error_text


def test_empty_required_list_is_not_triggered():
    # list_files has no required fields today; the empty call is valid
    res = tool_args.check(_schema("list_files"), {})
    assert res.ok


def test_the_executor_mode_leaves_value_limits_to_the_tools():
    """unknown="hint" (the executor's mode): an enum miss, an out-of-range
    number and an extra field pass -- the tools clamp or answer with their
    own message (schedule_wakeup clamps delay_seconds, subagent lists its
    types). Types and required fields are still checked."""
    schema = {"type": "object", "required": ["n"], "properties": {
        "n": {"type": "integer", "maximum": 10},
        "kind": {"type": "string", "enum": ["a", "b"]}}}
    ok = tool_args.check(schema, {"n": 99, "kind": "zz", "extra": 1}, unknown="hint")
    assert ok.ok, ok.error_text
    assert not tool_args.check(schema, {"n": "x"}, unknown="hint").ok
    missing = tool_args.check(schema, {"m": 1}, unknown="hint")
    assert not missing.ok and "missing required field 'n'" in missing.error_text
