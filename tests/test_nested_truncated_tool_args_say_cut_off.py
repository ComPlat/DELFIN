# SYNTHETIC test fixture written for DELFIN -- not ORCA output
"""A nested JSON string cut off by the output limit is named as such.

Field report (dashboard, Code mode, 2026-11-02): a call whose `options`
argument arrived as a JSON string that stopped mid-array was answered with
`options: expected array (array of object), got string: ...` -- a type error
that invites resending the same over-long call. The root-level cut-off
diagnosis ("arguments were cut off -- likely the output limit") already
exists; it must also apply to a nested field, because the cut lands wherever
the limit bites, not only at the end of the whole blob.
"""
import json

from delfin.agent import tool_args


_SCHEMA = {
    "type": "object",
    "properties": {
        "question": {"type": "string"},
        "options": {
            "type": "array",
            "minItems": 2,
            "maxItems": 6,
            "items": {
                "type": "object",
                "properties": {
                    "label": {"type": "string"},
                    "description": {"type": "string"},
                },
            },
        },
    },
    "required": ["question", "options"],
}


def test_nested_truncated_array_reports_cut_off_not_type_error():
    # A long options string whose JSON dies mid-array: how an output limit
    # cuts a call, nested inside the arguments object.
    blob = json.dumps([{"label": "Option %d" % i,
                        "description": "x" * 40} for i in range(60)])
    cut = blob[: len(blob) // 2]          # stops mid-array, no closing bracket
    assert tool_args._looks_truncated(cut)
    res = tool_args.check(_SCHEMA, {"question": "Which?",
                                    "options": cut})
    assert not res.ok
    assert "cut off" in (res.error_text or ""), res.error_text
    assert "expected array" not in (res.error_text or ""), res.error_text


def test_nested_truncated_object_reports_cut_off_not_type_error():
    obj = json.dumps({"question": "Which?", "payload": "x" * 400})
    cut = obj[: len(obj) // 2]
    schema = {"type": "object",
              "properties": {"payload": {"type": "object"}},
              "required": ["payload"]}
    res = tool_args.check(schema, {"payload": cut})
    assert not res.ok
    assert "cut off" in (res.error_text or ""), res.error_text


def test_well_formed_json_string_still_repairs_to_array():
    arr = json.dumps([{"label": "A", "description": "first"},
                      {"label": "B", "description": "second"}])
    res = tool_args.check(_SCHEMA, {"question": "Which?", "options": arr})
    assert res.ok, res.error_text
    assert isinstance(res.args["options"], list)


def test_malformed_non_truncated_string_keeps_type_error():
    # A parse error far from the end is malformed JSON, not a cut-off: the
    # existing type error stays (regression guard).
    res = tool_args.check(_SCHEMA, {"question": "Which?",
                                    "options": "[not json at all]"})
    assert not res.ok
    assert "cut off" not in (res.error_text or "")
