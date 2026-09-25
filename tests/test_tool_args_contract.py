"""The argument checker must hold against the REAL tool schemas.

Every tool advertised in ``_DOC_TOOLS_OPENAI`` is checked twice:
a minimal valid call (built from the schema itself) must pass, and
removing any single required field must produce an error naming that
field. Building the sample from the schema — rather than hand-writing
one call per tool — keeps this contract true when tools are added.
"""

import pytest

from delfin.agent import api_client as A
from delfin.agent import tool_args


def _schemas():
    for tool in A._DOC_TOOLS_OPENAI:
        fn = tool.get("function") or {}
        name = fn.get("name")
        params = fn.get("parameters") or {}
        if name:
            yield name, params


def _sample(schema):
    """The smallest value that satisfies this property schema."""
    if "enum" in schema and schema["enum"]:
        return schema["enum"][0]
    t = schema.get("type")
    if isinstance(t, list):
        t = t[0] if t else None
    if t == "array" or "items" in schema:
        items = schema.get("items") or {}
        n = max(schema.get("minItems", 1), 1)
        return [_sample(items) for _ in range(n)]
    if t == "object" or "properties" in schema:
        props = schema.get("properties") or {}
        req = schema.get("required") or []
        return {k: _sample(props[k]) for k in req if k in props}
    if t == "integer":
        return schema.get("minimum", 1)
    if t == "number":
        return float(schema.get("minimum", 1))
    if t == "boolean":
        return True
    return "x"


def _minimal_call(params):
    props = params.get("properties") or {}
    return {k: _sample(props[k]) for k in params.get("required") or []
            if k in props}


@pytest.mark.parametrize("name,params", list(_schemas()))
def test_minimal_valid_call_passes(name, params):
    res = tool_args.check(params, _minimal_call(params))
    assert res.ok, f"{name}: {res.error_text}"


@pytest.mark.parametrize("name,params", list(_schemas()))
def test_missing_each_required_field_is_named(name, params):
    required = params.get("required") or []
    if not required:
        pytest.skip("no required fields")
    call = _minimal_call(params)
    for field in required:
        broken = {k: v for k, v in call.items() if k != field}
        res = tool_args.check(params, broken)
        assert not res.ok, f"{name}: dropping {field} went through"
        assert field in (res.error_text or ""), (
            f"{name}: error does not name {field}: {res.error_text}")


def test_every_advertised_tool_has_a_schema_entry():
    names = [n for n, _ in _schemas()]
    assert len(names) == len(set(names)), "duplicate tool names"
    assert len(names) > 50, "catalogue unexpectedly small"
