"""A nested ``required:`` output_schema is not refused at load.

The ``output_schema`` self-check in ``subagents.validate_preset`` builds a
conforming sample and validates it, so that a preset whose own contract the
framework cannot fill is rejected before it becomes a live type. The sample
builder used to fill an object with an empty ``{}`` for every ``object``
property. That is a conforming value for a flat schema, but it is NOT a
conforming value for a nested object that carries its own ``required:`` —
and a nested block is exactly the shape a verdict contract needs (a
``live_run`` that must always say whether the run happened and where), so
the check refused the schema it exists to allow.

Input: a preset definition whose ``output_schema`` is a nested object with a
nested ``required:``. Output: ``validate_preset``'s problem list for it.

Semantics pinned here:

* a nested ``required:`` object is sampled to its own required properties,
  so the self-check accepts it;
* the builder and the validator stay on the SAME subset — for a battery of
  flat and nested schemas the sample ``validate_json_schema`` produces is
  itself accepted by ``validate_json_schema`` (if one keyword were added to
  one without the other, a self-check turns red);
* the never-worse half: a top-level schema that is not an object is still
  refused, because the return extraction keeps the first JSON OBJECT in the
  message and a non-object contract cannot be honoured;
* a flat schema (no nesting) validates exactly as before the builder
  learned to nest.

The red side of this file is measured at the baseline commit by a control
run: at ``origin/main`` the first test below fails with
``output_schema: does not validate a conforming sample`` and passes here.
"""
from __future__ import annotations

import pytest

from delfin.agent import subagents as S


# A nested ``required:`` — the shape a verdict contract needs (a live_run
# block that must always say whether the run happened and where).
NESTED_REQUIRED = (
    '---\n'
    'name: nested-required\n'
    'description: a preset with a nested-required schema\n'
    'mode: plan\n'
    'tools: [read_file]\n'
    'output_schema: {"type": "object", "properties": '
    '{"verdict": {"type": "string"}, "live_run": {"type": "object", '
    '"properties": {"ran": {"type": "boolean"}, '
    '"where_ran": {"type": "string"}}, '
    '"required": ["ran", "where_ran"]}}, '
    '"required": ["verdict", "live_run"]}\n'
    '---\n'
    'body\n'
)


def _frontmatter_schema(text: str) -> dict:
    """The parsed output_schema of a preset definition string."""
    from delfin.agent.skills import _parse_frontmatter
    meta, _body = _parse_frontmatter(text)
    return S._parse_preset_schema(meta.get("output_schema")) or {}


def test_a_nested_required_schema_is_not_refused_at_load():
    """The contract a verdict needs must be allowed, not refused.

    At the baseline the sample builder filled ``live_run`` with ``{}``, so
    the nested ``required: ["ran", "where_ran"]`` failed and the whole
    schema was refused -- a check that rejects the schema it exists to allow.
    """
    schema = _frontmatter_schema(NESTED_REQUIRED)
    assert schema, "the test reads no schema; the definition is malformed"
    problems = S.validate_preset(NESTED_REQUIRED)
    assert problems == [], (
        f"a nested required schema was refused at load: {problems}")


def test_the_builder_and_validator_agree_on_the_subset():
    """For a battery of shapes, the sample validates against its own schema.

    This is the lockstep the builder claims: if a new keyword or type were
    added to the validator without the builder, one of these would stop
    validating and the self-check would turn red.
    """
    shapes = [
        {"type": "string"},
        {"type": "integer"},
        {"type": "number"},
        {"type": "boolean"},
        {"type": "string", "enum": ["a", "b"]},
        {"type": "array", "items": {"type": "string"}},
        {"type": "object",
         "properties": {"a": {"type": "string"}, "b": {"type": "integer"}},
         "required": ["a", "b"]},
        {"type": "object",
         "properties": {"outer": {"type": "object",
                                  "properties": {"inner": {"type": "boolean"}},
                                  "required": ["inner"]}},
         "required": ["outer"]},
        {"type": "object",
         "properties": {"deep": {"type": "object",
                                 "properties": {"list":
                                                {"type": "array",
                                                 "items": {"type": "object",
                                                           "properties":
                                                           {"k": {"type": "string"}},
                                                           "required": ["k"]}}}}},
         "required": ["deep"]},
    ]
    for schema in shapes:
        sample = S._conforming_sample(schema)
        errs = S.validate_json_schema(sample, schema)
        assert errs == [], f"{schema!r} -> sample {sample!r} not accepted: {errs}"


def test_a_required_key_without_a_type_reads_as_a_string():
    """An unknown shape must never fabricate a failure.

    A ``required:`` entry the schema gives no type for is sampled as a
    string (the only value that is never the wrong type in a way the
    builder can predict), so the check stays silent rather than inventing a
    refusal for a schema it cannot see through.
    """
    schema = {"type": "object",
              "properties": {},
              "required": ["mystery"]}
    sample = S._conforming_sample(schema)
    assert sample == {"mystery": "x"}
    assert S.validate_json_schema(sample, schema) == []


def test_a_top_level_non_object_schema_is_still_refused():
    """The never-worse half: a non-object top schema cannot be honoured.

    The return extraction keeps the first JSON OBJECT in the delegate's
    message, so a contract that says "return an array" can never be met. The
    check must keep refusing it -- the builder learning to nest must not
    have silently started accepting what the extractor cannot deliver.
    """
    arr = ('---\nname: not-an-object\ndescription: d\nmode: plan\n'
           'tools: [read_file]\n'
           'output_schema: {"type": "array", "items": {"type": "string"}}\n'
           '---\nbody\n')
    problems = S.validate_preset(arr)
    assert any("does not validate a conforming sample" in p for p in problems), (
        f"a top-level array schema was accepted: {problems}")


def test_a_flat_schema_still_passes():
    """No nesting: the builder's new behaviour must not disturb the flat case
    it already handled -- the same pass it has always given."""
    flat = ('---\nname: flat\ndescription: d\nmode: plan\ntools: [read_file]\n'
            'output_schema: {"type": "object", "properties": '
            '{"verdict": {"type": "string"}, "evidence": {"type": "array", '
            '"items": {"type": "string"}}}, "required": ["verdict"]}\n'
            '---\nbody\n')
    assert S.validate_preset(flat) == []


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(pytest.main([__file__, "-q"]))
