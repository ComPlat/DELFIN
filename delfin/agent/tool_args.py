"""Validate tool-call arguments against the tool's JSON schema.

A model calling a tool with bad arguments (missing required field,
wrong name, wrong type, unknown extra field, arguments cut off by the
output limit) used to fail wherever the executor happened to notice --
usually as a Python-flavoured message the model cannot learn from.
This module answers those calls with a short, schema-derived error that
names the field, its type, and the correct shape.

Design rules (from the task that introduced it):

- Pure: no I/O, no imports beyond the standard library. jsonschema is
  NOT a dependency of this project, so the supported subset of JSON
  Schema (type, enum, required, properties, items, minItems/maxItems,
  minimum/maximum -- exactly what the tool schemas in
  ``delfin.agent.api_client`` use) is validated here by hand.
- Safe repairs only, each recorded in the result: a JSON *string* that
  holds an object/array where the schema demands one, and a quoted
  integer ("3") where the schema says integer. Nothing else is guessed.
- Confusions are a HINT, never a silent rename: "unknown field
  ``file_path`` -- did you mean ``path``?" The suggestion comes from
  the schema's own field names (difflib similarity), not from a
  hand-maintained list, so it stays correct when tools change.
- Truncated JSON is named as such: a string that starts like an object
  or array but cannot parse, failing at its very end, is reported as
  cut off by the output limit -- with the advice to send a smaller
  edit, which is the only thing that helps in that situation.
"""

from __future__ import annotations

import difflib
import json
from dataclasses import dataclass, field
from typing import Any, Optional

__all__ = ["Result", "check"]


@dataclass
class Result:
    """Outcome of checking one call's arguments.

    ``ok`` True means ``args`` (after any safe repairs) satisfies the
    schema. ``error_text`` is None then. ``repairs`` lists every value
    that was changed, as "field: JSON string -> object", so a caller
    can log what the checker did -- a repair is a change to what the
    model sent, and changes must be visible.
    """

    ok: bool
    args: Any
    error_text: Optional[str] = None
    repairs: list[str] = field(default_factory=list)


_TYPE_NAMES = {
    str: "string", bool: "boolean", int: "integer",
    float: "number", list: "array", dict: "object",
}


def _type_name(value: Any) -> str:
    # bool is a subclass of int in Python; name it before int sees it.
    return _TYPE_NAMES.get(type(value), type(value).__name__)


def _matches_type(value: Any, want: str) -> bool:
    if want == "string":
        return isinstance(value, str)
    if want == "integer":
        return isinstance(value, int) and not isinstance(value, bool)
    if want == "number":
        return isinstance(value, (int, float)) and not isinstance(value, bool)
    if want == "boolean":
        return isinstance(value, bool)
    if want == "array":
        return isinstance(value, list)
    if want == "object":
        return isinstance(value, dict)
    return True  # unknown keyword: nothing to check


def _describe(schema: dict) -> str:
    """One-line shape of a property schema, for error messages."""
    parts = []
    types = schema.get("type")
    if isinstance(types, list):
        parts.append("|".join(types))
    elif types:
        parts.append(str(types))
    if "enum" in schema and schema["enum"]:
        parts.append("one of " + json.dumps(schema["enum"]))
    if schema.get("type") == "array" and isinstance(schema.get("items"), dict):
        parts.append("of " + _describe(schema["items"]))
    return " ".join(parts)


def _suggestion(name: str, known: list[str]) -> Optional[str]:
    """The closest schema field name, from similarity alone."""
    if not known:
        return None
    best = max(known, key=lambda k: difflib.SequenceMatcher(
        None, name.lower(), k.lower()).ratio())
    ratio = difflib.SequenceMatcher(
        None, name.lower(), best.lower()).ratio()
    return best if ratio >= 0.6 else None


def _looks_truncated(text: str) -> bool:
    """A string that begins an object/array but dies mid-way.

    Heuristic, stated: only a JSONDecodeError whose position sits at
    (or one past) the end of the stripped text counts -- a parse that
    fails long before the end is malformed JSON, not cut off. The cut
    by an output limit always lands at the very end.
    """
    t = text.strip()
    if not (t.startswith("{") or t.startswith("[")):
        return False
    try:
        json.loads(t)
    except json.JSONDecodeError as exc:
        # An unterminated string always means the text simply stops;
        # other errors count only when the parser ran out of text at
        # the very end.
        if exc.msg.startswith("Unterminated string"):
            return True
        return exc.pos >= max(len(t) - 2, 0)
    return False


def _repair(value: Any, schema: dict, path: str, repairs: list[str]) -> Any:
    """Apply the two safe repairs, in place on the caller's structure.

    Repairs are conservative on purpose: a JSON string that parses to
    the object/array the schema demands, and a quoted integer. A quoted
    float, a "true"/"false" string, a number-as-string where the schema
    says string -- none of these are repaired, because each of them is
    also a plausible intentional value and guessing wrong here turns a
    readable error into a wrong execution.
    """
    want = schema.get("type")
    if isinstance(want, list):
        want = want[0] if want else None
    if isinstance(value, str):
        if want in ("object", "array") or (
                not want and ("properties" in schema or "items" in schema)):
            try:
                parsed = json.loads(value)
            except json.JSONDecodeError:
                return value
            if _matches_type(parsed, want or (
                    "object" if "properties" in schema else "array")):
                repairs.append(f"{path}: JSON string -> {want}")
                return parsed
        elif want == "integer" and value.strip().lstrip("-").isdigit():
            repairs.append(f"{path}: '{value}' -> integer")
            return int(value)
    return value


def _check_object(obj: dict, schema: dict, path: str,
                  repairs: list[str]) -> Optional[str]:
    props = schema.get("properties") or {}
    required = schema.get("required") or []
    # Unknown fields FIRST: when a required field is missing and an
    # unknown one is present, the confusion hint is the message that
    # teaches the model something -- "missing path" alone hides the
    # actual mistake (it sent `file_path`).
    for key, value in obj.items():
        here = f"{path}.{key}" if path else key
        if key not in props:
            hint = _suggestion(key, list(props))
            return (f"unknown field '{key}'"
                    + (f" -- did you mean '{hint}'?" if hint else "")
                    + f"; known fields: " + ", ".join(sorted(props)))
        # repair in place so the caller's dict carries the fixed value
        obj[key] = repaired = _repair(value, props[key], here, repairs)
        err = _check_value(repaired, props[key], here, repairs)
        if err:
            return err
    for key in required:
        if key not in obj:
            shape = _describe(props[key]) if key in props else "required"
            where = f"{path}: " if path else ""
            return (f"{where}missing required field '{key}'"
                    + (f" ({shape})" if shape != "required" else "")
                    + "; required: " + ", ".join(
                        f"{k}" + (f" ({_describe(props[k])})"
                                  if k in props else "")
                        for k in required))
    return None


def _check_array(arr: list, schema: dict, path: str,
                 repairs: list) -> Optional[str]:
    items = schema.get("items") or {}
    if "minItems" in schema and len(arr) < schema["minItems"]:
        return (f"{path}: needs at least {schema['minItems']} items, "
                f"got {len(arr)}")
    if "maxItems" in schema and len(arr) > schema["maxItems"]:
        return (f"{path}: allows at most {schema['maxItems']} items, "
                f"got {len(arr)}")
    for i, item in enumerate(arr):
        # repair array items too: a JSON string holding the object the
        # items schema demands is as safe here as at the property level
        item = _repair(item, items, f"{path}[{i}]", repairs)
        arr[i] = item
        err = _check_value(item, items, f"{path}[{i}]", repairs)
        if err:
            return err
    return None


def _check_value(value: Any, schema: dict, path: str,
                 repairs: list[str]) -> Optional[str]:
    if isinstance(schema.get("enum"), list) and schema["enum"]:
        if value not in schema["enum"]:
            return (f"{path}: must be one of "
                    f"{json.dumps(schema['enum'])}, got "
                    f"{json.dumps(value, default=str)}")
    want = schema.get("type")
    if isinstance(want, list):
        if not any(_matches_type(value, w) for w in want):
            return (f"{path}: expected {'|'.join(want)} "
                    f"({_describe(schema)}), got {_type_name(value)}")
    elif want and not _matches_type(value, want):
        return (f"{path}: expected {want} ({_describe(schema)}), "
                f"got {_type_name(value)}: "
                f"{json.dumps(value, default=str)[:120]}")
    if isinstance(value, dict):
        return _check_object(value, schema, path, repairs)
    if isinstance(value, list):
        return _check_array(value, schema, path, repairs)
    for bound, op in (("minimum", ">="), ("maximum", "<=")):
        if bound in schema and isinstance(value, (int, float)) \
                and not isinstance(value, bool):
            if op == ">=" and value < schema[bound]:
                return f"{path}: must be >= {schema[bound]}, got {value}"
            if op == "<=" and value > schema[bound]:
                return f"{path}: must be <= {schema[bound]}, got {value}"
    return None


def check(schema: dict, args: Any) -> Result:
    """Check ``args`` against a tool's parameter ``schema``.

    ``args`` is the arguments object as the caller received it -- a
    dict for every tool here, but a non-dict is reported, not crashed
    on. The returned ``args`` carries every applied repair; when ``ok``
    is False the caller must NOT execute the tool but answer with
    ``error_text``.
    """
    if isinstance(args, str):
        # The whole argument blob arrived as one JSON string: repair
        # (object expected) or, when it will not parse and looks cut
        # off, say so -- that call lost its tail to the output limit.
        if _looks_truncated(args):
            return Result(
                False, args,
                "arguments were cut off -- likely the output limit; "
                "send a smaller edit")
        try:
            parsed = json.loads(args)
        except json.JSONDecodeError:
            return Result(False, args,
                          "arguments must be a JSON object, not a bare "
                          "string; send the fields as a JSON object")
        args = parsed
        if not isinstance(args, dict):
            return Result(False, args,
                          "arguments must be a JSON object of fields, "
                          f"got {_type_name(args)}")
    repairs: list[str] = []
    if not isinstance(args, dict):
        return Result(False, args,
                      "arguments must be a JSON object of fields, got "
                      f"{_type_name(args)}")
    root = schema if schema.get("type") == "object" or "properties" in schema \
        else {"type": "object", "properties": {}}
    err = _check_object(args, root, "", repairs)
    if err:
        return Result(False, args, err)
    return Result(True, args, None, repairs)

