"""Schema-validated model output as a reusable harness service.

Historically the only schema-validated model path lived hardcoded inside
``subagents.run_subagent`` (the ``output_schema`` mechanism). This module
hoists the building blocks out so ANY caller — memory distillation, job
diagnosis, future tools — can request validated JSON without reimplementing
the extract/validate/retry loop:

  - ``validate_json_schema``: dependency-free JSON-Schema subset validator
  - ``extract_json_object``: first JSON object embedded in prose/fences
  - ``schema_instruction``: the system-prompt clause enforcing the contract
  - ``request_structured``: one-shot "give me a dict matching this schema"
    call against any client exposing ``stream_message``

``subagents`` imports the validator/extractor back from here (and keeps
re-exporting them), so existing import paths stay valid.

Supported schema subset (deliberately small; anything else in a schema is
IGNORED, never enforced):
  - "type": object / array / string / number / integer / boolean
  - "required" (on objects)
  - "properties" (on objects; validated only for keys that are present)
  - "enum"
  - "items" (single-schema form, applied to every array element)
"""

from __future__ import annotations

import json
from typing import Any


_SCHEMA_TYPES: dict[str, Any] = {
    "object": dict,
    "array": list,
    "string": str,
    "boolean": bool,
    "integer": int,
    "number": (int, float),
}


def validate_json_schema(value, schema: dict, path: str = "$") -> list[str]:
    """Validate ``value`` against the supported JSON-Schema subset.

    Returns a list of human-readable error strings; empty list = valid.
    Unknown/unsupported schema keywords are ignored (subset documented
    above). Never raises on malformed schemas — they just validate less.
    """
    errors: list[str] = []
    if not isinstance(schema, dict):
        return errors
    t = schema.get("type")
    if isinstance(t, str) and t:
        expected = _SCHEMA_TYPES.get(t)
        if expected is not None:
            ok = isinstance(value, expected)
            # bool is an int subclass in Python — a boolean must not
            # satisfy integer/number.
            if t in ("integer", "number") and isinstance(value, bool):
                ok = False
            # ...and an integer/number must not satisfy boolean.
            if t == "boolean" and not isinstance(value, bool):
                ok = False
            if not ok:
                errors.append(
                    f"{path}: expected {t}, got {type(value).__name__}")
                return errors  # wrong type — descending would only cascade
    if "enum" in schema and isinstance(schema.get("enum"), list):
        if value not in schema["enum"]:
            errors.append(
                f"{path}: value {value!r} not in enum {schema['enum']!r}")
    if isinstance(value, dict):
        req = schema.get("required")
        if isinstance(req, list):
            for k in req:
                if k not in value:
                    errors.append(f"{path}: missing required property {k!r}")
        props = schema.get("properties")
        if isinstance(props, dict):
            for k, sub in props.items():
                if k in value and isinstance(sub, dict):
                    errors.extend(
                        validate_json_schema(value[k], sub, f"{path}.{k}"))
    if isinstance(value, list):
        items = schema.get("items")
        if isinstance(items, dict):
            for i, item in enumerate(value):
                errors.extend(
                    validate_json_schema(item, items, f"{path}[{i}]"))
    return errors


def extract_json_object(text: str) -> dict | None:
    """First JSON *object* embedded in ``text``, or None.

    Tolerates markdown fences and surrounding prose: scans for ``{``
    candidates and raw-decodes from each until one parses to a dict.
    """
    if not isinstance(text, str) or not text.strip():
        return None
    dec = json.JSONDecoder()
    i = text.find("{")
    while i != -1:
        try:
            obj, _end = dec.raw_decode(text, i)
        except json.JSONDecodeError:
            i = text.find("{", i + 1)
            continue
        if isinstance(obj, dict):
            return obj
        i = text.find("{", i + 1)
    return None


def schema_instruction(schema: dict) -> str:
    """System-prompt clause enforcing a schema-shaped final message."""
    try:
        compact = json.dumps(schema, separators=(",", ":"),
                             ensure_ascii=False)
    except (TypeError, ValueError):
        compact = "{}"
    return (
        "\n\nStructured output contract: your FINAL message must be exactly "
        "one JSON object that validates against the JSON Schema below. No "
        "prose before or after it (a ```json fence around it is tolerated).\n"
        "Schema: " + compact
    )


def _collect_text(client, prompt: str, system: str,
                  max_tokens: int) -> tuple[str, str]:
    """One text-only model round. Returns ``(text, error)``; never raises."""
    parts: list[str] = []
    err = ""
    try:
        for event in client.stream_message(
            messages=[{"role": "user", "content": prompt}],
            system=system,
            max_tokens=max_tokens,
        ):
            if getattr(event, "type", "") == "text_delta" \
                    and getattr(event, "text", ""):
                parts.append(event.text)
    except Exception as exc:
        err = f"model call raised: {exc}"
    return "".join(parts).strip(), err


def request_structured(
    client,
    *,
    prompt: str,
    schema: dict,
    system: str = "",
    retries: int = 1,
    max_tokens: int = 800,
) -> dict:
    """Request a schema-conforming JSON object from any streaming client.

    ``client`` only needs a ``stream_message(messages=..., system=...,
    max_tokens=...)`` generator yielding events with ``type``/``text``
    attributes — the OpenAI client and test fakes both qualify.

    The schema contract (``schema_instruction``) is appended to ``system``;
    the reply is scanned for its first JSON object and validated against
    the subset validator. On mismatch the call is retried up to ``retries``
    more times with the validation errors appended to the prompt, so the
    model can self-correct. A transport failure stops retrying — a broken
    client will not get healthier mid-loop.

    Returns ``{"data": dict | None, "error": str, "raw": str,
    "attempts": int}``; ``data`` is the validated object or None with
    ``error`` explaining the last failure. ``raw`` always keeps the last
    reply text. Never raises.
    """
    result: dict = {"data": None, "error": "", "raw": "", "attempts": 0}
    if not isinstance(schema, dict) or not schema:
        result["error"] = "schema must be a non-empty dict"
        return result
    sys_prompt = (system or "") + schema_instruction(schema)
    base_prompt = str(prompt or "")
    current_prompt = base_prompt
    try:
        tries = 1 + max(0, int(retries))
    except (TypeError, ValueError):
        tries = 2
    last_errs: list[str] = []
    for _ in range(tries):
        result["attempts"] += 1
        text, err = _collect_text(client, current_prompt, sys_prompt,
                                  max_tokens)
        result["raw"] = text
        if err:
            # Transport-level failure: retrying the same broken client is
            # pointless; report and stop.
            result["error"] = err
            return result
        obj = extract_json_object(text)
        errs = (validate_json_schema(obj, schema)
                if obj is not None
                else ["no JSON object found in the reply"])
        if not errs:
            result["data"] = obj
            result["error"] = ""
            return result
        last_errs = errs
        current_prompt = (
            base_prompt
            + "\n\nYour previous reply did not satisfy the required output "
            + "schema: " + "; ".join(errs)[:600]
            + "\nReply again with ONLY one JSON object matching the schema "
            + "from your instructions — no prose, no explanation."
        )
    result["error"] = "; ".join(last_errs)[:600] or "no valid reply"
    return result


__all__ = [
    "validate_json_schema",
    "extract_json_object",
    "schema_instruction",
    "request_structured",
]
