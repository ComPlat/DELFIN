"""Repair helpers for ACTION-protocol tool calls.

Some roles operate the host UI through a text protocol: they emit
``ACTION: /command`` lines in plain assistant text, and the host parses and
executes those lines (with its own safety tiers and confirmation gates).
Models occasionally emit that protocol as a structured TOOL CALL instead —
a tool named ``ACTION`` (optionally ``mcp__*__``-namespaced), or the slash
command itself used as the tool name. Dispatching such a call can only
produce an unknown-tool or role-denial error, and a model that re-issues
the same call trips the consecutive-identical-error abort, turning a
recoverable protocol confusion into a hard turn failure.

The helpers here are pure and side-effect free:

- :func:`is_action_style_call` detects an ACTION-style tool call.
- :func:`extract_slash_command` recovers the intended slash command from
  the call's name and/or arguments (liberal in what it accepts).
- :func:`canonical_action_line` renders the canonical text-protocol line.
- :func:`build_repair_result` builds a structured, NON-error tool result
  that states the repaired interpretation and steers the model back to the
  plain-text protocol. The result never starts with ``{"error"`` — the
  identical-error abort machinery keys on that prefix, so a repaired call
  can never count as a failure round — and it varies with the occurrence
  counter, so repeats are never byte-identical either.

The caller decides whether to additionally register the repaired command on
the text channel (downstream ACTION parsers read accumulated assistant
text), and passes ``registered=True`` so the result says so. No slash
command list is hardcoded here: only the command SHAPE is validated, and
the downstream dispatcher remains the authority on which commands execute.
"""

from __future__ import annotations

import json
import re

__all__ = [
    "ACTION_PROTOCOL_ROLES",
    "role_uses_action_protocol",
    "strip_tool_namespace",
    "is_action_style_call",
    "extract_slash_command",
    "canonical_action_line",
    "build_repair_result",
]

# Roles that operate the ACTION text protocol. Mirrors the role roster in
# the execution allow-list (api_client) — extend both when a new
# text-protocol role is introduced.
ACTION_PROTOCOL_ROLES: frozenset[str] = frozenset({"dashboard_agent"})

_ACTION_KEYWORD = "ACTION"

# Slash-command shape: "/" + a bare first token (letters, digits, "_",
# "-"), then either end-of-string or whitespace before arguments.
# Deliberately rejects filesystem-path lookalikes ("/home/user/x" has no
# whitespace after a bare first token) and bare "/".
_SLASH_SHAPE_RE = re.compile(r"^/[A-Za-z][A-Za-z0-9_-]*(?:\s|$)")

# A bare first word that can be promoted to a slash command by prefixing
# "/" (used only for imperative argument keys, never for prose carriers).
_BARE_WORD_RE = re.compile(r"^[A-Za-z][A-Za-z0-9_-]*(?:\s|$)")

# Leading protocol keyword to strip from candidate strings ("ACTION: /x",
# "action /x", "Action./x") — mirrors the tolerant downstream line parser.
_ACTION_PREFIX_RE = re.compile(r"^\s*ACTION\s*[:.]?\s*", re.IGNORECASE)

# Argument keys likely to carry the slash command, in scan order.
_COMMAND_KEYS: tuple[str, ...] = (
    "command", "cmd", "action", "slash_command", "input", "line",
    "text", "value", "args", "arguments", "content", "query",
)

# Keys whose value is imperative enough to promote a bare "tab calc" into
# "/tab calc". Kept narrower than _COMMAND_KEYS so prose-carrying keys
# (query/content/...) never fabricate a command.
_IMPERATIVE_KEYS: frozenset[str] = frozenset(
    {"command", "cmd", "action", "slash_command", "input", "line"})

# Keys whose value may be appended as the argument tail when the tool NAME
# is already the slash command ("/orca" + {"command": "set ..."}).
_TAIL_KEYS: frozenset[str] = frozenset(
    {"command", "cmd", "args", "arguments", "input", "text", "value",
     "line"})


def role_uses_action_protocol(role: str) -> bool:
    """True when *role* drives the UI via ACTION text lines."""
    return (role or "") in ACTION_PROTOCOL_ROLES


def strip_tool_namespace(name: str) -> str:
    """Strip an ``mcp__server__`` namespace prefix from a tool name."""
    n = (name or "").strip()
    if n.startswith("mcp__"):
        return n.rsplit("__", 1)[-1]
    return n


def _clean(text: str) -> str:
    """Normalize a candidate string: fences, quotes, ACTION keyword."""
    s = (text or "").strip()
    s = s.strip("`").strip()
    if len(s) >= 2 and s[0] == s[-1] and s[0] in ("'", '"'):
        s = s[1:-1].strip()
    s = _ACTION_PREFIX_RE.sub("", s, count=1).strip()
    return s


def _is_slash_shaped(text: str) -> bool:
    return bool(_SLASH_SHAPE_RE.match(text or ""))


def _base_name(name: str) -> str:
    """De-namespaced tool name with a trailing colon stripped."""
    return strip_tool_namespace(name).rstrip(":").strip()


def is_action_style_call(name: str, arguments: object = None) -> bool:
    """True when the tool call is really a text-protocol invocation.

    Matches the tool name ``ACTION`` (case-insensitive, with or without an
    ``mcp__*__`` namespace prefix or a trailing colon) and tool names that
    are themselves slash-command shaped (``/tab``, ``/orca set ...``,
    including an ``ACTION:``-fused form). ``arguments`` is accepted for
    signature stability; detection is name-based.
    """
    del arguments  # detection is name-based
    base = _base_name(name)
    if base.upper() == _ACTION_KEYWORD:
        return True
    return _is_slash_shaped(_clean(base))


def _coerce_arguments(arguments: object) -> tuple[dict, str]:
    """Return a (dict view, raw-string view) of an arguments payload.

    Accepts a dict, a JSON-encoded object/string, or plain text. Exactly
    one of the two views is populated.
    """
    if isinstance(arguments, dict):
        return arguments, ""
    if isinstance(arguments, str):
        try:
            parsed = json.loads(arguments)
        except (json.JSONDecodeError, ValueError):
            return {}, arguments
        if isinstance(parsed, dict):
            return parsed, ""
        if isinstance(parsed, str):
            return {}, parsed
        return {}, arguments
    return {}, ""


def _string_items(arguments: object) -> list[tuple[str, str]]:
    """(key, cleaned string) pairs from the arguments payload.

    Known command keys come first (in ``_COMMAND_KEYS`` order), then the
    remaining string values in insertion order. Dict values are descended
    one level (``{"input": {"command": "/x"}}`` yields ``("command", "/x")``).
    """
    args, raw = _coerce_arguments(arguments)
    items: list[tuple[str, str]] = []

    def _add(key: str, value: object) -> None:
        if isinstance(value, str):
            cleaned = _clean(value)
            if cleaned:
                items.append((key, cleaned))
        elif isinstance(value, dict):
            for inner_key, inner_value in value.items():
                if isinstance(inner_value, str):
                    cleaned = _clean(inner_value)
                    if cleaned:
                        items.append((str(inner_key), cleaned))

    for key in _COMMAND_KEYS:
        if key in args:
            _add(key, args[key])
    for key, value in args.items():
        if key not in _COMMAND_KEYS:
            _add(str(key), value)
    if raw:
        cleaned = _clean(raw)
        if cleaned:
            items.append(("", cleaned))
    return items


def extract_slash_command(name: str, arguments: object = None) -> str:
    """Best-effort recovery of the intended slash command.

    Sources, in priority order:
    1. A slash-shaped string in the arguments (known command keys first,
       then any string value; one nesting level; JSON-string payloads and
       plain-text payloads accepted; ``ACTION:`` prefixes stripped).
    2. The tool name itself when slash-shaped — argument values under
       tail-capable keys are appended as the command's argument tail.
    3. A bare-word value under an imperative key, promoted with "/".

    Returns "" when nothing command-shaped can be recovered. The result is
    shape-validated only — the downstream dispatcher decides whether the
    command exists and whether it may run.
    """
    base = _clean(_base_name(name))
    items = _string_items(arguments)
    slash_values = [value for _key, value in items if _is_slash_shaped(value)]

    if _is_slash_shaped(base):
        if slash_values:
            # Prefer the most complete command string available.
            best = max(slash_values, key=len)
            return best if len(best) >= len(base) else base
        for key, value in items:
            if key in _TAIL_KEYS and not value.startswith("/"):
                return f"{base} {value}"
        return base

    if slash_values:
        return slash_values[0]
    for key, value in items:
        if key in _IMPERATIVE_KEYS and _BARE_WORD_RE.match(value):
            return "/" + value
    return ""


def canonical_action_line(command: str) -> str:
    """Render the canonical single-line form: ``ACTION: /command ...``.

    Real newlines are encoded as literal ``\\n`` — the downstream line
    parser decodes that form, so multiline payloads survive the
    single-line protocol.
    """
    cmd = (command or "").strip()
    cmd = cmd.replace("\r\n", "\n").replace("\r", "\n").replace("\n", "\\n")
    return f"{_ACTION_KEYWORD}: {cmd}"


def build_repair_result(
    command: str,
    *,
    role: str = "",
    attempt: int = 1,
    registered: bool = False,
    tool_name: str = "",
) -> str:
    """Structured tool result for a repaired ACTION-style call.

    Returns a JSON string whose FIRST key is ``"repaired"`` — never
    ``"error"`` — so the identical-error abort machinery (which keys on
    the ``{"error"`` prefix) can never classify a repaired round as a
    failure round. The ``attempt`` counter is embedded in both a field and
    the guidance text, so repeat occurrences are also never byte-identical.

    ``registered=True`` means the caller injected the canonical ACTION
    line on the text channel; the guidance then tells the model the
    command is already on its way and must not be re-issued. Without a
    recoverable command the result explains the expected plain-text form.
    """
    attempt = max(1, int(attempt or 1))
    line = canonical_action_line(command) if command else ""
    parts: list[str]
    if command:
        head = f"This tool call was understood as the text-protocol line `{line}`"
        if tool_name:
            head += f" (repaired from tool call '{tool_name}')"
        parts = [head + "."]
        if registered and attempt == 1:
            parts.append(
                "The command has been registered for execution as if that "
                "line had appeared in your reply — do not re-issue it.")
        elif registered:
            parts.append(
                f"Occurrence {attempt} of the same call — the command was "
                "already registered on the first occurrence; do not repeat "
                "it.")
        parts.append(
            "ACTION commands are plain text, not tool calls: write the "
            f"line `{line}` directly in your reply. Reserve tool calls for "
            "real tools such as search_docs.")
    else:
        parts = [
            "This tool call looks like an ACTION-protocol invocation, but "
            "no slash command could be recovered from its name or "
            "arguments.",
            "Write the command as a plain text line in your reply, e.g. "
            "`ACTION: /tab calc` — do not wrap it in a tool call.",
        ]
        if attempt > 1:
            parts.append(
                f"Occurrence {attempt} of this call shape — change approach "
                "now and emit the ACTION line as plain text.")
    payload: dict[str, object] = {
        "repaired": bool(command),
        "understood_as": line or None,
        "registered": bool(registered and command),
        "attempt": attempt,
        "guidance": " ".join(parts),
    }
    if role:
        payload["role"] = role
    return json.dumps(payload)
