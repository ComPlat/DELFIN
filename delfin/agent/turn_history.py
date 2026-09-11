"""What a turn did, kept for the turns after it.

The OpenAI-style tool loop builds its rounds -- the assistant row that
carries ``tool_calls``, the ``tool`` rows that answer them, the user's
mid-run steers -- in a list local to one request. The engine's history
received only the final text. Every later turn, and every correction
turn a guard forces, was therefore answered by a model that could not
see what it had read or run: told two cited paths did not exist, GLM
wrote "ich habe die Dateien nicht gelesen" after seven tool calls;
asked to say its English answer in German, DeepSeek added "eigentlich
hätte ich die Dateien zunächst selbst lesen sollen" after reading them
(both 2026-09-11). Keeping the rounds is how the assistant this project
measures itself against works.

Two functions. ``compact_turn_rows`` turns the loop's rows into history
rows: results and arguments cut to a head, incomplete rounds dropped
(an assistant row whose calls have no answers makes the next request
unsendable), only the last rounds of a long turn kept, with a note about
the rest. ``repair_tool_pairing`` keeps a history sendable after
compaction has cut through a round.
"""
from __future__ import annotations

import json
from typing import Any

_API_KEYS = ("role", "content", "tool_calls", "tool_call_id", "name")


def _as_text(value: Any) -> str:
    if value is None:
        return ""
    if isinstance(value, str):
        return value
    try:
        return json.dumps(value, ensure_ascii=False)
    except Exception:
        return str(value)


def _cut(text: str, limit: int, what: str) -> str:
    if len(text) <= limit:
        return text
    return (text[:limit]
            + f"\n… [{len(text) - limit} more chars of this {what} are not "
              "kept in the history]")


def _api_row(row: dict) -> dict:
    out = {k: row[k] for k in _API_KEYS if k in row and row[k] is not None}
    out["content"] = _as_text(out.get("content"))
    return out


def _call_ids(row: dict) -> list[str]:
    calls = row.get("tool_calls") or []
    ids = []
    for call in calls:
        try:
            ids.append(str(call.get("id") or ""))
        except Exception:
            ids.append("")
    return ids


def _cut_arguments(row: dict, args_chars: int) -> dict:
    calls = []
    for call in row.get("tool_calls") or []:
        try:
            call = dict(call)
            fn = dict(call.get("function") or {})
            args = _as_text(fn.get("arguments"))
            if len(args) > args_chars:
                # Still JSON: a cut string would not be, and an endpoint
                # that parses old arguments would refuse the request.
                fn["arguments"] = json.dumps({
                    "_not_kept_in_history": f"{len(args)} chars of arguments"})
            call["function"] = fn
            calls.append(call)
        except Exception:
            calls.append(call)
    row["tool_calls"] = calls
    return row


def _groups(rows: list[dict]) -> list[list[dict]]:
    """Split rows into units: a tool round (assistant + its tool rows) or
    a single other row. Tool rows that answer no round are dropped."""
    units: list[list[dict]] = []
    i = 0
    while i < len(rows):
        row = rows[i]
        if row.get("role") == "assistant" and row.get("tool_calls"):
            unit = [row]
            j = i + 1
            while j < len(rows) and rows[j].get("role") == "tool":
                unit.append(rows[j])
                j += 1
            units.append(unit)
            i = j
        elif row.get("role") == "tool":
            i += 1                      # an answer to nothing
        else:
            units.append([row])
            i += 1
    return units


def _round_is_complete(unit: list[dict]) -> bool:
    if not (unit and unit[0].get("tool_calls")):
        return True
    wanted = set(_call_ids(unit[0]))
    have = {str(r.get("tool_call_id") or "") for r in unit[1:]}
    return bool(wanted) and wanted <= have


def compact_turn_rows(rows: list[dict], *, keep_rounds: int = 16,
                      result_chars: int = 1200, args_chars: int = 2000,
                      text_chars: int = 4000) -> list[dict]:
    """History rows for the rows a tool loop added in one turn."""
    if not rows:
        return []
    cleaned: list[dict] = []
    for raw in rows:
        if not isinstance(raw, dict):
            continue
        row = _api_row(raw)
        role = row.get("role")
        if role == "tool":
            row["content"] = _cut(row["content"], result_chars, "result")
        elif role == "assistant" and row.get("tool_calls"):
            row = _cut_arguments(row, args_chars)
            row["content"] = _cut(row["content"], text_chars, "text")
        elif role in ("assistant", "user"):
            row["content"] = _cut(row["content"], text_chars, "text")
        else:
            continue                    # system/developer rows never belong here
        cleaned.append(row)
    units = [u for u in _groups(cleaned) if _round_is_complete(u)]
    rounds = [u for u in units if u[0].get("tool_calls")]
    dropped = 0
    if len(rounds) > keep_rounds:
        drop = set(id(u) for u in rounds[:len(rounds) - keep_rounds])
        dropped = len(drop)
        names: dict[str, int] = {}
        for u in rounds[:len(rounds) - keep_rounds]:
            for call in u[0].get("tool_calls") or []:
                try:
                    n = str((call.get("function") or {}).get("name") or "?")
                except Exception:
                    n = "?"
                names[n] = names.get(n, 0) + 1
        units = [u for u in units if id(u) not in drop]
        first = next((u for u in units if u[0].get("tool_calls")), None)
        if first is not None:
            summary = ", ".join(f"{k} ×{v}" if v > 1 else k
                                for k, v in sorted(names.items()))
            note = (f"[{dropped} earlier tool round(s) of this turn are not "
                    f"kept in the history: {summary}]")
            first[0]["content"] = (note + ("\n" + first[0]["content"]
                                           if first[0]["content"] else ""))
    return [row for u in units for row in u]


def repair_tool_pairing(messages: list[dict]) -> list[dict]:
    """The same history without rounds a cut has broken.

    An assistant row whose ``tool_calls`` are not all answered by the
    ``tool`` rows right after it, and a ``tool`` row that answers no
    such row, make the request unsendable. Compaction keeps the last N
    rows and summarises the rest, and N can fall in the middle of a
    round. Everything else is returned as it was, same objects.
    """
    if not messages:
        return list(messages or [])
    out: list[dict] = []
    for unit in _groups([m for m in messages if isinstance(m, dict)]):
        if _round_is_complete(unit):
            out.extend(unit)
    return out


def strip_tool_rows(messages: list[dict]) -> list[dict]:
    """The history without its tool rounds, for a backend that cannot send them.

    A session is saved with the rounds an OpenAI-shaped client kept and may
    be restored into another client. That client's conversion knows only
    user and assistant text; a ``tool`` row or a ``tool_calls`` key would
    reach the endpoint as-is and be refused. Assistant rows that carried
    only calls are dropped with them; the answer they led to stays.
    """
    out: list[dict] = []
    for m in messages or []:
        if not isinstance(m, dict):
            continue
        if m.get("role") == "tool":
            continue
        if m.get("tool_calls"):
            text = m.get("content")
            if not (isinstance(text, str) and text.strip()):
                continue
            m = {k: v for k, v in m.items() if k != "tool_calls"}
        out.append(m)
    return out
