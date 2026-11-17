"""Help lookup: the single entry point for "what is X" questions.

Reads the two sources of truth -- ``CONTROL_FIELD_SPECS`` (keys) and
``_TABS_HELP`` (dashboard tabs) -- and answers name/help/allowed queries.
No second copy of any explanation lives here: this module only reads.

Used by the agent's operating rule (knowledge/operating/help-lookup) and
available for anything else that wants grounded answers about DELFIN's
surface.
"""

from __future__ import annotations

from typing import Any


def _key_entries() -> list[dict[str, str]]:
    from delfin.common.control_validator import CONTROL_FIELD_SPECS

    entries = []
    for spec in CONTROL_FIELD_SPECS:
        entries.append({
            "kind": "key",
            "name": spec.name,
            "help": (spec.help or "").strip(),
            "allowed": (getattr(spec, "allowed", "") or "").strip(),
            "default": "" if spec.default is None else str(spec.default),
            "required": bool(spec.required),
        })
    return entries


def _tab_entries() -> list[dict[str, str]]:
    from delfin.dashboard.tab_help import _TABS_HELP

    entries = []
    for title, description in _TABS_HELP:
        entries.append({
            "kind": "tab",
            "name": title,
            "help": description,
            "allowed": "",
            "default": "",
            "required": False,
        })
    return entries


def _score(entry: dict[str, str], q: str) -> int:
    """Match strength: exact name > name contains > text contains."""
    name = entry["name"].lower()
    if q == name:
        return 3
    if q in name:
        return 2
    if q in entry["help"].lower() or (entry["allowed"] and q in entry["allowed"].lower()):
        return 1
    return 0


def lookup(query: str, *, kinds: tuple[str, ...] = ("key", "tab"),
           limit: int = 8) -> dict[str, Any]:
    """Look up a CONTROL key or dashboard tab by name or content.

    Returns ``{"found": bool, "entries": [...], "query": query}``.  Each
    entry is ``{kind, name, help, allowed, default, required}``; a miss is
    reported explicitly with ``found=False`` and an empty list -- callers
    must say so rather than guess.
    """
    q = (query or "").strip().lower()
    if not q:
        return {"found": False, "entries": [], "query": query or ""}
    scored = []
    for entry in _key_entries() + _tab_entries():
        if entry["kind"] not in kinds:
            continue
        s = _score(entry, q)
        if s:
            scored.append((s, entry))
    scored.sort(key=lambda p: (-p[0], p[1]["name"].lower()))
    entries = [e for _s, e in scored[:limit]]
    return {"found": bool(entries), "entries": entries, "query": query}


def format_answer(result: dict[str, Any]) -> str:
    """Render a lookup result as the grounded answer text to give a user.

    Explicitly states "no entry" on a miss -- that is the signal for the
    caller to derive the answer from the consuming code instead and
    propose a spec (the gap loop in the operating rule).
    """
    if not result["found"]:
        return (f"No entry for '{result['query']}' in CONTROL_FIELD_SPECS "
                "or the tab guide. Deriving from the consuming code is "
                "allowed; state the file:line and propose a FieldSpec or "
                "_TABS_HELP entry for the gap.")
    lines = []
    for e in result["entries"]:
        if e["kind"] == "key":
            head = (f"CONTROL key `{e['name']}`"
                    + (" (required)" if e["required"] else "")
                    + f" -- default: {e['default'] or '(empty)'}")
        else:
            head = f"Dashboard tab **{e['name']}**"
        body = e["help"] or "(no description)"
        tail = f"  Allowed: {e['allowed']}" if e["allowed"] else ""
        lines.append(f"{head}\n  {body}{tail}")
    return "\n".join(lines)
