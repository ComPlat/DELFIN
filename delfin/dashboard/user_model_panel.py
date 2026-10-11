"""The dashboard's "User model" panel (package L4, phase 2).

The typed user picture (role, expertise, preferences, working style,
corrections) is visible, editable, and offers "forget this entry" and
"forget everything". Everything here is a pure function over the
``delfin.agent.user_model.UserModel`` store: a click saves for real only
because these actions call the store's own ``forget`` / ``forget_all`` /
``update`` methods, and nothing is remembered here that the store refuses.

Deliberately no ipywidgets import at module scope, so the logic is testable
headless. ``delfin.agent.user_model`` is imported lazily inside each call --
a tree that lacks the module shows "could not be read" rather than breaking
the Agent tab.
"""

from __future__ import annotations

from html import escape
from typing import Any, Optional


def _user_model(path: Optional[str] = None):
    """Lazy handle to the user-model store (import errors surface in text)."""
    from delfin.agent import user_model
    return user_model.UserModel(path=path) if path else user_model.UserModel()


def collect_rows(path: Optional[str] = None) -> list[dict[str, Any]]:
    """One row per stored entry: field, value, source role, source msg, date.

    ``path`` lets a caller point the panel at an isolated store (tests);
    the default is the user's own ``~/.delfin/user_model.json``. If the store
    cannot be read, returns a single row carrying the error so the panel is
    never blank-silent.
    """
    try:
        model = _user_model(path)
        entries = model.entries()
    except Exception as exc:  # pragma: no cover - read-back failure path
        return [{
            "field": "error", "value": f"could not read user model: {exc}",
            "source_role": "", "source_message": "", "date": "",
        }]
    rows: list[dict[str, Any]] = []
    for field, bucket in entries.items():
        for index, entry in enumerate(bucket):
            rows.append({
                "field": field,
                "index": index,
                "value": str(entry.get("value", "")),
                "source_role": str(entry.get("source_role", "")),
                "source_message": str(entry.get("source_message", "")),
                "date": str(entry.get("date", "")),
            })
    return rows


def render(path: Optional[str] = None) -> str:
    """Text rendering of the whole model for display (HTML-escaped)."""
    rows = collect_rows(path)
    if not rows or (len(rows) == 1 and rows[0]["field"] == "error"):
        return rows[0]["value"] if rows else "(no user model yet)"
    lines: list[str] = ["--- user model ---"]
    for row in rows:
        lines.append(
            f"- {escape(row['field'])}: {escape(row['value'])}  "
            f"[source={escape(row['source_role'])}, {escape(row['date'])}]"
        )
    return "\n".join(lines)


def forget_entry(field: str, index: int, path: Optional[str] = None) -> dict:
    """"Forget this": remove entry ``index`` of ``field``.

    Returns ``{"removed": True}`` on success, ``{"removed": False, "error":
    ...}`` otherwise. Deleting is the store's own ``forget``, so it also
    removes the entry from what the next session loads.
    """
    try:
        model = _user_model(path)
        ok = model.forget(field, index)
    except Exception as exc:  # pragma: no cover - defensive
        return {"removed": False, "error": str(exc)}
    if not ok:
        return {"removed": False, "error": (
            f"no entry {index} in field '{field}'")}
    return {"removed": True}


def forget_everything(path: Optional[str] = None) -> dict:
    """"Forget everything": clear every field.

    Returns ``{"removed": N}`` where N is how many entries were cleared.
    """
    try:
        model = _user_model(path)
        removed = model.forget_all()
    except Exception as exc:  # pragma: no cover - defensive
        return {"removed": 0, "error": str(exc)}
    return {"removed": removed}
