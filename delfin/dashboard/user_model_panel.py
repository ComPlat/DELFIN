"""The dashboard's "User model" panel.

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


def forget_entry(field: str, index: int, path: Optional[str] = None,
                 value: Optional[str] = None) -> dict:
    """"Forget this": remove entry ``index`` of ``field``.

    With ``value`` (what the panel showed), the entry is found by that
    value instead: a line typed in the chat since the panel was drawn
    moves entries (said again -> newest place; past twelve -> oldest
    dropped), and the shown index then names a different entry. An
    entry no longer present is reported, not substituted.

    Returns ``{"removed": True}`` on success, ``{"removed": False, "error":
    ...}`` otherwise. Deleting is the store's own ``forget``, so it also
    removes the entry from what the next session loads.
    """
    try:
        model = _user_model(path)
        if value is not None:
            want = str(value).casefold()
            index = next((i for i, e in enumerate(
                model.entries().get(field, []))
                if str(e.get("value", "")).casefold() == want), -1)
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


def add_entry(text: str, path: Optional[str] = None) -> dict:
    """"Add": one line the user types into the panel about themselves.

    It goes through the store's own intake, so the panel can add nothing
    a typed chat message could not: one short line, a recognised field,
    nothing permission-shaped.
    """
    try:
        model = _user_model(path)
        before = model.entries()
        added = model.update_from_user_message(str(text or "")) != before
    except Exception as exc:  # pragma: no cover - defensive
        return {"added": False, "error": str(exc)}
    if not added:
        return {"added": False, "error": (
            "not stored: one short line about yourself (role, expertise, "
            "preference, working style, standing correction), and nothing "
            "that reads like a permission, a path or a mode")}
    return {"added": True}


def build_widget(path: Optional[str] = None):
    """The panel as a widget: the entries, add, forget one, forget all."""
    import ipywidgets as widgets

    view = widgets.HTML()
    status = widgets.HTML()
    pick = widgets.Dropdown(
        options=[], description="Entry:",
        layout=widgets.Layout(width="520px"))
    new_line = widgets.Text(
        placeholder="e.g. I am a computational chemist",
        layout=widgets.Layout(width="520px"))
    add_btn = widgets.Button(description="Add")
    forget_btn = widgets.Button(description="Forget this")
    all_btn = widgets.Button(
        description="Forget everything", button_style="warning")
    refresh_btn = widgets.Button(description="Refresh")

    def _refresh(*_):
        rows = [r for r in collect_rows(path) if r["field"] != "error"]
        view.value = (
            "<pre style='white-space:pre-wrap;margin:0'>"
            + render(path) + "</pre>")
        pick.options = [
            (f"{r['field']}: {r['value'][:60]}",
             (r["field"], r["index"], r["value"]))
            for r in rows]

    def _say(result: dict, ok: str):
        err = result.get("error")
        status.value = (
            f"<span style='color:#d32f2f;'>{escape(str(err))}</span>"
            if err else f"<span style='color:#2e7d32;'>{escape(ok)}</span>")

    def _add(_):
        result = add_entry(new_line.value, path)
        if result.get("added"):
            new_line.value = ""
        _say(result, "Stored.")
        _refresh()

    def _forget(_):
        if pick.value is None:
            return
        field, index, value = pick.value
        _say(forget_entry(field, index, path, value=value), "Forgotten.")
        _refresh()

    def _forget_all(_):
        result = forget_everything(path)
        _say(result, f"Forgot {result.get('removed', 0)} entries.")
        _refresh()

    add_btn.on_click(_add)
    forget_btn.on_click(_forget)
    all_btn.on_click(_forget_all)
    refresh_btn.on_click(_refresh)
    _refresh()
    return widgets.VBox([
        widgets.HTML(
            "<b>User model</b> &mdash; what you said about yourself, "
            "loaded into every new session. It is a description: it "
            "grants the agent nothing. Stored only in your own "
            "<code>~/.delfin/user_model.json</code>."),
        view,
        widgets.HBox([new_line, add_btn]),
        widgets.HBox([pick, forget_btn, all_btn, refresh_btn]),
        status,
    ])
