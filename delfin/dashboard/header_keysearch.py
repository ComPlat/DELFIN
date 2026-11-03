"""Global CONTROL-key search in the dashboard header.

A small search box left of the HOME-usage chip: type part of a key name or
help text, a dropdown of suggestions appears, picking one shows an
explanation card (what it is, default, allowed values) directly under the
bar.  Data comes from ``CONTROL_FIELD_SPECS`` -- the same table the
validator and the Help tab read, so it cannot drift from the code.
"""

from __future__ import annotations

from typing import Any

import ipywidgets as W


def _specs() -> list[dict[str, str]]:
    """One entry per CONTROL key: name, help, allowed, default, required."""
    from delfin.common.control_validator import CONTROL_FIELD_SPECS

    entries = []
    for spec in CONTROL_FIELD_SPECS:
        entries.append({
            "name": spec.name,
            "help": (spec.help or "").strip(),
            "allowed": (getattr(spec, "allowed", "") or "").strip(),
            "default": "" if spec.default is None else str(spec.default),
            "required": "required" if spec.required else "",
        })
    return sorted(entries, key=lambda e: e["name"].lower())


def build_key_search(max_suggestions: int = 8) -> W.DOMWidget:
    """Build the header search box with suggestion dropdown + explanation."""
    entries = _specs()
    names = [e["name"] for e in entries]

    search = W.Text(
        placeholder="🔍 search keys…",
        layout=W.Layout(width="180px"),
    )
    dropdown = W.VBox(layout=W.Layout(
        display="none", position="absolute", z_index="1000",
        max_height="300px", overflow_y="auto", background="white",
        border="1px solid #ccc", border_radius="4px", padding="2px",
        box_shadow="0 2px 8px rgba(0,0,0,.15)",
    ))
    card = W.HTML("", layout=W.Layout(display="none"))
    wrapper = W.VBox(
        [W.HBox([search, dropdown]), card],
        layout=W.Layout(position="relative"),
    )
    wrapper.add_class("delfin-keysearch")

    def _close(_b=None) -> None:
        dropdown.layout.display = "none"

    def _pick(name: str) -> None:
        """Show the explanation card for a key."""
        entry = entries[names.index(name)]
        allowed = entry["allowed"] or "— (free-form; see the help text)"
        req = entry["required"]
        card.value = (
            "<div class='delfin-keycard'>"
            f"<div><code style='font-weight:700'>{name}</code>"
            f"{' — ' + req if req else ''}"
            f" &nbsp; default: <code>{entry['default'] or '—'}</code></div>"
            f"<div style='margin-top:2px'>{entry['help']}</div>"
            "<div style='margin-top:2px; color:#555'>"
            f"<b>Allowed:</b> {allowed}</div>"
            "<div style='margin-top:2px; color:#888; font-size:12px'>"
            "In CONTROL.txt write <code>KEY=?</code> for this text; "
            "full reference in the <b>Help</b> tab.</div>"
            "</div>"
        )
        card.layout.display = ""
        _close()
        search.value = ""

    def _refresh() -> None:
        q = (search.value or "").strip().lower()
        if not q:
            _close()
            return
        scored = [
            e for e in entries
            if q in e["name"].lower()
            or q in e["help"].lower()
            or q in e["allowed"].lower()
        ][:max_suggestions]
        if not scored:
            dropdown.children = (
                W.HTML("<div style='padding:4px 8px; color:#999'>no key "
                       "matches</div>"),
            )
        else:
            buttons = []
            for e in scored:
                snippet = e["help"][:80] + ("…" if len(e["help"]) > 80 else "")
                btn = W.Button(
                    description=e["name"], button_style="",
                    tooltip=snippet,
                    layout=W.Layout(
                        width="100%", border="none", text_align="left",
                        padding="2px 8px", background="none",
                    ),
                )

                def _on(b, _name=e["name"]):
                    _pick(_name)

                btn.on_click(_on)
                buttons.append(btn)
            dropdown.children = tuple(buttons)
        dropdown.layout.display = ""

    search.observe(lambda ch: _refresh(), names="value")
    # No blur-close: ipywidgets Text has no reliable blur trait, and a timer
    # hack races the suggestion click.  The dropdown closes when a key is
    # picked (or the query is cleared), which is the only way forward anyway.

    return wrapper


def build_and_register() -> W.DOMWidget | None:
    """Build the widget; returns None where specs are unavailable."""
    try:
        return build_key_search()
    except Exception:
        return None
