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


def _tab_entries() -> list[dict[str, str]]:
    """One entry per dashboard tab, read from _TABS_HELP (single source)."""
    from delfin.dashboard.tab_help import _TABS_HELP

    entries = []
    for item in _TABS_HELP:
        name = item[0]
        purpose = (item[1] or "").strip()
        entries.append({
            "name": name,
            "purpose": purpose,
        })
    return entries


def build_key_search(
    max_suggestions: int = 8,
    on_open_key=None,
    on_open_tab=None,
) -> W.DOMWidget:
    """Build the header search box with suggestion dropdown + explanation.

    ``on_open_key(name)`` is optional: when given, picking a key suggestion
    also switches the dashboard to the Help tab pre-filtered to that key
    (in addition to the inline card).  ``on_open_tab(name)`` is the same for
    tab-name suggestions: it switches to that tab itself, not to Help.  The
    dashboard passes both when the Help tab is registered; without them
    picking just shows the card.
    """
    entries = _specs()
    names = [e["name"] for e in entries]
    tab_entries = _tab_entries()
    tab_names = [e["name"] for e in tab_entries]

    search = W.Text(
        placeholder="🔍 search keys…",
        layout=W.Layout(width="180px"),
    )
    # Dropdown and card are fixed overlays via CSS CLASS, not Layout: this
    # ipywidgets (8.1.8) silently drops Layout(position=...) -- there is no
    # position trait -- which is why earlier Layout-based fixes never took
    # effect and the page kept stretching.  The CSS block below sets
    # position:fixed on the classes; here only visibility is controlled.
    _overlay_css = W.HTML(
        "<style>"
        ".delfin-keysearch { position: relative; }"
        ".delfin-keysearch-dropdown { position: absolute; top: 100%; right: 0;"
        " width: 460px; max-height: 300px; overflow-y: auto; z-index: 10000;"
        " background: var(--jp-layout-color1, #fff);"
        " border: 1px solid #ccc; border-radius: 4px; padding: 2px;"
        " box-shadow: 0 2px 8px rgba(0,0,0,.15); }"
        ".delfin-keysearch-card { position: absolute; top: 100%; right: 0;"
        " width: 500px; z-index: 10001; margin-top: 4px;"
        " background: var(--jp-layout-color1, #fff);"
        " border: 1px solid #bbb; border-radius: 6px; padding: 10px 12px;"
        " box-shadow: 0 4px 16px rgba(0,0,0,.2); }"
        "</style>"
    )
    dropdown = W.VBox(layout=W.Layout(
        display="none", max_height="300px", overflow_y="auto", width="460px",
    ))
    card = W.VBox(
        layout=W.Layout(display="none", width="500px"),
    )
    wrapper = W.VBox(
        [_overlay_css, search, dropdown, card],
        layout=W.Layout(width="180px"),
    )
    wrapper.add_class("delfin-keysearch")
    dropdown.add_class("delfin-keysearch-dropdown")
    card.add_class("delfin-keysearch-card")

    close_btn = W.Button(
        description="\u2715", button_style="",
        tooltip="Hide this explanation",
        layout=W.Layout(width="26px", height="22px", padding="0",
                        margin="0 0 4px 0"),
    )

    def _close(_b=None) -> None:
        dropdown.layout.display = "none"

    def _pick(name: str) -> None:
        """Show the explanation card for a key."""
        entry = entries[names.index(name)]
        allowed = entry["allowed"] or "— (free-form; see the help text)"
        req = entry["required"]
        card.children = (close_btn, W.HTML(
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
        ))
        card.layout.display = ""
        _close()
        if on_open_key is not None:
            # Route to the Help tab pre-filtered to this key as well; the
            # inline card stays open so the explanation is visible either way.
            try:
                on_open_key(name)
            except Exception:
                pass
        search.value = ""

    def _pick_tab(name: str) -> None:
        """Switch to the named dashboard tab itself (not Help)."""
        if on_open_tab is not None:
            try:
                on_open_tab(name)
            except Exception:
                pass
        _close()
        search.value = ""

    def _close_card(_b=None) -> None:
        card.layout.display = "none"

    close_btn.on_click(_close_card)

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
        ]
        # Tab-name matches come first (they route directly); a match on a
        # key's help text or name counts like before.  Tabs are visually
        # distinct: a 📁 prefix and a "Tab" tag.
        tab_hits = [
            e for e in tab_entries
            if q in e["name"].lower()
            or q in e["purpose"].lower()
        ]
        scored = (tab_hits + scored)[:max_suggestions]
        if not scored:
            dropdown.children = (
                W.HTML("<div style='padding:4px 8px; color:#999'>no key "
                       "or tab matches</div>"),
            )
        else:
            buttons = []
            for e in scored:
                if e.get("purpose") is not None and e in tab_entries:
                    is_tab = True
                else:
                    is_tab = False
                if is_tab:
                    snippet = e["purpose"][:80] + (
                        "…" if len(e["purpose"]) > 80 else "")
                    label = "📁 Tab: " + e["name"]
                    def _on(b, _name=e["name"]):
                        _pick_tab(_name)
                else:
                    snippet = e["help"][:80] + ("…" if len(e["help"]) > 80 else "")
                    label = e["name"]
                    def _on(b, _name=e["name"]):
                        _pick(_name)
                btn = W.Button(
                    description=label, button_style="",
                    tooltip=snippet,
                    layout=W.Layout(
                        width="100%", border="none", text_align="left",
                        padding="2px 8px", background="none",
                    ),
                )
                btn.on_click(_on)
                buttons.append(btn)
            dropdown.children = tuple(buttons)
        dropdown.layout.display = ""

    search.observe(lambda ch: _refresh(), names="value")
    # No blur-close: ipywidgets Text has no reliable blur trait, and a timer
    # hack races the suggestion click.  The dropdown closes when a key is
    # picked (or the query is cleared), which is the only way forward anyway.

    return wrapper


def build_and_register(
    on_open_key=None,
    on_open_tab=None,
) -> W.DOMWidget | None:
    """Build the widget; returns None where specs are unavailable."""
    try:
        return build_key_search(
            on_open_key=on_open_key, on_open_tab=on_open_tab)
    except Exception:
        return None
