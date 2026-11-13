"""Which dashboard tabs show, and in which order, from the settings.

Input: the tab specs (``id``, ``default_order``, ``default_visible``,
``fixed``) and the person's ``ui.tabs`` settings -- ``order``, the
arrangement Settings > Dashboard Tabs saves (every tab of its day), and
``hidden``, the tabs switched off there.

Output: :func:`resolve_tab_preferences` -> ``(order, hidden)`` over the
movable tabs; the tab bar shows ``order`` minus ``hidden``, then the fixed
ones (Settings).

The defaults (:data:`FIRST_TABS`, applied by :func:`apply_first_tabs`)
decide every tab the saved arrangement does not name: all of them for a
person who never arranged the tabs -- a new installation first of all --
and a tab added after the arrangement was saved. A saved arrangement is
kept whole.

Rejected: telling a new installation from an old one by the settings file
(the installer writes it before the dashboard ever starts), and a
per-installation marker written at the first start (a person who never
arranged the tabs would see a different bar depending on when they first
opened the dashboard).
"""
from __future__ import annotations

#: The tabs a person sees before arranging any, in this order (2026-11-13).
#: Every other tab stays available: Settings > Dashboard Tabs switches it on.
FIRST_TABS = ("submit", "job_status", "agent", "home", "calculations")


def apply_first_tabs(specs: list) -> list:
    """Set ``default_visible`` and the leading ``default_order`` of every
    movable spec from FIRST_TABS (in place; returns ``specs``). Fixed tabs
    (Settings) are left as they are: always shown."""
    for spec in specs:
        if spec.get("fixed"):
            continue
        tab_id = spec.get("id")
        spec["default_visible"] = tab_id in FIRST_TABS
        if tab_id in FIRST_TABS:
            spec["default_order"] = FIRST_TABS.index(tab_id) + 1
    return specs


def resolve_tab_preferences(specs: list, prefs: dict | None) -> tuple[list, set]:
    """``(order, hidden)`` of the movable tabs for ``prefs`` (``ui.tabs``).

    The saved order first (unknown and fixed ids dropped), then every other
    tab by ``default_order``. A tab the saved order names is hidden when the
    saved hidden list says so; any other tab when its default says so."""
    prefs = prefs or {}
    spec_map = {spec.get("id"): spec for spec in specs or []}
    raw_order = [str(item) for item in (prefs.get("order", []) or []) if str(item).strip()]
    raw_hidden = {str(item) for item in (prefs.get("hidden", []) or []) if str(item).strip()}
    order = []
    seen = set()
    for tab_id in raw_order:
        spec = spec_map.get(tab_id)
        if not spec or spec.get("fixed") or tab_id in seen:
            continue
        seen.add(tab_id)
        order.append(tab_id)
    for spec in sorted(specs or [], key=lambda item: item.get("default_order", 10_000)):
        tab_id = spec.get("id")
        if not tab_id or spec.get("fixed") or tab_id in seen:
            continue
        seen.add(tab_id)
        order.append(tab_id)
    hidden = set()
    saved_ids = set(raw_order)
    for tab_id in order:
        spec = spec_map.get(tab_id) or {}
        if tab_id in saved_ids:
            if tab_id in raw_hidden:
                hidden.add(tab_id)
        elif not bool(spec.get("default_visible", True)):
            hidden.add(tab_id)
    return order, hidden
