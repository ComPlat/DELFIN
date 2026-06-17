"""Additive registry for dashboard tabs.

Lets code register extra Voila/dashboard tabs **without editing**
``create_dashboard``.  A registered tab provides a factory with the same
contract as the built-in ``tab_*.py`` modules — ``factory(ctx) -> (widget,
refs)`` (or just ``widget``) — plus display metadata.  ``create_dashboard``
appends the registered tabs *after* its hardcoded tabs, so the built-in tabs and
their behaviour are completely unchanged; with nothing registered the hook is a
no-op.

    from delfin.dashboard.tab_registry import register_tab

    def create_tab(ctx):
        return my_widget, {"my_panel": panel}

    register_tab("my_tab", "My Tab", create_tab, order=8200)
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Callable, Dict, List


@dataclass
class TabSpec:
    """A dynamically-registered dashboard tab."""

    id: str
    title: str
    factory: Callable[[Any], Any]   # ctx -> (widget, refs) | widget
    order: int = 8000               # placed after built-ins (which use 10..10000)
    visible: bool = True


_REGISTRY: Dict[str, TabSpec] = {}
_DISCOVERED = False

# Built-in dynamic tabs that register themselves on import.
_BUILTIN_DYNAMIC_TABS = ("delfin.dashboard.tab_application",)


def register_tab(
    tab_id: str,
    title: str,
    factory: Callable[[Any], Any],
    *,
    order: int = 8000,
    visible: bool = True,
) -> None:
    """Register a dashboard tab (idempotent by ``tab_id``)."""
    _REGISTRY[tab_id] = TabSpec(id=tab_id, title=title, factory=factory,
                                order=order, visible=visible)


def _discover() -> None:
    global _DISCOVERED
    if _DISCOVERED:
        return
    _DISCOVERED = True
    for mod in _BUILTIN_DYNAMIC_TABS:
        try:
            __import__(mod)
        except Exception:
            pass


def registered_tabs() -> List[TabSpec]:
    """All registered tabs, ordered (triggers lazy discovery)."""
    _discover()
    return sorted(_REGISTRY.values(), key=lambda t: (t.order, t.id))


def registered_tab_specs(ctx: Any) -> List[dict]:
    """Build dashboard ``tab_specs`` dicts for the registered tabs.

    Each factory is invoked defensively: a tab that raises is included but marked
    unavailable (``widget=None``) so it can never break dashboard startup.
    """
    specs: List[dict] = []
    for tab in registered_tabs():
        widget = None
        reason = ""
        try:
            result = tab.factory(ctx)
            widget = result[0] if isinstance(result, tuple) else result
        except Exception as exc:  # noqa: BLE001
            reason = f"tab '{tab.id}' failed to build: {exc}"
        specs.append({
            "id": tab.id,
            "title": tab.title,
            "widget": widget,
            "default_order": tab.order,
            "default_visible": tab.visible,
            "available": widget is not None,
            "fixed": False,
            "reason": reason,
        })
    return specs


__all__ = ["TabSpec", "register_tab", "registered_tabs", "registered_tab_specs"]
