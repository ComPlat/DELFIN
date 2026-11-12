"""Every tab shows its whole name, every control has one, and a phone fits.

An accessibility check (axe, WCAG 2 AA) and a look at the page found:
- the tab bar clipped its labels ("Submit Jo", "Calculatic") and the
  widgets give a tab no tooltip to read the rest from;
- an inactive tab's label measured below 4.5:1 contrast on the bar's grey;
- the "Resume a saved session…" select, the three A-Z sort selects and the
  structure editor's sliders had no accessible name;
- at phone width the agent's session list kept its 236 px column and left
  the chat a 50 px strip.
These pin the fixes at the level of the page's CSS and scripts.
"""
from __future__ import annotations

import inspect
import re

from delfin.dashboard import helpers


def _css(widget) -> str:
    return " ".join(widget.value.split())


def test_the_dashboard_tab_bar_wraps_instead_of_clipping():
    css = _css(helpers.create_page_css())
    bar = f".jupyter-widgets.jupyter-widget-tab.{helpers.DASHBOARD_TABS_CLASS}"
    assert f"{bar} > .lm-TabBar > .lm-TabBar-content {{ flex-wrap: wrap; }}" in css
    assert f"{bar} > .lm-TabBar .lm-TabBar-tab {{ flex: 0 0 auto; }}" in css
    assert f"{bar} > .lm-TabBar .lm-TabBar-tabLabel {{ overflow: visible; }}" in css


def test_an_inactive_tab_label_has_the_primary_text_colour():
    css = _css(helpers.create_page_css())
    assert (".lm-TabBar-tab:not(.lm-mod-current) { "
            "color: var(--jp-ui-font-color1); }") in css


def test_the_dashboard_marks_its_own_tab_bar_and_names_its_controls():
    from delfin import dashboard
    source = inspect.getsource(dashboard)
    assert "ctx.tabs_widget.add_class(DASHBOARD_TABS_CLASS)" in source
    assert "name_page_controls(ctx)" in source
    # Registered before the startup script is assembled, so it goes out
    # with the page and not after it.
    assert source.index("name_page_controls(ctx)") < source.index(
        "*ctx.init_js_parts,")


class _Ctx:
    def __init__(self):
        self.init_js_parts = []

    def add_init_js(self, script):
        self.init_js_parts.append(script)


def _naming_script() -> str:
    ctx = _Ctx()
    helpers.name_page_controls(ctx)
    assert len(ctx.init_js_parts) == 1
    return ctx.init_js_parts[0]


def test_every_tab_gets_its_full_name_as_a_tooltip():
    script = _naming_script()
    assert "'.lm-TabBar-tab'" in script
    assert "tab.setAttribute('title', name)" in script
    # Renamed in place: the title follows when the pointer or focus arrives.
    assert "'pointerover'" in script and "'focusin'" in script


def test_unnamed_selects_and_slider_handles_get_an_accessible_name():
    script = _naming_script()
    assert "select:not([aria-label])" in script
    assert ".noUi-handle:not([aria-label])" in script
    assert "el.setAttribute('aria-label', name)" in script
    # The words come from the widget itself: description, else tooltip.
    assert "widget-label" in script and "title" in script


def test_naming_adds_no_timer():
    script = _naming_script()
    for timer in ("setTimeout", "setInterval", "requestAnimationFrame"):
        assert timer not in script, timer
    assert script.count("new MutationObserver") == 1


def test_the_selects_carry_a_name_of_their_own():
    """A tooltip on a widget without a description is the ``title`` of its
    ``<select>`` (ipywidgets), a name even before the page script runs."""
    from delfin.dashboard import (agent_sessions, tab_calculations_browser,
                                  tab_literature, tab_remote_archive)
    assert 'tooltip="Resume a saved session"' in inspect.getsource(agent_sessions)
    for module in (tab_calculations_browser, tab_literature, tab_remote_archive):
        source = inspect.getsource(module)
        sort = re.search(r"options=\[\(['\"]A-Z['\"].{0,200}", source, re.S)
        assert sort and "tooltip=" in sort.group(0) and "Sort files" in sort.group(0), \
            module.__name__


def test_the_session_list_goes_above_the_chat_on_a_narrow_window():
    from delfin.dashboard import agent_sessions as AS
    css = " ".join(AS._SIDEBAR_CSS.split())
    start = css.index("@media (max-width: 800px)")
    block = css[start:css.index("} }", start) + 3]
    assert "flex-direction: column !important" in block
    assert ".delfin-sessions { position: static; flex: 0 0 auto !important;" in block
    assert "width: 100% !important" in block
    assert ".delfin-splitter-host { display: none !important; }" in block
