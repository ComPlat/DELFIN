"""The dashboard's tab bar: one row, whole names, a calm look.

Input: the page CSS (helpers.create_page_css) and the tab strip script
(helpers.keep_dashboard_tab_in_view). Output pinned here:

* one row: the strip never wraps (on main a second row held "Settings"
  alone at 1366 px); tabs are as wide as their names, never clipped
  ("Submit Jo" before 2026-11-12); where the names do not fit, the strip
  scrolls sideways -- wheel, trackpad or keyboard -- with a fade on the
  side that has more, and the current tab is scrolled into view;
* flat tabs: no grey boxes; inactive labels in the theme's secondary text
  colour (4.6:1 on the white bar, AA), the current one in the primary
  colour with an accent underline, the hover a light background;
* the padding follows the window (clamp), so the eighteen tabs fit one row
  at 1366 px and do not crowd at 1920 px;
* the height stays that of the old single row (28 px): the tab bodies are
  calc(100vh - 145px) tall.

Measured in the real page with scripts/page_check.py (rows, clipped labels,
scroll) at 1100, 1366 and 1920 px; this file pins the decisions in Python.
"""
from __future__ import annotations

import inspect

from delfin.dashboard import helpers

BAR = f".jupyter-widgets.jupyter-widget-tab.{helpers.DASHBOARD_TABS_CLASS} > .lm-TabBar"


def _css() -> str:
    return " ".join(helpers.create_page_css().value.split())


def _rule(css: str, selector: str) -> str:
    start = css.index(selector + " {")
    return css[start:css.index("}", start) + 1]


def test_the_strip_is_one_row_that_scrolls_sideways():
    css = _css()
    strip = _rule(css, f"{BAR} > .lm-TabBar-content")
    assert "flex-wrap: nowrap" in strip
    assert "overflow-x: auto" in strip and "overflow-y: hidden" in strip
    assert "flex-wrap: wrap" not in css


def test_a_tab_is_as_wide_as_its_name_and_never_clipped():
    css = _css()
    tab = _rule(css, f"{BAR} .lm-TabBar-tab")
    assert "flex: 0 0 auto" in tab
    assert "padding: 0 clamp(" in tab
    assert "height: 28px" in tab
    assert "background: transparent" in tab and "border: none" in tab
    label = _rule(css, f"{BAR} .lm-TabBar-tabLabel")
    assert "overflow: visible" in label


def test_inactive_current_and_hover_are_told_apart_calmly():
    css = _css()
    inactive = _rule(css, f"{BAR} .lm-TabBar-tab:not(.lm-mod-current)")
    assert "color: var(--jp-ui-font-color2)" in inactive
    current = _rule(css, f"{BAR} .lm-TabBar-tab.lm-mod-current")
    assert "color: var(--jp-ui-font-color0)" in current
    assert "border-bottom-color: var(--delfin-accent" in current
    assert "transform: none" in current
    hover = _rule(css, f"{BAR} .lm-TabBar-tab:hover:not(.lm-mod-current)")
    assert "background: var(--jp-layout-color2)" in hover
    assert f"{BAR} .lm-TabBar-tab.lm-mod-current::before {{ display: none; }}" in css


def test_the_side_with_more_tabs_fades():
    css = _css()
    for side in ("more-left", "more-right"):
        assert f"{BAR}.delfin-tabs-{side} > .lm-TabBar-content {{" in css
    assert "mask-image" in css


def _strip_script() -> str:
    class Ctx:
        def __init__(self):
            self.init_js_parts = []

        def add_init_js(self, script):
            self.init_js_parts.append(script)
    ctx = Ctx()
    helpers.keep_dashboard_tab_in_view(ctx)
    assert len(ctx.init_js_parts) == 1
    return ctx.init_js_parts[0]


def test_the_strip_script_marks_edges_reveals_the_current_tab_and_scrolls_by_wheel():
    script = _strip_script()
    assert helpers.DASHBOARD_TABS_CLASS in script
    assert "delfin-tabs-more-left" in script and "delfin-tabs-more-right" in script
    assert "lm-mod-current" in script
    assert "'wheel'" in script and "passive: false" in script
    # Only when the strip overflows and the wheel is vertical: a page that
    # fits keeps its own wheel.
    assert "scrollWidth <= c.clientWidth" in script
    for timer in ("setTimeout", "setInterval"):
        assert timer not in script, timer


def test_the_dashboard_registers_the_strip_script():
    from delfin import dashboard
    source = inspect.getsource(dashboard)
    assert "keep_dashboard_tab_in_view(ctx)" in source
    assert source.index("keep_dashboard_tab_in_view(ctx)") < source.index("*ctx.init_js_parts,")
