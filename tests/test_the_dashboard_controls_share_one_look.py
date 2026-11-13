"""Text fields, dropdowns, text areas and buttons share one look.

Input: helpers.DASHBOARD_FORM_CSS, part of the page CSS (create_page_css).
Output pinned here: the widgets' own controls (square inputs and selects,
2 px buttons, no focus mark beyond a border colour) take the header's
shape -- 6 px corners, the accent border and a soft ring on focus, a
visible keyboard focus on buttons -- across every tab, from one place
instead of per-widget styles.

The layer must lose to every DELFIN rule written for a particular control
(the composer, the chips, the session list): its selectors are the
widgets' own (same specificity) and carry no !important. The page CSS
sits in the header cell, before every tab's styles, so at equal
specificity a tab's rule wins; ipywidgets' styles live in the document
head, before it, so this layer wins over them.
"""
from __future__ import annotations

import re

from delfin.dashboard import helpers


def _css() -> str:
    return " ".join(helpers.DASHBOARD_FORM_CSS.split())


def test_the_page_css_carries_the_form_layer():
    page = " ".join(helpers.create_page_css().value.split())
    assert _css() in page


def test_fields_dropdowns_and_buttons_are_rounded_like_the_header():
    css = _css()
    for selector in (".widget-text input[type='text']", ".widget-text input[type='number']",
                     ".widget-text input[type='password']", ".widget-textarea textarea",
                     ".widget-dropdown > select"):
        assert selector in css, selector
    assert "border-radius: 6px" in css
    assert re.search(r"\.jupyter-button \{ border-radius: 6px; \}", css)


def test_focus_is_shown_without_raising_specificity():
    css = _css()
    assert ":where(:focus)" in css
    assert "box-shadow: 0 0 0 2px" in css
    assert ".jupyter-button:where(:focus-visible)" in css


def test_the_layer_never_overrides_a_rule_written_for_one_control():
    css = _css()
    assert "!important" not in css
    # No id, no tab or DELFIN class: only the widgets' own selectors.
    assert "#" not in css.replace("#000", "")
    assert "delfin-" not in css.replace("--delfin-accent", "")
