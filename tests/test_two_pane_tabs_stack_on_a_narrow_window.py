"""Two-pane tabs put their panes one above the other on a narrow window.

Input: Submit Job (form | molecule preview) and the file browsers
(Calculations, Archive, Office, Home: list | content). Measured in the real
page (page_check driver):

* 390x844 (a phone): Submit Job kept both panes side by side at about
  170 px each; the preview's control panel covered the Job Name field and
  the Time Limit buttons (three controls painted over), labels wrapped
  into one word per line. The Calculations browser gave the list its
  fixed width and left the content pane no width at all: a selected
  file could not be read.
* 900x1400: Submit Job's left pane (about 430 px) cut its third Time
  Limit button column and the Job Name field at the pane edge.

Output pinned: a media rule per tab -- Submit Job stacks at 960 px and
below, the browsers at 700 px and below (list on top, 45 % of the height,
content below; the splitter, which only drags sideways, is hidden).
Above those widths nothing changes.
"""
from __future__ import annotations

from pathlib import Path

ROOT = Path(__file__).resolve().parents[1] / "delfin" / "dashboard"


def _flat(name: str) -> str:
    return " ".join((ROOT / name).read_text(encoding="utf-8").split())


def test_submit_job_stacks_its_panes_at_960px():
    src = _flat("tab_submit.py")
    start = src.index("@media (max-width: 960px) {")
    block = src[start:src.index("} }", start) + 3]
    assert ".submit-split-root.widget-hbox { flex-direction: column !important; }" in block
    assert ".submit-split-root > .submit-split-pane { flex: 0 0 auto !important; width: 100% !important;" in block


def test_the_file_browsers_stack_their_panes_at_700px():
    src = _flat("tab_calculations_browser.py")
    assert "add_class('calc-main-row')" in src and "add_class('calc-splitter-host')" in src
    start = src.index("@media (max-width: 700px) {")
    block = src[start:src.index("}'", src.index("calc-splitter-host", start)) + 2]
    assert ".calc-tab .calc-main-row { flex-direction:column !important; }" in block
    assert ".calc-tab .calc-main-row > .calc-left { flex:0 0 45% !important;" in block
    assert ".calc-tab .calc-main-row > .calc-right { flex:1 1 0 !important;" in block
    assert ".calc-tab .calc-main-row > .calc-splitter-host { display:none !important; }" in block


def test_the_time_limit_buttons_fit_a_phone():
    """The Time Limit toggle is 550 px wide (its layout) with 148 px
    buttons; stacked on a 390 px phone its third column ran past the pane
    and was cut ('36h', 'Custom'). Below 600 px each button takes a third
    of the row and the toggle never exceeds its pane."""
    src = _flat("helpers.py")
    assert ".delfin-time-limit-toggle.widget-toggle-buttons { max-width: 100% !important; }" in src
    start = src.index("@media (max-width: 600px) {")
    block = src[start:src.index("} }", start) + 3]
    assert (".delfin-time-limit-toggle.widget-toggle-buttons .widget-toggle-button "
            "{ width: calc(33.333% - 4px) !important; min-width: 0 !important; }") in block
