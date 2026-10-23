"""The agent tab's "Export as notebook" handler (button now hidden).

Source-level test: the dashboard button row is UI glue, so this pins the
wiring, not pixel behaviour — the button exists next to the Markdown
export, its handler exports via ``session_export`` into the tab's
exports dir, and errors surface as a system message instead of killing
the tab.
"""
from __future__ import annotations

from pathlib import Path

SRC = Path(__file__).resolve().parents[1].joinpath(
    "delfin", "dashboard", "tab_agent.py").read_text()


def test_both_exports_left_the_visible_row():
    """Taken out of the top row (2026-10-09): the git buttons took their
    place, and a sent bug report downloads a copy instead. The widgets
    stay wired, as the other retired buttons do."""
    i = SRC.index("for _hidden in (")
    hidden = SRC[i:i + 400]
    assert "export_btn" in hidden and "nb_export_btn" in hidden
    j = SRC.index("_controls_hbox = widgets.HBox(")
    row = SRC[j:j + 300]
    assert "export_btn" not in row and "git_group" in row


def test_the_handler_exports_via_session_export():
    i = SRC.index("def _on_export_notebook(button):")
    body = SRC[i:i + 2500]
    assert "session_export" in body
    assert "exports" in body


def test_the_button_is_wired_to_the_handler():
    assert "nb_export_btn.on_click(_on_export_notebook)" in SRC


def test_export_failures_surface_as_system_messages():
    i = SRC.index("def _on_export_notebook(button):")
    body = SRC[i:i + 2500]
    assert "_append_system_message" in body
