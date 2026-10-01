"""The agent tab has an \"Export as notebook\" button.

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


def test_the_button_sits_next_to_the_markdown_export():
    i = SRC.index("export_btn = widgets.Button(")
    j = SRC.index("nb_export_btn = widgets.Button(")
    assert abs(i - j) < 1500, "notebook button should share the export row"
    k = SRC.index("nb_export_btn,")
    # the row layout lists both buttons together
    assert "export_btn," in SRC[k - 400:k]


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
