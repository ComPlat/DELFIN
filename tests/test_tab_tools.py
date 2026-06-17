"""Tests for the dashboard Tools & Platform tab (renderers are display-free)."""

from __future__ import annotations

import pytest

from delfin.dashboard import tab_tools


def test_renderers_return_expected_html():
    assert "orca_sp" in tab_tools.render_capabilities_html()
    assert "opt_freq_energy" in tab_tools.render_applications_html()
    assert "functional" in tab_tools.render_keys_html()
    assert "capabilities ready" in tab_tools.render_environment_html()
    assert isinstance(tab_tools.render_runs_html(), str)
    assert isinstance(tab_tools.render_all_html(), str)


def test_tool_status_table_lists_tools_with_install_state():
    html = tab_tools.render_tool_status_html()
    assert "Tool status" in html
    assert "orca" in html                                   # tracked tools are shown
    assert ("✓ installed" in html) or ("✗ missing" in html)  # each has a state


def test_environment_renders_manual_orca_guidance(monkeypatch):
    """ORCA must be shown as manual/license-restricted with source + hint, never auto."""
    from delfin.tools import platform

    monkeypatch.setattr(platform, "probe", lambda: [])
    monkeypatch.setattr(platform, "install_plan", lambda: {
        "auto_binaries": ["xtb"],
        "manual": [{
            "tool": "orca",
            "source": "https://orcaforum.kit.edu",
            "hint": "ORCA is free for academic use; DELFIN cannot install it.",
        }],
        "python": [],
        "installer": "",
    })
    html = tab_tools.render_environment_html()
    assert "orca" in html.lower()
    assert "orcaforum.kit.edu" in html
    assert "does <u>not</u> install" in html
    assert "xtb" in html  # the auto-installable one is listed separately


def test_renderers_never_raise(monkeypatch):
    from delfin.tools import platform

    def _boom(*a, **k):
        raise RuntimeError("boom")

    monkeypatch.setattr(platform, "probe", _boom)
    monkeypatch.setattr(platform, "list_runs", _boom)
    # render_all_html wraps each section, so a failure shows a message, not a crash
    out = tab_tools.render_all_html()
    assert "unavailable" in out


def test_create_tab_builds_widget():
    if not tab_tools.HAS_WIDGETS:
        pytest.skip("ipywidgets not installed")
    widget, refs = tab_tools.create_tab(ctx=None)
    assert widget is not None
    assert "tools_panel" in refs
