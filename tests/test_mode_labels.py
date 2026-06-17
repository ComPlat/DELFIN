"""The internal 'solo' id must surface to the user as 'Code' everywhere — the
mode was renamed (dropdown shows Dashboard + Code) but the badges leaked the
old internal name (Jerome report 20260617-152728: 'noch paar Stellen mit Solo
nicht Code')."""

from __future__ import annotations


def test_solo_role_and_mode_render_as_code():
    from delfin.dashboard.tab_agent import _format_role_label, _mode_label
    assert _format_role_label("solo_agent") == "Code Agent"
    assert _format_role_label("solo") == "Code Agent"
    assert _mode_label("solo") == "Code"
    # other modes/roles unchanged
    assert _mode_label("dashboard") == "Dashboard"
    assert _format_role_label("dashboard_agent") == "Dashboard Agent"
    assert _format_role_label("") == "Agent"


def test_no_user_facing_solo_label_leaks():
    """The status badge + the 'switch the dropdown to ...' hints must say Code,
    not the internal 'solo'."""
    from pathlib import Path
    src = (Path(__file__).resolve().parent.parent / "delfin" / "dashboard"
           / "tab_agent.py").read_text(encoding="utf-8")
    assert "Switch the Mode dropdown to **solo**" not in src
    assert "In **solo mode**" not in src
    assert 'mode-badge">{_html.escape(_mode_label(mode))}' in src
