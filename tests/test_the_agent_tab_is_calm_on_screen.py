"""The agent tab shows what matters and nothing that alarms for no reason.

Measured in a real page at 1920x1080 with a filled chemistry session
(2026-10-22):
- the control row wrapped its last button onto a second line, and
  "Perms:" was cut to "Per...";
- the git field and four git buttons sat greyed out in a folder that is
  not a repository;
- the closed Details header read "⚠ Containment 0 blocked / 2 flagged",
  the last tool's internal name and token counts;
- an idle Stop was pale red;
- the input placeholder carried a two-line KIT tip;
- plots were drawn at matplotlib's defaults and shown at 576 px in a
  1500 px chat, a 300 dpi energy profile filled the whole chat, and a
  filename ending in .svg came back as ".svg.png".
"""

from __future__ import annotations

import re


def _tab(tmp_path, monkeypatch, workspace):
    from delfin.agent import scheduler as S
    from delfin.dashboard import tab_agent
    from delfin.dashboard.context import DashboardContext

    S._GLOBAL = S.Scheduler(path=tmp_path / "cron.json")
    for name in ("calc", "archive", "office"):
        (tmp_path / name).mkdir(exist_ok=True)
    monkeypatch.setenv("DELFIN_AGENT_WORKSPACE", str(workspace))
    monkeypatch.setenv("DELFIN_LAUNCH_CWD", str(workspace))
    ctx = DashboardContext(calc_dir=tmp_path / "calc",
                           archive_dir=tmp_path / "archive",
                           office_dir=tmp_path / "office")
    ctx.run_js = lambda script, **kw: None
    root, refs = tab_agent.create_tab(ctx)
    return root, refs


def _walk(n):
    stack = [n]
    while stack:
        x = stack.pop()
        yield x
        stack.extend(getattr(x, "children", ()) or ())


def test_git_controls_are_hidden_outside_a_repository(tmp_path, monkeypatch):
    import ipywidgets as widgets
    ws = tmp_path / "plain-folder"
    ws.mkdir()
    root, refs = _tab(tmp_path, monkeypatch, ws)
    try:
        base = next(w for w in _walk(root) if isinstance(w, widgets.Combobox)
                    and w.description == "Base:")
        assert base.layout.display == "none"
        for label in ("Merge", "Branch", "Push", "PR"):
            hits = [w for w in _walk(root) if isinstance(w, widgets.Button)
                    and label in (w.description or "")
                    and w.layout.display != "none"]
            assert not hits, label
    finally:
        refs["shutdown"]()


def test_the_permission_label_is_not_cut_and_stop_is_grey_when_idle(
        tmp_path, monkeypatch):
    import ipywidgets as widgets
    ws = tmp_path / "ws"
    ws.mkdir()
    root, refs = _tab(tmp_path, monkeypatch, ws)
    try:
        perms = next(w for w in _walk(root) if isinstance(w, widgets.Dropdown)
                     and w.description == "Perms:")
        assert int(perms.style.description_width.rstrip("px")) >= 48
        stop = next(w for w in _walk(root) if isinstance(w, widgets.Button)
                    and w.description == "Stop")
        assert stop.disabled and stop.button_style == ""
        textarea = next(w for w in _walk(root)
                        if isinstance(w, widgets.Textarea)
                        and "agent-input" in " ".join(w._dom_classes))
        assert "\n" not in textarea.placeholder
        assert "KIT tip" not in textarea.placeholder
    finally:
        refs["shutdown"]()


def test_the_details_header_warns_only_about_what_was_blocked(
        tmp_path, monkeypatch):
    ws = tmp_path / "ws"
    ws.mkdir()
    root, refs = _tab(tmp_path, monkeypatch, ws)
    try:
        state, title = refs["state"], refs["fold_title"]
        state["_fold_ctx"] = 48.0
        state["_fold_tools"] = 10
        state["_fold_last_tool"] = "mcp__delfin-docs__show_molecule"
        state["_fold_security"] = (2, 0)          # flagged, nothing blocked
        state["_fold_tokens"] = (1200, 120, 0.0)
        assert title() == "Context 48% · 10 tool calls"
        state["_fold_security"] = (3, 1)
        assert title() == "Context 48% · 10 tool calls · ⚠ 1 blocked"
        state.update(_fold_ctx=None, _fold_tools=0, _fold_security=(0, 0))
        assert title() == "Details"
    finally:
        refs["shutdown"]()


def test_a_line_plot_marks_its_points_and_takes_its_format_from_the_name(
        tmp_path, monkeypatch):
    from delfin.agent import chat_plots
    spec = {"kind": "line", "x": [0, 15, 30], "y": [5.0, 4.3, 2.8],
            "xlabel": "angle", "ylabel": "E", "units": "kcal/mol",
            "filename": "scan.svg"}
    marked = chat_plots.plot(dict(spec), out_dir=str(tmp_path / "a"))
    monkeypatch.setattr(chat_plots, "_MARKERS_UP_TO", 0)
    plain = chat_plots.plot(dict(spec), out_dir=str(tmp_path / "b"))
    assert marked.path.name == "scan.svg"
    uses = lambda p: p.read_text().count("<use ")          # noqa: E731
    # The same axes and ticks; the difference is one marker per point.
    assert uses(marked.path) - uses(plain.path) == 3
    pes = {"kind": "energy_profile", "filename": "pes.svg", "ylabel": "ΔE",
           "series": [{"label": "a", "states": [
               {"name": "A", "value": 0.0}, {"name": "B", "value": 1.0}]}],
           "from_user": True}
    assert chat_plots.resolve_out_path(pes, tmp_path).name == "pes.svg"
    pes.pop("filename")
    # Without a name the group's 300 dpi PNG stays the default.
    assert chat_plots.resolve_out_path(pes, tmp_path).suffix == ".png"


def test_figures_are_shown_at_reading_size(tmp_path):
    from delfin.dashboard import tab_agent
    svg = tmp_path / "f.svg"
    svg.write_text('<svg xmlns="http://www.w3.org/2000/svg" width="10" '
                   'height="10"></svg>')
    html = tab_agent._render_artifact_body(svg)
    assert "width:min(900px,100%)" in html and "margin:0 auto" in html
    png = tmp_path / "f.png"
    png.write_bytes(bytes.fromhex(
        "89504e470d0a1a0a0000000d4948445200000001000000010806000000"
        "1f15c4890000000d49444154789c6360000002000154a24f5d0000000049454e44ae426082"))
    html = tab_agent._render_artifact_body(png)
    assert "max-width:min(900px,100%)" in html
