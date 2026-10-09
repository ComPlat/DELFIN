"""The agent tab fits the window, and the approvals sit above the input.

User, 2026-10-09: approvals were missed because they sat below the input
and the task list, and the page had to be scrolled to find them; the tab
should fit the window, panels should not push it longer, suggestions
belong in the row of the "/" button, and subagent buttons side by side.
"""

from __future__ import annotations


def _tab(tmp_path):
    from delfin.agent import scheduler as S
    from delfin.dashboard import tab_agent
    from delfin.dashboard.context import DashboardContext

    S._GLOBAL = S.Scheduler(path=tmp_path / "cron.json")
    for name in ("calc", "archive", "office"):
        (tmp_path / name).mkdir(exist_ok=True)
    ctx = DashboardContext(calc_dir=tmp_path / "calc",
                           archive_dir=tmp_path / "archive",
                           office_dir=tmp_path / "office")
    ctx.run_js = lambda script, **kw: None
    tab = tab_agent.create_tab(ctx)
    return tab[0] if isinstance(tab, tuple) else tab


def _classes(w):
    return set(getattr(w, "_dom_classes", ()) or ())


def _index_of(children, pred):
    for i, c in enumerate(children):
        if pred(c):
            return i
    raise AssertionError("not found")


def test_the_root_fits_the_window_and_the_chat_takes_the_rest(tmp_path):
    from delfin.dashboard.tab_agent import _AGENT_CSS
    root = _tab(tmp_path)
    assert "delfin-agent-root" in _classes(root)
    assert "--delfin-agent-top" in _AGENT_CSS
    assert ".delfin-agent-root > .delfin-agent-chat-host" in _AGENT_CSS


def test_the_approvals_sit_directly_above_the_input(tmp_path):
    import ipywidgets as widgets
    root = _tab(tmp_path)
    kids = list(root.children)
    dock = _index_of(kids, lambda c: "delfin-agent-dock" in _classes(c))
    chat = _index_of(kids, lambda c: "delfin-agent-chat-host" in _classes(c))

    def _has_input(c):
        stack = [c]
        while stack:
            n = stack.pop()
            if isinstance(n, widgets.Textarea):
                return True
            stack.extend(getattr(n, "children", ()) or ())
        return False

    typing = _index_of(kids, _has_input)
    assert chat < dock < typing
    # Nothing between the dock and the input but the "/" row.
    assert typing - dock <= 3


def test_the_panels_below_are_capped(tmp_path):
    from delfin.dashboard.tab_agent import _AGENT_CSS
    root = _tab(tmp_path)
    assert "delfin-agent-below" in _classes(root.children[-1])
    assert "max-height: 30vh" in _AGENT_CSS


def test_suggestions_share_the_row_of_the_slash_button(tmp_path):
    import ipywidgets as widgets
    root = _tab(tmp_path)
    rows = [n for n in root.children if isinstance(n, widgets.HBox)
            and any(isinstance(b, widgets.Button) and b.description == "/"
                    for b in n.children)]
    assert rows, "no row with the / button"
    assert len(rows[0].children) == 3          # /, filter, suggestions
    assert rows[0].layout.flex_flow == "row wrap"


def test_subagent_buttons_wrap_side_by_side():
    from pathlib import Path
    src = Path(__file__).resolve().parents[1].joinpath(
        "delfin", "dashboard", "tab_agent.py").read_text()
    i = src.index("agent_view_chips = widgets.")
    assert src[i:i + 200].startswith("agent_view_chips = widgets.HBox(")
    assert 'flex_flow="row wrap"' in src[i:i + 200]
