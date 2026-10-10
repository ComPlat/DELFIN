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


def _find(w, pred):
    stack = [w]
    while stack:
        n = stack.pop()
        if pred(n):
            return n
        stack.extend(getattr(n, "children", ()) or ())
    return None


def _has_input(c):
    import ipywidgets as widgets
    return _find(c, lambda n: isinstance(n, widgets.Textarea)) is not None


def _frame(root):
    kids = list(root.children)
    return kids[_index_of(
        kids, lambda c: "delfin-agent-chat-frame" in _classes(c))]


def test_the_root_fits_the_window_and_the_chat_takes_the_rest(tmp_path):
    from delfin.dashboard.tab_agent import _AGENT_CSS
    root = _tab(tmp_path)
    assert "delfin-agent-root" in _classes(root)
    assert "--delfin-agent-top" in _AGENT_CSS
    assert ".delfin-agent-root > .delfin-agent-chat-frame" in _AGENT_CSS
    assert ".delfin-agent-chat-frame > .delfin-agent-chat-host" in _AGENT_CSS
    frame = _frame(root)
    assert "delfin-agent-chat-host" in _classes(frame.children[0])


def test_the_approvals_sit_at_the_foot_of_the_chat_above_the_input(tmp_path):
    root = _tab(tmp_path)
    kids = list(root.children)
    frame_i = _index_of(
        kids, lambda c: "delfin-agent-chat-frame" in _classes(c))
    frame = kids[frame_i]
    # The dock closes the chat frame, and the input follows the frame.
    assert "delfin-agent-dock" in _classes(frame.children[-1])
    assert _index_of(kids, _has_input) == frame_i + 1


def test_the_panels_below_are_folded(tmp_path):
    import ipywidgets as widgets
    root = _tab(tmp_path)
    below = root.children[-1]
    assert "delfin-agent-below" in _classes(below)
    assert isinstance(below, widgets.Accordion)
    assert below.selected_index is None


def test_the_slash_button_sits_beside_the_input(tmp_path):
    import ipywidgets as widgets
    root = _tab(tmp_path)
    row = _find(root, lambda n: isinstance(n, widgets.HBox) and any(
        isinstance(b, widgets.Button) and b.description == "/"
        for b in n.children))
    assert row is not None, "no row with the / button"
    assert _has_input(row)
    # The suggestions sit at the foot of the chat, in the palette row.
    assert _find(_frame(root),
                 lambda n: "delfin-next-steps" in _classes(n)) is not None


def test_the_task_line_unfolds_to_the_list_and_lives_only_there(tmp_path):
    import ipywidgets as widgets
    root = _tab(tmp_path)
    strip = _find(_frame(root),
                  lambda n: "delfin-agent-task-strip" in _classes(n))
    assert isinstance(strip, widgets.Accordion)
    assert strip.selected_index is None
    assert strip.layout.display == "none"          # no tasks, no line
    ticker = strip.children[0]
    below = root.children[-1]
    assert _find(below, lambda n: n is ticker) is None


def test_subagent_buttons_wrap_side_by_side():
    from pathlib import Path
    src = Path(__file__).resolve().parents[1].joinpath(
        "delfin", "dashboard", "tab_agent.py").read_text()
    i = src.index("agent_view_chips = widgets.")
    assert src[i:i + 200].startswith("agent_view_chips = widgets.HBox(")
    assert 'flex_flow="row wrap"' in src[i:i + 200]


def test_an_open_request_does_not_hide_the_task_line():
    from delfin.dashboard.tab_agent import _AGENT_CSS
    rules = [r.split("{")[0] for r in _AGENT_CSS.split("}")
             if ".delfin-agent-request:not(" in r.split("{")[0]
             and ".delfin-agent-chat-foot" in r.split("{")[0]]
    assert rules
    assert all(":not(.delfin-agent-task-strip)" in r for r in rules)


def test_an_open_request_makes_the_chat_give_way():
    """The confirm buttons were clipped under the chat frame's bottom edge;
    while a request is open the transcript may shrink and the dock scrolls
    in itself. The chat's inline min-height must be overridden."""
    from delfin.dashboard.tab_agent import _AGENT_CSS
    css = " ".join(_AGENT_CSS.split())
    assert ("> .delfin-agent-chat-host { /* !important: the widget carries "
            "min-height 200px inline. */ min-height: 60px !important;") in css
    assert (".delfin-agent-chat-frame > .delfin-agent-dock { flex: 0 1 auto "
            "!important; min-height: 0; overflow-y: auto;") in css
    # A KIT confirmation takes the task line's place; a question does not.
    assert ("kit-confirm > :not([style*=\"display: none\"])) > "
            ".delfin-agent-task-strip { display: none !important;") in css


def test_open_tasks_are_clicked_in_the_list_and_suggestions_offer_the_rest():
    from pathlib import Path
    src = Path(__file__).resolve().parents[1].joinpath(
        "delfin", "dashboard", "tab_agent.py").read_text()
    i = src.index("def _task_buttons(rows):")
    body = src[i:src.index("def _refresh_task_ticker", i)]
    assert 'b.on_click(_fill_input(r["subject"]))' in body
    assert "disabled=not open_" in body
    # The buttons beside the box: only the agent's offer, never a task.
    j = src.index("def _refresh_next_steps")
    nxt = src[j:src.index("def _", j + 30)]
    assert "shown = [offer] if offer and _norm(offer) not in" in nxt
    assert "for step in shown:" in nxt


def test_tab_takes_the_suggestion_behind_its_own_guard():
    """Two keyboard scripts share __delfinAgentKeys; the one that ran
    second returned early, so the Tab listener was never installed."""
    from pathlib import Path
    src = Path(__file__).resolve().parents[1].joinpath(
        "delfin", "dashboard", "tab_agent.py").read_text()
    i = src.index("window.__delfinTabTakesSuggestion = true;")
    guard = src.index("if (window.__delfinAgentKeys) return;", i - 2000)
    assert i < guard, "the Tab listener must be installed before the guard"
    assert "document.addEventListener('keydown'" in src[i:guard]
    assert ", true);" in src[i:guard]
