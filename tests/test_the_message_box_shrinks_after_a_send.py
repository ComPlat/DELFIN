"""The message box returns to its resting size when a send empties it.

Input: the agent tab's message box after a multi-line message was sent.
Before: the box grew in the browser on 'input' events and was
re-measured once, 60 ms after a click in the send row. The send empties
the box from Python, after a kernel round trip and without an 'input'
event, so on a slow connection the box kept its grown height until the
next keystroke. Measured in the real page with 300 ms latency: after
Enter on a six-line message the box was empty at 2 s and still 138 px.
After: the box carries the class ``delfin-agent-input-empty`` whenever
its value is empty, set by Python in the same message that empties it,
and the stylesheet holds an empty box at its resting 80 px. No timing in
the browser is involved, and the test needs no JavaScript runtime.
"""
from __future__ import annotations

import re


def _tab(tmp_path, monkeypatch):
    from delfin.agent import scheduler as S
    from delfin.dashboard import tab_agent
    from delfin.dashboard.context import DashboardContext

    S._GLOBAL = S.Scheduler(path=tmp_path / "cron.json")
    for name in ("calc", "archive", "office", "ws"):
        (tmp_path / name).mkdir(exist_ok=True)
    monkeypatch.setenv("DELFIN_AGENT_WORKSPACE", str(tmp_path / "ws"))
    monkeypatch.setenv("DELFIN_LAUNCH_CWD", str(tmp_path / "ws"))
    ctx = DashboardContext(calc_dir=tmp_path / "calc",
                           archive_dir=tmp_path / "archive",
                           office_dir=tmp_path / "office")
    ctx.run_js = lambda script, **kw: None
    return tab_agent.create_tab(ctx)


def _box(root):
    import ipywidgets as widgets
    stack = [root]
    while stack:
        w = stack.pop()
        if (isinstance(w, widgets.Textarea)
                and "delfin-agent-input" in w._dom_classes):
            return w
        stack.extend(getattr(w, "children", ()) or ())
    raise AssertionError("no message box")


def test_an_emptied_box_is_marked_in_the_same_message(tmp_path, monkeypatch):
    root, refs = _tab(tmp_path, monkeypatch)
    try:
        box = _box(root)
        assert "delfin-agent-input-empty" in box._dom_classes
        box.value = "line one\nline two\nline three"
        assert "delfin-agent-input-empty" not in box._dom_classes
        box.value = ""                      # what every send path does
        assert "delfin-agent-input-empty" in box._dom_classes
    finally:
        refs["shutdown"]()


def test_the_stylesheet_holds_an_empty_box_at_its_resting_height():
    from delfin.dashboard import tab_agent
    css = " ".join(tab_agent._AGENT_CSS.split())
    rule = re.search(r"\.delfin-agent-input-empty textarea \{([^}]*)\}", css)
    assert rule, "no rule for the empty box"
    assert "height: 80px !important" in rule.group(1)
