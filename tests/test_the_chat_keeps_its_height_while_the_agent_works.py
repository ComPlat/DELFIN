"""The chat keeps its height while the agent works.

Input: the chat frame, a column of the transcript, the working line (one
HTML widget that also carries the slow-phase hint), the task line, the
Background box and the docks. Measured in the real page (2026-10-22,
1920x1080): the transcript's viewport was 737 px idle, 691 px once the
working line appeared, 657 px when the slow-phase hint joined it after
4 s, 691 px again at the first token and 737 px at the end. Each change
moved the whole transcript of a reader at the end by that much.

After: the frame is a grid whose first cell holds the transcript and,
over its lower edge, the working line; the task line, the Background box
and the docks follow in rows below. The transcript keeps a fixed space
below its last message, so a reader at the end sees the working line
there and nothing under it. The "Newest" control is lifted above the
working line while it is shown.
"""
from __future__ import annotations

import re


def _css() -> str:
    from delfin.dashboard import tab_agent
    return " ".join(re.sub(r"/\*.*?\*/", "", tab_agent._AGENT_CSS,
                           flags=re.S).split())


def _rule(selector: str) -> str:
    out = ""
    for head, body in re.findall(r"([^{}]+)\{([^{}]*)\}", _css()):
        if selector in [s.strip() for s in head.split(",")]:
            out += body + ";"
    return out


def test_the_working_line_overlays_the_transcript_in_one_cell():
    frame = _rule(".delfin-agent-root .delfin-agent-chat-frame")
    assert "display: grid !important" in frame
    assert "grid-template-rows: minmax(0, 1fr)" in frame
    host = _rule(".delfin-agent-chat-frame > .delfin-agent-chat-host")
    assert "grid-area: 1 / 1" in host
    slot = _rule(".delfin-agent-chat-frame > .delfin-agent-working-slot")
    assert "grid-area: 1 / 1" in slot and "align-self: end" in slot


def test_the_transcript_keeps_room_below_its_last_message():
    room = _rule(".delfin-agent-chat-host .delfin-agent-chat")
    assert re.search(r"padding-bottom: \d+px", room), room


def test_the_working_line_widget_carries_the_slot_class(tmp_path,
                                                        monkeypatch):
    import ipywidgets as widgets
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
    root, refs = tab_agent.create_tab(ctx)
    try:
        stack, frame = [root], None
        while stack:
            w = stack.pop()
            if "delfin-agent-chat-frame" in getattr(w, "_dom_classes", ()):
                frame = w
                break
            stack.extend(getattr(w, "children", ()) or ())
        assert frame is not None
        assert isinstance(frame.children[1], widgets.HTML)
        assert "delfin-agent-working-slot" in frame.children[1]._dom_classes
    finally:
        refs["shutdown"]()


def test_the_newest_control_is_lifted_over_the_working_line():
    lifted = [body for head, body in
              re.findall(r"([^{}]+)\{([^{}]*)\}", _css())
              if ".delfin-chat-jump" in head
              and "delfin-agent-working-slot" in head]
    assert any("bottom:" in b for b in lifted), lifted
