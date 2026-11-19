"""Enter in the agent tab sends what is in the box; Ctrl+L clears nothing.

Measured 11.11. in a real browser: typing fast (or pasting) and pressing
Enter at once sent "g" of "gamma three" -- the click on Send read the
kernel's copy of the box, which lags behind ipywidgets' throttled value
updates. And Ctrl+L sent /clear (a fresh conversation, the draft lost) from
any focus. Convention: a chat sends what the box shows; neither the TUI nor
Claude Code clears a session on Ctrl+L, and the browser keeps its own.
"""
from pathlib import Path

SRC = (Path(__file__).resolve().parent.parent / "delfin" / "dashboard"
       / "tab_agent.py").read_text(encoding="utf-8")


def test_enter_carries_the_text_with_the_request():
    assert "window.__delfinAgentSubmit" in SRC
    assert "submit_channel.add_class(\"delfin-agent-submit\")" in SRC
    assert "submit_channel.observe(_on_submit_channel, names=\"value\")" in SRC
    # both key handlers use it, neither clicks Send on Enter any more
    assert SRC.count("window.__delfinAgentSubmit(e.target);") == 2


def test_ctrl_l_no_longer_clears_the_chat():
    assert "e.key === 'l'" not in SRC
    assert "Ctrl+L           \\u2014 Clear chat" not in SRC
