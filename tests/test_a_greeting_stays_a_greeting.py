"""A first "hallo" in the dashboard cost 37,304 prompt tokens where the
same word through the engine alone cost 11,101 (measured 2026-09-11).
The dashboard glued its session primer in front of the user's text, so
the engine's bare-greeting rule -- which withholds the tool schemas from
a message that is nothing but a greeting -- never saw a bare greeting.
The primer now travels as memory context, which the engine places in
the system prompt; the user's text stays the user's text.
"""

from __future__ import annotations

import pathlib

_TAB = (pathlib.Path(__file__).resolve().parents[1]
        / "delfin" / "dashboard" / "tab_agent.py")


def test_the_engine_judges_the_bare_text_and_a_primer_would_hide_it():
    from delfin.agent.engine import AgentEngine
    assert AgentEngine.is_bare_greeting("hallo") is True
    assert AgentEngine.is_bare_greeting("# Session boot\nrecent outcomes: none\n\nhallo") is False


def test_the_primer_is_no_longer_glued_to_the_users_text():
    text = _TAB.read_text(encoding="utf-8")
    assert 'current_msg = f"{_boot}\\n\\n{current_msg}"' not in text
    i = text.index("_boot = _build_dashboard_session_boot()")
    window = text[i:i + 1200]
    assert window.count('state["_pending_boot_brief"] = _boot') == 2, (
        "both the dashboard and the solo primer travel as memory context")


def test_the_primer_reaches_the_engine_as_memory_context():
    text = _TAB.read_text(encoding="utf-8")
    i = text.index('_boot_brief = str(state.pop("_pending_boot_brief", "") or "")')
    after = text[i:i + 400]
    assert "_memory = (_boot_brief" in after
    # and the memory block is what every stream_response call hands over
    assert text.count("memory_context=_memory,") >= 4
