"""/retry works after a turn that never started.

2026-09-16, GLM behind a stalling gateway: the turn ended after 260 s with
"The request never started streaming ... /retry sends the same message
again", and /retry answered "Nothing to retry." The engine had dropped the
user message (no text came, so nothing was appended), and the command
refused on an empty engine history although the chat still held the
message it re-sends from.
"""
import inspect
from pathlib import Path

SRC = Path(inspect.getfile(__import__("delfin.dashboard.tab_agent", fromlist=["x"]))).read_text()


def test_an_empty_engine_history_is_not_nothing_to_retry():
    i = SRC.index("def _retry_last(")
    body = SRC[i:i + 2500]
    assert "if not engine:" in body
    assert "if not engine or not engine.messages:" not in body
    assert 'state["chat_messages"][-1]["role"] == "user"' in body


def test_the_notice_still_promises_retry():
    assert "/retry sends the same message" in SRC
