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


def test_the_gateways_own_timeout_is_named_not_pasted():
    """The same failed turn ended with "Error: Error code: 400 - {'detail':
    'Open WebUI: Server Connection Error'}" -- the gateway's wait for the
    model ran out. The chat says that, and what helps."""
    from delfin.dashboard.tab_agent import _gateway_gave_up
    assert _gateway_gave_up("Error code: 400 - {'detail': 'Open WebUI: Server Connection Error'}")
    assert not _gateway_gave_up("Error code: 401 - invalid key")
    assert not _gateway_gave_up("")
    i = SRC.index("elif _gateway_gave_up(error_text):")
    body = SRC[i:i + 1600]
    assert "could not reach" in body and "`/retry`" in body and "kit.deepseek-v4-flash" in body
