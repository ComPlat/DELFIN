"""The agent knows when its session is connected to the user's phone.

Field report (Pi): after /sync the user wrote "schreib mir hallo auf
simplex"; the model searched calculation folders, offered image generation
and finally said it could not reach SimpleX. Nothing in its context said
the session was connected.

Output, through the real engine and prompt loader (a client that records
what it is sent):

* while the session is enabled for the messenger, the turn's tail carries
  one short, stable line: connected to the user's phone through SimpleX,
  replies appear there, offered files can be fetched there;
* with /sync off the line is gone;
* the line never enters the system prompt: its bytes are the first turn's
  across on, off and on again (the prefix cache keeps the history).
"""
from __future__ import annotations

import copy
from unittest.mock import MagicMock, patch

import pytest

from delfin.agent import api_client as A
from delfin.agent import engine as E
from delfin.agent.api_client import StreamEvent
from delfin.messenger import activation


@pytest.fixture
def eng(tmp_path, monkeypatch):
    monkeypatch.setattr(activation, "_DIR", tmp_path / "act")
    ws = tmp_path / "ws"
    ws.mkdir()
    with patch("delfin.agent.engine.create_client", return_value=MagicMock()):
        e = E.AgentEngine(repo_dir=ws, backend="api", provider="kit",
                          model="kit.glm-5.3", mode="solo")
    e.client._permissions = A.KitToolPermissions(workspace=ws)
    e.loader.workspace_root = ws
    e.session_id = "s-phone"
    seen = []

    def _gen(**kw):
        seen.append({"system": kw.get("system"), "messages": copy.deepcopy(kw.get("messages"))})
        yield StreamEvent(type="text_delta", text="ok")
        yield StreamEvent(type="message_delta", input_tokens=5, output_tokens=7, cost_usd=0.0)

    e.client.stream_message = _gen
    e._seen = seen
    return e


def _tail(request):
    last = next(m for m in reversed(request["messages"]) if m.get("role") == "user")
    return str(last["content"])


def test_the_line_is_there_while_on_and_gone_when_off(eng):
    eng.stream_response("Hallo")
    activation.configure("s-phone")
    eng.stream_response("schreib mir hallo auf simplex")
    activation.configure("s-phone", "off")
    eng.stream_response("Weiter")
    activation.configure("s-phone")
    eng.stream_response("Und jetzt?")
    first, on, off, again = eng._seen
    assert E.MESSENGER_LINE not in _tail(first)
    assert E.MESSENGER_LINE in _tail(on) and E.MESSENGER_LINE in _tail(again)
    assert E.MESSENGER_LINE not in _tail(off)
    for req in (on, off, again):
        assert req["system"] == first["system"]
        assert E.MESSENGER_LINE not in req["system"]


def test_the_line_says_what_the_agent_needs():
    line = E.MESSENGER_LINE
    assert "SimpleX" in line and "phone" in line
    assert "appear" in line and "file" in line
    assert len(line) < 400 and "\n" not in line
