"""A phone's prompt reaches the session exactly as the person wrote it.

Field report (TUI): the bridge appended a block of German instructions to
EVERY prompt from the phone ("[Messenger: Dateien/Bilder nur auf ...
forum.txt ... max. 1400 UTF-8-Bytes ...]"); the TUI showed it as part of
the person's message and queued it with it.

Output, each with its control:

* through the real gateway and bridge, the text the agent service gets is
  the phone's text, byte for byte;
* the worker keeps that text as the line and announces it with
  ``via: "phone"`` (a badge in the TUI and the dashboard), never a suffix;
  a window's line carries no ``via``;
* the file-offer and forum rules are said once per turn in the turn's
  volatile messenger context line, in English, while the session is on
  the phone -- with the session's own outbox path; the person's message
  in the request is the text alone;
* the dashboard's live view shows the badge, not a suffix.
"""
from __future__ import annotations

import copy
import threading
from types import SimpleNamespace
from unittest.mock import MagicMock, patch

import pytest

from _messenger_people_rig import Account
from delfin.messenger import activation, store

PHONE_TEXT = "rechne das: E(S1) – bitte  mit 2 Kernen"


@pytest.fixture
def acct(tmp_path, monkeypatch):
    account = Account(tmp_path, monkeypatch)
    store.save("21", person="anna", label="anna phone")
    account.enable("anna", "s-anna")
    store.put_room("#a1", "s-anna", "Butan", person="anna", ready=True)
    return account


def test_the_service_gets_the_phone_text_byte_for_byte(acct):
    acct.say("21", "#a1", PHONE_TEXT)
    prompts = [p["text"] for p in acct.services["anna"].calls("session.prompt")]
    assert prompts == [PHONE_TEXT]


def _host():
    events = []
    agent = SimpleNamespace(cond=threading.Condition(), closing=threading.Event(), busy=False, idle=True,
                            inbox=[], queued=[], wake=threading.Event(), _turn_active=threading.Event(),
                            post=lambda fn: fn(), steer_text=lambda text: False)
    host = SimpleNamespace(agent=agent, emitter=SimpleNamespace(
        session_id="s7", emit=lambda kind, **fields: events.append((kind, fields))))
    return host, events


def test_the_worker_marks_a_phone_line_by_its_origin_only():
    from delfin.agent.service import worker
    host, events = _host()
    worker.Worker._m_session_prompt(host, {"text": PHONE_TEXT, "via": "messenger"})
    assert host.agent.inbox == [PHONE_TEXT]
    kind, fields = events[-1]
    assert kind == "user.message" and fields["text"] == PHONE_TEXT and fields.get("via") == "phone"
    control, events2 = _host()
    worker.Worker._m_session_prompt(control, {"text": "from a window"})
    assert "via" not in events2[-1][1]


def test_the_protocol_names_the_origin():
    from delfin.agent.service import protocol
    assert protocol.EVENTS["user.message"].get("via") == "string?"


@pytest.fixture
def eng(tmp_path, monkeypatch):
    from delfin.agent import api_client as A
    from delfin.agent import engine as E
    from delfin.agent.api_client import StreamEvent
    monkeypatch.setattr(activation, "_DIR", tmp_path / "act")
    ws = tmp_path / "ws"
    ws.mkdir()
    with patch("delfin.agent.engine.create_client", return_value=MagicMock()):
        e = E.AgentEngine(repo_dir=ws, backend="api", provider="kit", model="kit.glm-5.3", mode="solo")
    e.client._permissions = A.KitToolPermissions(workspace=ws)
    e.loader.workspace_root = ws
    e.session_id = "s-phone"
    seen = []

    def _gen(**kw):
        seen.append(copy.deepcopy(kw.get("messages")))
        yield StreamEvent(type="text_delta", text="ok")
        yield StreamEvent(type="message_delta", input_tokens=5, output_tokens=7, cost_usd=0.0)

    e.client.stream_message = _gen
    e._seen, e._ws = seen, ws
    return e


def test_the_rules_ride_once_in_the_messenger_line(eng):
    from delfin.messenger import media
    activation.configure("s-phone")
    eng.stream_response(PHONE_TEXT)
    last = next(m for m in reversed(eng._seen[-1]) if m.get("role") == "user")
    content = str(last["content"])
    outbox = str(eng._ws / media.ROOT / media.session_hash("s-phone") / "outbox")
    assert content.count("forum.txt") == 1 and content.count(outbox) == 1
    assert "Dateien" not in content and "Keine" not in content
    assert content.startswith(PHONE_TEXT)


def test_the_live_view_shows_a_badge_not_a_suffix():
    from delfin.dashboard import live_tab
    T = SimpleNamespace(_md_to_html=lambda text: text)
    html = live_tab._render_entry(T, {"role": "user", "text": PHONE_TEXT, "via": "phone"})
    assert PHONE_TEXT in html and "via phone" in html
    assert "via phone" not in live_tab._render_entry(T, {"role": "user", "text": "x"})
