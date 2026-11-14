"""A phone's text is never a window's command: no /sync, no /pair, no slash or shell line.

Review finding 3: ``@NAME /sync off`` and ``@NAME /pair`` in the forum
reached the service worker as a prompt, and the worker ran its /sync and
/pair handlers (and any REPL builtin) for whoever sent it. A phone could
switch a session or start a pairing.

Output, each with its control:

* the gateway refuses slash-led text after ``@NAME`` (control: plain text
  reaches the session);
* the agent service marks a prompt that comes from the messenger bridge
  (``via: messenger``, set by the service from the client's name and
  stripped from anyone else's request);
* the worker refuses a slash- or shell-led prompt that came via the
  messenger -- /pair, /sync, any builtin, ``!cmd`` -- except the read-only
  ``/sync status`` the bridge shows the session; the same lines from a
  window work (control);
* the dashboard adapter types no slash or shell line from the phone;
* through the real gateway and bridge, a verified phone never starts or
  confirms a pairing: no request appears in the person's pairing directory.
"""
from __future__ import annotations

import asyncio
from types import SimpleNamespace

import pytest

from _messenger_people_rig import Account
from delfin.messenger import activation, pairing, store


@pytest.fixture
def acct(tmp_path, monkeypatch):
    account = Account(tmp_path, monkeypatch)
    store.save("21", person="anna", label="anna phone")
    account.enable("anna", "s-anna")
    store.put_room("#a1", "s-anna", "Butan", person="anna", ready=True)
    store.put_room("#af", None, "Forum", person="anna", ready=True)
    return account


def _prompts(acct):
    return [p["text"] for p in acct.services["anna"].calls("session.prompt")]


@pytest.mark.parametrize("text", ["@Butan /sync off", "@Butan /pair", "@Butan /pair 123456",
                                  "@Butan /mode bypass", "@Butan !ls"])
def test_the_gateway_refuses_a_command_after_at_name(acct, text):
    acct.say("21", "#af", text)
    assert not _prompts(acct)
    acct.say("21", "#af", "@Butan hallo")             # control
    assert _prompts(acct) and _prompts(acct)[-1].startswith("hallo")


@pytest.mark.parametrize("text", ["/pair", "/pair 123456", "pair", "@Butan /pair"])
def test_a_verified_phone_never_starts_or_confirms_a_pairing(acct, text):
    acct.say("21", "#a1" if not text.startswith("@") else "#af", text)
    acct.say("21", "@21", text)
    assert not list(acct.pairing_dir("anna").glob("req-*.json")) if acct.pairing_dir("anna").exists() else True
    assert not any(t.lstrip().startswith("/pair") for t in _prompts(acct))


def test_the_service_marks_prompts_from_the_messenger_bridge():
    from delfin.agent.service import supervisor as S
    sup = S.Supervisor.__new__(S.Supervisor)
    sup._waking, sup.dormant = {}, {}
    seen = []

    class Session:
        async def request(self, method, params, timeout=None):
            seen.append(params)
            return {"accepted": True}

    async def awake(client, sid):
        return Session()

    sup._awake = awake
    sup.sessions = {"s1": Session()}
    bridge = SimpleNamespace(name="delfin-messenger-bridge")
    tui = SimpleNamespace(name="delfin-tui")
    asyncio.run(sup._forward(bridge, "session.prompt", {"sessionId": "s1", "text": "/pair"}))
    asyncio.run(sup._forward(tui, "session.prompt", {"sessionId": "s1", "text": "/pair",
                                                     "via": "messenger"}))
    asyncio.run(sup._forward(tui, "session.prompt", {"sessionId": "s1", "text": "x",
                                                     "via": "window"}))
    assert seen[0]["via"] == "messenger"
    assert "via" not in seen[1] and "via" not in seen[2]


def _worker(tmp_path, monkeypatch, sid="s7"):
    from delfin.agent.service import worker
    monkeypatch.setattr(activation, "_DIR", tmp_path / "act")
    monkeypatch.setattr(pairing, "_DIR", tmp_path / "pairing")
    monkeypatch.setattr(pairing, "_spawn", lambda fn: None)
    monkeypatch.setattr(pairing, "_door_dir", lambda directory=None: tmp_path)
    monkeypatch.setattr(pairing, "_unavailable", lambda: "")   # a messenger runs on this host
    notices = []
    host = type("W", (), {})()
    host.emitter = type("E", (), {"session_id": sid,
                                  "emit": lambda self, kind, text=None, **k: notices.append(text),
                                  "emit_live": lambda self, text: notices.append(text)})()
    host._m_session_setTitle = lambda params: {"title": params["title"]}
    return worker, host, notices


@pytest.mark.parametrize("text", ["/pair", "/pair 123456", "/sync off", "/sync", "/sync name X",
                                  "/mode bypass", "!ls"])
def test_the_worker_refuses_a_command_that_came_via_the_messenger(tmp_path, monkeypatch, text):
    worker, host, _ = _worker(tmp_path, monkeypatch)
    activation.configure("s7")
    with pytest.raises(worker._Refusal):
        worker.Worker._m_session_prompt(host, {"text": text, "via": "messenger"})
    assert activation.enabled("s7")
    assert not list((tmp_path / "pairing").glob("req-*.json")) if (tmp_path / "pairing").exists() else True


def test_the_read_only_status_from_the_bridge_is_shown(tmp_path, monkeypatch):
    worker, host, notices = _worker(tmp_path, monkeypatch)
    activation.configure("s7")
    assert worker.Worker._m_session_prompt(host, {"text": "/sync status", "via": "messenger"}) == {"accepted": True}
    assert "on" in notices[-1] and activation.enabled("s7")


def test_the_same_lines_from_a_window_work_control(tmp_path, monkeypatch):
    worker, host, _ = _worker(tmp_path, monkeypatch)
    worker.Worker._m_session_prompt(host, {"text": "/sync"})
    assert activation.enabled("s7")
    worker.Worker._m_session_prompt(host, {"text": "/pair"})
    assert list((tmp_path / "pairing").glob("req-*.json"))
    pairing.forget("session:s7")


@pytest.mark.parametrize("text", ["/pair", "/sync off", "!ls"])
def test_the_dashboard_adapter_types_no_command_from_the_phone(text):
    from delfin.messenger.dashboard import DashboardSession
    typed = []
    ui = DashboardSession("one", {"active_session_id": "one"}, send=lambda t: typed.append(t) or True,
                          stop=lambda: None, answer=lambda t: False, plan=lambda yes: None)
    reply = ui({"method": "session.prompt", "text": text})
    assert not typed and not reply["accepted"]
    ui({"method": "session.prompt", "text": "rechne das"})       # control
    assert typed == ["rechne das"]
