"""A person brings their running session to the phone by number or name.

Input: the phone commands ``sessions``, ``hol N|NAME`` (``get``), ``weg``
(``off``) in a session's chat, ``name NAME``, ``hilfe``; in the session
``/sync NAME``, ``/pair``. Through the REAL gateway and bridge
(tests/_messenger_people_rig.py). Before this, a session reached the phone
only by ``/sync`` typed in the session.

Output, each refusal with its control:

* ``sessions`` lists the person's OWN running sessions, numbered, named;
* ``hol N`` / ``hol NAME`` switches the session on for the messenger (the
  record says by the phone, which one), invites to its chat, tells the
  session and is audited;
* ``hol`` refuses another person's session, a session switched off in the
  session (``/sync off``), and anything from an unverified phone;
* ``weg`` takes it off again; then ``hol`` may bring it back;
* the German and English words are both understood (lang.COMMAND_WORDS);
* ``name NAME`` renames the session (the service's title) and its group;
* ``/sync NAME`` in a service session enables it and sets its title;
* ``/pair`` and ``delfin-agent messenger`` refuse an agent's process; the
  pairing directory is masked from agent file tools and shells.
"""
from __future__ import annotations

import fnmatch

import pytest

from _messenger_people_rig import Account
from delfin.agent import security_events
from delfin.messenger import activation, lang, pairing, store


@pytest.fixture
def acct(tmp_path, monkeypatch):
    security_events.clear()
    account = Account(tmp_path, monkeypatch)
    store.save("7")
    store.save("21", person="anna", label="anna phone")
    store.save("31", person="ben", label="ben phone")
    account.enable("ben", "s-ben")
    yield account
    security_events.clear()


def _short_door_dir(monkeypatch):
    """The pairing door in a short private directory of the test's own, not
    the user's runtime directory (Unix socket path limit)."""
    import shutil
    import tempfile
    path = tempfile.mkdtemp(prefix="dpair-", dir="/tmp")
    monkeypatch.setattr(pairing, "_door_dir", lambda directory=None: path)
    monkeypatch.setattr(pairing, "_OPEN", {})
    import atexit
    atexit.register(shutil.rmtree, path, True)


def _act(acct, person, sid):
    return activation.record(sid, directory=acct.activation_dir(person))


def test_sessions_lists_the_persons_own_running_sessions(acct):
    acct.say("21", "@21", "sessions")
    text = acct.daemon.to("@21")[-1]
    assert "1. Butan" in text and "2. Ethan" in text
    assert "Benzol" not in text and "Eigene" not in text


@pytest.mark.parametrize("words", [("sitzungen", "hol 2"), ("sessions", "get 2"), ("sessions", "hol Ethan"),
                                   ("sessions", "get ethan")])
def test_hol_brings_a_session_to_the_phone(acct, words):
    for text in words:
        acct.say("21", "@21", text)
    record = _act(acct, "anna", "s-anna2")
    assert record["enabled"] is True and record["by"] == "phone" and record["phone"] == "anna phone"
    assert not _act(acct, "anna", "s-anna").get("enabled")
    chat = store.room_of_session("s-anna2", "anna")
    assert chat and (chat, "21") in acct.daemon.invited
    assert {"sessionId": "s-anna2", "text": "/sync status"} in acct.services["anna"].calls("session.prompt")
    assert "messenger_session_attached" in [e.kind for e in security_events.recent()]
    assert "Ethan" in acct.daemon.to("@21")[-1]


def test_hol_refuses_another_persons_session(acct):
    acct.say("21", "@21", "hol Benzol")
    acct.say("21", "@21", "hol s-ben")
    assert "Keine deiner laufenden Sitzungen" in acct.daemon.to("@21")[-1]
    assert not acct.services["ben"].asked
    assert not store.room_of_session("s-ben", "anna") and not store.room_of_session("s-ben", "ben")


def test_hol_refuses_a_session_switched_off_in_the_session(acct):
    acct.enable("anna", "s-anna")
    acct.enable("anna", "s-anna", "off")                 # /sync off typed in the session
    acct.say("21", "@21", "hol Butan")
    assert "/sync off" in acct.daemon.to("@21")[-1]
    assert _act(acct, "anna", "s-anna")["enabled"] is False
    # the control: a session never switched brings fine
    acct.say("21", "@21", "hol Ethan")
    assert _act(acct, "anna", "s-anna2")["enabled"] is True


def test_hol_from_an_unverified_phone_does_nothing(acct):
    store.bind("25", person="anna", label="new", request="r" * 32, code="123456")
    acct.say("25", "@25", "hol Butan")
    assert not acct.services["anna"].asked
    assert not _act(acct, "anna", "s-anna")


def test_weg_takes_the_session_off_and_hol_brings_it_back(acct):
    acct.say("21", "@21", "hol Butan")
    chat = store.room_of_session("s-anna", "anna")
    acct.say("21", chat, "weg")
    record = _act(acct, "anna", "s-anna")
    assert record["enabled"] is False and record["by"] == "phone"
    assert "messenger_session_detached" in [e.kind for e in security_events.recent()]
    acct.say("21", "@21", "rechne", mid="x1")
    acct.say("21", chat, "noch da?")
    assert not [p for p in acct.services["anna"].calls("session.prompt")
                if p["text"].startswith("noch da")]
    acct.say("21", "@21", "hol Butan")                   # the phone's off: may come back
    assert _act(acct, "anna", "s-anna")["enabled"] is True


def test_name_renames_the_session_title_and_its_group(acct):
    acct.say("21", "@21", "hol Butan")
    chat = store.room_of_session("s-anna", "anna")
    acct.say("21", chat, "name Butan-Konformere")
    assert {"sessionId": "s-anna", "title": "Butan-Konformere"} in acct.services["anna"].calls("session.setTitle")
    assert acct.daemon.renamed[-1] == (chat, "DELFIN · Butan-Konformere")
    assert store.rooms()[chat]["title"] == "Butan-Konformere"


def test_a_name_is_trimmed_and_cut_to_fifty_characters():
    assert activation.clean_name("  a\x1b[31mb  ") == "a b"
    assert len(activation.clean_name("x" * 80)) == activation.NAME_MAX == 50


def test_both_languages_name_the_new_commands():
    pairs = {"hol": "get", "get": "get", "weg": "off", "off": "off",
             "entkoppeln": "revoke", "revoke": "revoke", "handys": "phones",
             "phones": "phones", "hilfe": "help", "help": "help",
             "sitzungen": "sessions", "sessions": "sessions"}
    assert {w: lang.command_word(w) for w in pairs} == pairs


def test_help_lists_the_commands_briefly(acct):
    acct.say("21", "@21", "hilfe")
    text = acct.daemon.to("@21")[-1]
    for word in ("sessions", "hol", "neu", "weg", "@NAME", "/sync NAME"):
        assert word in text


def test_the_phone_unpairs_itself(acct):
    store.save("22", person="anna", label="tablet")
    acct.say("21", "@21", "entkoppeln")
    assert store.contact("21") is None
    assert any("Entkoppelt" in t for t in acct.daemon.to("@21"))
    assert any("anna phone" in t for t in acct.daemon.to("@22"))
    acct.say("21", "@21", "sessions", mid="late")       # nothing reaches it any more
    assert not any("Butan" in t for t in acct.daemon.to("@21"))


def test_sync_name_in_a_service_session_sets_its_title(tmp_path, monkeypatch):
    from delfin.agent.service import worker
    monkeypatch.setattr(activation, "_DIR", tmp_path / "act")
    titles, notices = [], []
    host = type("W", (), {})()
    host.emitter = type("E", (), {"session_id": "s5",
                                  "emit": lambda self, kind, text=None, **k: notices.append(text)})()
    host._m_session_setTitle = lambda params: titles.append(params["title"]) or {"title": params["title"]}
    reply = worker.Worker._m_session_prompt(host, {"text": "/sync name Butan Konformere"})
    assert reply == {"accepted": True}
    assert activation.enabled("s5") and titles == ["Butan Konformere"]
    assert "Butan Konformere" in notices[-1]


def test_pair_in_a_service_session_goes_to_a_live_notice(tmp_path, monkeypatch):
    from delfin.agent.service import worker
    monkeypatch.setattr(pairing, "_DIR", tmp_path / "pairing")
    _short_door_dir(monkeypatch)
    monkeypatch.setattr(pairing, "_spawn", lambda fn: None)
    durable, live = [], []
    host = type("W", (), {})()
    host.emitter = type("E", (), {"session_id": "s6",
                                  "emit": lambda self, kind, text=None, **k: durable.append(text),
                                  "emit_live": lambda self, text: live.append(text)})()
    worker.Worker._m_session_prompt(host, {"text": "/pair"})
    assert pairing.open_request("session:s6")
    assert durable and "invitation" in durable[-1]
    assert list((tmp_path / "pairing").glob("req-*.json"))


def test_an_agent_process_cannot_pair(monkeypatch, tmp_path):
    from delfin.agent.approval_guard import ENV_MARKER
    monkeypatch.setattr(pairing, "_DIR", tmp_path / "pairing")
    _short_door_dir(monkeypatch)
    said = []
    monkeypatch.setenv(ENV_MARKER, "1")
    assert pairing.cli_main(["pair"], out=said.append, timeout=0) == 2
    assert "Refused" in said[-1] and not (tmp_path / "pairing").exists()
    monkeypatch.delenv(ENV_MARKER)
    assert pairing.cli_main(["pair"], out=said.append, timeout=0) == 1     # control: asks
    assert list((tmp_path / "pairing").glob("req-*.json"))


def test_the_pairing_directory_is_hidden_from_agents():
    from delfin.agent import api_client, boundaries, sandbox, state_paths
    assert ".delfin/messenger_pairing" in sandbox._HOME_SECRET_DIRS
    assert ".delfin/messenger_pairing" in boundaries.MASKED_HOME_ENTRIES
    assert any(fnmatch.fnmatch("messenger_pairing/req-x.json", g)
               for g in api_client._DEFAULT_PATH_DENY_GLOBS)
    assert ("delfin.messenger.pairing", "_DIR", "messenger_pairing") in state_paths.USER_STATE_SINKS
