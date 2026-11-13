"""The one-time invitation link is never written to any file.

Input: a person's /pair through the REAL bridge and gateway
(tests/_messenger_people_rig.py), shown in the dashboard's pairing panel
(delfin.messenger.dashboard.PairPanel). Before this, the link stood in
plain text in the person's 0600 answer file while the invitation was open,
and the dashboard showed it in a chat line that a save could keep.

Output:

* during the invitation window and after verification, no file under the
  test's tree -- answer files, gateway store, audit log, saved session,
  session export, outcome history -- holds the link or any part of it
  longer than the 12-character tag; neither does the SimpleX security code;
* the link and the security code never reach a file at all: the bridge
  hands them to the asking process over that process's own pairing door (a
  Unix socket in its private runtime directory, peer uid and agent-process
  checked on both ends), and they stay in that process's memory; the
  answer file holds the state only;
* only the asking process gets them: another requester reads no link; a
  request naming a door outside the person's runtime directory, or another
  person's, gets nothing delivered and its invitation is withdrawn; a
  forged delivery to a door (wrong request id) is ignored;
* the dashboard shows the link in a panel outside the chat transcript and
  the saved session, and removes it on success, expiry and cancel; the
  chat gets only lines without it.
"""
from __future__ import annotations

import os

import pytest

from _messenger_people_rig import Account
from delfin.agent import audit_log, outcome_tracker, session_store
from delfin.messenger import pairing, store
from delfin.messenger.dashboard import PairPanel


@pytest.fixture(autouse=True)
def _german_phone(monkeypatch):
    """These tests read the German phone texts (lang.DE); English is the
    default (test_the_phone_speaks_the_persons_language)."""
    monkeypatch.setenv("DELFIN_MESSENGER_LANG", "de")

SECURITY = "12345 67890 11111"


@pytest.fixture
def acct(tmp_path, monkeypatch):
    monkeypatch.setattr(audit_log, "_default_log_path", lambda: tmp_path / "audit.log")
    monkeypatch.setattr(session_store, "_SESSIONS_DIR", tmp_path / "sessions")
    monkeypatch.setattr(outcome_tracker, "_DEFAULT_PATH", tmp_path / "outcome_history.jsonl")
    monkeypatch.setattr(pairing, "_unavailable", lambda: "")   # a messenger runs on this host
    account = Account(tmp_path, monkeypatch)
    account.tmp = tmp_path
    return account


def _leaks(root, secret):
    """Files under ``root`` holding ``secret`` or any 13-character part of it."""
    parts = {secret[i:i + 13] for i in range(max(1, len(secret) - 12))}
    found = []
    for base, _dirs, files in os.walk(root):
        for name in files:
            path = os.path.join(base, name)
            try:
                data = open(path, "rb").read().decode("utf-8", "replace")
            except OSError:
                continue
            if any(p in data for p in parts):
                found.append(path)
    return found


def _all_leaks(acct, secret):
    found = _leaks(acct.tmp, secret)
    for root in acct.runtime.values():
        found += _leaks(root, secret)
    return found


def _code(acct, cid):
    return acct.daemon.to("@" + cid)[0].split("Code ")[1].split(".")[0]


def test_no_file_holds_the_link_during_the_window_or_after(acct):
    rid = pairing.begin("session:p1", directory=acct.pairing_dir("anna"))
    acct.round()
    link = acct.daemon.links[0]
    assert pairing.read(rid, "session:p1", directory=acct.pairing_dir("anna"))["link"] == link
    assert _all_leaks(acct, link) == []
    acct.connect("21", "901")
    shown = pairing.read(rid, "session:p1", directory=acct.pairing_dir("anna"))
    assert shown["state"] == "connected" and shown["securityCode"] == SECURITY
    assert _all_leaks(acct, SECURITY.replace(" ", "")) == [] and _all_leaks(acct, SECURITY) == []
    pairing.confirm("session:p1", _code(acct, "21"), directory=acct.pairing_dir("anna"))
    acct.round()
    assert store.verified_person("21") == "anna"
    assert _all_leaks(acct, link) == []


def test_only_the_asking_process_gets_the_link(acct):
    rid = pairing.begin("session:p1", directory=acct.pairing_dir("anna"))
    acct.round()
    raw = pairing.result(rid, directory=acct.pairing_dir("anna"))
    assert raw["state"] == "offered" and set(raw) <= {"kind", "state", "expires"}
    assert pairing.read(rid, "session:other", directory=acct.pairing_dir("anna")).get("link") is None
    assert pairing.read(rid, "session:p1", directory=acct.pairing_dir("anna"))["link"] == acct.daemon.links[0]


@pytest.mark.parametrize("where", ["elsewhere", "other person"])
def test_a_request_naming_a_door_outside_the_persons_runtime_dir_gets_nothing(acct, tmp_path, where):
    import socket
    import tempfile
    root = (tempfile.mkdtemp(prefix="dpairx-", dir="/tmp") if where == "elsewhere"
            else acct.runtime["ben"])
    path = os.path.join(str(root), "pair-forged.sock")
    listener = socket.socket(socket.AF_UNIX, socket.SOCK_STREAM)
    listener.bind(path)
    listener.listen(1)
    listener.settimeout(0.5)
    rid = pairing.submit("pair", {"requester": "x", "nonce_sha": "a" * 64, "door": path},
                         directory=acct.pairing_dir("anna"))
    acct.round()
    try:
        listener.accept()
        reached = True
    except OSError:
        reached = False
    listener.close()
    assert not reached
    res = pairing.result(rid, directory=acct.pairing_dir("anna"))
    assert res["state"] == "refused" and "link" not in res
    assert acct.daemon.dropped                       # the invitation withdrawn
    assert not store.invitations()


def test_a_forged_delivery_to_a_door_is_ignored(acct):
    import json
    import socket
    rid = pairing.begin("session:p1", directory=acct.pairing_dir("anna"))
    door = pairing._OPEN["session:p1"][2]
    with socket.socket(socket.AF_UNIX, socket.SOCK_STREAM) as sock:
        sock.connect(door.path)
        sock.sendall((json.dumps({"id": "f" * 32, "link": "simplex:/forged"}) + "\n").encode())
        assert b'"ok": false' in sock.recv(100)
    assert door.values == {}
    assert (pairing.read(rid, "session:p1", directory=acct.pairing_dir("anna")) or {}).get("link") is None


def test_a_door_refuses_an_agent_process(acct, monkeypatch):
    rid = pairing.begin("session:p1", directory=acct.pairing_dir("anna"))
    monkeypatch.setattr(pairing, "_peer_is_agent", lambda sock: True)
    acct.round()
    assert pairing.read(rid, "session:p1", directory=acct.pairing_dir("anna")).get("link") is None
    assert pairing.result(rid, directory=acct.pairing_dir("anna"))["state"] == "refused"


def test_the_door_is_private_and_closed_when_the_pairing_ends(acct):
    import stat
    rid = pairing.begin("session:p1", directory=acct.pairing_dir("anna"))
    door = pairing._OPEN["session:p1"][2]
    assert stat.S_IMODE(os.stat(door.path).st_mode) == 0o600
    assert os.path.dirname(door.path) == str(acct.runtime["anna"])
    pairing.forget("session:p1")
    assert not os.path.exists(door.path)
    del rid


def test_the_dashboard_panel_shows_the_link_outside_the_transcript(acct):
    rendered, chat = [], []
    panel = PairPanel(render=rendered.append, say=chat.append)
    reply = pairing.session_command("", requester="dashboard:d1", show=panel.show,
                                    spawn=lambda fn: None, directory=acct.pairing_dir("anna"))
    chat.append(reply)
    acct.round()
    rid = pairing.open_request("dashboard:d1")
    panel.show(pairing.offer_text(pairing.read(rid, "dashboard:d1", directory=acct.pairing_dir("anna"))))
    link = acct.daemon.links[0]
    import html
    assert link in html.unescape(rendered[-1])
    acct.connect("21", "901")
    panel.show(pairing.state_text(pairing.read(rid, "dashboard:d1", directory=acct.pairing_dir("anna"))))
    assert SECURITY in rendered[-1]
    pairing.confirm("dashboard:d1", _code(acct, "21"), directory=acct.pairing_dir("anna"))
    acct.round()
    panel.show(pairing.state_text(pairing.read(rid, "dashboard:d1", directory=acct.pairing_dir("anna"))))
    assert rendered[-1] == ""                                  # removed on success
    assert all(link[:13] not in line and SECURITY not in line for line in chat)
    # what the dashboard saves and exports is its chat
    session_store.save_session("d1", chat_messages=[{"role": "system", "content": c} for c in chat])
    from delfin.agent.session_export import export_session
    exported = str(export_session(session_store.load_session("d1")))
    assert link[:13] not in exported
    assert _all_leaks(acct, link) == []


@pytest.mark.parametrize("end", ["expired", "cancel"])
def test_the_panel_is_removed_on_expiry_and_cancel(acct, end):
    rendered, chat = [], []
    panel = PairPanel(render=rendered.append, say=chat.append, ttl_s=0.0)
    panel.show(pairing.offer_text({"link": "simplex:/invitation#/?v=2&smp=SECRETKEYMATERIAL0009",
                                   "expires": 0}))
    assert "SECRETKEY" in rendered[-1]
    if end == "cancel":
        panel.clear()
    else:
        panel.expire_if_due()
    assert rendered[-1] == ""

