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
* the answer file carries the link sealed to a key that lives only in the
  memory of the process that asked (finite-field Diffie-Hellman, RFC 3526
  group 14); that process reads it, another cannot;
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

SECURITY = "12345 67890 11111"


@pytest.fixture
def acct(tmp_path, monkeypatch):
    monkeypatch.setattr(audit_log, "_default_log_path", lambda: tmp_path / "audit.log")
    monkeypatch.setattr(session_store, "_SESSIONS_DIR", tmp_path / "sessions")
    monkeypatch.setattr(outcome_tracker, "_DEFAULT_PATH", tmp_path / "outcome_history.jsonl")
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


def _code(acct, cid):
    return acct.daemon.to("@" + cid)[0].split("Code ")[1].split(".")[0]


def test_no_file_holds_the_link_during_the_window_or_after(acct):
    rid = pairing.begin("session:p1", directory=acct.pairing_dir("anna"))
    acct.round()
    link = acct.daemon.links[0]
    assert pairing.read(rid, "session:p1", directory=acct.pairing_dir("anna"))["link"] == link
    assert _leaks(acct.tmp, link) == []
    acct.connect("21", "901")
    shown = pairing.read(rid, "session:p1", directory=acct.pairing_dir("anna"))
    assert shown["state"] == "connected" and shown["securityCode"] == SECURITY
    assert _leaks(acct.tmp, SECURITY.replace(" ", "")) == [] and _leaks(acct.tmp, SECURITY) == []
    pairing.confirm("session:p1", _code(acct, "21"), directory=acct.pairing_dir("anna"))
    acct.round()
    assert store.verified_person("21") == "anna"
    assert _leaks(acct.tmp, link) == []


def test_only_the_asking_process_can_open_the_link(acct):
    rid = pairing.begin("session:p1", directory=acct.pairing_dir("anna"))
    acct.round()
    raw = pairing.result(rid, directory=acct.pairing_dir("anna"))
    assert raw["state"] == "offered" and "link" not in raw and raw.get("sealed")
    assert pairing.read(rid, "session:other", directory=acct.pairing_dir("anna")).get("link") is None
    tampered = dict(raw["sealed"], ct=raw["sealed"]["ct"][:-4] + "AAAA")
    assert pairing.unseal_with(pairing._OPEN["session:p1"][2], tampered) is None


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
    assert _leaks(acct.tmp, link) == []


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


def test_the_group_is_rfc_3526_group_14():
    """A mistyped prime makes every seal fail (measured: one doubled digit
    raised OverflowError on the first answer)."""
    import hashlib
    p = pairing._MODP
    assert p.bit_length() == 2048
    assert hashlib.sha256(format(p, "X").encode()).hexdigest()[:16] == "dcd8538e629d7b8b"
    q = (p - 1) // 2                       # a safe prime: q prime too; 2 generates q's subgroup
    assert pow(3, p - 1, p) == 1 and pow(3, q - 1, q) == 1 and pow(2, q, p) == 1
