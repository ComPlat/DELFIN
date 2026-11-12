"""Every person pairs their own phone, without admin rights, and verifies it.

Input: a person's ``/pair`` request in their own ``~/.delfin/messenger_pairing``
(delfin.messenger.pairing), collected by the REAL bridge, answered by the
REAL gateway (Rooms.pairing_round / contact_connected) with a fake SimpleX
daemon (tests/_messenger_people_rig.py). Before this, only an admin could
pair, and only one phone per installation.

Output, each attack with its control:

* the invitation is one-time, 0600, bound to the person and the request;
  the link never reaches the audit log, the panel or stderr (a 12-character
  tag does);
* a connected phone is bound UNVERIFIED and is told only to verify first --
  a message in that window commands nothing;
* only the process that started /pair can confirm (the nonce); a verify
  forged with another nonce is refused;
* the right code verifies; the person's other phones are told;
* a replayed invitation (a second phone through the same link), an expired
  one and a replayed request id bind nothing and are audited;
* the phone limit and the invitation rate refuse with a reason;
* every event kind the messenger records is registered with the panel.
"""
from __future__ import annotations

import asyncio
import os
import shutil
import stat
import time

import pytest

from _messenger_people_rig import Account
from delfin.agent import audit_log, security_events
from delfin.messenger import audit, pairing, store


@pytest.fixture(autouse=True)
def _german_phone(monkeypatch):
    """These tests read the German phone texts (lang.DE); English is the
    default (test_the_phone_speaks_the_persons_language)."""
    monkeypatch.setenv("DELFIN_MESSENGER_LANG", "de")


@pytest.fixture
def acct(tmp_path, monkeypatch):
    security_events.clear()
    log = tmp_path / "audit.log"
    monkeypatch.setattr(audit_log, "_default_log_path", lambda: log)
    account = Account(tmp_path, monkeypatch)
    account.log = log
    yield account
    security_events.clear()


def _pair(acct, person="anna", requester="session:s-anna"):
    rid = pairing.begin(requester, directory=acct.pairing_dir(person))
    acct.round()
    return rid, pairing.read(rid, requester, directory=acct.pairing_dir(person))


def _code(acct, cid):
    text = acct.daemon.to("@" + cid)[0]
    return text.split("Code ")[1].split(".")[0]


def _kinds():
    return [e.kind for e in security_events.recent(200)]


def test_a_person_gets_a_one_time_invitation_and_the_link_stays_out_of_the_logs(acct, capsys):
    rid, res = _pair(acct)
    assert res["state"] == "offered" and res["link"] == acct.daemon.links[0]
    path = acct.pairing_dir("anna") / f"res-{rid}.json"
    assert stat.S_IMODE(os.stat(path).st_mode) == 0o600
    assert stat.S_IMODE(os.stat(acct.pairing_dir("anna")).st_mode) == 0o700
    inv = store.open_invitation("anna", rid)
    assert inv and inv[1]["expires"] - time.time() <= pairing.INVITATION_TTL_S
    assert "messenger_invitation_issued" in _kinds()
    secret = "SECRETKEYMATERIAL"
    assert secret not in acct.log.read_text()
    assert pairing.link_tag(res["link"]) in acct.log.read_text()
    assert all(secret not in e.detail for e in security_events.recent(200))
    assert secret not in capsys.readouterr().err
    # The owner's pairing directory got nothing.
    assert not (acct.pairing_dir("") / f"res-{rid}.json").exists()


def test_a_connected_phone_may_do_nothing_until_its_code_is_confirmed(acct):
    rid, _ = _pair(acct)
    acct.connect("21", "901", name="Anna Pixel")
    rec = store.contact("21")
    assert rec["person"] == "anna" and rec["state"] == "unverified"
    assert pairing.result(rid, directory=acct.pairing_dir("anna"))["state"] == "connected"
    code = _code(acct, "21")
    assert len(code) == pairing.CODE_DIGITS and code.isdigit()
    # A message in the verification window: told to verify, nothing done.
    before = len(acct.client.commands)
    acct.say("21", "@21", "neu Angriff")
    assert len(acct.client.commands) == before
    assert not acct.services["anna"].calls("session.create")
    assert "Erst bestätigen" in acct.daemon.to("@21")[-1]
    assert "messenger_unverified_message" in _kinds()
    # The control: the same message after verification starts the session.
    vid = pairing.confirm("session:s-anna", code, directory=acct.pairing_dir("anna"))
    acct.round()
    assert pairing.result(vid, directory=acct.pairing_dir("anna"))["state"] == "verified"
    assert pairing.result(rid, directory=acct.pairing_dir("anna"))["state"] == "verified"
    assert store.verified_person("21") == "anna"
    acct.say("21", "@21", "neu Angriff")
    assert acct.services["anna"].calls("session.create")
    assert not acct.services["ben"].calls("session.create")


def test_only_the_session_that_started_pairing_can_confirm_it(acct):
    rid, _ = _pair(acct)
    acct.connect("21", "901")
    code = _code(acct, "21")
    # another session of anna holds no nonce for this pairing
    assert pairing.confirm("session:other", code, directory=acct.pairing_dir("anna")) is None
    # a verify forged with the right code but without the nonce
    forged = pairing.submit("verify", {"pair": rid, "nonce": "0" * 32, "code": code},
                            directory=acct.pairing_dir("anna"))
    acct.round()
    assert pairing.result(forged, directory=acct.pairing_dir("anna"))["state"] == "refused"
    assert store.contact("21")["state"] == "unverified"
    # the control
    pairing.confirm("session:s-anna", code, directory=acct.pairing_dir("anna"))
    acct.round()
    assert store.verified_person("21") == "anna"


def test_wrong_codes_then_the_phone_is_removed(acct):
    _pair(acct)
    acct.connect("21", "901")
    for i in range(pairing.MAX_ATTEMPTS):
        vid = pairing.confirm("session:s-anna", "000000", directory=acct.pairing_dir("anna"))
        acct.round()
    assert pairing.result(vid, directory=acct.pairing_dir("anna"))["state"] == "failed"
    assert store.contact("21") is None and "21" in acct.daemon.deleted
    assert "messenger_verification_failed" in _kinds()


def test_a_second_phone_through_the_same_invitation_is_refused(acct):
    _pair(acct)
    acct.connect("21", "901")
    acct.connect("22", "901")                       # the same link, again
    assert store.contact("22") is None and "22" in acct.daemon.deleted
    assert "messenger_invitation_replayed" in _kinds()
    assert store.contact("21")["person"] == "anna"   # the first one stays


def test_an_expired_invitation_binds_nothing(acct, monkeypatch):
    rid, _ = _pair(acct)
    later = time.time() + pairing.INVITATION_TTL_S + 30
    monkeypatch.setattr(time, "time", lambda: later)
    acct.round()
    assert pairing.result(rid, directory=acct.pairing_dir("anna"))["state"] == "expired"
    assert "901" in acct.daemon.dropped
    acct.connect("21", "901")
    assert store.contact("21") is None and "21" in acct.daemon.deleted
    assert "messenger_invitation_expired" in _kinds()


def test_a_replayed_request_id_gets_no_second_invitation(acct):
    rid, _ = _pair(acct)
    acct.connect("21", "901")
    pairing.confirm("session:s-anna", _code(acct, "21"), directory=acct.pairing_dir("anna"))
    acct.round()
    # the finished request file put back, as if replayed
    pairing._write(acct.pairing_dir("anna") / f"req-{rid}.json",
                   {"id": rid, "kind": "pair", "created": time.time(),
                    "requester": "session:s-anna", "nonce_sha": "a" * 64})
    acct.round()
    assert len(acct.daemon.links) == 1


def test_the_other_phones_hear_of_a_new_one_and_of_a_revocation(acct):
    store.save("20", person="anna", label="old phone")
    _pair(acct)
    acct.connect("21", "901", name="new phone")
    pairing.confirm("session:s-anna", _code(acct, "21"), directory=acct.pairing_dir("anna"))
    acct.round()
    assert any("Neues Handy verbunden: new phone" in t for t in acct.daemon.to("@20"))
    assert any("new phone" in n["text"] for n in pairing.notes(acct.pairing_dir("anna")))
    # anna unpairs the new phone from her session
    rid = pairing.submit("revoke", {"label": "new phone"}, directory=acct.pairing_dir("anna"))
    acct.round()
    assert pairing.result(rid, directory=acct.pairing_dir("anna"))["revoked"] == ["new phone"]
    assert store.contact("21") is None
    assert any("entkoppelt: new phone" in t for t in acct.daemon.to("@20"))
    assert "messenger_phone_revoked" in _kinds()
    # ben's phones heard nothing of anna's
    assert not any(c.startswith("@3") for c, _ in acct.daemon.sent)


def test_the_phone_limit_and_the_rate_limit_refuse_with_a_reason(acct):
    for i in range(pairing.MAX_PHONES):
        store.save(str(40 + i), person="anna")
    _, res = _pair(acct)
    assert res == {"kind": "pair", "state": "refused", "reason": "limit"}
    assert "messenger_invitation_refused" in _kinds()
    # ben: three invitations in the window, the fourth refused
    states = [_pair(acct, "ben", f"session:b{i}")[1]["state"]
              for i in range(pairing.INVITES_PER_WINDOW + 1)]
    assert states == ["offered"] * pairing.INVITES_PER_WINDOW + ["refused"]


def test_a_new_pair_request_supersedes_the_open_invitation(acct):
    first, _ = _pair(acct, requester="session:a")
    second, res = _pair(acct, requester="session:b")
    assert res["state"] == "offered"
    assert pairing.result(first, directory=acct.pairing_dir("anna"))["state"] == "cancelled"
    acct.connect("21", "901")                     # the first link: no longer valid
    assert store.contact("21") is None


def test_a_request_file_of_another_owner_or_mode_is_ignored(acct):
    path = acct.pairing_dir("anna")
    path.mkdir(parents=True, exist_ok=True)
    rid = "b" * 32
    (path / f"req-{rid}.json").write_text(
        '{"id": "%s", "kind": "pair", "created": %f, "nonce_sha": "%s"}' % (rid, time.time(), "c" * 64))
    os.chmod(path / f"req-{rid}.json", 0o644)      # readable by others: refused
    acct.round()
    assert not acct.daemon.links and not (path / f"req-{rid}.json").exists()


def test_a_concurrent_connect_for_one_invitation_binds_one_phone(acct):
    _pair(acct)

    async def both():
        event = lambda cid: {"type": "contactConnected",       # noqa: E731
                             "contact": {"contactId": cid, "activeConn": {"connId": 901}}}
        await asyncio.gather(acct.rooms.contact_connected(event(21)),
                             acct.rooms.contact_connected(event(22)))
    asyncio.run(both())
    bound = [c for c in ("21", "22") if store.contact(c)]
    assert len(bound) == 1


def test_the_session_command_shows_the_invitation_live_and_confirms(acct, monkeypatch):
    monkeypatch.setattr(shutil, "which", lambda name: None)        # no qrencode: link only
    shown = []
    out = pairing.session_command("", requester="session:s9", show=shown.append,
                                  spawn=lambda fn: None, directory=acct.pairing_dir("anna"))
    assert "invitation" in out
    rid = pairing.open_request("session:s9")
    acct.round()
    text = pairing.offer_text(pairing.read(rid, "session:s9", directory=acct.pairing_dir("anna")))
    assert "1. /pair" in text and "2. scan" in text and "3. type" in text
    assert acct.daemon.links[0] in text
    acct.connect("21", "901")
    code = _code(acct, "21")
    reply = pairing.session_command(code, requester="session:s9", show=shown.append,
                                    spawn=lambda fn: fn(), directory=acct.pairing_dir("anna"),
                                    timeout=0.0)
    assert "Code sent" in reply
    acct.round()
    assert store.verified_person("21") == "anna"


def test_every_messenger_event_kind_is_registered_with_the_panel():
    assert set(audit.KINDS) <= security_events.known_kinds()
    import ast
    from pathlib import Path
    root = Path(__file__).resolve().parents[1] / "delfin" / "messenger"
    used = set()
    for src in root.glob("*.py"):
        for node in ast.walk(ast.parse(src.read_text())):
            if (isinstance(node, ast.Call) and isinstance(node.func, ast.Attribute)
                    and node.func.attr == "event" and node.args
                    and isinstance(node.args[0], ast.Constant)):
                used.add(node.args[0].value)
            if (isinstance(node, ast.IfExp) and isinstance(node.body, ast.Constant)
                    and str(node.body.value).startswith("messenger_")):
                used |= {node.body.value, node.orelse.value}
    assert used and used <= set(audit.KINDS), used - set(audit.KINDS)


def test_a_revoked_phone_is_removed_from_its_groups(acct):
    _pair(acct)
    acct.connect("21", "901")
    pairing.confirm("session:s-anna", _code(acct, "21"), directory=acct.pairing_dir("anna"))
    acct.round()
    forum = store.forum_of("anna")
    assert forum and "21" in store.rooms()[forum]["members"]
    pairing.submit("revoke", {"label": "all"}, directory=acct.pairing_dir("anna"))
    acct.round()
    assert (forum, "21") in acct.daemon.removed and "21" in acct.daemon.deleted
    assert not store.uncleaned()


class _AdminDaemon:
    invitation_id = 77

    def __init__(self, *args, **kwargs):
        self.sent = []

    async def __aenter__(self):
        return self

    async def __aexit__(self, *args):
        pass

    async def subscribe(self):
        pass

    async def pending_invitation(self):
        return "simplex:/invitation#/?v=2&smp=ADMINLINK"

    async def events(self):
        yield {"type": "contactConnected",
               "contact": {"contactId": 5, "activeConn": {"connId": 77}}}

    async def contact_code(self, cid):
        return None

    async def send_text(self, cid, text):
        self.sent.append(text)

    async def delete_contact(self, cid):
        pass


def test_the_admin_pairing_of_the_owner_needs_the_code_too(acct, monkeypatch):
    from delfin.messenger import __main__ as admin
    monkeypatch.setattr(admin.simplex, "Daemon", _AdminDaemon)
    monkeypatch.setattr(admin.pairing, "new_code", lambda: "135790")
    monkeypatch.setattr(admin, "_ask_code", lambda: None)        # no terminal answer
    assert admin._pair() == 0
    assert store.contact("5")["state"] == "unverified" and store.paired_contact() is None
    assert admin.main(["verify", "000000"]) == 1                  # a wrong code
    assert store.paired_contact() is None
    assert admin.main(["verify", "135790"]) == 0
    assert store.paired_contact() == "5" and store.verified_person("5") == ""
    assert admin.main(["people"]) == 0
