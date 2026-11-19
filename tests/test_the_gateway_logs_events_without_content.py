"""The gateway's journal says what happened to each message, never what
it said.

Field debugging on the Pi was impossible: the journal had only the
"notifications:" lines. Now one content-free line per message handled
(delfin/messenger/eventlog.py): direction, kind (prompt, a command's name,
approval answer, question answer, file, voice), the person, the room and
the session as short tags (never their ids), the result (delivered,
refused: reason, no session, failed) and the latency in ms. On by default,
off with DELFIN_MESSENGER_EVENTS=off; a repeated line is collapsed.

Output, through the real gateway and bridge (the rig), each with its
control: the lines are there with their fields; no message text, link,
code or file name ever appears in the log.
"""
from __future__ import annotations

import re

import pytest

from _messenger_people_rig import CONFIRM, Account
from delfin.messenger import store

SECRET = "the-secret-prompt-text"
LINK = "simplex:/invitation#/?v=2&smp=SECRETKEYMATERIAL0001"


@pytest.fixture
def acct(tmp_path, monkeypatch):
    from delfin.messenger import eventlog
    monkeypatch.setattr(eventlog, "_last", {})       # no collapsing across tests
    monkeypatch.delenv("DELFIN_MESSENGER_EVENTS", raising=False)
    account = Account(tmp_path, monkeypatch)
    store.save("21", person="anna", label="anna phone")
    account.enable("anna", "s-anna")
    store.put_room("#a1", "s-anna", "Butan", person="anna", ready=True)
    return account


def events(err):
    return [line for line in err.splitlines() if line.startswith("[gateway] event ")]


def test_a_prompt_is_logged_without_its_text(acct, capsys):
    acct.say("21", "#a1", f"{SECRET} {LINK} 123456")
    err = capsys.readouterr().err
    lines = events(err)
    assert lines, err
    line = next(x for x in lines if x.startswith("[gateway] event in "))
    assert " in " in line and "kind=prompt" in line and "person=anna" in line
    assert "result=delivered" in line and re.search(r"ms=\d+", line)
    assert re.search(r"room=r[0-9a-f]{4}\b", line) and re.search(r"session=s[0-9a-f]{4}\b", line)
    for leak in (SECRET, "SECRETKEYMATERIAL", "simplex:/", "123456", "#a1", "s-anna", "Butan"):
        assert leak not in err, leak


def test_commands_and_answers_are_named_by_kind(acct, capsys):
    acct.services["anna"].pending["s-anna"] = [{**CONFIRM, "requestId": "r1"}]
    acct.say("21", "#a1", "pending")
    mid = next(m for (c, m) in acct.rooms._reply_targets if c == "#a1")
    acct.say("21", "#a1", "yes", quote=mid)
    kinds = [re.search(r"kind=(\S+)", line).group(1) for line in events(capsys.readouterr().err)]
    assert "pending" in kinds and "approval_answer" in kinds


def test_a_refused_message_says_why(acct, capsys):
    acct.say("99", "#a1", SECRET)                       # not a paired phone
    store.bind("25", person="anna", label="x", request="r" * 32, code="123456")
    acct.say("25", "#a1", SECRET)                       # not verified yet
    lines = events(capsys.readouterr().err)
    assert any("result=refused:unknown" in line for line in lines)
    assert any("result=refused:unverified" in line for line in lines)


def test_repeats_are_collapsed(acct, capsys):
    for i in range(20):
        acct.say("99", "#a1", f"{SECRET} {i}")
    lines = events(capsys.readouterr().err)
    assert 1 <= len(lines) <= 3, lines


def test_off_means_off(acct, capsys, monkeypatch):
    monkeypatch.setenv("DELFIN_MESSENGER_EVENTS", "off")
    acct.say("21", "#a1", SECRET)
    assert not events(capsys.readouterr().err)


def test_an_outgoing_message_is_logged_without_its_text(acct, capsys):
    acct.say("21", "@21", "help")
    lines = [line for line in events(capsys.readouterr().err) if " out " in line]
    assert lines and all("Try:" not in line and "sessions" not in line for line in lines)
