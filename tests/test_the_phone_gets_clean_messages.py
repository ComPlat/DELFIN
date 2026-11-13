"""What the phone is shown is short, clean and in SimpleX markdown.

Field report (Pi, owner session "Office"): every message in a session's
chat started with "[Office · 0a56d860]", approvals were a page of text
with every option and rule, the run's end was the bare state word "idle".

Output, rendered through the real gateway and bridge (the rig), each with
its control:

* a session chat gets no "[name · id]" prefix (the group name says which
  session); a forum post starts with the session's name in bold;
* no raw ids or hashes in any outgoing text;
* an approval is one compact card: what (tool and a one-line summary),
  why (the agent's reason), the command in a code block when it is short
  (at most 6 lines), and "Reply yes / no"; a long one is shown on ``show``;
* one "working" message per turn, edited when the turn ends;
* the final answer is sent as the model wrote it, without the state
  word; a long one is cut with a ``more`` hint, and ``more`` sends the
  rest;
* the first look at an idle session sends nothing.
"""
from __future__ import annotations

import asyncio
import re

import pytest

from _messenger_people_rig import Account
from delfin.messenger import store
from delfin.messenger.gateway import _PERSON

HEXISH = re.compile(r"\b[0-9a-f]{8,}\b")


@pytest.fixture
def acct(tmp_path, monkeypatch):
    monkeypatch.setenv("DELFIN_MESSENGER_LANG", "en")
    account = Account(tmp_path, monkeypatch)
    store.save("21", person="anna", label="anna phone")
    account.enable("anna", "s-anna")
    store.put_room("#a1", "s-anna", "Butan", person="anna", ready=True)
    store.put_room("#f1", None, "Forum", person="anna", ready=True)
    account.service = account.services["anna"]
    account.row = account.service.sessions[0]
    return account


def refresh(acct):
    async def go():
        token = _PERSON.set("anna")
        try:
            await acct.rooms.refresh()
        finally:
            _PERSON.reset(token)
    asyncio.run(go())


def clean(text):
    assert not text.startswith("["), text
    assert "s-anna" not in text and not HEXISH.search(text), text


def test_a_short_approval_is_one_compact_card(acct):
    acct.service.pending["s-anna"] = [{
        "requestId": "r1", "kind": "confirm", "tool": "bash",
        "preview": "$ ls -l results", "purpose": "list the results",
        "flags": {"shell": True}, "options": [{"key": "y", "label": "Once"},
                                               {"key": "n", "label": "No"}]}]
    refresh(acct)
    card = acct.daemon.to("#a1")[-1]
    clean(card)
    lines = card.splitlines()
    assert lines[0].startswith("🔐 *run a command*"), card
    assert "list the results" in card
    assert "```\n$ ls -l results\n```" in card
    assert re.fullmatch(r"Tap /'yes [0-9a-f]{4}' or /'no [0-9a-f]{4}', or reply \*yes\* / \*no\* on the card\.", lines[-1]), card
    assert len(lines) <= 8 and len(card) <= 400, card


def test_a_write_card_names_the_file_on_its_own_line(acct):
    preview = "--- a/notiz.txt\n+++ b/notiz.txt\n@@ -0,0 +1 @@\n+hallo\n"
    acct.service.pending["s-anna"] = [{
        "requestId": "r4", "kind": "confirm", "tool": "write_file",
        "preview": preview, "flags": {"write": True},
        "options": [{"key": "y", "label": "Once"}, {"key": "n", "label": "No"}]}]
    refresh(acct)
    card = acct.daemon.to("#a1")[-1]
    clean(card)
    lines = card.splitlines()
    # The file name is its own line, above the code block.
    assert any(line == "Datei: notiz.txt" or line == "File: notiz.txt"
               for line in lines), card
    # The diff scaffolding (headers and hunk line) is not shown; the
    # content stays in the code block, on its own lines.
    assert "```\n+hallo\n```" in card, card
    # The raw diff headers are not shown as content lines.
    assert not any(line.startswith(("--- a/", "+++ b/")) for line in lines), card


def test_a_short_write_card_shows_file_and_content(acct):
    preview = "--- a/notiz.txt\n+++ b/notiz.txt\n@@ -0,0 +1 @@\n+hallo\n"
    acct.service.pending["s-anna"] = [{
        "requestId": "r5", "kind": "confirm", "tool": "write_file",
        "preview": preview, "flags": {"write": True},
        "options": [{"key": "y", "label": "Once"}, {"key": "n", "label": "No"}]}]
    refresh(acct)
    card = acct.daemon.to("#a1")[-1]
    lines = card.splitlines()
    assert any(line == "File: notiz.txt" for line in lines), card
    assert "```\n+hallo\n```" in card, card


def test_the_show_hint_names_the_command_in_the_phones_language(acct, monkeypatch):
    preview = "\n".join(f"+line {i}" for i in range(20))
    acct.service.pending["s-anna"] = [{
        "requestId": "r6", "kind": "confirm", "tool": "write_file",
        "preview": "notiz.txt\n" + preview, "flags": {"write": True},
        "options": [{"key": "y", "label": "Once"}, {"key": "n", "label": "No"}]}]
    refresh(acct)
    card_en = acct.daemon.to("#a1")[-1]
    assert "reply `show`" in card_en, card_en
    monkeypatch.setenv("DELFIN_MESSENGER_LANG", "de")
    from delfin.messenger import lang, language
    language.set_contact("21", "de")
    lang.CURRENT.set(None)
    # A fresh request: unchanged pending items are not re-sent.
    acct.service.pending["s-anna"][0]["requestId"] = "r6-de"
    refresh(acct)
    card_de = acct.daemon.to("#a1")[-1]
    assert "mit `zeige` antworten" in card_de, card_de
    assert "`show`" not in card_de, card_de
    # Both spellings still work: 'zeige' on the German card answers it.
    mid = next(m for (c, m) in acct.rooms._reply_targets if c == "#a1")
    lang.CURRENT.set(None)
    monkeypatch.setenv("DELFIN_MESSENGER_LANG", "en")
    acct.say("21", "#a1", "zeige", quote=mid)
    full = acct.daemon.to("#a1")[-1]
    assert "line 19" in full, full


def test_a_card_without_a_preview_offers_no_whole_action():
    # Nothing to show: no 'show' and no 'more' hint ('more' continues a
    # cut message, it never shows an action).
    from delfin.messenger.render import approval_card
    card = approval_card({"kind": "confirm", "tool": "bash", "text": "",
                          "phoneApproval": True})
    assert "Whole action" not in card and "`more`" not in card, card


def test_a_long_approval_shows_its_preview_on_show(acct):
    preview = "\n".join(f"+line {i}" for i in range(20))
    acct.service.pending["s-anna"] = [{
        "requestId": "r2", "kind": "confirm", "tool": "write_file",
        "preview": "notes.txt\n" + preview, "flags": {"write": True},
        "options": [{"key": "y", "label": "Once"}, {"key": "n", "label": "No"}]}]
    refresh(acct)
    card = acct.daemon.to("#a1")[-1]
    clean(card)
    assert card.startswith("🔐 *write a file*: notes.txt"), card
    assert "+line 19" not in card and "`show`" in card
    mid = next(m for (c, m) in acct.rooms._reply_targets if c == "#a1")
    acct.say("21", "#a1", "show", quote=mid)
    full = acct.daemon.to("#a1")[-1]
    assert "+line 19" in full and full.startswith("```"), full
    assert acct.service.pending["s-anna"], "show never answers"


def test_a_question_card_has_no_prefix_and_no_id(acct):
    acct.service.pending["s-anna"] = [{
        "requestId": "q1", "kind": "ask", "question": "Which solvent?",
        "choices": ["water", "toluene"], "multiSelect": False}]
    refresh(acct)
    card = acct.daemon.to("#a1")[-1]
    clean(card)
    assert card.startswith("❓ *Which solvent?*"), card


def test_the_forum_names_the_session_in_bold(acct):
    acct.service.pending["s-anna"] = [{
        "requestId": "r3", "kind": "confirm", "tool": "bash", "preview": "$ ls",
        "flags": {"shell": True}, "options": [{"key": "y", "label": "Once"},
                                               {"key": "n", "label": "No"}]}]
    refresh(acct)
    post = acct.daemon.to("#f1")[-1]
    clean(post)
    assert post.startswith("*Butan*"), post


def test_one_working_message_per_turn_and_a_clean_answer(acct):
    refresh(acct)
    assert acct.daemon.to("#a1") == [], "the first look at an idle session sends nothing"
    acct.say("21", "#a1", "compute the barrier")     # a turn started from the phone
    acct.row["state"] = "running"
    refresh(acct)
    refresh(acct)
    working = [t for t in acct.daemon.to("#a1") if "Working" in t]
    assert len(working) == 1, acct.daemon.to("#a1")
    acct.row.update(state="idle", answer="The *barrier* is 12.3 kcal/mol.")
    refresh(acct)
    answer = acct.daemon.to("#a1")[-1]
    assert answer == "The *barrier* is 12.3 kcal/mol.", answer
    assert acct.daemon.updated and acct.daemon.updated[-1][0] == "#a1"
    forum = acct.daemon.to("#f1")[-1]
    assert forum.startswith("*Butan*\n") and "12.3" in forum and "idle" not in forum


def test_a_long_answer_is_cut_and_more_sends_the_rest(acct):
    refresh(acct)
    long = " ".join(f"Sentence {i} of the long answer." for i in range(120))
    acct.row.update(state="running")
    refresh(acct)
    acct.row.update(state="idle", answer=long)
    refresh(acct)
    first = acct.daemon.to("#a1")[-1]
    assert len(first) <= 1000 and first.endswith("/'more'"), first[-80:]
    acct.say("21", "#a1", "more")
    rest = acct.daemon.to("#a1")[-1]
    assert "Sentence 119" in rest or "/'more'" in rest
    assert rest[:30] not in first


def test_a_reply_in_the_session_chat_has_no_prefix(acct):
    acct.say("21", "#a1", "calculate the barrier")
    out = acct.daemon.to("#a1")[-1]
    clean(out)
