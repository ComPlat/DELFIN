"""A phone's display name is shown only cleaned: printable, no escapes, 50 characters.

Review finding 8: the SimpleX display name a phone chooses went unchanged
into the store, the person's session ("phone connected: <label>"), the
pairing notes and the other phones' messages -- terminal escapes and
control characters included.

Output: the label the gateway binds (the person's pairing path and the
admin path) is activation.clean_name of the display name: printable only,
no ESC or other control character, at most 50 characters; what the other
phones and the notes get is that label. Control: a plain name is kept.
"""
from __future__ import annotations


import pytest

from _messenger_people_rig import Account
from delfin.messenger import pairing, store

EVIL = "Evil\x1b[2J\x1b]0;owned\x07Phone\r\n" + "x" * 80


def _clean(text):
    return all(c.isprintable() for c in text) and "\x1b" not in text and len(text) <= 50


@pytest.fixture
def acct(tmp_path, monkeypatch):
    return Account(tmp_path, monkeypatch)


def _pair_and_verify(acct, name):
    pairing.begin("session:s", directory=acct.pairing_dir("anna"))
    acct.round()
    acct.connect("21", "901", name=name)
    code = acct.daemon.to("@21")[0].split("Code ")[1].split(".")[0]
    pairing.confirm("session:s", code, directory=acct.pairing_dir("anna"))
    acct.round()


def test_the_bound_label_and_every_message_are_clean(acct):
    store.save("20", person="anna", label="old phone")
    _pair_and_verify(acct, EVIL)
    label = store.contact("21")["label"]
    assert _clean(label) and label.startswith("Evil")
    for text in acct.daemon.to("@20"):
        assert "\x1b" not in text and "\x07" not in text
    for note in pairing.notes(acct.pairing_dir("anna")):
        assert _clean(note["text"][:60]) or "\x1b" not in note["text"]
    pairing.forget("session:s")


def test_a_label_written_to_the_store_directly_is_cleaned_too(acct):
    store.bind("22", person="anna", label=EVIL, request="r" * 32, code="123456")
    assert _clean(store.contact("22")["label"])


def test_a_plain_name_is_kept_control(acct):
    _pair_and_verify(acct, "Anna Pixel 8")
    assert store.contact("21")["label"] == "Anna Pixel 8"
    pairing.forget("session:s")
