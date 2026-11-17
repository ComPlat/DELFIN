"""``/sync WORD`` -- one unknown word -- switches on and names the session.

Seen 2026-11-12: ``/sync OFFICE`` answered with the usage text and the
session kept the title ``/sync OFFICE``; the person meant "on, named
OFFICE". One word that is not a keyword (on/off/status/policy/name) is
now a name: it switches the session on and names it. Everything that is
not one clean word -- several words, a bare ``name``, unknown ``policy``
arguments -- still changes nothing (the guard from review finding 2
stays for those).
"""
from __future__ import annotations

import pytest

from delfin.messenger import activation


@pytest.fixture(autouse=True)
def _records(tmp_path, monkeypatch):
    monkeypatch.setattr(activation, "_DIR", tmp_path / "act")


@pytest.mark.parametrize("word", ["OFFICE", "office", "Butan", "Lab", "Office-2", "Büro"])
def test_a_single_word_switches_on_and_names(word):
    assert activation.parse_sync(word) == ("on", word)
    activation.configure("s1", word)
    assert activation.enabled("s1")


@pytest.mark.parametrize("args", ["aus bitte", "on off",
                                  "Bitte benenne die Session Büro",
                                  "name", "OFFICE 2"])
def test_what_is_not_one_word_still_changes_nothing(args):
    assert activation.parse_sync(args)[0] == "usage"
    activation.configure("s1", args)
    assert activation.record("s1") == {}


def test_the_worker_switches_on_and_names_by_one_word():
    from delfin.agent.service import worker
    titles, notices = [], []
    host = type("W", (), {})()
    host.emitter = type("E", (), {"session_id": "s9",
                                  "emit": lambda self, kind, text=None, **k: notices.append(text)})()
    host._m_session_setTitle = lambda params: titles.append(params["title"]) or {"title": params["title"]}
    worker.Worker._m_session_prompt(host, {"text": "/sync OFFICE"})
    assert activation.enabled("s9") and titles == ["OFFICE"]
    assert "OFFICE" in str(notices[-1])


@pytest.mark.parametrize("word", [
    "nein", "no", "of", "ende", "stopp", "halt", "exit", "help", "hilfe",
    "ja", "yes", "off!", "status?", "x", "0", "1", "\"Off"])
def test_a_stop_answer_or_typo_never_switches_the_phone_on(word):
    # '/sync nein' or a typo of 'off' switching the phone ON would send
    # the session out against the person's intent: such a word changes
    # nothing.
    assert activation.parse_sync(word) == ("usage", "")
