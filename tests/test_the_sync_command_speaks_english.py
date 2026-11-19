"""/sync answers in English, in the TUI and the dashboard alike.

The answers of ``/sync`` and its line in the command list were German
("SimpleX für diese Sitzung on|off|status", "Noch keine Sitzung. Zuerst
eine Sitzung starten.") in an interface that is English everywhere else.
What the phone is told stays with the messenger's own language setting
(delfin.messenger.lang); this is what the person at the screen reads.
"""
from __future__ import annotations

import re

import pytest

from delfin.agent import repl_commands
from delfin.messenger import activation

_GERMAN = re.compile(
    r"[äöüÄÖÜß]|\b(?:Sitzung|nicht|bitte|Zuerst|Verwendung|aktiviert|"
    r"ausgeschaltet|für|Gilt|Noch|keine)\b")


def _english(text: str) -> bool:
    return bool(text) and not _GERMAN.search(text)


@pytest.fixture
def store(tmp_path, monkeypatch):
    monkeypatch.setattr(activation, "_DIR", tmp_path)


@pytest.mark.parametrize("args", ["", "on", "status", "off", "bogus"])
def test_every_answer_of_sync_is_english(store, args):
    answer = activation.configure("s1", args, endpoint="")
    assert _english(answer), answer


def test_sync_without_a_session_says_so_in_english(store):
    answer = activation.configure("", "on")
    assert _english(answer), answer
    assert "session" in answer.lower()


def test_status_names_the_state(store):
    activation.configure("s1", "on", endpoint="")
    # "on." and, with no messenger bridge on this machine, the notice that
    # nothing reaches the phone yet
    assert activation.configure("s1", "status").startswith("SimpleX for this session: on.")
    activation.configure("s1", "off", endpoint="")
    assert activation.configure("s1", "status").endswith("off.")


def test_the_command_list_describes_sync_in_english():
    summary = repl_commands.BUILTINS["/sync"].summary
    assert _english(summary), summary
    assert "on|off|status" in summary
