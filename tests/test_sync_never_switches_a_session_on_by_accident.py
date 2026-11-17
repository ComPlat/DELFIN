"""``/sync`` switches on only when asked to; anything unknown changes nothing.

Review finding 2: ``/sync <anything>`` switched a session ON and named it
``<anything>`` -- ``/sync aus`` enabled the session and called it "aus".

Output, through activation.configure and the service worker's handler:

* off: off, aus, ausschalten, disable, stop; on: (nothing), on, an, ein,
  enable; status: status;
* a name only with ``/sync name NAME`` or a quoted ``/sync "NAME"``;
* anything else: the usage text, the record unchanged, no title.
"""
from __future__ import annotations

import pytest

from delfin.messenger import activation


@pytest.fixture(autouse=True)
def _records(tmp_path, monkeypatch):
    monkeypatch.setattr(activation, "_DIR", tmp_path / "act")


@pytest.mark.parametrize("word", ["off", "aus", "ausschalten", "disable", "stop", "AUS"])
def test_the_off_words_switch_off(word):
    activation.configure("s1")
    assert activation.enabled("s1")
    activation.configure("s1", word)
    assert not activation.enabled("s1")


@pytest.mark.parametrize("word", ["", "on", "an", "ein", "enable"])
def test_the_on_words_switch_on(word):
    activation.configure("s1", word)
    assert activation.enabled("s1")


# 2026-11-12: one clean word now names the session (see
# test_sync_word_names_the_session.py); what is NOT one word still
# changes nothing.
@pytest.mark.parametrize("args", ["aus bitte", "on off", "OFFICE 2",
                                   "Bitte benenne die Session B\u00fcro"])
def test_anything_else_changes_nothing(args):
    out = activation.configure("s1", args)
    assert "Usage" in out
    assert activation.record("s1") == {}
    assert activation.parse_sync(args)[0] == "usage"


@pytest.mark.parametrize("args,name", [("name Butan", "Butan"),
                                       ("name Butan Konformere", "Butan Konformere"),
                                       ('"Butan Konformere"', "Butan Konformere"),
                                       ("'Ethan'", "Ethan")])
def test_a_name_only_by_name_or_quotes(args, name):
    assert activation.parse_sync(args) == ("on", name)
    activation.configure("s1", args)
    assert activation.enabled("s1")


def _worker(sid="s5"):
    from delfin.agent.service import worker
    titles, notices = [], []
    host = type("W", (), {})()
    host.emitter = type("E", (), {"session_id": sid,
                                  "emit": lambda self, kind, text=None, **k: notices.append(text)})()
    host._m_session_setTitle = lambda params: titles.append(params["title"]) or {"title": params["title"]}
    return worker, host, titles, notices


def test_the_worker_names_only_by_name():
    worker, host, titles, notices = _worker()
    worker.Worker._m_session_prompt(host, {"text": "/sync aus"})
    assert not activation.enabled("s5") and titles == []
    worker.Worker._m_session_prompt(host, {"text": "/sync aus bitte"})
    assert not activation.enabled("s5") and titles == [] and "Usage" in notices[-1]
    worker.Worker._m_session_prompt(host, {"text": "/sync name Butan"})
    assert activation.enabled("s5") and titles == ["Butan"]
