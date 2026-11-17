"""A session record carries an optional person id, as data only.

Each person gets their own workspace (a later step of the plan); this
step only lets the record carry a ``person`` string so that later steps
can group and separate sessions by person. Nothing reads it to make a
decision yet — the field is a label, not a boundary.

Rules this step pins:

* ``save_session`` takes an optional ``person`` and stores it; the
  default is ``""``, so today's single-user saves change nothing.
* ``list_sessions`` returns it; ``load_session`` keeps it in the record.
* A record written before the field existed (no ``person`` key) reads
  back as ``""`` and everything about it keeps working.
* The terminal resume row (``cli_resume.list_sessions``) carries it too.
* The agent service declares it on ``session.info`` as an optional
  field, so a session carries it along where it already carries the
  session id; an entry without a person stays valid.
"""

from __future__ import annotations

import json
import time

import pytest

from delfin.agent import session_store as ss
from delfin.agent.cli_resume import list_sessions as cli_list_sessions
from delfin.agent.service import protocol


@pytest.fixture
def sessions_dir(tmp_path):
    d = tmp_path / "agent_sessions"
    d.mkdir()
    with pytest.MonkeyPatch.context() as mp:
        mp.setattr(ss, "_SESSIONS_DIR", d)
        yield d


@pytest.fixture
def store(tmp_path):
    d = tmp_path / "agent_sessions"
    d.mkdir()
    return d


def _write_old_record(d, session_id, **fields):
    # A record predating the person field: no ``person`` key at all.
    data = {"session_id": session_id, "updated_at": time.time(),
            "mode": "solo", "title": "before person existed"}
    data.update(fields)
    (d / f"{session_id}.json").write_text(json.dumps(data))


def test_save_session_stores_an_optional_person(sessions_dir):
    ss.save_session("s1", person="alice",
                    chat_messages=[{"role": "user", "content": "hi"}])
    ss.release_session_lock("s1")
    row = next(s for s in ss.list_sessions() if s["session_id"] == "s1")
    assert row["person"] == "alice"
    assert ss.load_session("s1")["person"] == "alice"


def test_default_person_is_empty(sessions_dir):
    ss.save_session("s2", chat_messages=[{"role": "user", "content": "hi"}])
    ss.release_session_lock("s2")
    row = next(s for s in ss.list_sessions() if s["session_id"] == "s2")
    assert row["person"] == ""
    assert ss.load_session("s2")["person"] == ""


def test_a_record_before_the_field_reads_back_empty(sessions_dir):
    _write_old_record(sessions_dir, "legacy")
    row = next(s for s in ss.list_sessions() if s["session_id"] == "legacy")
    assert row["person"] == ""
    # Everything about the legacy record still works: it lists and loads.
    assert ss.load_session("legacy")["session_id"] == "legacy"
    assert [s["session_id"] for s in ss.list_sessions()] == ["legacy"]


def test_the_terminal_resume_row_carries_person(store):
    import json as _json
    data = {"session_id": "res1", "created_at": 1000.0, "updated_at": 2000.0,
            "model": "kit", "workspace": "/tmp", "person": "bob",
            "chat_messages": [{"role": "user", "content": "work"},
                              {"role": "assistant", "content": "done"}]}
    (store / "res1.json").write_text(_json.dumps(data))
    rows = cli_list_sessions(sessions_dir=store)
    assert rows[0]["person"] == "bob"


def test_a_terminal_row_without_person_reads_empty(store):
    import json as _json
    data = {"session_id": "res2", "created_at": 1000.0, "updated_at": 2000.0,
            "model": "kit", "workspace": "/tmp",
            "chat_messages": [{"role": "user", "content": "work"}]}
    (store / "res2.json").write_text(_json.dumps(data))
    rows = cli_list_sessions(sessions_dir=store)
    assert rows[0]["person"] == ""


def test_the_service_declares_person_on_session_info():
    assert protocol.EVENTS["session.info"]["person"] == "string?"
