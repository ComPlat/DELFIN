"""A schedule's confirmation is signed with the account key.

The content hash alone is public: anything that can write the store could
compute it and mark its own entry confirmed. The confirmation is therefore
an HMAC under the account key (``trust_store``), which no process carrying
the agent-session marker obtains and every cage masks. A plain hash is
taken only where no account key can be had at all (a refused store), so a
broken key file does not stop the user's schedules.
"""

from __future__ import annotations

import json

import pytest

from delfin.agent import approval_guard, schedules, trust_store


@pytest.fixture
def ws(tmp_path, monkeypatch):
    monkeypatch.setattr(schedules, "_DIR", tmp_path / "user_schedules")
    monkeypatch.setattr(schedules, "_hostname", lambda: "here")
    trust_store.forget_cached_keys()
    folder = tmp_path / "ws"
    folder.mkdir()
    yield folder
    trust_store.forget_cached_keys()


def _forge(entry_id: int, tag: str) -> None:
    data = json.loads(schedules._file().read_text())
    data["entries"][str(entry_id)]["confirmed_sha"] = tag
    schedules._file().write_text(json.dumps(data))


def test_the_users_confirmation_is_keyed(ws):
    entry = schedules.add("alle 1h x", workspace=str(ws), created_by="t", confirm=True)
    tag = json.loads(schedules._file().read_text())["entries"][str(entry["id"])]["confirmed_sha"]
    assert tag and tag != schedules.fingerprint(entry)
    assert schedules.get(entry["id"])["confirmed"]


@pytest.mark.parametrize("forged", ["plain", "sha"])
def test_a_written_hash_does_not_confirm(ws, forged):
    entry = schedules.add("alle 1h x", workspace=str(ws), created_by="agent", confirm=False)
    fp = schedules.fingerprint(schedules.get(entry["id"]))
    _forge(entry["id"], fp if forged == "plain" else "sha:" + fp)
    assert not schedules.get(entry["id"])["confirmed"]
    assert schedules.take_due("here", now=entry["next_fire_at"] + 1) == []


def test_inside_an_agent_session_nothing_is_confirmed(ws, monkeypatch):
    entry = schedules.add("alle 1h x", workspace=str(ws), created_by="agent", confirm=False)
    monkeypatch.setenv(approval_guard.ENV_MARKER, "1")
    trust_store.forget_cached_keys()
    with pytest.raises(schedules.ScheduleError):
        schedules.confirm(entry["id"])
    assert not schedules.get(entry["id"])["confirmed"]


def test_without_any_account_key_the_plain_hash_still_works(ws, monkeypatch):
    monkeypatch.setattr(trust_store, "account_key", lambda **_: None)
    trust_store.forget_cached_keys()
    entry = schedules.add("alle 1h x", workspace=str(ws), created_by="t", confirm=True)
    assert schedules.get(entry["id"])["confirmed"]
