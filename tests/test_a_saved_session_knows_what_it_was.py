"""A saved session says what it is — and whether it still lives.

Report 20260919: of the 52 session records then on disk, 52 carried an
empty ``model`` and ``provider`` (half of that was fixed on 2026-09-18)
and none carried a heartbeat, a pid or a host. Whoever ran ``--resume``
saw a list that could not say which model a session ran on, nor whether
its process was still alive — exactly the question that matters with
several parallel sessions on several nodes.

The record now carries model, provider, host and a heartbeat (pid +
process-start fingerprint + timestamp, refreshed while the session
works). Liveness is judged through ``proc_identity`` — a pid alone
names nothing, because the system hands the number out again.
"""

from __future__ import annotations

import json
import time

import pytest

from delfin.agent import session_store as ss
from delfin.agent.session_store import (
    list_sessions,
    load_session,
    save_session,
)


def _touch_heartbeat(session_id):
    # Imported lazily: with the heartbeat work absent, THIS name is the
    # only thing that may fail; every other case must fail on its own
    # assertion, so the control run shows what each one measures.
    from delfin.agent.session_store import touch_heartbeat
    return touch_heartbeat(session_id)


@pytest.fixture
def sessions_dir(tmp_path):
    d = tmp_path / "agent_sessions"
    d.mkdir()
    with pytest.MonkeyPatch.context() as mp:
        mp.setattr(ss, "_SESSIONS_DIR", d)
        yield d


def _write(d, session_id, **fields):
    data = {"session_id": session_id, "updated_at": time.time(),
            "mode": "solo", "title": "t"}
    data.update(fields)
    (d / f"{session_id}.json").write_text(json.dumps(data))


def test_a_saved_record_names_itself(sessions_dir):
    save_session("hb-1", mode="solo", model="kit", provider="glm",
                 chat_messages=[{"role": "user", "content": "hi"}])
    data = load_session("hb-1")
    assert data["model"] == "kit"
    assert data["provider"] == "glm"
    # heartbeat: this process, its start fingerprint, and when
    assert data["pid"] == pytest.importorskip("os").getpid()
    assert data["proc_start"]            # set, not guessed at creation
    assert data["host"]
    assert data["heartbeat_at"]


def test_touch_heartbeat_refreshes_without_touching_the_rest(sessions_dir):
    save_session("hb-2", mode="solo", cost_usd=1.25,
                 chat_messages=[{"role": "user", "content": "hi"}])
    before = load_session("hb-2")
    time.sleep(0.01)
    assert _touch_heartbeat("hb-2") is True
    after = load_session("hb-2")
    assert after["heartbeat_at"] > before["heartbeat_at"]
    assert after["cost_usd"] == 1.25          # the conversation is untouched
    assert after["chat_messages"] == before["chat_messages"]


def test_touch_heartbeat_on_unknown_session_is_false(sessions_dir):
    assert _touch_heartbeat("no-such-session") is False


def test_listing_says_whether_a_session_still_lives(sessions_dir):
    # A live one: this very process.
    save_session("hb-live", mode="solo", model="kit",
                 chat_messages=[{"role": "user", "content": "hi"}])
    # A dead one: a pid nothing answers to, with a fingerprint.
    _write(sessions_dir, "hb-dead", model="kit",
           pid=999_999_999, proc_start="12345", host="",
           heartbeat_at=time.time())
    # An old record: predates the fields entirely.
    _write(sessions_dir, "hb-old", model="")

    rows = {r["session_id"]: r for r in list_sessions(limit=10)}
    assert rows["hb-live"]["alive"] is True
    assert rows["hb-live"]["model"] == "kit"
    assert rows["hb-dead"]["alive"] is False
    # Old records must not fall over — they say "unknown" instead.
    assert rows["hb-old"]["alive"] is None
    assert rows["hb-old"]["model"] == ""
    assert rows["hb-live"]["heartbeat_at"] > 0


def test_liveness_of_a_reused_pid_is_false(sessions_dir):
    """A pid alone answers nothing: the number is handed out again.

    A stranger that took the pid has a DIFFERENT start fingerprint, so
    the record's answer must be False even though ``os.kill(pid, 0)``
    succeeds.
    """
    import os
    _write(sessions_dir, "hb-stranger",
           pid=os.getpid(), proc_start="1",   # wrong start: a stranger
           host="", heartbeat_at=time.time())
    rows = {r["session_id"]: r for r in list_sessions(limit=10)}
    assert rows["hb-stranger"]["alive"] is False


def test_a_foreign_host_is_not_judged_here(sessions_dir):
    """A pid from another machine names nothing here; guessing could
    cost somebody their running work, so the answer is None."""
    _write(sessions_dir, "hb-elsewhere",
           pid=1234, proc_start="77", host="some-other-node",
           heartbeat_at=time.time())
    rows = {r["session_id"]: r for r in list_sessions(limit=10)}
    assert rows["hb-elsewhere"]["alive"] is None


def test_liveness_works_without_proc(sessions_dir, monkeypatch):
    """Everything must work off Linux too: when /proc is gone, the
    ps fallback answers the same question. Proven by switching the
    /proc reader off, not by finding a machine without it."""
    from delfin.agent import proc_identity
    import os
    monkeypatch.setattr(proc_identity, "start_ticks", lambda pid: None)
    # ps must still give a fingerprint for THIS process...
    mine = proc_identity.process_start(os.getpid())
    assert mine, "no /proc and no ps: nothing left to ask"
    # ...and liveness must trust it: same fingerprint -> alive,
    # different one -> a stranger.
    _write(sessions_dir, "hb-noproc-live", pid=os.getpid(),
           proc_start=mine, host="", heartbeat_at=time.time())
    _write(sessions_dir, "hb-noproc-stranger", pid=os.getpid(),
           proc_start="not-the-real-one", host="", heartbeat_at=time.time())
    rows = {r["session_id"]: r for r in list_sessions(limit=10)}
    assert rows["hb-noproc-live"]["alive"] is True
    assert rows["hb-noproc-stranger"]["alive"] is False


def test_the_record_carries_nothing_private(sessions_dir):
    """The repository is public: no token, no home path, no account
    name may land in a session record's identity fields."""
    save_session("hb-priv", mode="solo",
                 chat_messages=[{"role": "user", "content": "hi"}])
    raw = (sessions_dir / "hb-priv.json").read_text()
    import os
    home = os.path.expanduser("~")
    assert home not in raw
    assert "proc_start" in raw and "pid" in raw and "host" in raw
