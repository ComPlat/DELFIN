"""A session working through one long turn stays addressable.

Three sessions worked together on a cluster and five messages between
them came back "no other open session '<key>'" -- for keys their peers
had reached minutes earlier in the same run. The messages were not
queued; they were refused and lost.

Two mechanisms, measured:

`announce` is called from the terminal's IDLE poll and from the
dashboard's refresh timer. A session in a long turn is not idle and may
have no window watching it, so it says nothing, and a session that says
nothing for _STALE_S counts as not open.

Then `_reap` DELETED its record. `_known_key` reads that file to decide
whether a message can be queued for a session that is closed or
mid-restart -- its own docstring says so -- so deleting it turned
"queue this for later" into a refusal. The sequence needs two steps and
both are ordinary: the session goes quiet, ANY session lists its peers
(the roster, or session_message with no `to`, both of which call
open_sessions), and the next message to that key is dropped.

So: quiet is not dead, and a tool call is a heartbeat.
"""

from __future__ import annotations

import json
import os
import socket
import time

import pytest

from delfin.agent import session_messages as M
from delfin.agent import session_presence as P


@pytest.fixture
def world(tmp_path, monkeypatch):
    """A presence directory and a mailbox of this test's own."""
    monkeypatch.setattr(P, "_DIR", tmp_path / "presence")
    monkeypatch.setattr(M, "_DIR", tmp_path / "mail")
    (tmp_path / "presence").mkdir()
    (tmp_path / "mail").mkdir()
    monkeypatch.setattr(P, "_last_written", {})
    monkeypatch.setattr(P, "_last_touched", {})
    return tmp_path


KEY = "3eb6c29f"


def _record(quiet_for: float = 0.0, pid: int | None = None,
            host: str | None = None, **extra) -> None:
    P._path(KEY).write_text(json.dumps({
        "key": KEY,
        "session_id": "c" * 32,
        "title": "Session C",
        "workspace": "/somewhere/worktree",
        "branch": "agent/c",
        "host": socket.gethostname() if host is None else host,
        "pid": os.getpid() if pid is None else pid,
        "updated_at": time.time() - quiet_for,
        **extra,
    }), encoding="utf-8")


def _survives() -> bool:
    P.open_sessions()
    return P._path(KEY).exists()


class TestQuietIsNotDead:
    def test_a_quiet_session_keeps_its_record(self, world):
        """The record is what lets a peer queue a message for it."""
        _record(quiet_for=P._STALE_S + 60)
        assert _survives(), "a session that is merely quiet is not gone"

    def test_it_is_still_not_on_the_roster(self, world):
        """Kept is not the same as open: the roster must not show a
        session nobody has heard from as available."""
        _record(quiet_for=P._STALE_S + 60)
        assert KEY not in [r.get("key") for r in P.open_sessions()]

    def test_a_dead_process_is_reaped_at_once(self, world):
        """The case the branch was written for: a crashed kernel never
        withdraws its record, and one file per crash was read whole on
        every refresh."""
        _record(quiet_for=P._STALE_S + 60, pid=2 ** 22 - 1)
        assert not _survives(), "proof of death is proof"

    def test_a_long_stale_record_is_still_reaped(self, world):
        """Even with a live pid. A record nobody has touched in an hour
        belongs to a session that is not coming back to it."""
        _record(quiet_for=P._REAP_AFTER_S + 60)
        assert not _survives()

    def test_a_pid_this_host_may_not_ask_about_is_left_alone(self, world,
                                                             monkeypatch):
        """Another user's process answers PermissionError, not
        ProcessLookupError. Guessing it dead costs its peers their
        messages, which is the expensive direction to be wrong in."""
        def _refuse(pid, sig):
            raise PermissionError(1, "Operation not permitted")

        monkeypatch.setattr(P.os, "kill", _refuse)
        _record(quiet_for=P._STALE_S + 60)
        assert _survives()

    def test_a_record_with_no_pid_is_not_proof_either(self, world):
        _record(quiet_for=P._STALE_S + 60, pid=0)
        assert _survives()


class TestTheMessageSurvivesTheSequence:
    def test_a_peer_can_still_leave_a_message_after_the_roster_was_read(
            self, world):
        """The whole defect, end to end, through the function
        `session_message` actually asks."""
        _record(quiet_for=P._STALE_S + 60)
        P.open_sessions()               # the roster, or any peer listing
        assert M.deliverable(KEY), (
            "the message is refused as 'no other open session' and lost, "
            "for a session that is running")


class TestAToolCallIsAHeartbeat:
    def test_touch_refreshes_the_clock_and_nothing_else(self, world):
        _record(quiet_for=P._STALE_S + 60)
        before = json.loads(P._path(KEY).read_text(encoding="utf-8"))
        assert P.touch(KEY) is True
        after = json.loads(P._path(KEY).read_text(encoding="utf-8"))
        assert after["updated_at"] > before["updated_at"]
        assert after["pid"] == os.getpid()
        for field in ("key", "session_id", "title", "workspace", "branch"):
            assert after[field] == before[field], (
                f"{field} changed; a heartbeat that rewrites the record "
                "blanks what the announcer knew and this side does not")
        assert KEY in [r.get("key") for r in P.open_sessions()], (
            "the point of the heartbeat is to be open again")

    def test_it_creates_nothing(self, world):
        """A session announces itself, with its title and its workspace.
        A record this made would be a nameless roster entry."""
        assert P.touch(KEY) is False
        assert not P._path(KEY).exists()
        assert P.open_sessions() == []

    def test_it_will_not_adopt_another_sessions_record(self, world):
        """A file whose own `key` disagrees is not this session's."""
        P._path(KEY).write_text(json.dumps({"key": "somebody-else"}),
                                encoding="utf-8")
        assert P.touch(KEY) is False

    def test_it_writes_at_most_once_per_heartbeat_interval(self, world):
        """Called once per tool call; a tool call is not a write."""
        _record(quiet_for=P._STALE_S + 60)
        assert P.touch(KEY) is True
        first = P._path(KEY).stat().st_mtime_ns
        for _ in range(20):
            assert P.touch(KEY) is True
        assert P._path(KEY).stat().st_mtime_ns == first, (
            "the rate limit is what makes this safe to call per tool call")

    def test_the_engine_touches_it_on_every_tool_call(self, world,
                                                      monkeypatch):
        """Driven through the engine's own path, not by calling touch:
        a heartbeat nobody calls is not a heartbeat."""
        from delfin.agent import engine as engine_mod

        _record(quiet_for=P._STALE_S + 60)
        touched: list[str] = []
        monkeypatch.setattr(P, "touch", lambda key: touched.append(key))

        eng = engine_mod.AgentEngine.__new__(engine_mod.AgentEngine)
        eng._trace_pending = []
        # kit_permissions is a read-only property reading the CLIENT's
        # permissions, so the key is supplied where the engine looks for it.
        eng.client = type("_C", (), {"_permissions": type(
            "_P", (), {"presence_key": KEY})()})()
        eng.trace_session = lambda: "s"
        eng._record_tool_trace("bash", "out")
        assert touched == [KEY]

    def test_a_session_with_no_presence_key_is_not_invented(self, world,
                                                            monkeypatch):
        from delfin.agent import engine as engine_mod

        called: list[str] = []
        monkeypatch.setattr(P, "touch", lambda key: called.append(key))
        eng = engine_mod.AgentEngine.__new__(engine_mod.AgentEngine)
        eng._trace_pending = []
        eng.client = type("_C", (), {"_permissions": type(
            "_P", (), {"presence_key": ""})()})()
        eng.trace_session = lambda: "s"
        eng._record_tool_trace("bash", "out")
        assert called == []

    def test_a_client_with_no_permissions_is_not_a_crash(self, world,
                                                         monkeypatch):
        """A non-KIT client has none; kit_permissions is then None."""
        from delfin.agent import engine as engine_mod

        called: list[str] = []
        monkeypatch.setattr(P, "touch", lambda key: called.append(key))
        eng = engine_mod.AgentEngine.__new__(engine_mod.AgentEngine)
        eng._trace_pending = []
        eng.client = type("_C", (), {})()
        eng.trace_session = lambda: "s"
        eng._record_tool_trace("bash", "out")      # must not raise
        assert called == []
