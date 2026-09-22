"""A published question outlives its session, and must not outlive it long.

Within minutes of the watch going live it reported eight questions
"waiting", the oldest for 3190 seconds. Every one of them belonged to a
process that no longer existed -- probes and unit drives whose brokers
were never answered, each leaving its record behind. A supervisor
reading that list learns nothing and acts on nothing, and the file this
test's own sibling already gives the reason for:

    a stale question read as open is worse than none at all

That sibling only handled the answered case. This is the other one: the
asker is gone and nobody will ever answer.

And it is the lesson session 4 was given, arriving in this corner too:
a pid is a number the system hands out again, so the record carries the
process START as well, and ``proc_identity`` -- the one reader of that
fact -- settles it. Where it cannot be asked (another machine, a kernel
that will not say) the question STAYS: guessing a session dead costs
somebody their running work, and that is the expensive direction.
"""

from __future__ import annotations

import json
import os

import pytest

from delfin.agent import terminal_confirm as tc


@pytest.fixture
def room(tmp_path, monkeypatch):
    monkeypatch.setattr(tc, "_PENDING_DIR", tmp_path / "pending")
    return tmp_path / "pending"


def _ask(session_key="runde2b-s1"):
    broker = tc.TerminalConfirmBroker(session_id="s-1", session_key=session_key)
    return broker, broker._enqueue(
        tc.ConfirmRequest(kind=tc.CONFIRM, tool="bash",
                          args={"command": "ls"}, preview="x"))


def _rewrite(room, **fields):
    path = next(room.glob("*.request.json"))
    record = json.loads(path.read_text(encoding="utf-8"))
    record.update(fields)
    path.write_text(json.dumps(record), encoding="utf-8")
    return path


class TestTheRecord:
    def test_it_carries_the_process_start(self, room):
        _ask()
        rec = tc.pending_at_terminals()[0]
        assert rec.get("pid") == os.getpid()
        assert rec.get("proc_start"), "a pid alone is a name a later process shares"


class TestWhatIsDropped:
    def test_a_live_asker_keeps_its_question(self, room):
        _ask()
        assert len(tc.pending_at_terminals()) == 1

    def test_a_dead_asker_loses_it(self, room):
        _ask()
        _rewrite(room, pid=999_999_999)
        assert tc.pending_at_terminals() == []

    def test_the_file_goes_too(self, room):
        _ask()
        _rewrite(room, pid=999_999_999)
        tc.pending_at_terminals()
        assert list(room.glob("*.request.json")) == [], "it would be read again"

    def test_a_recycled_pid_does_not_revive_it(self, room):
        # The number is in use -- by this very process -- but by a
        # different start than the record names.
        _ask()
        _rewrite(room, proc_start="1")
        assert tc.pending_at_terminals() == []

    def test_a_question_from_another_machine_stays(self, room):
        # It cannot be asked there, and guessing a session dead costs
        # somebody their running work.
        _ask()
        _rewrite(room, host="some-other-node", pid=999_999_999)
        assert len(tc.pending_at_terminals()) == 1

    def test_a_record_without_a_pid_stays(self, room):
        _ask()
        _rewrite(room, pid=0)
        assert len(tc.pending_at_terminals()) == 1


class TestItStillNeverRaises:
    def test_a_broken_record_is_skipped_not_fatal(self, room):
        _ask()
        next(room.glob("*.request.json")).write_text("{not json",
                                                     encoding="utf-8")
        assert tc.pending_at_terminals() == []

    def test_an_unreadable_room_reads_as_empty(self, tmp_path, monkeypatch):
        monkeypatch.setattr(tc, "_PENDING_DIR", tmp_path / "no" / "deeper")
        assert tc.pending_at_terminals() == []
