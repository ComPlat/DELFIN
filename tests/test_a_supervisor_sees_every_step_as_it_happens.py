"""A supervisor sees every step a session takes, as it takes it.

Supervising five sessions on 2026-09-20 meant reading the audit log with
grep and a hand-written python one-liner, once per question. Two things
were wanted and neither was at hand: what a session is doing right now,
and -- the one that matters -- what it TRIED to do and was refused.

The audit log already holds both. It is the only cross-process record
there is: security_events lives in the process that recorded it, so one
session cannot read another's. This reads the log the way a supervisor
reads it.

  every step, newest last     the order things happened in
  a refusal is marked         "denied" is the line to stop on, and the
                              reason is carried, not summarised away
  by session                  five worktrees produce five streams into
                              one file
  an offset, not a count      so a follower prints each record once,
                              however many arrive between two looks
"""

from __future__ import annotations

import json

import pytest

from delfin.agent import audit_log as AL


@pytest.fixture
def log(tmp_path):
    return tmp_path / "audit.log"


def _write(path, records):
    with path.open("a", encoding="utf-8") as fh:
        for r in records:
            fh.write(json.dumps(r) + "\n")


def _rec(session="s1", decision="ok", tool="bash", command="ls", **kw):
    out = {"ts": "2026-09-20T10:00:00Z", "session_id": session,
           "tool": tool, "decision": decision, "command": command}
    out.update(kw)
    return out


class TestTheStream:
    def test_it_reads_what_is_there(self, log):
        _write(log, [_rec(), _rec(command="pwd")])
        rows, _ = AL.steps_since(0, log_path=log)
        assert [r["command"] for r in rows] == ["ls", "pwd"]

    def test_the_offset_makes_each_record_arrive_once(self, log):
        _write(log, [_rec(command="one")])
        first, offset = AL.steps_since(0, log_path=log)
        _write(log, [_rec(command="two")])
        second, _ = AL.steps_since(offset, log_path=log)
        assert [r["command"] for r in first] == ["one"]
        assert [r["command"] for r in second] == ["two"]

    def test_nothing_new_is_nothing(self, log):
        _write(log, [_rec()])
        _, offset = AL.steps_since(0, log_path=log)
        rows, again = AL.steps_since(offset, log_path=log)
        assert rows == [] and again == offset

    def test_a_missing_log_is_empty_not_an_error(self, tmp_path):
        rows, offset = AL.steps_since(0, log_path=tmp_path / "nope.log")
        assert rows == [] and offset == 0

    def test_a_half_written_line_is_skipped_not_fatal(self, log):
        log.write_text('{"ts": "x", "tool": "bash"}\n{"broken\n', encoding="utf-8")
        rows, _ = AL.steps_since(0, log_path=log)
        assert len(rows) == 1

    def test_one_session_can_be_singled_out(self, log):
        _write(log, [_rec(session="aaa"), _rec(session="bbb")])
        rows, _ = AL.steps_since(0, log_path=log, session="bbb")
        assert [r["session_id"] for r in rows] == ["bbb"]

    def test_only_the_refusals(self, log):
        _write(log, [_rec(), _rec(decision="denied", reason="outside roots"),
                     _rec(decision="error")])
        rows, _ = AL.steps_since(0, log_path=log, denied_only=True)
        assert len(rows) == 1
        assert rows[0]["reason"] == "outside roots"


class TestTheRendering:
    def test_a_refusal_is_marked_and_keeps_its_reason(self):
        line = AL.render_step(_rec(decision="denied",
                                   reason="path is outside the allowed roots",
                                   command="cat /etc/shadow"))
        assert "DENIED" in line
        assert "outside the allowed roots" in line

    def test_an_ordinary_step_is_quiet(self):
        line = AL.render_step(_rec())
        assert "DENIED" not in line
        assert "ls" in line

    def test_a_failing_command_is_not_a_refusal(self):
        # exit code 1 from grep is not the gate saying no, and reading it
        # as one is how a first draft of this reported 46 false alarms.
        line = AL.render_step(_rec(decision="error", reason="exited with code 1"))
        assert "DENIED" not in line

    def test_the_session_is_named(self):
        assert "abcd1234"[:8] in AL.render_step(_rec(session="abcd1234ef"))

    def test_a_write_names_its_path(self):
        line = AL.render_step(_rec(tool="write_file", command=None,
                                   path="/w/x.py"))
        assert "/w/x.py" in line

    def test_rendering_never_raises(self):
        assert isinstance(AL.render_step({}), str)
