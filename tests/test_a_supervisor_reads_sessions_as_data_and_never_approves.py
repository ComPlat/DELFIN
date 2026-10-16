"""A supervisor reads running sessions as data, and never answers an approval.

``delfin-agent ops`` gives a supervisor (a person with a script, or an
operator agent running a wave) events and state instead of terminal screens,
and three actions that cannot widen anything: a message shown as from the
operator, an escalation to the person, a pause. There is no approve.

Everything here runs on a redirected home with the real recorders writing:
``session_presence.announce``, ``turn_metrics`` lines, the audit log, the
approval room and ``session_messages``.
"""
from __future__ import annotations

import json
import subprocess
import sys
import time
from pathlib import Path

import pytest

from delfin.agent import (file_confirm, ops_interface,
                          session_messages, session_presence, session_store,
                          terminal_confirm, turn_metrics)


@pytest.fixture
def home(tmp_path, monkeypatch):
    h = tmp_path / "home"
    (h / ".delfin").mkdir(parents=True)
    monkeypatch.setattr(Path, "home", lambda: h)
    monkeypatch.setattr(session_messages, "_DIR", h / ".delfin" / "session_inbox")
    monkeypatch.setattr(session_presence, "_DIR", h / ".delfin" / "session_presence")
    monkeypatch.setattr(session_store, "_SESSIONS_DIR", h / ".delfin" / "agent_sessions")
    monkeypatch.setattr(terminal_confirm, "_PENDING_DIR", h / ".delfin" / "terminal_confirmations")
    monkeypatch.setattr(turn_metrics, "_DIR", h / ".delfin" / "turn_metrics")
    from delfin.agent import session_pause
    monkeypatch.setattr(session_pause, "_DIR", h / ".delfin" / "session_pause")
    rooms = h / ".delfin" / "confirmations"
    rooms.mkdir()
    monkeypatch.setattr(file_confirm, "requests_dir", lambda: rooms)
    monkeypatch.delenv("DELFIN_AGENT_SESSION", raising=False)
    try:
        from delfin.agent import approval_guard
        monkeypatch.delenv(approval_guard.ENV_MARKER, raising=False)
    except Exception:
        pass
    ws = tmp_path / "ws"
    ws.mkdir()
    subprocess.run(["git", "init", "-q", str(ws)], check=True)
    return h, ws, rooms


def _open(key, sid, ws):
    session_presence.announce(key, session_id=sid, title=f"title {key}",
                              workspace=str(ws))


def _turn(h, sid, **fields):
    d = h / ".delfin" / "turn_metrics"
    d.mkdir(parents=True, exist_ok=True)
    rec = {"ts": time.time(), "tool_calls": 3, "stopped": False, "error": ""}
    rec.update(fields)
    with open(d / f"{sid}.jsonl", "a", encoding="utf-8") as fh:
        fh.write(json.dumps(rec) + "\n")


def _audit(h, **fields):
    rec = {"ts": "2026-10-16T10:00:00Z", "tool": "bash", "decision": "ok",
           "mode": "default", "exit_code": 0}
    rec.update(fields)
    with open(h / ".delfin" / "audit.log", "a", encoding="utf-8") as fh:
        fh.write(json.dumps(rec) + "\n")


def _request(rooms, rid, sid, command):
    (rooms / f"{rid}.request.json").write_text(json.dumps({
        "id": rid, "session_id": sid, "tool": "bash", "command": command,
        "asked_at": time.time(), "timeout_s": 900}), encoding="utf-8")


def _kinds(events):
    return [e["kind"] for e in events]


# --------------------------------------------------------------------------
# state
# --------------------------------------------------------------------------

def test_an_open_session_is_listed_with_its_state(home):
    h, ws, rooms = home
    _open("wave-s11", "a" * 32, ws)
    _turn(h, "a" * 32, tool_calls=7)
    _request(rooms, "r1", "a" * 32, "git clean -f x")
    rows = ops_interface.sessions()
    assert [r["name"] for r in rows] == ["wave-s11"]
    row = rows[0]
    assert row["session_id"] == "a" * 32
    assert row["open_approvals"] == ["r1"]
    assert row["last_turn"]["tool_calls"] == 7
    assert row["idle_s"] is not None and row["idle_s"] < 60
    assert ops_interface.session("wave-s11")["name"] == "wave-s11"
    assert ops_interface.session("aaaaaaaa")["name"] == "wave-s11"


# --------------------------------------------------------------------------
# events
# --------------------------------------------------------------------------

def test_a_new_cursor_starts_now_and_then_reports_what_happens(home):
    h, ws, rooms = home
    sid = "b" * 32
    _open("wave-s12", sid, ws)
    _turn(h, sid)                                   # history, not an event
    state = {}
    assert ops_interface.poll(state) == []

    _turn(h, sid, stopped=True)
    _request(rooms, "r2", sid, "find data -delete")
    session_messages.send("operator", "phase 2 committed", from_key="wave-s12")
    _audit(h, session_id=sid, decision="denied", command="git push origin x")
    _audit(h, session_id=sid, command="git commit -F .gate/m.txt -- tests/t.py",
           cwd=str(ws))
    ev = ops_interface.poll(state)
    kinds = _kinds(ev)
    assert "turn_ended" in kinds and "approval_open" in kinds
    assert "message" in kinds and "denied" in kinds and "commit" in kinds
    turn = next(e for e in ev if e["kind"] == "turn_ended")
    assert turn["session"] == "wave-s12" and turn["end"] == "stopped"
    msg = next(e for e in ev if e["kind"] == "message")
    assert msg["to"] == "operator" and msg["text"] == "phase 2 committed"
    denied = next(e for e in ev if e["kind"] == "denied")
    assert denied["session"] == "wave-s12" and "git push" in denied["detail"]
    assert all(e["schema"] == ops_interface.SCHEMA for e in ev)
    seqs = [e["seq"] for e in ev]
    assert seqs == sorted(seqs) and len(set(seqs)) == len(seqs)

    (rooms / "r2.answer.json").write_text("{}", encoding="utf-8")
    assert _kinds(ops_interface.poll(state)) == ["approval_closed"]


def test_a_reader_that_restarts_neither_repeats_nor_drops(home):
    h, ws, _ = home
    sid = "c" * 32
    _open("wave-s13", sid, ws)
    ops_interface.save_cursor(_first(), "r")
    for i in range(5):
        _turn(h, sid, tool_calls=i)
    first = _drain("r")
    for i in range(3):
        _turn(h, sid, tool_calls=10 + i)
    second = _drain("r")                            # a new process, same cursor
    tc = [e["tool_calls"] for e in first + second if e["kind"] == "turn_ended"]
    assert tc == [0, 1, 2, 3, 4, 10, 11, 12]
    assert _drain("r") == []
    seqs = [e["seq"] for e in first + second]
    assert seqs == list(range(seqs[0], seqs[0] + len(seqs)))


def _first():
    state = {}
    ops_interface.poll(state)
    return state


def _drain(cursor):
    state = ops_interface.load_cursor(cursor)
    ev = ops_interface.poll(state)
    ops_interface.save_cursor(state, cursor)
    return ev


def test_a_session_idle_long_enough_is_reported_once(home):
    h, ws, _ = home
    sid = "d" * 32
    _open("wave-s14", sid, ws)
    _turn(h, sid, ts=time.time() - 3600)
    state = {}
    ops_interface.poll(state, idle_s=600)
    ev = ops_interface.poll(state, idle_s=600)
    assert [e["kind"] for e in ev] == ["session_idle"]
    assert ev[0]["idle_s"] >= 3600 - 5
    assert ops_interface.poll(state, idle_s=600) == []
    _turn(h, sid)                                   # works again
    assert _kinds(ops_interface.poll(state, idle_s=600)) == ["turn_ended"]


def test_a_session_waiting_for_a_person_is_not_called_idle(home):
    h, ws, rooms = home
    sid = "e" * 32
    _open("wave-s15", sid, ws)
    _turn(h, sid, ts=time.time() - 3600)
    _request(rooms, "r3", sid, "git clean -f x")
    state = {}
    ops_interface.poll(state, idle_s=600)
    assert "session_idle" not in _kinds(ops_interface.poll(state, idle_s=600))


@pytest.mark.parametrize("spec,hit", [
    ("kind=approval_open", True),
    ("kind=approval_open,session=wave-s11", True),
    ("kind=approval_open,session=wave-s99", False),
    ("idle_s>=1800", True),
    ("idle_s>=4000", False),
    ("kind=", False),
    ("nonsense", False),
])
def test_the_wait_filter(spec, hit):
    ev = {"kind": "approval_open", "session": "wave-s11", "idle_s": 1900}
    assert ops_interface.matches(ev, spec) is hit


def test_wait_returns_the_first_matching_event_and_times_out(home):
    h, ws, rooms = home
    sid = "f" * 32
    _open("wave-s16", sid, ws)
    ops_interface.save_cursor(_first(), "w")
    _request(rooms, "r4", sid, "git clean -f x")
    ev = ops_interface.wait("kind=approval_open", cursor="w", timeout=0.2,
                            interval=0.05)
    assert ev and ev["id"] == "r4"
    assert ops_interface.wait("kind=commit", cursor="w", timeout=0.2,
                              interval=0.05) is None


# --------------------------------------------------------------------------
# actions
# --------------------------------------------------------------------------

def test_a_message_from_the_operator_is_never_the_user(home):
    h, ws, _ = home
    _open("wave-s17", "g" * 32, ws)
    msg = ops_interface.send("wave-s17", "please report your hash")
    taken = session_messages.take("wave-s17")
    assert [m["id"] for m in taken] == [msg["id"]]
    shown = session_messages.render(taken[0])
    assert "not from the user" in shown
    assert "operator" in shown and "please report your hash" in shown


def test_a_secret_in_a_message_does_not_travel(home):
    h, ws, _ = home
    token = "ghp_" + "A1b2C3d4E5f6G7h8I9j0K1l2M3n4O5p6Q7r8"
    msg = ops_interface.send("wave-s18", f"use {token} to push")
    assert token not in json.dumps(msg)


def test_escalate_puts_a_waiting_approval_before_the_person_and_leaves_it_open(
        home):
    h, ws, rooms = home
    sid = "h" * 32
    _request(rooms, "r5", sid, "git clean -f x")
    eid = ops_interface.escalate("r5", "deletes a file; your call")
    from delfin.agent import attention
    assert any(e.get("id") == eid for e in attention.list_pending())
    assert [a["id"] for a in ops_interface.pending_approvals()] == ["r5"]
    with pytest.raises(ops_interface.OpsRefused):
        ops_interface.escalate("nope")


def test_pause_marks_the_session(home):
    from delfin.agent import session_pause
    ops_interface.pause("wave-s19", "stalled for an hour")
    assert session_pause.is_paused("wave-s19")


@pytest.mark.parametrize("call", [
    lambda: ops_interface.send("wave-s20", "x"),
    lambda: ops_interface.pause("wave-s20", "x"),
    lambda: ops_interface.escalate("r6"),
])
def test_an_agent_session_cannot_act_as_the_operator(home, monkeypatch, call):
    from delfin.agent import approval_guard
    monkeypatch.setenv(approval_guard.ENV_MARKER, "nacht-s20")
    with pytest.raises(ops_interface.OpsRefused):
        call()


# --------------------------------------------------------------------------
# no way to a yes
# --------------------------------------------------------------------------

def test_no_ops_command_can_answer_an_approval(home, monkeypatch):
    """Every ops subcommand, run through the real CLI parser, with the two
    functions that write an answer replaced by tripwires."""
    h, ws, rooms = home
    sid = "i" * 32
    _open("wave-s21", sid, ws)
    _request(rooms, "r7", sid, "git clean -f x")
    called = []
    monkeypatch.setattr(file_confirm, "answer",
                        lambda *a, **k: called.append(("file", a)) or True)
    monkeypatch.setattr(terminal_confirm, "answer_waiting",
                        lambda *a, **k: called.append(("term", a)) or True)
    from delfin.agent import cli
    parser = cli.build_parser()
    argvs = [
        ["ops", "sessions"], ["ops", "session", "wave-s21"],
        ["ops", "approvals"], ["ops", "events"],
        ["ops", "wait", "--until", "kind=commit", "--timeout", "0.1",
         "--interval", "0.05"],
        ["ops", "send", "wave-s21", "--text", "yes"],
        ["ops", "escalate", "r7", "--note", "approve please"],
        ["ops", "pause", "wave-s21", "--reason", "yes"],
    ]
    for argv in argvs:
        args = parser.parse_args(argv)
        rc = args.func(args)
        assert rc in (0, 2, 3), argv
    assert called == []
    assert [a["id"] for a in ops_interface.pending_approvals()] == ["r7"]
    assert not (rooms / "r7.answer.json").exists()


def test_the_ops_commands_are_exactly_these():
    """A new ops subcommand is a decision: list it here, and never one that
    answers."""
    from delfin.agent import cli
    import argparse
    parser = cli.build_parser()
    sub = next(a for a in parser._actions
               if isinstance(a, argparse._SubParsersAction))
    ops = sub.choices["ops"]
    ops_sub = next(a for a in ops._actions
                   if isinstance(a, argparse._SubParsersAction))
    names = set(ops_sub.choices)
    assert names == {"status", "sessions", "session", "approvals", "events",
                     "wait", "send", "escalate", "pause"}
    assert not names & {"approve", "answer", "allow", "grant", "accept"}


# --------------------------------------------------------------------------
# shapes seen in the three-session acceptance run
# --------------------------------------------------------------------------

def test_a_commit_at_the_end_of_a_command_chain_is_reported(home):
    h, ws, rooms = home
    sid = "c" * 32
    _open("probe-a", sid, ws)
    state = {}
    ops_interface.poll(state)
    _audit(h, session_id=sid, cwd=str(ws), command=(
        'git add util.py test_util.py && git commit -m "Add mean() with test"'))
    _audit(h, session_id=sid, cwd=str(ws), command="echo git commit")
    _audit(h, session_id=sid, cwd=str(ws), exit_code=1,
           command="git status; git commit -m nothing")
    commits = [e for e in ops_interface.poll(state) if e["kind"] == "commit"]
    assert len(commits) == 1 and commits[0]["session"] == "probe-a"


def test_a_message_taken_before_the_next_poll_is_still_reported(home):
    h, ws, rooms = home
    _open("probe-c", "d" * 32, ws)
    state = {}
    ops_interface.poll(state)
    ops_interface.send("probe-c", "please add a docstring")
    assert session_messages.take("probe-c")          # gone from the inbox
    session_messages.send("operator", "done", from_key="probe-c")
    msgs = [e for e in ops_interface.poll(state) if e["kind"] == "message"]
    assert [(m["from"], m["to"]) for m in msgs] == [
        ("operator", "probe-c"), ("probe-c", "operator")]
    assert msgs[0]["taken_at"] and msgs[1]["text"] == "done"
    assert [e for e in ops_interface.poll(state) if e["kind"] == "message"] == []


def test_a_session_whose_heartbeat_pauses_is_not_reported_as_exited(home):
    h, ws, rooms = home
    _open("probe-b", "e" * 32, ws)
    state = {}
    ops_interface.poll(state)
    rec_path = next((h / ".delfin" / "session_presence").glob("*.json"))
    rec = json.loads(rec_path.read_text())
    rec["updated_at"] = time.time() - 4 * 3600       # waits at a dialog
    rec_path.write_text(json.dumps(rec))
    assert _kinds(ops_interface.poll(state)) == []
    gone = subprocess.Popen([sys.executable, "-c", "pass"])
    gone.wait()
    rec["pid"] = gone.pid                            # the process has ended
    rec_path.write_text(json.dumps(rec))
    assert _kinds(ops_interface.poll(state)) == ["session_exited"]


def test_a_closed_approval_names_its_session_and_tool(home):
    h, ws, rooms = home
    sid = "f" * 32
    _open("probe-b", sid, ws)
    state = {}
    ops_interface.poll(state)
    _request(rooms, "r9", sid, "rm notes.txt")
    assert _kinds(ops_interface.poll(state)) == ["approval_open"]
    (rooms / "r9.request.json").rename(rooms / "r9.answered.json")
    closed = [e for e in ops_interface.poll(state) if e["kind"] == "approval_closed"]
    assert closed and closed[0]["session"] == "probe-b" and closed[0]["tool"] == "bash"
    # the cursor written by the earlier version (a plain id list) still loads
    state["approvals"] = ["old"]
    assert [e["id"] for e in ops_interface.poll(state)
            if e["kind"] == "approval_closed"] == ["old"]


def test_a_session_that_just_took_a_message_is_not_reported_idle(home):
    h, ws, rooms = home
    sid = "9" * 32
    _open("probe-c", sid, ws)
    _turn(h, sid, ts=time.time() - 3600)             # last turn an hour ago
    state = {}
    ops_interface.poll(state, idle_s=600)
    ops_interface.send("probe-c", "next step")
    session_messages.take("probe-c")                 # a new turn begins
    assert "session_idle" not in _kinds(ops_interface.poll(state, idle_s=600))
    _audit(h, session_id=sid, command="ls")          # and it works on
    assert "session_idle" not in _kinds(ops_interface.poll(state, idle_s=600))
