"""`where --attach` and `where --ssh-config`: the way back, in one command.

The note already knows the node, the port, the tmux session and the token.
Telling the user to type three things by hand is making them re-derive what
the machine already wrote down. Two flags close that gap:

  --attach       walk into the tmux session the note names -- on this node
                 directly, from another node over ssh with the port forward
                 filled in from the note
  --ssh-config   print (or, with --write, append) one ssh config entry:
                 HostName, RemoteCommand tmux new -A -s, LocalForward

Both are universal: a host without tmux or ssh says so in one line and
exits 1 -- no guess, no hostname baked in. The token lives in the note so
the note alone is the way back; it must never reach the ssh config, which
is a file people copy around, and never reach the attach command, which a
terminal logs.
"""

from __future__ import annotations

import json
import socket
import sys
from pathlib import Path

from delfin.agent import where


def _note(tmp_path, **over):
    record = {
        "host": socket.gethostname(),
        "pid": 1,
        "port": 8867,
        "token": "SECRET-TOKEN-7f3a",
        "tmux": "delfin",
    }
    record.update(over)
    note = tmp_path / "dashboard_here.json"
    note.write_text(json.dumps(record), encoding="utf-8")
    monkey_record = json.loads(note.read_text())
    return monkey_record


# -- attach ------------------------------------------------------------------

def test_attach_inside_the_named_session_says_so(tmp_path, monkeypatch,
                                                 capsys):
    record = _note(tmp_path)
    monkeypatch.setattr(where, "_tmux_session", lambda: "delfin")

    rc = where.run_attach(record)

    assert rc == 0
    out = capsys.readouterr().out
    assert "already" in out.lower()
    assert "token" not in out and "SECRET" not in out


def test_attach_on_this_node_runs_tmux_directly(tmp_path, monkeypatch):
    record = _note(tmp_path)
    monkeypatch.setattr(where, "_tmux_session", lambda: "")
    monkeypatch.setattr(where.shutil, "which", lambda n: "/usr/bin/" + n)
    seen = {}
    monkeypatch.setattr(where.os, "execvp",
                        lambda name, argv: seen.update(name=name, argv=argv))

    where.run_attach(record)

    assert seen["argv"] == ["tmux", "attach", "-t", "delfin"]


def test_attach_from_another_node_carries_the_forward(tmp_path, monkeypatch):
    record = _note(tmp_path, host="login-node-2")
    monkeypatch.setattr(where, "_tmux_session", lambda: "")
    monkeypatch.setattr(where.shutil, "which", lambda n: "/usr/bin/" + n)
    seen = {}
    monkeypatch.setattr(where.os, "execvp",
                        lambda name, argv: seen.update(name=name, argv=argv))

    where.run_attach(record)

    assert seen["argv"] == ["ssh", "-t", "-L", "8867:localhost:8867",
                            "login-node-2", "tmux", "attach", "-t", "delfin"]


def test_attach_without_the_tool_says_so_in_one_line(tmp_path, monkeypatch,
                                                     capsys):
    record = _note(tmp_path)
    monkeypatch.setattr(where, "_tmux_session", lambda: "")
    monkeypatch.setattr(where.shutil, "which", lambda n: None)
    ran = []
    monkeypatch.setattr(where.os, "execvp", lambda *a: ran.append(a))

    rc = where.run_attach(record)

    assert rc == 1
    assert not ran
    assert "tmux" in capsys.readouterr().out


def test_attach_without_a_note_says_so(tmp_path, capsys):
    rc = where.run_attach({})
    assert rc == 1
    assert "no" in capsys.readouterr().out.lower()


def test_attach_with_no_tmux_in_the_note_refuses(tmp_path, capsys):
    rc = where.run_attach(_note(tmp_path, tmux=""))
    assert rc == 1


# -- ssh-config ----------------------------------------------------------------

def test_the_entry_names_the_node_the_session_and_the_port(tmp_path):
    entry = where.ssh_config_entry(_note(tmp_path, host="uc3n991"))
    assert "HostName uc3n991" in entry
    assert "RemoteCommand tmux new -A -s delfin" in entry
    assert "LocalForward 8867 localhost:8867" in entry
    assert "RequestTTY" in entry


def test_the_entry_never_carries_the_token(tmp_path):
    entry = where.ssh_config_entry(_note(tmp_path))
    assert "SECRET" not in entry, "the token reached a file people copy around"
    assert "token" not in entry.lower()


def test_write_appends_when_the_alias_is_new(tmp_path, monkeypatch, capsys):
    home = tmp_path / "home"
    (home / ".ssh").mkdir(parents=True)
    (home / ".ssh" / "config").write_text("Host other\n    HostName x\n")
    monkeypatch.setattr(Path, "home", staticmethod(lambda: home))

    rc = where.write_ssh_config(_note(tmp_path, host="uc3n991"))

    assert rc == 0
    text = (home / ".ssh" / "config").read_text()
    assert "Host other" in text, "the existing config was clobbered"
    assert "HostName uc3n991" in text


def test_write_refuses_an_existing_entry_and_says_what_would_change(
        tmp_path, monkeypatch, capsys):
    home = tmp_path / "home"
    (home / ".ssh").mkdir(parents=True)
    (home / ".ssh" / "config").write_text(
        "Host delfin\n    HostName old-node\n")
    monkeypatch.setattr(Path, "home", staticmethod(lambda: home))

    rc = where.write_ssh_config(_note(tmp_path, host="uc3n991"))

    assert rc == 1
    text = (home / ".ssh" / "config").read_text()
    assert "old-node" in text, "an existing entry was changed"
    out = capsys.readouterr().out
    assert "old-node" in out and "uc3n991" in out, (
        "the refusal does not show what would change")


def test_write_creates_the_config_when_there_is_none(tmp_path, monkeypatch):
    home = tmp_path / "home"
    home.mkdir()
    monkeypatch.setattr(Path, "home", staticmethod(lambda: home))

    rc = where.write_ssh_config(_note(tmp_path))

    assert rc == 0
    assert (home / ".ssh" / "config").exists()
    mode = (home / ".ssh" / "config").stat().st_mode & 0o777
    assert mode == 0o600


# -- liveness: the contract, pinned -------------------------------------------

def test_liveness_on_another_host_is_cannot_tell(tmp_path, monkeypatch):
    monkeypatch.setattr(where, "record_path",
                        lambda: _note(tmp_path, host="another-node",
                                      pid=12345).get("_path",
                                                     tmp_path / "x"))
    # write the note properly and ask
    note = tmp_path / "dashboard_here.json"
    note.write_text(json.dumps({
        "host": "another-node", "pid": 12345, "port": 8867,
    }), encoding="utf-8")
    monkeypatch.setattr(where, "record_path", lambda: note)

    record = where.dashboard()
    assert record["running"] is None, "a pid on another node was guessed at"
    out = where.format_text(record, [])
    assert "cannot tell" in out


def test_liveness_gone_when_the_pid_is_dead_here(tmp_path, monkeypatch):
    note = tmp_path / "dashboard_here.json"
    note.write_text(json.dumps({
        "host": socket.gethostname(), "pid": -1, "port": 8867,
    }), encoding="utf-8")
    monkeypatch.setattr(where, "record_path", lambda: note)

    record = where.dashboard()
    assert record["running"] is False
    assert "gone" in where.format_text(record, [])


# -- the CLI wires the flags through -------------------------------------------

def test_where_main_dispatches_attach(monkeypatch):
    called = {}
    def fake_attach(r):
        called["r"] = r
        return 0
    monkeypatch.setattr(where, "run_attach", fake_attach)
    monkeypatch.setattr(where, "dashboard", lambda: {"host": "h"})

    rc = where.main(["--attach"])
    assert rc == 0
    assert called["r"] == {"host": "h"}


def test_where_main_prints_the_config_entry(monkeypatch, capsys):
    monkeypatch.setattr(where, "dashboard",
                        lambda: {"host": "h", "port": 8867, "tmux": "delfin"})
    rc = where.main(["--ssh-config"])
    assert rc == 0
    assert "RemoteCommand tmux new -A -s delfin" in capsys.readouterr().out
