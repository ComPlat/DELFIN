"""`delfin-voila --keep` is the whole recipe in one flag.

Surviving a closed terminal took three things done by hand, in order,
every time: start tmux, start the dashboard inside it, then walk into
each window and switch "Keep session" on. Miss the third and a session
that pauses between turns is gone; miss the first and a dropped
connection takes everything. A recipe that must be followed exactly is
one that will be got wrong, and was.

  the flag starts tmux         once, reusing a session of that name
  already inside, carry on     no second wrapping
  no tmux on this host         say so; do not pretend it will survive
  the sessions arm themselves  including the record the server reads
  the switch still owns it     turning it off disarms, as always
"""

from __future__ import annotations

import os

import pytest

from delfin import cli_voila
from delfin.dashboard import session as S


def test_inside_tmux_it_does_not_wrap_again(monkeypatch):
    monkeypatch.setenv("TMUX", "/tmp/tmux-1000/default,1,0")
    assert cli_voila._reexec_under_tmux([]) is None


def test_a_second_pass_does_not_wrap_again(monkeypatch):
    monkeypatch.delenv("TMUX", raising=False)
    monkeypatch.setenv(cli_voila._TMUX_MARKER, "1")
    assert cli_voila._reexec_under_tmux([]) is None


def test_without_tmux_it_says_so_and_runs_here(monkeypatch, capsys):
    monkeypatch.delenv("TMUX", raising=False)
    monkeypatch.delenv(cli_voila._TMUX_MARKER, raising=False)
    monkeypatch.setattr(cli_voila.shutil, "which", lambda name: None,
                        raising=False)
    import shutil

    monkeypatch.setattr(shutil, "which", lambda name: None)
    assert cli_voila._reexec_under_tmux([]) is None
    said = capsys.readouterr().err
    assert "tmux is not installed" in said
    assert "closing this terminal ends the dashboard" in said


def test_it_starts_the_command_inside_tmux(monkeypatch):
    monkeypatch.delenv("TMUX", raising=False)
    monkeypatch.delenv(cli_voila._TMUX_MARKER, raising=False)
    import shutil
    import subprocess

    monkeypatch.setattr(shutil, "which", lambda name: "/usr/bin/tmux")
    seen = {}

    class _Done:
        returncode = 0

    def _run(cmd, env=None, **kw):
        seen["cmd"] = cmd
        seen["marked"] = (env or {}).get(cli_voila._TMUX_MARKER)
        return _Done()

    monkeypatch.setattr(subprocess, "run", _run)
    assert cli_voila._reexec_under_tmux(["--port", "8899"]) == 0
    assert seen["cmd"][:5] == ["tmux", "new", "-A", "-s",
                               cli_voila.TMUX_SESSION]
    assert "--port" in seen["cmd"] and "8899" in seen["cmd"]
    assert seen["marked"] == "1", "the child must not wrap itself again"


# -- the half that is easiest to forget -----------------------------------

def test_the_sessions_of_a_keep_dashboard_start_armed(monkeypatch):
    monkeypatch.setenv("DELFIN_KEEP_SESSIONS", "1")
    assert S._keep_by_default() is True


@pytest.mark.parametrize("value", ["", "0", "no", "off"])
def test_otherwise_they_start_as_they_always_did(monkeypatch, value):
    monkeypatch.setenv("DELFIN_KEEP_SESSIONS", value)
    assert S._keep_by_default() is False


def test_an_armed_start_writes_the_record_the_server_reads(monkeypatch,
                                                           tmp_path):
    """The observer fires on a CHANGE; a control that starts armed never
    changed, so without this the switch would look on and the session
    would still die with its window."""
    pytest.importorskip("ipywidgets")
    monkeypatch.setenv("DELFIN_KEEP_SESSIONS", "1")
    monkeypatch.setattr(S, "RECORD_DIR", str(tmp_path))
    monkeypatch.setattr(S, "kernel_id", lambda: "k-armed")
    monkeypatch.setattr(S, "session_name", lambda: "armed-one")
    # A kernel announces itself on the SERVER's terminal, by writing to its
    # parent's file descriptor -- which pytest does not capture. Without
    # this the run prints "Session armed-one is kept" into whatever
    # terminal started the suite.
    monkeypatch.setattr(S, "_server_stdout", lambda: None)

    strip = S.build_status_strip()
    assert strip is not None
    assert list(tmp_path.glob("*.json")), "the record is on disk"
    assert S.is_kept_alive() is True
