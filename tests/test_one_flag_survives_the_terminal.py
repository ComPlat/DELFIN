"""`delfin-voila --keep`: the old recipe, now the default plus one thing.

Surviving a closed terminal took tmux, the dashboard inside it, and the
Keep switch in every window. ``delfin-voila`` alone now survives the
terminal (it runs detached: test_the_dashboard_outlives_its_terminal.py),
so ``--keep`` is that default and no longer starts tmux. What it still
adds is the third step: every window the dashboard opens starts with
Keep dashboard on.

  no tmux                    the flag starts no tmux and needs none
  said once                  one line says what the flag does now
  the sessions arm themselves  including the record the server reads
  the switch still owns it   turning it off disarms, as always
"""

from __future__ import annotations

import sys

import pytest

from delfin import cli_voila
from delfin.dashboard import session as S


class _Tty:
    def isatty(self):
        return True


def test_keep_starts_no_tmux_and_says_what_it_does_now(monkeypatch, capsys):
    for name in ("TMUX", "STY", "SLURM_JOB_ID", "DELFIN_VOILA_FOREGROUND",
                 "DELFIN_VOILA_DETACHED"):
        monkeypatch.delenv(name, raising=False)
    for name in ("DELFIN_KEEP_SESSIONS", "DELFIN_LAUNCH_CWD"):
        monkeypatch.setenv(name, "-")
        monkeypatch.delenv(name)
    monkeypatch.setattr(sys, "stdin", _Tty())
    monkeypatch.setattr(cli_voila, "_voila_is_available", lambda: True)
    from delfin.agent import process_guard
    monkeypatch.setattr(process_guard, "protect",
                        lambda *a: (_ for _ in ()).throw(RuntimeError("off")))
    monkeypatch.setattr(cli_voila.subprocess, "run",
                        lambda cmd, *a, **k: pytest.fail(f"started {cmd}"))
    started = []
    monkeypatch.setattr(cli_voila, "_start_detached",
                        lambda args, argv: started.append(list(argv)) or 0)
    assert cli_voila.main(["--keep", "--no-browser"]) == 0
    assert started == [["--keep", "--no-browser"]]
    err = capsys.readouterr().err
    assert err.count("--keep") == 1
    assert "tmux" in err and "Keep dashboard" in err
    import os
    assert os.environ.get("DELFIN_KEEP_SESSIONS") == "1"


def test_there_is_no_tmux_re_execution_left():
    assert not hasattr(cli_voila, "_reexec_under_tmux")


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
