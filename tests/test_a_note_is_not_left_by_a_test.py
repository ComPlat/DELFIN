"""The note that leads back to a running dashboard is not a test's to touch.

``delfin-voila`` leaves one file behind that says which node it runs on,
on which port, with which token -- the whole way back, on a cluster where
the next login lands somewhere else. Seven tests drive the launcher, and
every one of them wrote that note with its own pid; the launcher's
``atexit`` hook then deleted it as pytest exited. Measured, not feared:
during a full run the real note read ``pid 498725``, which was the pytest
process, and after the run there was no note at all. A dashboard serving
in tmux at that moment would have been reported as "not running" by
``delfin-agent where``, and its address would have been gone.

Two things were wrong, and both are fixed here:

  the path was frozen at import   ~/.delfin/... resolved once, before any
                                  redirection could reach it
  the note had no owner           whoever ran last deleted it, however
                                  little it had to do with them

The suite already redirects every user-wide sink the product names, so
the note joins that table rather than growing a mechanism of its own --
and this file asserts the redirect from the PRODUCT's table, because the
copy of it that stood in conftest.py silently shadowed the import.
"""

from __future__ import annotations

import importlib
import inspect
import json
import os
from pathlib import Path

import pytest

from delfin.agent import state_paths, where


# -- the redirect covers what the product declares --------------------------

def _resolve(mod, attr, tmp_path):
    fn = getattr(mod, attr)
    try:
        takes_root = len(inspect.signature(fn).parameters) >= 1
    except (TypeError, ValueError):       # a builtin or a mock
        takes_root = False
    return Path(fn(tmp_path) if takes_root else fn())


@pytest.mark.parametrize("entry", state_paths.USER_STATE_RESOLVERS,
                         ids=lambda e: f"{e[0].split('.')[-1]}.{e[1]}")
def test_every_resolver_the_product_names_is_redirected(entry, tmp_path):
    """Not "the ones conftest happens to list" -- the product's own table.

    A second, literal copy of this table used to sit in conftest.py and
    overwrite the imported one on the very next statement. The copies
    agreed, so nothing failed; an entry added to the product table would
    simply never have reached the suite.
    """
    mod_name, attr, _rel = entry
    try:
        mod = importlib.import_module(mod_name)
    except Exception:                      # an optional dependency is absent
        pytest.skip(f"{mod_name} does not import here")
    if not callable(getattr(mod, attr, None)):
        pytest.skip(f"{mod_name}.{attr} is not a resolver here")

    real_home = Path.home()
    got = _resolve(mod, attr, tmp_path)
    assert real_home not in got.parents, (
        f"{mod_name}.{attr} still resolves inside the user's home ({got}); "
        "the suite's redirect does not reach it")


def test_the_note_is_one_of_them():
    """The entry that this file exists for, named explicitly, so that
    removing it from the table fails here and not somewhere quiet."""
    assert any(e[0] == "delfin.agent.where" and e[1] == "record_path"
               for e in state_paths.USER_STATE_RESOLVERS)


def test_the_path_is_decided_per_call_not_at_import(monkeypatch, tmp_path):
    """A module constant freezes ``$HOME`` at import time.

    Moving ``Path.home`` is also how the suite's redirect steps aside, so
    what answers here is the resolver the product ships -- and it follows
    the home it is given, which a constant could not.
    """
    assert not hasattr(where, "RECORD_PATH"), (
        "a module-level path is resolved once, before anything can "
        "redirect it")
    home = tmp_path / "elsewhere"
    home.mkdir()
    monkeypatch.setattr(Path, "home", staticmethod(lambda: home))
    assert home in where.record_path().parents


# -- a note belongs to the process that wrote it ----------------------------

def test_a_note_someone_else_wrote_is_not_removed(monkeypatch, tmp_path):
    note = tmp_path / "dashboard_here.json"
    monkeypatch.setattr(where, "record_path", lambda: note)
    note.write_text(json.dumps({
        "host": where.socket.gethostname(),
        "pid": os.getpid() + 1,          # not this process
        "proc_start": "1",
        "port": 8867,
    }), encoding="utf-8")

    where.withdraw_dashboard()

    assert note.exists(), (
        "a process removed a note it did not write -- that is the address "
        "of a dashboard that is still serving")


def test_a_process_removes_the_note_it_wrote_itself(monkeypatch, tmp_path):
    note = tmp_path / "dashboard_here.json"
    monkeypatch.setattr(where, "record_path", lambda: note)

    assert where.announce_dashboard(port=8867, token="t") == str(note)
    assert note.exists()

    where.withdraw_dashboard()
    assert not note.exists()


def test_the_same_pid_in_another_life_is_not_the_dashboard(monkeypatch,
                                                           tmp_path):
    """A pid is handed out again. The note says it is running; the number
    belongs to something else entirely by then."""
    note = tmp_path / "dashboard_here.json"
    monkeypatch.setattr(where, "record_path", lambda: note)
    where.announce_dashboard(port=8867, token="t")

    record = json.loads(note.read_text(encoding="utf-8"))
    assert record["proc_start"], "a fingerprint is recorded at all"
    assert where._still_running(record) is True

    record["proc_start"] = str(int(record["proc_start"] or 0) + 99) \
        if record["proc_start"].isdigit() else "a different life"
    assert where._still_running(record) is False

    # ... and that note is not this process's to delete, either.
    note.write_text(json.dumps(record), encoding="utf-8")
    where.withdraw_dashboard()
    assert note.exists()


def test_a_note_from_before_the_fingerprint_still_works(monkeypatch,
                                                        tmp_path):
    """An older note carries no ``proc_start``. It must still be readable,
    and its own process must still be able to take it back."""
    note = tmp_path / "dashboard_here.json"
    monkeypatch.setattr(where, "record_path", lambda: note)
    note.write_text(json.dumps({
        "host": where.socket.gethostname(),
        "pid": os.getpid(),
        "port": 8867,
    }), encoding="utf-8")

    assert where.dashboard().get("running") is True
    where.withdraw_dashboard()
    assert not note.exists()
