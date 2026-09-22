"""A run with no session id must not share one undo journal with all of them.

The undo journal is per session by design -- ``~/.delfin/undo/<sid>/``.
The sanitiser that makes a directory name out of the id ends with

    str(session_id or "") or "session"

so every run that has NO id lands in one bucket literally called
"session". Measured on this machine 2026-09-19: that bucket held 528
records from **163 different path roots** -- other people's pytest runs,
scratch checkouts, a ``/tmp/dgit/repo`` from two days earlier.

That is exactly the shape a headless agent has. A session driven with
``-p`` and no id would write there and, worse, READ there: ``undo_changes``
answers with the entries it could not use, naming each one --

    "skipped": [{"path": "/tmp/dgit/repo/f.txt", "reason": "outside workspace"}]

-- so a sandboxed run learns file paths from workspaces it can never
reach. The revert itself was correctly bounded; the disclosure was not.

Two changes, because one of them alone would still leak:

* A run without an id gets a bucket of its OWN, named after this process
  and the moment it started -- the same fingerprint ``proc_identity``
  gives everything else that has to tell one life of a pid from the next.
* A path outside the workspace is counted, never named. The agent needs
  to know that something was not reverted; it does not need the path, and
  the path is the part that belongs to somebody else.
"""

from __future__ import annotations

import json
import os

import pytest

from delfin.agent import change_journal as CJ


def test_two_processes_without_an_id_do_not_share_a_bucket():
    """The name carries this process. Another process gets another name."""
    mine = CJ._safe_session_id("")
    assert mine != "session", (
        "every id-less run in the world writes to one directory called "
        "'session' -- 528 records from 163 roots when this was written")
    assert str(os.getpid()) in mine, mine


def test_the_same_process_keeps_the_same_bucket():
    """Undo has to reach back over several calls of one run."""
    assert CJ._safe_session_id("") == CJ._safe_session_id(None)


def test_a_real_id_is_untouched():
    assert CJ._safe_session_id("abc-123") == "abc-123"


def test_traversal_is_still_impossible():
    assert "/" not in CJ._safe_session_id("../../etc/passwd")
    assert ".." not in CJ._safe_session_id("../../etc/passwd")


def test_a_skipped_path_outside_the_workspace_is_not_named(tmp_path,
                                                           monkeypatch):
    """It is counted, so the agent knows. It is not named, because the
    name belongs to a workspace this run cannot see."""
    import inspect
    src = inspect.getsource(CJ)
    marker = '{"path": path_str, "reason": "outside workspace"}'
    assert marker not in src, (
        "the path of a file outside the workspace is handed back to the "
        "model; it is the one field that discloses another run's work")
    assert "outside workspace" in src, "the fact must still be reported"
