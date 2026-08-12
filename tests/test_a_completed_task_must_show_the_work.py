"""Marking a task completed is a claim, and it is now checked.

THE INCIDENT. ``TaskStore.update`` guarded the ``in_progress``
transition and had no branch for ``completed`` at all: any task, of any
shape, became done for the asking. The one advisory check ran AFTER the
write (the store already said completed on disk while the note was
composed), and it asked whether an extension appeared as a SUBSTRING of
the session's observed-files ledger joined into one blob — a ledger that
counts files merely READ. So "Create the PDF report" was satisfied by
having opened any unrelated .pdf earlier in the session.

Worse, the ledger it read (``_observed_files_session``) lives on the
CLIENT and the check runs on the document EXECUTOR, which never has that
attribute. ``getattr(self, ..., None) or ()`` therefore evaluated to the
empty tuple on every call, so the check could only ever fire: every
completed task whose subject said Bericht / Brief / Tabelle / report /
PDF was told its artefact was missing, whether or not it existed. A
check that accuses honest work is worse than no check — it teaches the
model to phrase subjects so nothing can key on them.

The repo's own end-to-end test proved the gap from the other side: it
created "Add mylib/optimizers/wrapper.py", marked it completed, and
asserted only that the status string round-tripped.

WHAT IS CHECKED NOW, per shape, from ledgers that record only what
happened: a named path must appear as a WRITE (the change journal, which
records mutations and nothing else) unless the subject only promised to
look at it; an artefact noun must match a written file per path; a test
task needs a green run since the task started; an edit task needs a
mutation. Anything else is recorded as ``verified: "unchecked"`` — an
honest unknown, distinguishable from a checked one.
"""

from __future__ import annotations

import json
import time
from pathlib import Path

import pytest

from delfin.agent.api_client import (
    KitToolPermissions,
    _doc_executor,
    _paths_in_text,
    _unmet_artifact,
    check_completion_claim,
)


# ---------------------------------------------------------------------------
# The blob substring, and the read/write merge that made it hollow
# ---------------------------------------------------------------------------

def test_an_extension_inside_another_filename_is_not_that_file():
    """"".pdf" in "notes.pdf.bak"" is true and neither is a PDF."""
    assert _unmet_artifact("Create the PDF report", ["notes.pdf.bak"]) == "pdf"
    assert _unmet_artifact("Create the PDF report", ["out/report.pdf"]) == ""


def test_a_pdf_that_was_only_read_does_not_satisfy_a_pdf_task():
    """The headline hole: the evidence set now holds writes only, so an
    unrelated PDF opened earlier in the session proves nothing."""
    res = check_completion_claim(
        "PDF-Bericht erstellen",
        changes=[],                       # nothing written
        observed=["/data/handbuch.pdf"],  # a PDF, read
    )
    assert res["verdict"] == "unmet"
    assert res["kind"] == "artifact"


def test_a_pdf_that_was_written_satisfies_it():
    res = check_completion_claim(
        "PDF-Bericht erstellen",
        changes=[{"path": "/proj/out/bericht.pdf", "ts": 10.0, "created": True}],
        observed=[],
    )
    assert res["verdict"] == "verified"


# ---------------------------------------------------------------------------
# A path in the subject
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("text,expected", [
    ("Add mylib/optimizers/wrapper.py", ["mylib/optimizers/wrapper.py"]),
    ("Refactor delfin/agent/engine.py and tests/test_x.py",
     ["delfin/agent/engine.py", "tests/test_x.py"]),
    ("Fix the 1.2 regression in v3.4", []),      # versions are not paths
    ("Rückfrage klären (z.B. wegen der Frist)", []),
    ("Update the README", []),                   # no extension, no claim
])
def test_only_real_paths_are_read_out_of_a_subject(text, expected):
    assert _paths_in_text(text) == expected


def test_a_named_file_that_was_never_written_is_unmet():
    res = check_completion_claim(
        "Add mylib/optimizers/wrapper.py",
        changes=[{"path": "/proj/other.py", "ts": 5.0, "created": False}],
        observed=["/proj/mylib/optimizers/wrapper.py"],   # only read
    )
    assert res["verdict"] == "unmet"
    assert res["kind"] == "path_unwritten"
    assert "wrapper.py" in res["note"]


def test_a_named_file_that_was_written_is_verified():
    res = check_completion_claim(
        "Add mylib/optimizers/wrapper.py",
        changes=[{"path": "/home/u/proj/mylib/optimizers/wrapper.py",
                  "ts": 5.0, "created": True}],
    )
    assert res["verdict"] == "verified"
    assert res["kind"] == "path_write"


def test_a_reading_task_is_satisfied_by_reading():
    """"analysiere core.py" is honestly done with no write at all. A check
    that demanded one would be training the model to stop naming files."""
    res = check_completion_claim(
        "analysiere mylib/core.py",
        changes=[],
        observed=["/proj/mylib/core.py"],
    )
    assert res["verdict"] == "verified"
    assert res["kind"] == "path_read"


def test_a_reading_task_without_a_read_ledger_is_unchecked_not_unmet():
    res = check_completion_claim(
        "analysiere mylib/core.py", changes=[], observed=None)
    assert res["verdict"] == "unchecked"


# ---------------------------------------------------------------------------
# Test tasks, edit tasks, and the honest unknown
# ---------------------------------------------------------------------------

def test_a_test_task_needs_a_green_run_since_it_started():
    green = [{"tool": "run_tests", "failed": 0, "status": "ok", "ts": 200.0}]
    assert check_completion_claim(
        "Regressionstests laufen lassen", tests=green, window_start=100.0,
    )["verdict"] == "verified"


def test_a_test_task_with_only_red_runs_is_unmet():
    red = [{"tool": "run_tests", "failed": 3, "status": "failed", "ts": 200.0}]
    res = check_completion_claim(
        "Run the test suite", tests=red, window_start=100.0)
    assert res["verdict"] == "unmet"
    assert res["kind"] == "tests_red"


def test_a_test_task_whose_only_run_predates_it_is_unmet():
    stale = [{"tool": "run_tests", "failed": 0, "status": "ok", "ts": 50.0}]
    res = check_completion_claim(
        "Run the test suite", tests=stale, window_start=100.0)
    assert res["verdict"] == "unmet"
    assert res["kind"] == "tests_stale"


def test_a_test_task_without_a_reachable_ledger_is_unchecked():
    assert check_completion_claim(
        "Run the test suite", tests=None)["verdict"] == "unchecked"


def test_an_edit_task_needs_a_mutation_somewhere():
    res = check_completion_claim("Refactor the optimizer wrapper", changes=[])
    assert res["verdict"] == "unmet"
    assert res["kind"] == "no_change"


def test_an_edit_task_whose_work_predates_the_status_flip_still_counts():
    """Editing first and flipping the status afterwards is doing the work,
    not faking it. The window is provenance, not a trap."""
    res = check_completion_claim(
        "Refactor the optimizer wrapper",
        changes=[{"path": "/p/w.py", "ts": 10.0, "created": False}],
        window_start=100.0,
    )
    assert res["verdict"] == "verified"
    assert res["kind"] == "journal_session"


@pytest.mark.parametrize("subject", [
    "Rückfrage mit Jerome klären",
    "Decide on the review order",
    "Mit dem Kunden telefonieren",
])
def test_a_subject_with_nothing_to_key_on_is_unchecked(subject):
    res = check_completion_claim(subject, changes=[], observed=[], tests=[])
    assert res["verdict"] == "unchecked"


def test_a_broken_input_does_not_raise():
    assert check_completion_claim(None, None)["verdict"] == "unchecked"
    assert check_completion_claim("x", changes=object())["verdict"] in (
        "unchecked", "unmet", "verified")


# ---------------------------------------------------------------------------
# End to end, through the tool the model actually calls
# ---------------------------------------------------------------------------

@pytest.fixture
def agent_home(monkeypatch, tmp_path):
    """Redirect the change journal (``~/.delfin/undo``) into the tmp dir."""
    monkeypatch.setattr(Path, "home", lambda: tmp_path / "home")
    (tmp_path / "home").mkdir()
    ws = tmp_path / "ws"
    ws.mkdir()
    from delfin.agent import agent_tasks
    agent_tasks._STORES.clear()
    return ws


def _perms(ws: Path, sid: str = "sess-verify") -> KitToolPermissions:
    return KitToolPermissions(
        workspace=ws, mode="acceptEdits", task_session_id=sid)


def _create_and_start(perms, subject: str) -> int:
    out = json.loads(_doc_executor.execute(
        "task_create", {"subject": subject}, perms))
    tid = out["task"]["id"]
    _doc_executor.execute(
        "task_update", {"task_id": tid, "status": "in_progress"}, perms)
    return tid


def test_completing_an_unwritten_file_task_is_flagged_and_recorded(agent_home):
    """The end-to-end case the old test asserted a bare status string for."""
    perms = _perms(agent_home)
    tid = _create_and_start(perms, "Add mylib/optimizers/wrapper.py")
    out = json.loads(_doc_executor.execute(
        "task_update", {"task_id": tid, "status": "completed"}, perms))
    assert "wrapper.py" in out["note"]
    assert out["task"]["verified"] == "unmet"


def test_completing_a_file_task_that_was_actually_written_says_nothing(
        agent_home):
    """The honest completion must be accepted, or the check is just noise
    the model learns to route around."""
    perms = _perms(agent_home)
    tid = _create_and_start(perms, "Add mylib/optimizers/wrapper.py")
    target = agent_home / "mylib" / "optimizers" / "wrapper.py"
    written = _doc_executor.execute("write_file", {
        "path": str(target), "content": "def wrap():\n    return 1\n",
    }, perms)
    assert target.exists(), written
    out = json.loads(_doc_executor.execute(
        "task_update", {"task_id": tid, "status": "completed"}, perms))
    assert "note" not in out, out
    assert out["task"]["verified"] == "verified"


def test_the_check_runs_before_the_store_says_completed(agent_home):
    """The verdict is written in the same update as the status, so the
    record never passes through a state where it claims done and knows
    nothing about it."""
    perms = _perms(agent_home)
    tid = _create_and_start(perms, "Rückfrage mit Jerome klären")
    out = json.loads(_doc_executor.execute(
        "task_update", {"task_id": tid, "status": "completed"}, perms))
    task = out["task"]
    assert task["status"] == "completed"
    assert task["verified"] == "unchecked"
    assert task["verify_note"]


def test_an_artifact_task_is_judged_on_files_the_session_wrote(agent_home):
    """The German subject the office runs actually use."""
    perms = _perms(agent_home)
    tid = _create_and_start(perms, "PDF-Bericht für Juni erstellen")
    # A docx exists, the PDF never got written — the original incident.
    _doc_executor.execute("write_file", {
        "path": str(agent_home / "bericht.docx"), "content": "x",
    }, perms)
    out = json.loads(_doc_executor.execute(
        "task_update", {"task_id": tid, "status": "completed"}, perms))
    assert "pdf" in out["note"].lower()
    assert out["task"]["verified"] == "unmet"


def test_the_journal_window_is_compared_in_the_same_timezone(agent_home):
    """The change journal writes LOCAL naive timestamps and the task store
    writes UTC ones. Compared raw, every change in the window looks hours
    out of it (or hours inside a window it predates)."""
    from delfin.agent.api_client import _journal_ts_epoch, _task_ts_epoch
    from delfin.agent import change_journal as _cj
    from delfin.agent.agent_tasks import _now_iso

    before = time.time()
    _cj.record_change("tz-check", tool="write_file", path="/p/x.py",
                      old_text=None, new_text="y")
    rec = _cj.list_changes("tz-check")[-1]
    after = time.time() + 1
    assert before - 1 <= _journal_ts_epoch(rec["ts"]) <= after
    assert before - 1 <= _task_ts_epoch(_now_iso()) <= after
