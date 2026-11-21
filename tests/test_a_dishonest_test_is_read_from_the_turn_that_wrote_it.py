"""The scan's record source reads what the turn wrote, through a real write.

Input: a workspace, the product's own write path
(``api_client._DocToolExecutor.execute`` — the path a model's call takes),
and the undo journal that path fills.

Output: ``honesty_scan.turn_records`` returns records in ``check``'s
``(path, start, end, source)`` order for exactly the changes made after the
mark, with paths relative to the workspace and an added-line range that
excludes what the file already contained.

What this file does NOT prove: that the scan's verdict reaches a caller at
turn end. That needs the engine hunk, which the Self-Modification Guard
refuses without a human approval; ``scratch/PR.md`` records it and
``scratch/engine.diff`` carries it. These tests pin the half that can be
proven here — that a real write becomes a correctly shaped record, and that
``scan_for_dishonest_tests`` actually fires on it. A record source that
silently returned nothing would pass a "the scan runs" test and measure
nothing, which is the failure mode this wave is about.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from delfin.agent import change_journal as CJ
from delfin.agent import honesty_scan
from delfin.agent import verify_guard as vg
from delfin.agent.api_client import KitToolPermissions, _doc_executor


class _Engine:
    """The two things the record source reads off the engine."""

    def __init__(self, ws: Path, session: str):
        self.repo_dir = ws
        self._session = session

    def _edit_signal_key(self) -> str:
        return self._session


@pytest.fixture
def box(tmp_path, monkeypatch):
    home = tmp_path / "home"
    home.mkdir()
    monkeypatch.setenv("HOME", str(home))
    monkeypatch.setattr(Path, "home", lambda: home)
    ws = tmp_path / "ws"
    (ws / "tests").mkdir(parents=True)
    perms = KitToolPermissions(mode="acceptEdits", workspace=str(ws))
    perms.task_session_id = f"hs-{tmp_path.name}"[:40]
    return ws, perms


def _write(perms, ws: Path, rel: str, content: str) -> None:
    """Write through the product's tool path, as a model's call does."""
    out = str(_doc_executor.execute(
        "write_file", {"path": str(ws / rel), "content": content}, perms))
    assert '"error"' not in out, out
    assert (ws / rel).read_text() == content


def test_a_real_write_becomes_a_record_for_the_scan(box):
    ws, perms = box
    engine = _Engine(ws, perms.task_session_id)
    assert honesty_scan.mark(engine) == 0        # nothing journalled yet

    _write(perms, ws, "tests/test_new.py", "def test_x():\n    assert 1\n")

    records = honesty_scan.turn_records(engine, since=0)
    assert [r[0] for r in records] == ["tests/test_new.py"]
    path, start, end, source = records[0]
    assert source.startswith("def test_x()")     # the file's NEW text
    assert (start, end) == (1, 2)                # created: all of it is new


def test_the_range_is_what_this_turn_added_not_the_whole_file(box):
    """A pre-existing ``except:`` earlier in the file must stay outside the
    range, or kind 1 would report old history as this turn's defect."""
    ws, perms = box
    engine = _Engine(ws, perms.task_session_id)
    old = ("import os\n"
           "\n"
           "def load(p):\n"
           "    try:\n"
           "        return open(p).read()\n"
           "    except:\n"
           "        pass\n")
    _write(perms, ws, "tests/test_old.py", old)      # turn 1 created it
    mark = honesty_scan.mark(engine)

    # Turn 2 appends an honest test; the bare except above is untouched.
    _write(perms, ws, "tests/test_old.py",
           old + "\ndef test_ok():\n    assert 1 == 1\n")

    records = honesty_scan.turn_records(engine, since=mark)
    assert len(records) == 1
    _path, start, _end, _source = records[0]
    assert start > 7, f"the range re-covers the old handler (start={start})"
    # And the scanner therefore says nothing about the old handler.
    assert vg.scan_for_dishonest_tests(records) == []


def test_a_dishonest_change_the_scan_does_name(box):
    """The record source is wired to the scan, not merely shaped like it."""
    ws, perms = box
    engine = _Engine(ws, perms.task_session_id)

    _write(perms, ws, "delfin/mod.py",
           "def load(p):\n"
           "    try:\n"
           "        return open(p).read()\n"
           "    except:\n"
           "        pass\n")

    findings = vg.scan_for_dishonest_tests(
        honesty_scan.turn_records(engine, since=0))
    assert findings, "the scan saw no text — the records are not reaching it"
    assert findings[0].path == "delfin/mod.py"
    assert "bare" in findings[0].reason
    # And the report the model would read names it.
    text = vg.test_honesty_feedback(findings)
    assert "delfin/mod.py" in text


def test_a_committed_probe_file_is_named_by_path(box):
    """Kind 2 is path-anchored at ``^tests/``: the record source must hand
    the scan a workspace-relative path, or this kind can never fire."""
    ws, perms = box
    engine = _Engine(ws, perms.task_session_id)

    _write(perms, ws, "tests/test_zz_probe.py", "def test_probe():\n    pass\n")

    findings = vg.scan_for_dishonest_tests(
        honesty_scan.turn_records(engine, since=0))
    assert [f.path for f in findings] == ["tests/test_zz_probe.py"]


def test_an_earlier_turns_change_is_not_re_reported(box):
    ws, perms = box
    engine = _Engine(ws, perms.task_session_id)
    _write(perms, ws, "tests/test_zz_probe.py", "def test_probe():\n    pass\n")
    mark = honesty_scan.mark(engine)
    assert mark >= 1

    _write(perms, ws, "tests/test_later.py", "def test_y():\n    assert 1\n")

    records = honesty_scan.turn_records(engine, since=mark)
    assert [r[0] for r in records] == ["tests/test_later.py"]


def test_a_broken_journal_yields_no_records_and_never_raises(box, monkeypatch):
    """The scan must degrade to quiet, not to a crash, on store trouble."""
    ws, perms = box
    engine = _Engine(ws, perms.task_session_id)
    _write(perms, ws, "tests/test_new.py", "def test_x():\n    assert 1\n")

    def _boom(*_a, **_k):
        raise OSError("journal unreadable")

    monkeypatch.setattr(CJ, "list_changes", _boom)
    assert honesty_scan.mark(engine) == 0
    assert honesty_scan.turn_records(engine, since=0) == []


def test_an_engine_with_no_session_yields_no_records(box):
    """A client without a session id has nothing to report on."""
    ws, _perms = box
    assert honesty_scan.turn_records(_Engine(ws, ""), since=0) == []
    assert honesty_scan.mark(_Engine(ws, "")) == 0
