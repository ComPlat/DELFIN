"""A damaged PROVENANCE.json is kept, not replaced, when a new one is written.

Provenance_record.py holds a run's ORIGIN: which files (sha256), which program
version, which DELFIN commit produced a number.  That is not reconstructable
from anything else in the run folder -- the numbers a result publishes stay on
disk while the record of how they came to be is gone.  That is why the writer
here must not overwrite a file it could not read (measured defect class of
2026-11-17: the agent memory store and the credential store).

Before this test existed, ``write_record`` did a plain ``write_text``: it
truncated the damaged file first, and a provenance record replaced by a new
one silently erased the previous record of how the previous numbers were made.
"""
from __future__ import annotations

import json

from delfin.provenance_record import PROVENANCE_FILE, write_record


def _run(tmp_path):
    run = tmp_path / "run"
    run.mkdir()
    (run / "CONTROL.txt").write_text(
        "# CONTROL\nsmiles=CCO\n", encoding="utf-8")
    return run


def _aside(run):
    return sorted(p for p in run.iterdir() if "unreadable" in p.name)


def test_writing_a_record_keeps_a_damaged_previous_one(tmp_path, monkeypatch):
    """The old bytes survive a write that lands on a file we cannot parse."""
    run = _run(tmp_path)
    damaged = run / PROVENANCE_FILE
    damaged.write_text('{"schema_version": "1", "delfin"', encoding="utf-8")

    def fake_state(_root):
        return {"commit": "abc123", "dirty": False}

    def fake_env(_exe):
        return {"error": "provenance_unavailable"}

    monkeypatch.setattr(
        "delfin.provenance_record.cbatch_code_state", fake_state)
    monkeypatch.setattr(
        "delfin.provenance_record.cbatch_probe_env", fake_env)

    written = write_record(
        {"schema_version": "1", "result": {"type": "calculation",
                                           "created": "now"},
         "step_id": "pv1:x", "chain": []},
        run,
    )
    assert written.name == PROVENANCE_FILE
    assert json.loads(written.read_text(encoding="utf-8"))["step_id"] == "pv1:x"

    kept = _aside(run)
    assert len(kept) == 1, "the damaged record is moved aside, not deleted"
    assert "delfin" in kept[0].read_text(encoding="utf-8"), (
        "what the damaged record held must still be there")


def test_writing_a_record_over_a_good_previous_one_replaces_it(tmp_path):
    """A readable record is a normal replace: no aside copy is left behind."""
    run = _run(tmp_path)
    (run / PROVENANCE_FILE).write_text(
        '{"schema_version": "1", "result": {"type": "calculation"}}',
        encoding="utf-8")
    write_record({"schema_version": "1", "result": {"type": "calculation",
                                                    "created": "now"},
                  "step_id": "pv1:x", "chain": []}, run)
    assert json.loads((run / PROVENANCE_FILE).read_text())["step_id"] == "pv1:x"
    assert not _aside(run), "a healthy file is not moved aside"
