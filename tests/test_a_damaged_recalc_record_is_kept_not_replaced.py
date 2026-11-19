"""A damaged recalc record is kept when a recalc edit records the old CONTROL.

``recalc_control`` holds which CONTROL.txt a completed run was computed with
(``.delfin_last_run.json``), which is what makes a recalc able to say what an
edit changes.  It is not reconstructable: the CONTROL the run actually used
is replaced the moment the Recalc tab saves a new one.

Measured defect shape: ``remember_before_edit`` answers "no record yet" from a
read that conflates "absent" and "damaged" (``previous_run``), so one truncated
byte in the record file makes it write a fresh record over it -- and the
CONTROL the finished run was computed with is gone, silently.
"""
from __future__ import annotations

import json

from delfin import recalc_control


def _record(job_dir):
    return job_dir / recalc_control.RECORD_NAME


def test_remember_before_edit_keeps_a_damaged_record(tmp_path, monkeypatch):
    job_dir = tmp_path
    rec = _record(job_dir)
    rec.write_text('{"control": "smiles=CCO\\nPAL 8", "inp', encoding="utf-8")

    warnings = []
    monkeypatch.setattr(
        recalc_control.logger, "warning",
        lambda msg, *a: warnings.append(msg % a if a else msg))

    recorded = recalc_control.remember_before_edit(job_dir, "smiles=CCC\nPAL 4")
    assert recorded is True, "the edit proceeds; only the file's fate differs"

    aside = sorted(p for p in job_dir.iterdir() if "unreadable" in p.name)
    assert len(aside) == 1, "the damaged record is moved aside, not overwritten"
    assert "smiles=CCO" in aside[0].read_text(encoding="utf-8"), (
        "what the damaged record held must still be there")

    new = json.loads(rec.read_text(encoding="utf-8"))
    assert new["control"] == "smiles=CCC\nPAL 4", "the new record is written"
    assert warnings and "delfin_last_run" in warnings[0], (
        "the warning names the file, so a person can act")


def test_remember_before_edit_over_a_readable_record_does_nothing(tmp_path):
    job_dir = tmp_path
    recalc_control.record_completed_run(job_dir, "smiles=CCO\nPAL 8")
    recorded = recalc_control.remember_before_edit(job_dir, "smiles=CCC")
    assert recorded is False, "a completed run already recorded its CONTROL"
    aside = [p for p in job_dir.iterdir() if "unreadable" in p.name]
    assert not aside, "a healthy record is never moved aside"
