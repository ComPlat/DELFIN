"""A run record survives a write that fails midway.

``tools/_runtime.py`` keeps one JSON file per submitted run (id, status,
outputs, events, slurm job id).  Its readers all answer ``None`` for an
unreadable file and every writer re-reads before saving, so the
silent-overwrite defect is absent here -- but ``RunStore.save`` used a plain
``write_text``, which truncates the file first.  A crash or a full disk
between the truncation and the last byte leaves a record the readers then
answer ``None`` to: the run's history (events, outputs, slurm ids) is gone
while the calculation may still be running on a node.
"""
from __future__ import annotations

import json
import os
import stat
from pathlib import Path

from delfin.tools._runtime import RunRecord, RunStore


def _record(run_id="r1"):
    return RunRecord(id=run_id, kind="application", name="orca_sp",
                     created_at="t0", work_dir="w",
                     events=[{"t": "t0", "event": "run_started"}])


def test_a_failing_write_keeps_the_previous_record(tmp_path):
    store = RunStore(base=tmp_path)
    rec = _record()
    store.save(rec)
    before = (tmp_path / "r1.json").read_text(encoding="utf-8")

    real_open = os.open

    def broken_open(path, flags, *a, **k):
        if Path(path).name == "r1.json.tmp":
            raise OSError(28, "No space left on device")
        return real_open(path, flags, *a, **k)

    saved_open = os.open
    os.open = broken_open
    try:
        import pytest
        with pytest.raises(OSError):
            store.save(rec)
    finally:
        os.open = saved_open

    assert (tmp_path / "r1.json").read_text(encoding="utf-8") == before, (
        "the previous record must still be there after a failed write")
    assert not list(tmp_path.glob("*.tmp")), "no half-written file is left"


def test_a_healthy_save_still_round_trips(tmp_path):
    store = RunStore(base=tmp_path)
    rec = _record()
    rec.status = "success"
    store.save(rec)
    back = store.get("r1")
    assert back is not None and back.status == "success"
    assert back.events[0]["event"] == "run_started"


def test_a_saved_record_is_not_world_readable(tmp_path):
    store = RunStore(base=tmp_path)
    store.save(_record())
    # 0o600 was the helper default for secret stores; run records hold work
    # dirs and metrics, and the store sits in ~/.delfin -- private is right.
    assert stat.S_IMODE((tmp_path / "r1.json").stat().st_mode) == 0o600
