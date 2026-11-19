"""A damaged lab run checkpoint is kept, not replaced, when a run resumes.

``lab/workflow.py`` keeps one lab run's whole state in
``checkpoint.json``: what steps are done, where every sample physically is,
what waits for approval.  A damaged file is not reconstructable -- the
samples' positions are the physical truth of the bench, and the step log
says what was already dispensed.  A resume that quietly rebuilds the
checkpoint from the workflow document would claim the run never happened
while the bench says otherwise.

Measured defect shape: ``load_checkpoint`` answers "absent" and "damaged"
alike with ``None``, and ``Coordinator.run`` treats that as "start fresh" --
so one truncated byte and a resume throw the run's history away.
"""
from __future__ import annotations

import json

from _lab_support import FakeClock, lab_fixture  # noqa: F401
from delfin.lab import workflow as wf


def _aside(run_dir):
    return sorted(p for p in run_dir.iterdir() if "unreadable" in p.name)


def _minimal_checkpoint(run_id):
    # A resume lands on the approval-waiting branch only if the module locks
    # were already taken; simpler: make the run look freshly started but
    # completed -- then run() returns immediately and must not have touched
    # the damaged file.
    return {"run_id": run_id, "status": "completed", "done": ["stir"],
            "samples": {}, "log": [{"step": "stir", "outcome": {"ok": True}}],
            "measurements": [], "deviations": [], "pending": None,
            "started": 0.0, "coordinator": "c1", "person": "alice"}


def test_run_over_a_damaged_checkpoint_keeps_it_and_reports(lab, monkeypatch):
    run_id = "run-damaged"
    run_dir = wf.run_dir(run_id)
    (run_dir / "checkpoint.json").write_text(
        '{"run_id": "%s", "status": "runn' % run_id, encoding="utf-8")

    # The coordinator is built with an empty manifest/clients: the damaged
    # file must be detected before any node is contacted, so nothing here
    # ever reaches an edge.
    class _M:
        version = "v1"
        checksum = "sum"
        devices = {}
        modules = {}

    clk = FakeClock(0.0)
    coord = wf.Coordinator(_M(), {}, person="alice", name="c1",
                           clock=clk, sleep=clk.sleep)

    # The run may not raise where it used to answer None: it reports.
    # record.write is imported inside the functions, so patch its module.
    from delfin.lab import record as _record_mod
    reports = []
    monkeypatch.setattr(_record_mod, "write", lambda cp, m: reports.append(cp))
    cp = coord.run({"id": "wf-x", "samples": {}, "steps": []}, run_id)
    assert cp["status"] in ("refused", "failed"), cp
    assert cp["deviations"], "the loss is visible in the run's deviations"
    assert any("checkpoint" in json.dumps(d) for d in cp["deviations"])
    assert _aside(run_dir), "the damaged file is moved aside, kept"
    assert any("status" in p.read_text(encoding="utf-8")
               for p in _aside(run_dir)), "the old bytes survive"


def test_run_over_a_completed_checkpoint_returns_it_unchanged(lab):
    run_id = "run-good"
    cp_ok = _minimal_checkpoint(run_id)
    run_dir = wf.run_dir(run_id)
    (run_dir / "checkpoint.json").write_text(
        json.dumps(cp_ok), encoding="utf-8")

    class _M:
        version = "v1"
        checksum = "sum"
        devices = {}
        modules = {}

    wf.Coordinator(_M(), {}, person="alice", name="c1",
                   clock=FakeClock(0.0), sleep=lambda s: None).run(
        {"id": "wf-x", "samples": {}, "steps": []}, run_id)
    now = json.loads((run_dir / "checkpoint.json").read_text())
    assert now["status"] == "completed"
    assert not _aside(run_dir), "a healthy checkpoint is never moved aside"
