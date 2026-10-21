"""Probe (test_zz_): why does the worker child write no metrics line?

The builder's in-process test greens (record written through repl.run_turn
hook repl.py:908-914); the same turn driven through the REAL worker child
(tests/tui_worker_client.py, private HOME + DELFIN_SCRATCH_STATE via
conftest.child_env) leaves agent_metrics.jsonl ABSENT under w.home. This
probe captures the child's stderr and the scratch/home layout to say where
the record went -- or why the write never happened.
"""
from __future__ import annotations

from tui_worker_client import WorkerClient


def test_probe_worker_metrics_location(tmp_path):
    with WorkerClient(tmp_path) as w:
        start = len(w.events)
        w.request("session.prompt", {"text": "hello"})
        w.wait_for("turn.ended", after=start)

    # Every place a line could be.
    for name, root in (("w.home", w.home), ("scratch", w.scratch)):
        d = root / ".delfin"
        if d.is_dir():
            listing = sorted(p.name for p in d.iterdir())
        else:
            listing = []
        print(f"[probe] {name}={root} .delfin entries={listing}")
        p = root / ".delfin" / "agent_metrics.jsonl"
        print(f"[probe] metrics file at {p}: exists={p.exists()}"
              + (f" bytes={p.stat().st_size}" if p.exists() else ""))
    print(f"[probe] child stderr:\n{w.stderr()[-2000:]}")
    # Deliberately always true: the probe prints what it measured.
    assert True
