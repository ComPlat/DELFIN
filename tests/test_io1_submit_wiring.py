"""IO1 wiring: cluster_tools.submit runs hpc_io_check on each planned
job's script-ish files BEFORE _ask; a "refuse" finding stops the submit
with the finding's text; "warn" findings ride along in the tool result.

Real call site: the product calls cluster_tools.submit(executor, perms,
cluster=..., jobs=...) -- these tests drive that function with a fake
executor/perms/conn, like the other cluster_tools tests do.
"""
from __future__ import annotations

import json
from pathlib import Path

import pytest

from delfin.agent import cluster_tools
from delfin import hpc_io_check


def _job_folder(tmp_path: Path, script_text: str) -> Path:
    folder = tmp_path / "job1"
    folder.mkdir(exist_ok=True)
    (folder / "CONTROL.txt").write_text("PAL 1\n", encoding="utf-8")
    (folder / "input.txt").write_text("positions\n", encoding="utf-8")
    (folder / "run.sh").write_text(script_text, encoding="utf-8")
    return folder


def _run_submit(monkeypatch, tmp_path, folder, *, lustre=True):
    """Call the real cluster_tools.submit with a faked cluster plumbing;
    everything except the check wiring is stubbed out."""
    calls = {}

    from types import SimpleNamespace
    job = SimpleNamespace(folder=folder)
    batch = SimpleNamespace(ok=True, install="", commit="abc",
                            core_hours=0.0, notes=[], problems=[],
                            jobs=[job],
                            preview=lambda self=None: "preview")

    class FakeConn:
        name = "testcluster"

        def __init__(self):
            pass

        def is_live(self):
            return True

    import delfin.remote.jobs as remote_jobs

    monkeypatch.setattr(remote_jobs, "plan",
                        lambda *a, **k: batch)
    monkeypatch.setattr(remote_jobs, "submit",
                        lambda *a, **k: calls.setdefault("submitted", []))
    monkeypatch.setattr(cluster_tools, "_connect",
                        lambda *a, **k: (FakeConn(), None))
    monkeypatch.setattr(cluster_tools, "_ask", lambda *a, **k: None)
    monkeypatch.setattr(cluster_tools, "_watch_submitted", lambda *a, **k: [])
    monkeypatch.setattr(cluster_tools, "_audit", lambda *a, **k: None)

    class FakePerms:
        workspace = tmp_path
        task_session_id = "t"

    return cluster_tools.submit(None, FakePerms(), cluster="testcluster",
                                jobs=[{"folder": str(folder), "kind": "delfin",
                                       "input": ""}]), calls


WALK = "#!/bin/bash\ndu -sh /pfs/data/workspace/\n"

LUSTRE_MOUNTS = ("mgmt@o2ib:fs /pfs lustre rw 0 0\n")


def test_refuse_finding_stops_submit(monkeypatch, tmp_path):
    folder = _job_folder(tmp_path, WALK)
    monkeypatch.setattr(hpc_io_check.MountsSource, "read",
                        lambda self: LUSTRE_MOUNTS)
    result, calls = _run_submit(monkeypatch, tmp_path, folder)
    parsed = json.loads(result)
    assert parsed.get("error") or parsed.get("code")
    assert "tree_walk_on_lustre" in result
    assert "submitted" not in calls


def test_warn_finding_rides_along(monkeypatch, tmp_path):
    folder = _job_folder(tmp_path,
                         "#!/bin/bash\nfor f in a b; do python /pfs/env/x.py; done\n")
    monkeypatch.setattr(hpc_io_check.MountsSource, "read",
                        lambda self: LUSTRE_MOUNTS)
    result, calls = _run_submit(monkeypatch, tmp_path, folder)
    parsed = json.loads(result)
    assert "lustre_interpreter_in_loop" in result
    assert parsed.get("cluster") == "testcluster"


def test_clean_script_submits_without_findings(monkeypatch, tmp_path):
    folder = _job_folder(tmp_path, "#!/bin/bash\nsleep 1\n")
    monkeypatch.setattr(hpc_io_check.MountsSource, "read", lambda self: "")
    result, calls = _run_submit(monkeypatch, tmp_path, folder)
    assert "tree_walk" not in result
