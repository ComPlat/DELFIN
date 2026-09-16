"""A background job's record cannot name somebody else's file as its output.

The registry is ``<workspace>/.delfin/bash_jobs.json``, inside the folder
the agent may write. Driven on 2026-09-16 in a throwaway home: a record
written there with ``stdout_path`` set to ``~/.ssh/id_ed25519`` put the
key into the next completion notice, which goes into the model's context,
and a second record dated eight days back had its "output", ``~/.bashrc``,
deleted by the prune.
"""
import json
import os
import time
from pathlib import Path

import pytest

from delfin.agent import background_view as BV
from delfin.agent import bash_jobs as BJ


@pytest.fixture
def ws(tmp_path):
    d = tmp_path / "ws"
    (d / ".delfin").mkdir(parents=True)
    return d


def _forge(ws, **jobs):
    (ws / ".delfin" / "bash_jobs.json").write_text(json.dumps({"jobs": jobs}))


def _rec(jid, path, **over):
    rec = {"job_id": jid, "command": "x", "pid": 0, "started_at": time.time() - 10,
           "finished_at": time.time(), "exit_code": 0,
           "stdout_path": str(path), "stderr_path": ""}
    rec.update(over)
    return rec


def test_a_forged_record_puts_no_foreign_file_into_the_notice(ws, tmp_path):
    key = tmp_path / "id_ed25519"
    key.write_text("-----BEGIN OPENSSH PRIVATE KEY-----\nPROBE\n")
    _forge(ws, **{"bg-read": _rec("bg-read", key)})
    events = BJ.drain_finished_events(str(ws))
    assert events and "PROBE" not in json.dumps(events)
    assert "PROBE" not in BV.peek(ws, "shells", "bg-read")


def test_pruning_a_forged_record_deletes_no_foreign_file(ws, tmp_path):
    victim = tmp_path / "bashrc"
    victim.write_text("# rc\n")
    _forge(ws, **{"bg-del": _rec("bg-del", victim, started_at=1.0, finished_at=2.0)})
    BJ.drain_finished_events(str(ws))
    assert victim.exists()


def test_a_link_named_like_an_output_is_not_one(ws, tmp_path):
    key = tmp_path / "id_ed25519"
    key.write_text("SECRET\n")
    sym = ws / "kit_bg_a.stdout"
    sym.symlink_to(key)
    hard = ws / "kit_bg_b.stdout"
    os.link(key, hard)
    assert not BJ.is_own_output_file(sym)
    assert not BJ.is_own_output_file(hard)
    _forge(ws, **{"s": _rec("s", sym), "h": _rec("h", hard)})
    assert "SECRET" not in json.dumps(BJ.drain_finished_events(str(ws)))


def test_a_real_jobs_output_is_still_read(ws):
    reg = BJ.get_registry()
    job = reg.start("echo hello-from-the-job", cwd=str(ws), workspace=str(ws))
    for _ in range(100):
        if job.poll() is not None:
            break
        time.sleep(0.05)
    assert BJ.is_own_output_file(job.stdout_path)
    assert "hello-from-the-job" in BJ.read_output(job)["stdout"]
