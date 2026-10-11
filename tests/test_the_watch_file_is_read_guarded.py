"""Every reader of a job_monitor watch file opens it bounded in time and size.

Input: a watch path (``<ws>/.delfin/agent_watched_jobs.json``, the shared
``watched_jobs.json``, the workspace index) that holds a FIFO, a symlink,
a file larger than ``job_monitor.WATCH_FILE_MAX_BYTES``, or JSON of the
wrong shape. Output: the empty view (``{"jobs": {}}``) within the join
timeout, nothing read through the symlink, nothing raised; a writer that
meets such an entry replaces it with a regular file and leaves a
symlink's target untouched.

The workspace is writable by the model, and ``check_agent_jobs`` runs at
the start of every turn: a FIFO there blocked the turn start forever on
6d6295af. Every read runs in a daemon thread joined with a timeout, and a
blocked FIFO reader is released by opening the FIFO for writing, so the
suite cannot hang.
"""

from __future__ import annotations

import json
import os
import threading

import pytest

from delfin.agent import api_client
from delfin.agent import job_monitor as jm

pytestmark = pytest.mark.skipif(not hasattr(os, "mkfifo"),
                                reason="needs os.mkfifo")

_JOBS = {"77": {"kind": "slurm", "description": "x", "folder": "",
                "added_at": 0.0, "last_state": ""}}


def _bounded(fn, fifo=None, timeout=5.0):
    """Run ``fn`` in a daemon thread; fail if it is still running after
    ``timeout`` seconds (releasing a reader blocked on ``fifo`` first)."""
    out: dict = {}

    def _run():
        try:
            out["value"] = fn()
        except BaseException as exc:          # noqa: BLE001 - reported below
            out["error"] = exc

    t = threading.Thread(target=_run, daemon=True)
    t.start()
    t.join(timeout)
    if t.is_alive():
        if fifo is not None:
            try:
                fd = os.open(str(fifo), os.O_WRONLY | os.O_NONBLOCK)
                os.close(fd)
            except OSError:
                pass
            t.join(timeout)
        pytest.fail(f"{getattr(fn, '__name__', fn)} blocked on the watch file")
    if "error" in out:
        raise out["error"]
    return out.get("value")


def _ws_with(tmp_path, shape):
    """A workspace whose agent watch path holds ``shape``; returns
    (workspace, watch path, symlink target or None)."""
    ws = tmp_path / "ws"
    (ws / ".delfin").mkdir(parents=True)
    p = jm._agent_watch_path(ws)
    target = None
    if shape == "fifo":
        os.mkfifo(p)
    elif shape == "symlink":
        target = tmp_path / "elsewhere.json"
        target.write_text(json.dumps({"jobs": dict(_JOBS)}))
        p.symlink_to(target)
    elif shape == "huge":
        jobs = dict(_JOBS)
        jobs.update({str(10_000_000 + i): {"kind": "slurm"}
                     for i in range(100_000)})
        p.write_text(json.dumps({"jobs": jobs}))
        assert p.stat().st_size > (1 << 20)
    elif shape == "list":
        p.write_text(json.dumps([1, 2, 3]))
    elif shape == "jobs_list":
        p.write_text(json.dumps({"jobs": [1, 2, 3]}))
    return ws, p, target


_SHAPES = ["fifo", "symlink", "huge", "list", "jobs_list"]


@pytest.mark.parametrize("shape", _SHAPES)
def test_load_watched_gives_the_empty_view(tmp_path, shape):
    ws, p, _ = _ws_with(tmp_path, shape)
    data = _bounded(lambda: jm.load_watched(p),
                    fifo=p if shape == "fifo" else None)
    assert data == {"jobs": {}}


@pytest.mark.parametrize("shape", _SHAPES)
def test_check_agent_jobs_at_a_turn_start_returns_nothing(tmp_path, shape):
    ws, p, target = _ws_with(tmp_path, shape)
    before = target.read_text() if target else None
    done = _bounded(lambda: jm.check_agent_jobs(
                        ws, run_fn=lambda cmd: "", session_id="S"),
                    fifo=p if shape == "fifo" else None)
    assert done == []
    if target is not None:
        assert target.read_text() == before


@pytest.mark.parametrize("shape", ["fifo", "symlink", "huge"])
def test_a_writer_replaces_the_entry_and_not_its_target(tmp_path, shape):
    ws, p, target = _ws_with(tmp_path, shape)
    before = target.read_text() if target else None
    _bounded(lambda: jm.register_agent_job(ws, "4242", "probe"),
             fifo=p if shape == "fifo" else None)
    assert p.is_file() and not p.is_symlink()
    assert set(json.loads(p.read_text())["jobs"]) == {"4242"}
    if target is not None:
        assert target.read_text() == before


@pytest.mark.parametrize("shape", ["fifo", "symlink", "huge"])
def test_the_workspace_index_is_read_guarded(tmp_path, monkeypatch, shape):
    idx = tmp_path / "agent_watch_index.json"
    monkeypatch.setattr(jm, "_AGENT_WATCH_INDEX_PATH", idx)
    ws = tmp_path / "ws_real"
    (ws / ".delfin").mkdir(parents=True)
    jm._agent_watch_path(ws).write_text(json.dumps({"jobs": {}}))
    target = None
    if shape == "fifo":
        os.mkfifo(idx)
    elif shape == "symlink":
        target = tmp_path / "idx_elsewhere.json"
        target.write_text(json.dumps({"workspaces": {str(ws): 9e12}}))
        idx.symlink_to(target)
    else:
        idx.write_text(json.dumps({"workspaces": {
            f"/nonexistent/{i:07d}": 9e12 for i in range(40_000)}}))
        assert idx.stat().st_size > (1 << 20)
    fifo = idx if shape == "fifo" else None
    assert _bounded(jm.agent_watch_workspaces, fifo=fifo) == []
    before = target.read_text() if target else None
    _bounded(lambda: jm._note_agent_watch_workspace(ws), fifo=fifo)
    assert idx.is_file() and not idx.is_symlink()
    if target is not None:
        assert target.read_text() == before


def test_a_regular_watch_file_still_reads(tmp_path):
    ws, p, _ = _ws_with(tmp_path, "none")
    p.write_text(json.dumps({"jobs": dict(_JOBS)}))
    assert jm.load_watched(p)["jobs"] == _JOBS
    assert api_client._read_watch_jobs(p) == _JOBS


@pytest.mark.parametrize("shape", _SHAPES)
def test_the_pending_note_reader_and_load_watched_agree(tmp_path, shape):
    """Two readers of one file must not drift: the pending-jobs note's
    reader is load_watched's ``jobs``."""
    ws, p, _ = _ws_with(tmp_path, shape)
    fifo = p if shape == "fifo" else None
    a = _bounded(lambda: api_client._read_watch_jobs(p), fifo=fifo)
    b = _bounded(lambda: jm.load_watched(p)["jobs"], fifo=fifo)
    assert a == b == {}


def test_one_size_bound_for_both_readers():
    assert not hasattr(api_client, "_PENDING_WATCH_MAX_BYTES")
    assert jm.WATCH_FILE_MAX_BYTES == 1 << 20
