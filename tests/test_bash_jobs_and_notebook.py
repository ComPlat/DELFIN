"""Tests for the background-bash registry and notebook tools.

Covers:
- bash_background returns a job_id, status flips after completion,
  output is readable while running and after exit.
- bash_kill terminates a long-running job (SIGTERM, then SIGKILL).
- The bash safety gate (deny-list, secret-scanner, sandbox) applies
  to bash_background just like to bash.
- notebook_read returns cell records with idx + type + source +
  output_summary; outputs are summarised, not dumped.
- notebook_edit supports replace / insert_before / insert_after /
  delete, and refuses out-of-range indices and missing read baseline.
"""

from __future__ import annotations

import json
import time
from pathlib import Path

import pytest

from delfin.agent.api_client import KitToolPermissions, _doc_executor


@pytest.fixture
def workspace(tmp_path) -> Path:
    ws = tmp_path / "ws"
    ws.mkdir()
    return ws


# ---------------------------------------------------------------------------
# Background bash
# ---------------------------------------------------------------------------


def test_bash_background_round_trip(workspace):
    perms = KitToolPermissions(workspace=workspace, mode="bypassPermissions")
    out = _doc_executor.execute("bash_background", {
        "command": "for i in 1 2 3; do echo line$i; done; echo DONE",
        "description": "demo",
    }, perms)
    payload = json.loads(out)
    assert payload["status"] == "started"
    job_id = payload["job_id"]
    assert len(job_id) == 8

    # Wait for completion (loop runs in <1s).
    for _ in range(50):
        s = json.loads(_doc_executor.execute(
            "bash_status", {"job_id": job_id}, perms))
        if not s["running"]:
            break
        time.sleep(0.1)
    assert s["running"] is False
    assert s["exit_code"] == 0

    o = json.loads(_doc_executor.execute(
        "bash_output", {"job_id": job_id}, perms))
    assert "line1" in o["stdout"]
    assert "DONE" in o["stdout"]
    assert o["exit_code"] == 0


def test_bash_background_output_while_running(workspace):
    perms = KitToolPermissions(workspace=workspace, mode="bypassPermissions")
    out = _doc_executor.execute("bash_background", {
        "command": "echo first; sleep 0.6; echo second",
        "description": "streaming",
    }, perms)
    job_id = json.loads(out)["job_id"]
    time.sleep(0.2)
    o = json.loads(_doc_executor.execute(
        "bash_output", {"job_id": job_id}, perms))
    assert "first" in o["stdout"]
    assert "second" not in o["stdout"]
    assert o["running"] is True


def test_bash_kill_terminates_long_job(workspace):
    perms = KitToolPermissions(workspace=workspace, mode="bypassPermissions")
    out = _doc_executor.execute("bash_background", {
        "command": "sleep 30; echo NEVER",
        "description": "long",
    }, perms)
    job_id = json.loads(out)["job_id"]
    time.sleep(0.2)
    k = json.loads(_doc_executor.execute(
        "bash_kill", {"job_id": job_id}, perms))
    assert k["status"] == "ok"
    time.sleep(0.3)
    s = json.loads(_doc_executor.execute(
        "bash_status", {"job_id": job_id}, perms))
    assert s["running"] is False
    assert s["exit_code"] != 0  # killed -> negative or non-zero


def test_bash_background_unknown_job_id(workspace):
    perms = KitToolPermissions(workspace=workspace, mode="default")
    s = json.loads(_doc_executor.execute(
        "bash_status", {"job_id": "ffffffff"}, perms))
    assert "error" in s
    assert "unknown" in s["error"].lower()


def test_bash_background_deny_list_applies(workspace):
    perms = KitToolPermissions(workspace=workspace, mode="bypassPermissions")
    out = json.loads(_doc_executor.execute("bash_background", {
        "command": "rm -rf /tmp/something",
        "description": "evil",
    }, perms))
    assert "error" in out
    assert "deny" in out["error"].lower()


def test_bash_background_secret_scanner_applies(workspace):
    perms = KitToolPermissions(workspace=workspace, mode="bypassPermissions")
    out = json.loads(_doc_executor.execute("bash_background", {
        "command": "cat /home/user/.ssh/id_rsa",
        "description": "leak",
    }, perms))
    assert "error" in out
    assert "secret" in out["error"].lower() or "deny" in out["error"].lower()


def test_bash_background_sandbox_cwd_check(workspace, tmp_path):
    foreign = tmp_path / "foreign"
    foreign.mkdir()
    perms = KitToolPermissions(workspace=workspace, mode="bypassPermissions")
    out = json.loads(_doc_executor.execute("bash_background", {
        "command": "ls",
        "description": "list",
        "cwd": str(foreign),
    }, perms))
    assert "error" in out


def test_bash_background_default_mode_needs_auto_allow(workspace):
    """Same gate as foreground bash: default mode without auto_allow rejects."""
    perms = KitToolPermissions(workspace=workspace, mode="default")
    out = json.loads(_doc_executor.execute("bash_background", {
        "command": "rare-custom-cmd --flag",
        "description": "x",
    }, perms))
    assert "error" in out
    assert "auto-allow" in out["error"].lower()


# ---------------------------------------------------------------------------
# Notebook read / edit
# ---------------------------------------------------------------------------


def _write_nb(path: Path, cells: list[dict]) -> None:
    nb = {
        "cells": cells,
        "metadata": {},
        "nbformat": 4,
        "nbformat_minor": 5,
    }
    path.write_text(json.dumps(nb))


def test_notebook_read_basic(workspace):
    nb_path = workspace / "test.ipynb"
    _write_nb(nb_path, [
        {"cell_type": "markdown", "source": "# Title\n", "metadata": {}},
        {"cell_type": "code", "source": "x = 1\n", "outputs": [],
         "execution_count": None, "metadata": {}},
    ])
    perms = KitToolPermissions(workspace=workspace, mode="default")
    out = json.loads(_doc_executor.execute(
        "notebook_read", {"path": str(nb_path)}, perms))
    assert out["cell_count"] == 2
    assert out["cells"][0]["cell_type"] == "markdown"
    assert out["cells"][1]["source"] == "x = 1\n"


def test_notebook_read_summarises_outputs(workspace):
    nb_path = workspace / "with_output.ipynb"
    _write_nb(nb_path, [
        {"cell_type": "code", "source": "print('hi')\n", "outputs": [
            {"output_type": "stream", "name": "stdout", "text": ["hi\n"]},
            {"output_type": "execute_result",
             "data": {"text/plain": ["3"], "image/png": ["BASE64..."]},
             "execution_count": 1, "metadata": {}},
        ], "execution_count": 1, "metadata": {}},
    ])
    perms = KitToolPermissions(workspace=workspace, mode="default")
    out = json.loads(_doc_executor.execute(
        "notebook_read", {"path": str(nb_path)}, perms))
    summary = out["cells"][0]["output_summary"]
    assert "stream" in summary
    assert "execute_result" in summary
    # Image base64 must NOT have leaked into the summary.
    assert "BASE64" not in summary


def test_notebook_read_rejects_non_ipynb(workspace):
    plain = workspace / "plain.txt"
    plain.write_text("hello")
    perms = KitToolPermissions(workspace=workspace, mode="default")
    out = json.loads(_doc_executor.execute(
        "notebook_read", {"path": str(plain)}, perms))
    assert "error" in out
    assert ".ipynb" in out["error"]


def test_notebook_edit_replace(workspace):
    nb_path = workspace / "rep.ipynb"
    _write_nb(nb_path, [
        {"cell_type": "code", "source": "x = 1\n", "outputs": [],
         "execution_count": None, "metadata": {}},
    ])
    perms = KitToolPermissions(workspace=workspace, mode="default")
    _doc_executor.execute("notebook_read", {"path": str(nb_path)}, perms)
    out = json.loads(_doc_executor.execute("notebook_edit", {
        "path": str(nb_path), "cell_idx": 0, "mode": "replace",
        "source": "x = 99\n", "cell_type": "code",
    }, perms))
    assert out["status"] == "ok"
    assert out["cells_before"] == out["cells_after"] == 1
    nb = json.loads(nb_path.read_text())
    src = nb["cells"][0]["source"]
    if isinstance(src, list):
        src = "".join(src)
    assert src == "x = 99\n"


def test_notebook_edit_insert_after(workspace):
    nb_path = workspace / "ins.ipynb"
    _write_nb(nb_path, [
        {"cell_type": "code", "source": "a", "outputs": [],
         "execution_count": None, "metadata": {}},
    ])
    perms = KitToolPermissions(workspace=workspace, mode="default")
    _doc_executor.execute("notebook_read", {"path": str(nb_path)}, perms)
    out = json.loads(_doc_executor.execute("notebook_edit", {
        "path": str(nb_path), "cell_idx": 0, "mode": "insert_after",
        "source": "b", "cell_type": "markdown",
    }, perms))
    assert out["cells_before"] == 1
    assert out["cells_after"] == 2
    nb = json.loads(nb_path.read_text())
    assert nb["cells"][1]["cell_type"] == "markdown"


def test_notebook_edit_delete(workspace):
    nb_path = workspace / "del.ipynb"
    _write_nb(nb_path, [
        {"cell_type": "code", "source": "a", "outputs": [],
         "execution_count": None, "metadata": {}},
        {"cell_type": "code", "source": "b", "outputs": [],
         "execution_count": None, "metadata": {}},
    ])
    perms = KitToolPermissions(workspace=workspace, mode="default")
    _doc_executor.execute("notebook_read", {"path": str(nb_path)}, perms)
    out = json.loads(_doc_executor.execute("notebook_edit", {
        "path": str(nb_path), "cell_idx": 0, "mode": "delete",
    }, perms))
    assert out["cells_before"] == 2
    assert out["cells_after"] == 1
    nb = json.loads(nb_path.read_text())
    src = nb["cells"][0]["source"]
    if isinstance(src, list):
        src = "".join(src)
    assert src == "b"


def test_notebook_edit_requires_read_baseline(workspace):
    nb_path = workspace / "noread.ipynb"
    _write_nb(nb_path, [
        {"cell_type": "code", "source": "a", "outputs": [],
         "execution_count": None, "metadata": {}},
    ])
    perms = KitToolPermissions(workspace=workspace, mode="default")
    out = json.loads(_doc_executor.execute("notebook_edit", {
        "path": str(nb_path), "cell_idx": 0, "mode": "replace",
        "source": "b", "cell_type": "code",
    }, perms))
    assert "error" in out
    assert "before editing" in out["error"]


def test_notebook_edit_out_of_range(workspace):
    nb_path = workspace / "oor.ipynb"
    _write_nb(nb_path, [
        {"cell_type": "code", "source": "a", "outputs": [],
         "execution_count": None, "metadata": {}},
    ])
    perms = KitToolPermissions(workspace=workspace, mode="default")
    _doc_executor.execute("notebook_read", {"path": str(nb_path)}, perms)
    out = json.loads(_doc_executor.execute("notebook_edit", {
        "path": str(nb_path), "cell_idx": 5, "mode": "replace",
        "source": "b", "cell_type": "code",
    }, perms))
    assert "error" in out
    assert "out of range" in out["error"]


def test_notebook_edit_outside_workspace_rejected(workspace, tmp_path):
    foreign = tmp_path / "foreign"
    foreign.mkdir()
    foreign_nb = foreign / "x.ipynb"
    _write_nb(foreign_nb, [{"cell_type": "code", "source": "a",
                             "outputs": [], "execution_count": None,
                             "metadata": {}}])
    perms = KitToolPermissions(workspace=workspace, mode="default")
    out = json.loads(_doc_executor.execute("notebook_edit", {
        "path": str(foreign_nb), "cell_idx": 0, "mode": "replace",
        "source": "b", "cell_type": "code",
    }, perms))
    assert "error" in out
