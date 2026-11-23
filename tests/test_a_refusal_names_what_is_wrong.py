"""Invalid calls through the product's own dispatch: the refusal names
what is wrong.

A tool that answers a bad call with a bare failure (or worse, a success
shape) leaves the model guessing; each refusal below names the job id,
the missing argument, or the gate that refused. These run through
``_doc_executor.execute`` -- the path the model's call takes -- with a
temporary workspace and no real job running.
"""

from __future__ import annotations

import json

import pytest

from delfin.agent.api_client import KitToolPermissions, _doc_executor


@pytest.fixture
def perms(tmp_path):
    return KitToolPermissions(workspace=str(tmp_path / "ws"),
                              mode="bypassPermissions")


def _exec(name, args, perms):
    return json.loads(_doc_executor.execute(name, args, perms))


def test_status_output_kill_name_the_unknown_job_id(perms):
    for name in ("bash_status", "bash_output", "bash_kill"):
        out = _exec(name, {"job_id": "ffffffff"}, perms)
        assert "error" in out, f"{name} must refuse, got {out}"
        assert "unknown" in out["error"].lower(), out["error"]


def test_missing_job_id_is_named_as_the_what_is_wrong(perms):
    for name in ("bash_status", "bash_output", "bash_kill", "watch_job"):
        out = _exec(name, {}, perms)
        assert "error" in out
        assert "job_id is required" in out["error"], out["error"]


def test_schedule_wakeup_without_a_delay_refuses_with_the_reason(perms):
    out = _exec("schedule_wakeup", {}, perms)
    assert "error" in out, out
    assert "delay_seconds" in out["error"], out["error"]


def test_cron_delete_without_an_id_refuses_with_the_reason(perms):
    out = _exec("cron_delete", {}, perms)
    assert "error" in out, out
    assert "id" in out["error"].lower(), out["error"]


def test_watch_job_unknown_id_does_not_promise_a_notice(perms):
    # A bash job id the registry does not know can never complete;
    # watching it would promise a notification that cannot arrive.
    out = _exec("watch_job", {"job_id": "ffffffff"}, perms)
    assert "error" in out, out
    assert "unknown" in out["error"].lower() or "not" in out["error"].lower()
