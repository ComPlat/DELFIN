"""bash_kill said {"status": "error"}; everything else says {"error": ...}.

Found by driving the background-job family end to end. Every other tool
in it reports a failure the same way -- bash_status, bash_output,
watch_job, and bash_kill's OWN "job_id is required" check two lines
above -- so a caller that tests for an "error" key is doing exactly what
the rest of the surface teaches it to do. On the one branch that reported
a failed kill, that test came back False and the failure read as a
success.

The job family is otherwise sound, which is why this was the only thing
to fix in it: start, block with wait_seconds, read while running, read
after, unknown ids named clearly, a huge wait_seconds capped rather than
hung, and killing a finished job answering "already finished (rc=0)".

The general rule this file pins is worth more than the one call: a tool
result is a data format with consumers, and one member of a family that
spells the same event differently is a bug in the format, not a detail of
the tool.
"""

from __future__ import annotations

import json
import tempfile

import pytest

import delfin.agent.api_client as A


@pytest.fixture
def perms():
    with tempfile.TemporaryDirectory() as tmp:
        yield A.KitToolPermissions(mode="bypassPermissions", workspace=tmp)


def _call(name, args, perms):
    return json.loads(A._doc_executor.execute(name, args, perms))


@pytest.mark.parametrize("tool", [
    "bash_status", "bash_output", "bash_kill", "watch_job",
])
def test_an_unknown_job_is_an_error_everywhere(tool, perms):
    args = {"job_id": "no-such-job"}
    if tool == "watch_job":
        args["description"] = "x"
    out = _call(tool, args, perms)
    assert "error" in out, f"{tool} -> {out}"
    assert "no-such-job" in out["error"]


def test_a_successful_kill_still_reports_status_ok(perms):
    started = _call("bash_background",
                    {"command": "sleep 30", "description": "a job to stop"},
                    perms)
    job = started["job_id"]
    out = _call("bash_kill", {"job_id": job}, perms)
    assert out.get("status") == "ok", out
    assert "error" not in out


def test_killing_a_finished_job_is_not_an_error(perms):
    """It reports what happened rather than failing: the caller's intent
    -- that the job be stopped -- is satisfied."""
    started = _call("bash_background",
                    {"command": "true", "description": "a job that ends"},
                    perms)
    job = started["job_id"]
    _call("bash_status", {"job_id": job, "wait_seconds": 10}, perms)
    out = _call("bash_kill", {"job_id": job}, perms)
    assert out.get("status") == "ok", out
    assert "already finished" in out.get("message", "")


def test_a_caller_testing_for_error_is_never_misled(perms):
    """The failure this file exists for, stated as the contract it broke."""
    failed = _call("bash_kill", {"job_id": "nope"}, perms)
    ok_start = _call("bash_background",
                     {"command": "sleep 20", "description": "d"}, perms)
    killed = _call("bash_kill", {"job_id": ok_start["job_id"]}, perms)
    assert ("error" in failed) is True
    assert ("error" in killed) is False
