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


# ---------------------------------------------------------------------------
# The same trap, the other way round: one tool spelled SUCCESS differently
# ---------------------------------------------------------------------------
#
# `status: "ok"` is not a convention of the whole surface, and should not
# become one -- 41 of 57 executors omit it, and for a read tool the key
# would be noise beside the data the caller actually wants.
#
# It IS unanimous among the tools whose entire result is "a file now
# exists": create_docx, create_pdf, merge_pdfs, split_pdf,
# fill_docx_template and fill_pdf_form all open with it. draft_email did
# not, and its payload -- recipients, byte count, a note -- gives a
# caller nothing that reads as a verdict. Asking that family's question
# of it returned a written draft as a failure. Found by driving it: the
# probe that caught it was written with the whole surface in view and
# still made the mistake.
#
# fill_series is deliberately outside this: it answers in prose, not
# JSON, because its result is a per-row report. A status key on a string
# would be a fiction.

_FILE_MAKERS = ("create_docx", "create_pdf", "merge_pdfs", "split_pdf",
                "fill_docx_template", "fill_pdf_form", "draft_email")


def test_a_written_draft_says_it_succeeded(tmp_path):
    perms = A.KitToolPermissions(workspace=str(tmp_path))
    perms.mode = "acceptEdits"
    perms.task_session_id = "draft-ok"
    out = _call("draft_email", {"path": "m.eml", "to": "max@example.org",
                                "subject": "Ergebnis", "body": "Anbei."},
                perms)
    assert out.get("status") == "ok", out
    assert "error" not in out
    assert (tmp_path / "m.eml").is_file()


def test_the_payload_a_caller_needs_is_still_there(tmp_path):
    """A status key must not become the whole answer: the recipients and
    the note that this is NOT sent are the point of the tool."""
    perms = A.KitToolPermissions(workspace=str(tmp_path))
    perms.mode = "acceptEdits"
    perms.task_session_id = "draft-payload"
    out = _call("draft_email", {"path": "m.eml", "to": "max@example.org",
                                "subject": "Ergebnis", "body": "Anbei."},
                perms)
    assert out["to"] == ["max@example.org"]
    assert "NOT sent" in out.get("note", "")
    assert out.get("bytes", 0) > 0


def test_a_refused_draft_still_answers_with_error_only(tmp_path):
    perms = A.KitToolPermissions(workspace=str(tmp_path))
    perms.mode = "acceptEdits"
    perms.task_session_id = "draft-bad"
    out = _call("draft_email", {"path": "m.eml", "to": "not-an-address",
                                "subject": "s", "body": "b"}, perms)
    assert "error" in out
    assert out.get("status") != "ok"


def test_the_family_rule_so_the_next_one_added_does_not_drift():
    """Every tool whose result is a written file opens with status:ok.

    Read off the source rather than by calling them: several need real
    inputs, and the point is the shape they are written to return.
    """
    import re
    from pathlib import Path

    src = Path(A.__file__).read_text(encoding="utf-8")
    for name in _FILE_MAKERS:
        m = re.search(rf"def _execute_{name}\b", src)
        assert m, f"{name} has no executor any more"
        nxt = re.search(r"\n    def _execute_", src[m.start() + 10:])
        body = (src[m.start(): m.start() + 10 + nxt.start()] if nxt
                else src[m.start():])
        assert re.search(r'"status":\s*"ok"|setdefault\("status", "ok"\)',
                         body), f"{name} no longer says status:ok on success"
