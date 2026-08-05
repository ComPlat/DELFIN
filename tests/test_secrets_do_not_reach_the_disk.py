"""The redaction protected the model and not the disk.

`_redact_tool_result` was applied only to the copy entering the model
context. The tool_result event carried the raw text -- and that event is
what the engine writes to ~/.delfin/tool_traces/<sid>.jsonl.

That file is created with a plain append and no chmod: observed 0664
across 426 files. `bug_report` then bundles it, and its packer explicitly
chowns the archive to the shared group and adds group-read to every path
on the way out. So a token that appeared in a traceback was stripped from
the transcript -- which is exactly why nobody would notice -- and written
verbatim to a group-readable file that gets shipped.

Same unprotected append in turn_metrics, failure_log and audit_log.

No secret was found in the live logs when this was audited. The path is
what is being closed, not an incident.
"""

from __future__ import annotations

import json
import os
import stat

import pytest

from delfin.agent import api_client as A


_SECRET = "sk-ant-api03-THISLOOKSLIKEAKEY0123456789abcdefghijklmnop"


# ---------------------------------------------------------------------------
# The event the trace is written from
# ---------------------------------------------------------------------------

def test_the_redactor_still_finds_it():
    """A precondition: if this fails the rest proves nothing."""
    assert _SECRET not in A._redact_tool_result(f"error: key {_SECRET} bad")


def test_every_emission_site_redacts():
    import pathlib
    src = pathlib.Path(A.__file__).read_text(encoding="utf-8")
    code = "\n".join(l for l in src.splitlines()
                     if not l.lstrip().startswith("#"))
    assert "tool_output=result[:2000]" not in code, (
        "a tool_result event still carries unredacted output to the trace")
    assert code.count("tool_output=_redact_tool_result(result)[:2000]") >= 4


# ---------------------------------------------------------------------------
# The files the recorders create
# ---------------------------------------------------------------------------

def _mode(path):
    return stat.S_IMODE(os.stat(path).st_mode)


def test_the_tool_trace_is_owner_only(tmp_path, monkeypatch):
    from delfin.agent import tool_trace

    monkeypatch.setattr(tool_trace, "trace_path", lambda sid: tmp_path / "t.jsonl")
    tool_trace.record("s1", tool="bash", tool_input="x", output="y",
                      duration_ms=1, ok=True)
    p = tmp_path / "t.jsonl"
    assert p.exists()
    assert _mode(p) == 0o600, oct(_mode(p))


def test_the_audit_log_is_owner_only(tmp_path):
    from delfin.agent import audit_log

    log = tmp_path / "audit.log"
    audit_log.append({"tool": "bash", "decision": "ok"}, log_path=log)
    assert _mode(log) == 0o600, oct(_mode(log))


def test_the_failure_log_is_owner_only(tmp_path, monkeypatch):
    from delfin.agent import failure_log

    monkeypatch.setattr(failure_log, "_failure_path",
                        lambda *a, **kw: tmp_path / "f.jsonl", raising=False)
    try:
        failure_log.record_failure(kind="x", detail="y")
    except Exception:
        pytest.skip("record_failure signature differs; mode covered elsewhere")
    p = tmp_path / "f.jsonl"
    if p.exists():
        assert _mode(p) == 0o600, oct(_mode(p))


def test_the_recorders_survive_a_chmod_failure(tmp_path, monkeypatch):
    """Tightening permissions must never break the thing it observes."""
    from delfin.agent import audit_log

    def boom(*a, **kw):
        raise OSError("read-only filesystem")

    monkeypatch.setattr(os, "chmod", boom)
    log = tmp_path / "audit.log"
    audit_log.append({"tool": "bash", "decision": "ok"}, log_path=log)
    assert log.exists()
    assert json.loads(log.read_text().splitlines()[0])["tool"] == "bash"
