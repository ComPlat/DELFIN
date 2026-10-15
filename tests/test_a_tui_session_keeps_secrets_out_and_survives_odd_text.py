"""What a TUI session shows and keeps is safe to show and keep.

A credential in a tool call is masked before it reaches a window or the
event log, as the terminal masks it; a file name that is not UTF-8 does not
take the session down; and a session that ends is saved for a resume.
"""

from __future__ import annotations

import json

from tui_worker_client import WorkerClient

TOKEN = "f4k3t0k3n" * 4


def _log_lines(w: WorkerClient) -> list[str]:
    path = w.log_path()
    assert path.exists(), "no event log written"
    return path.read_text().splitlines()


def test_a_credential_in_a_tool_call_is_masked_live_and_in_the_log(tmp_path):
    with WorkerClient(tmp_path) as w:
        start, end = w.turn("secret please")
        started = next(e for e in w.events[start:end] if e["type"] == "tool.started")
        ended = next(e for e in w.events[start:end] if e["type"] == "tool.ended")
        assert TOKEN not in json.dumps(started) and TOKEN not in json.dumps(ended)
        assert "redacted" in started["input"]["command"]
        assert "redacted" in ended["output"]
        # The identifiers are left alone.
        assert started["callId"] == ended["callId"]
        assert started["session"] == w.init_result["session"]["sessionId"]
    assert not any(TOKEN in line for line in _log_lines(w))
    assert not any(TOKEN.encode() in line for line in w.raw_lines)


def test_a_file_name_that_is_not_utf8_does_not_end_the_session(tmp_path):
    with WorkerClient(tmp_path) as w:
        start, end = w.turn("latin please")
        ended = next(e for e in w.events[start:end] if e["type"] == "tool.ended")
        assert ended["output"].startswith("caf")
        assert w.events[end].get("error") is None
        _, end2 = w.turn("hello")              # still alive, still talking
        assert w.events[end2].get("error") is None
    assert w.proc.returncode == 0


def test_a_session_that_ends_is_saved_for_a_resume(tmp_path):
    with WorkerClient(tmp_path, "--fake-save") as w:
        w.turn("hello")
        session = w.init_result["session"]["sessionId"]
    saved = [p for root in (w.home, tmp_path / "child_scratch")
             for p in root.rglob("*.json") if session in p.read_text(errors="replace")]
    assert saved, "the ended session was not saved"
