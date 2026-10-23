"""Every model request is timed on its own.

Report 20260915-132613 carried one entry for a 2.9-hour turn: 64 tool calls,
3.2M input tokens, 11% cached, and no way to tell the cold requests (200+ s
to the first token on kit.glm-5.3) from the warm ones (7-12 s). Reports
filed during a turn carried no timing at all, because the turn entry is
written when the turn ends.
"""

from __future__ import annotations

import pytest

from delfin.agent import bug_report
from delfin.agent import turn_metrics as tm


@pytest.fixture(autouse=True)
def _metrics_dir(tmp_path, monkeypatch):
    monkeypatch.setattr(tm, "_DIR", tmp_path / "turn_metrics")


def test_a_cold_and_a_warm_request_are_told_apart():
    tm.record_request("s1", model="kit.glm-5.3", ttft_ms=212_000,
                      total_ms=230_000, input_tokens=100_000,
                      cached_tokens=0, output_tokens=400, tool_calls=1,
                      finish_reason="tool_calls")
    tm.record_request("s1", model="kit.glm-5.3", ttft_ms=9_000,
                      total_ms=15_000, input_tokens=100_000,
                      cached_tokens=90_000, output_tokens=300,
                      finish_reason="stop")

    entries = tm.read_requests("s1")
    assert [e["ttft_ms"] for e in entries] == [212_000, 9_000]

    summary = tm.format_requests(entries)
    head, cold, warm = summary.splitlines()
    assert "2 request(s)" in head and "cold=1" in head and "cached=45%" in head
    assert "⚠ cold" in cold and "cached=0%" in cold
    assert "⚠ cold" not in warm and "cached=90%" in warm


def test_a_report_carries_the_requests():
    tm.record_request("s1", ttft_ms=212_000, total_ms=230_000,
                      input_tokens=100_000)
    md = bug_report._render_markdown(
        {}, [], request_metrics=tm.read_requests("s1"))
    assert "## Requests" in md and "⚠ cold" in md


def _streaming_client(session):
    """A client whose one streamed request answers and reports usage."""
    import threading
    import types

    from delfin.agent import api_client as ac

    text = types.SimpleNamespace(
        usage=None, choices=[types.SimpleNamespace(
            delta=types.SimpleNamespace(content="ok", tool_calls=None,
                                        reasoning_content=None),
            finish_reason="stop")])
    usage = types.SimpleNamespace(
        usage=types.SimpleNamespace(
            prompt_tokens=1000, completion_tokens=5,
            prompt_tokens_details=types.SimpleNamespace(cached_tokens=900)),
        choices=[])

    class _Stream:
        def __iter__(self):
            return iter((text, usage))

        def close(self):
            pass

    class _Stub:
        class chat:
            class completions:
                @staticmethod
                def create(**kw):
                    return _Stream()

    c = ac.OpenAIClient.__new__(ac.OpenAIClient)
    c.client = _Stub()
    c.model = "kit.glm-5.3"
    c._provider = "kit"
    c._base_url = "https://ki-toolbox.scc.kit.edu/api/v1"
    c._api_key = "x"
    c.effort = ""
    c._permissions = None
    c.on_model_switched = None
    c._steer_lock = threading.Lock()
    c._steer_queue = []
    c._run_notes = []
    c._stop_flag = False
    c.metrics_session = session
    return c


def test_the_client_writes_a_line_as_each_request_ends():
    client = _streaming_client("live")
    try:
        for ev in client.stream_message(
                messages=[{"role": "user", "content": "hi"}],
                system="t", max_tokens=32):
            # Written before the turn goes on: a report filed now has it.
            if getattr(ev, "type", "") == "message_start":
                break
    except Exception:
        pass
    (entry,) = tm.read_requests("live")
    assert (entry["input_tokens"], entry["cached_tokens"],
            entry["output_tokens"]) == (1000, 900, 5)
    assert entry["ttft_ms"] is not None and entry["finish_reason"] == "stop"


def test_nothing_recorded_means_no_section():
    assert tm.read_requests("never") == []
    assert tm.format_requests([]) == ""
    assert "## Requests" not in bug_report._render_markdown({}, [])
