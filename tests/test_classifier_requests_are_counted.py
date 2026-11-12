"""The action classifier's verdict calls are requests against the endpoint's quota.

KIT counts every request of a key (2,000 a day, 2026-11-11). The turn's own
model calls were written to turn_metrics' request log; the classifier's were
not, so on 2026-11-11 about 860 of 2,000 counted requests had no trace here.
Each verdict call is now one line marked ``purpose: classifier``; the request
summary counts them on their own and keeps its per-request cache analysis to
the turn's own calls.
"""
from __future__ import annotations

from delfin.agent import turn_metrics as TM


def _isolate(tmp_path, monkeypatch):
    monkeypatch.setattr(TM, "_DIR", tmp_path)


def test_a_classifier_request_is_written_with_its_purpose(tmp_path, monkeypatch):
    _isolate(tmp_path, monkeypatch)
    TM.record_request("s1", model="kit.deepseek-v4-flash", purpose="classifier",
                      total_ms=900, finish_reason="verdict")
    TM.record_request("s1", model="kit.glm-5.3-flash", ttft_ms=1200, total_ms=4000,
                      input_tokens=9000, cached_tokens=8000, tool_calls=1,
                      finish_reason="tool_calls", prefix_hash="p1", system_hash="s1")
    rows = TM.read_requests("s1")
    assert [r.get("purpose") for r in rows] == ["classifier", None]
    assert rows[0]["model"] == "kit.deepseek-v4-flash"


def test_the_summary_counts_classifier_requests_apart(tmp_path, monkeypatch):
    _isolate(tmp_path, monkeypatch)
    TM.record_request("s2", model="kit.glm-5.3-flash", ttft_ms=1200, total_ms=4000,
                      input_tokens=9000, cached_tokens=8000, tool_calls=1)
    for _ in range(3):
        TM.record_request("s2", model="kit.deepseek-v4-flash", purpose="classifier",
                          total_ms=800, finish_reason="verdict")
    text = TM.format_requests(TM.read_requests("s2"))
    head = text.splitlines()[0]
    assert head.startswith("1 request(s)")
    assert "+3 classifier request(s)" in head
    assert len(text.splitlines()) == 2              # one per-request line, the turn's own


def test_a_plain_request_has_no_purpose_field(tmp_path, monkeypatch):
    _isolate(tmp_path, monkeypatch)
    TM.record_request("s3", model="kit.glm-5.3-flash", total_ms=10)
    assert "purpose" not in TM.read_requests("s3")[0]


def test_the_clients_verdict_call_is_written_to_its_session(tmp_path, monkeypatch):
    """Through the client's own runner (_classifier_complete), success and failure."""
    from types import SimpleNamespace
    from delfin.agent import api_client as A
    _isolate(tmp_path, monkeypatch)

    class _Completions:
        def __init__(self, fail):
            self.fail = fail

        def create(self, **kw):
            if self.fail:
                raise ConnectionError("quota")
            return SimpleNamespace(choices=[SimpleNamespace(
                message=SimpleNamespace(content="ALLOW: routine"))])

    def _client(fail):
        c = A.OpenAIClient.__new__(A.OpenAIClient)
        c.model = "kit.mistral-small-4-119b-a8b"
        c._provider = "kit"
        c.metrics_session = "s4"
        chat = SimpleNamespace(completions=_Completions(fail))
        c.client = SimpleNamespace(with_options=lambda **kw: SimpleNamespace(chat=chat))
        return c

    assert _client(False)._classifier_complete("judge this", "kit.deepseek-v4-flash")
    try:
        _client(True)._classifier_complete("judge this", "kit.glm-5.3-flash")
    except ConnectionError:
        pass
    rows = TM.read_requests("s4")
    assert [(r["purpose"], r["model"], r["finish_reason"]) for r in rows] == [
        ("classifier", "kit.deepseek-v4-flash", "verdict"),
        ("classifier", "kit.glm-5.3-flash", "error")]
