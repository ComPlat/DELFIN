"""Each request records a fingerprint of the prompt it sent, so a cold request
can be told apart: the prompt changed, or the endpoint lost its cache.

2026-09-16: three GLM sessions building `delfin doctor` served 7 of 14, 6 of
11 and 20 of 41 requests cold (60 s or more to the first token), several of
them right after a 99 %-cached request. The report could not say why.
"""
import inspect

from delfin.agent import api_client as A
from delfin.agent import turn_metrics as TM


def test_the_fingerprints_follow_the_prompt():
    msgs = [{"role": "system", "content": "S"}, {"role": "user", "content": "q"},
            {"role": "assistant", "content": "a"}]
    one = A._prompt_fingerprints(msgs, 3)
    assert one == A._prompt_fingerprints(list(msgs), 3)
    changed = A._prompt_fingerprints([{"role": "system", "content": "S2"}] + msgs[1:], 3)
    assert changed["system_hash"] != one["system_hash"]
    assert changed["prefix_hash"] != one["prefix_hash"]
    assert A._prompt_fingerprints([], 0) == {}


def test_the_report_says_why_a_request_was_cold(tmp_path, monkeypatch):
    monkeypatch.setattr(TM, "_DIR", tmp_path)
    TM.record_request("s", model="m", ttft_ms=2000, total_ms=3000, input_tokens=100,
                      cached_tokens=99, system_hash="a", prefix_hash="p1")
    TM.record_request("s", model="m", ttft_ms=200000, total_ms=210000, input_tokens=100,
                      cached_tokens=0, system_hash="a", prefix_hash="p1")
    TM.record_request("s", model="m", ttft_ms=200000, total_ms=210000, input_tokens=100,
                      cached_tokens=0, system_hash="b", prefix_hash="p2")
    text = TM.format_requests(TM.read_requests("s"))
    assert "(same prefix: endpoint cache lost)" in text
    assert "(system prompt changed)" in text


def test_the_client_hands_the_fingerprints_over():
    src = inspect.getsource(A.OpenAIClient.stream_message)
    i = src.index("_tm_req.record_request(")
    assert "_prompt_fingerprints(" in src[i:i + 900]
