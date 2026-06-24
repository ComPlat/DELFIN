"""Transient-API-error retry on the streaming call.

A flaky shared-proxy hiccup (timeout / 5xx / rate-limit) on a long KIT run
should be retried with backoff instead of killing the whole turn — and the
retry must NOT duplicate already-streamed output, nor fire on a deterministic
error (400/401) that would just fail again. No live model: create() is faked.
"""

from __future__ import annotations

import pytest

from delfin.agent import mcp_client as M
from delfin.agent import model_capabilities as mc


# --- minimal OpenAI-stream fakes -------------------------------------------
class _Delta:
    def __init__(self, content=None, tool_calls=None):
        self.content = content
        self.tool_calls = tool_calls


class _Choice:
    def __init__(self, delta, finish=None):
        self.delta = delta
        self.finish_reason = finish


class _Usage:
    prompt_tokens = 5
    completion_tokens = 3


class _Chunk:
    def __init__(self, choices, usage=None):
        self.choices = choices
        self.usage = usage


class _Stream:
    def __init__(self, chunks):
        self._chunks = chunks

    def __iter__(self):
        return iter(self._chunks)

    def close(self):
        pass


def _final_text(txt="done"):
    return _Stream([_Chunk([_Choice(_Delta(content=txt), finish="stop")],
                           usage=_Usage())])


class _FakeRegistry:
    def discover_all(self):
        return []

    def discover_resources(self):
        return []

    def discover_prompts(self):
        return []


class APITimeoutError(Exception):
    """Class name matches the transient classifier's 'timeout' hint."""


class BadRequestError(Exception):
    def __init__(self, msg):
        super().__init__(msg)
        self.status_code = 400


@pytest.fixture
def client(monkeypatch, tmp_path):
    caps = mc.ModelCapabilities(model="m", provider="kit",
                                context_window=200_000, supports_tools=True)
    monkeypatch.setattr(mc, "resolve", lambda *a, **k: caps)
    monkeypatch.setattr(M, "get_registry", lambda *a, **k: _FakeRegistry())
    import delfin.agent.api_client as A
    monkeypatch.setattr(A.time, "sleep", lambda *a, **k: None)  # skip real backoff
    return A.create_client(backend="api", provider="kit", api_key="",
                           model="kit.qwen3.5-397b-A17b", cwd=str(tmp_path))


def _drive(client):
    return [ev.text for ev in client.stream_message(
        "sys", [{"role": "user", "content": "go"}], max_tokens=100)
        if ev.type == "text_delta" and ev.text]


def test_transient_error_is_retried_then_succeeds(client):
    calls = {"n": 0}

    def _create(**kwargs):
        calls["n"] += 1
        if calls["n"] == 1:
            raise APITimeoutError("upstream timed out")
        return _final_text("done")

    client.client.chat.completions.create = _create
    joined = "".join(_drive(client))
    assert calls["n"] == 2                  # retried exactly once
    assert joined.count("done") == 1        # answer present once — no duplication
    assert "retrying" in joined.lower()     # the retry was surfaced to the user


def test_nontransient_error_is_not_retried(client):
    calls = {"n": 0}

    def _create(**kwargs):
        calls["n"] += 1
        raise BadRequestError("bad request")

    client.client.chat.completions.create = _create
    with pytest.raises(BadRequestError):
        _drive(client)
    assert calls["n"] == 1                   # a deterministic 400 is not retried
