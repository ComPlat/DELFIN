"""A verdict cut off before it was written is retried, not latched.

Measured on a compute node 2026-11-20, on a session whose provider is not
KIT: ``_classifier_complete`` asks the session's OWN model, and that model
carries a reasoning channel. At ``max_tokens=120`` one call in three ended
``finish_reason="length"`` with ``content=""`` -- the reasoning had spent
the whole budget before a verdict was written. The empty string was read as
"no verdict", which latched the gate into reporting the watcher as down,
and harmless partner messages then fell to the person as approvals. At 600
tokens three of three answered in 1.3-1.7 s with completions of 98-169
tokens.

Input: a client whose completion returns an empty message. Output: whether
a second request is made and what the caller receives. Semantics: the retry
fires for exactly one shape -- no content AND the model ran out of room. An
empty answer for any other reason, a refusal and an error are answers, and
asking twice for them would double the request count against the
endpoint's daily limit for nothing.
"""
from __future__ import annotations

from types import SimpleNamespace

import pytest

from delfin.agent import api_client


def _choice(content: str, finish_reason: str) -> SimpleNamespace:
    return SimpleNamespace(message=SimpleNamespace(content=content),
                           finish_reason=finish_reason)


class _Client:
    """Records every request and replies from a scripted list."""

    def __init__(self, replies):
        self.replies = list(replies)
        self.budgets: list[int] = []

    def with_options(self, **_kw):
        return self

    @property
    def chat(self):
        return SimpleNamespace(completions=SimpleNamespace(create=self._create))

    def _create(self, *, model, messages, max_tokens, stream):
        self.budgets.append(max_tokens)
        reply = self.replies.pop(0) if self.replies else _choice("", "stop")
        return SimpleNamespace(choices=[reply])


def _session(client) -> api_client.OpenAIClient:
    """A client object with only what _classifier_complete reads."""
    obj = object.__new__(api_client.OpenAIClient)
    obj.client = client
    obj.model = "deepseek-v4.1-flash-at-hzdr"
    obj._provider = "blablador"
    obj.metrics_session = ""
    return obj


def test_the_first_request_has_room_for_a_reasoning_channel():
    """120 was the measured cause; the budget is the fix."""
    assert api_client._CLASSIFIER_MAX_TOKENS >= 600, (
        "the first request's budget is back below the measured floor: at 120 "
        "one call in three was truncated before a verdict")
    assert (api_client._CLASSIFIER_RETRY_TOKENS
            > api_client._CLASSIFIER_MAX_TOKENS), (
        "a retry at the same budget asks the same question twice")


def test_a_truncated_reply_is_retried_at_a_larger_budget():
    client = _Client([_choice("", "length"), _choice("allow", "stop")])
    out = _session(client)._classifier_complete("p", "kit.glm-5.3-flash")
    assert out == "allow", (
        "a verdict truncated before it was written reached the gate as an "
        "empty reply, which is what latched it into 'watcher is down'")
    assert client.budgets == [api_client._CLASSIFIER_MAX_TOKENS,
                              api_client._CLASSIFIER_RETRY_TOKENS], (
        f"expected one retry at a larger budget, got {client.budgets}")


def test_an_empty_answer_that_was_not_truncated_is_not_retried():
    """Asking twice for an answer the model finished giving doubles the
    request count against a daily limit for nothing."""
    client = _Client([_choice("", "stop")])
    out = _session(client)._classifier_complete("p", "kit.glm-5.3-flash")
    assert out == ""
    assert len(client.budgets) == 1, (
        f"retried a finished empty answer: {client.budgets}")


def test_a_refusal_is_an_answer_and_is_not_retried():
    client = _Client([_choice("refuse: writes outside the workspace", "stop")])
    out = _session(client)._classifier_complete("p", "kit.glm-5.3-flash")
    assert out.startswith("refuse")
    assert len(client.budgets) == 1


def test_a_truncation_twice_gives_up_rather_than_asking_again():
    """One retry only: the budget is not raised until it works."""
    client = _Client([_choice("", "length"), _choice("", "length")])
    out = _session(client)._classifier_complete("p", "kit.glm-5.3-flash")
    assert out == ""
    assert len(client.budgets) == 2, (
        f"asked more than twice for one verdict: {client.budgets}")


def test_a_session_off_the_kit_proxy_asks_its_own_model():
    """The swap is what exposes the reasoning channel: the configured
    classifier model lives on KIT, so a session elsewhere asks its own."""
    client = _Client([_choice("allow", "stop")])
    sess = _session(client)
    asked: list[str] = []
    orig = client._create

    def spy(*, model, **kw):
        asked.append(model)
        return orig(model=model, **kw)

    client._create = spy
    sess._classifier_complete("p", "kit.glm-5.3-flash")
    assert asked == ["deepseek-v4.1-flash-at-hzdr"], asked


def test_an_error_is_not_turned_into_a_retry():
    """An exception must reach the caller: a transport failure is not a
    truncation, and retrying it here would hide it from the gate."""

    class _Boom(_Client):
        def _create(self, **kw):
            self.budgets.append(kw["max_tokens"])
            raise RuntimeError("connection reset")

    client = _Boom([])
    with pytest.raises(RuntimeError):
        _session(client)._classifier_complete("p", "kit.glm-5.3-flash")
    assert len(client.budgets) == 1
