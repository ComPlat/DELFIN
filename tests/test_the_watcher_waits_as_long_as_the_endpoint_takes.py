"""The second opinion waits as long as the endpoint currently takes.

Measured 2026-10-19 on the KIT proxy: median first-token time 18.5 s, 45 %
of 137 requests past 20 s. The watcher, on the same endpoint with a fixed
20 s, failed again and again; under Bypass each failure became a question to
the user. Without a setting the wait now follows the endpoint (1.5 x the
75th percentile of the last first-token times), never below the default and
never above _CLASSIFIER_ADAPTIVE_MAX_S. A setting still decides alone.
"""

from __future__ import annotations

import pytest

from delfin.agent import api_client as A


@pytest.fixture(autouse=True)
def _fresh(monkeypatch):
    monkeypatch.setattr(A, "_ENDPOINT_TTFT_S",
                        A.collections.deque(maxlen=20))
    monkeypatch.setattr("delfin.user_settings.load_settings", lambda: {})


def test_without_samples_the_default_holds():
    assert A._classifier_timeout_s() == A._CLASSIFIER_TIMEOUT_S


def test_a_slow_endpoint_lengthens_the_wait():
    for t in (18.5, 40.7, 12.0, 27.0, 45.0, 22.0, 30.0):
        A._note_endpoint_latency(t)
    wait = A._classifier_timeout_s()
    assert 40.0 < wait <= A._CLASSIFIER_ADAPTIVE_MAX_S


def test_a_fast_endpoint_never_shortens_it():
    for _ in range(10):
        A._note_endpoint_latency(0.5)
    assert A._classifier_timeout_s() == A._CLASSIFIER_TIMEOUT_S


def test_a_stalled_endpoint_is_capped():
    for _ in range(10):
        A._note_endpoint_latency(300.0)
    assert A._classifier_timeout_s() == A._CLASSIFIER_ADAPTIVE_MAX_S


def test_the_setting_decides_alone(monkeypatch):
    for _ in range(10):
        A._note_endpoint_latency(60.0)
    monkeypatch.setattr("delfin.user_settings.load_settings",
                        lambda: {"agent": {"classifier": {"timeout_s": 15}}})
    assert A._classifier_timeout_s() == 15.0
