"""A KIT outage of several minutes must not end an unattended turn.

Measured 2026-11-06/07: bursts of 503/504 lasting minutes ended dozens of
sessions, because transient errors were retried 3 times over 10.5 s in
total. Production now retries up to 7 times, backing off to minutes.
"""
import importlib
import os

import pytest


@pytest.fixture(autouse=True)
def _api_client_restored():
    """Put api_client back exactly as it was after each test.

    importlib.reload re-runs the module in place: every class and function
    in it becomes a new object. A test module that imported one by name
    (``from delfin.agent.api_client import KitToolPermissions``) keeps the
    old object, whose isinstance checks then meet the new classes -- in a
    shared run a later test failed for no reason of its own
    (test_workspace_bias, 'delfin qm_check' no longer auto-allowed). The
    module's namespace is restored, so the old objects are the live ones
    again.
    """
    import delfin.agent.api_client as ac
    saved = dict(ac.__dict__)
    yield
    ac.__dict__.clear()
    ac.__dict__.update(saved)


def _fresh(monkeypatch, **env):
    for k in ("DELFIN_STREAM_RETRY_MAX", "DELFIN_STREAM_RETRY_DELAY_SCALE"):
        monkeypatch.delenv(k, raising=False)
    for k, v in env.items():
        monkeypatch.setenv(k, v)
    import delfin.agent.api_client as ac
    return importlib.reload(ac)


def test_production_rides_out_several_minutes(monkeypatch):
    ac = _fresh(monkeypatch)
    try:
        assert ac._STREAM_RETRY_MAX >= 7
        total = sum(ac._transient_retry_delay(i) for i in range(1, ac._STREAM_RETRY_MAX + 1))
        assert total >= 6 * 60          # several minutes, not ten seconds
        assert ac._transient_retry_delay(1) <= 2.0   # a hiccup is still retried fast
        assert ac._transient_retry_delay(50) <= 240 * 1.2
    finally:
        _fresh(monkeypatch, DELFIN_STREAM_RETRY_MAX="3", DELFIN_STREAM_RETRY_DELAY_SCALE="0")


def test_the_suite_does_not_sleep(monkeypatch):
    ac = _fresh(monkeypatch, DELFIN_STREAM_RETRY_MAX="3", DELFIN_STREAM_RETRY_DELAY_SCALE="0")
    assert ac._STREAM_RETRY_MAX == 3
    assert ac._transient_retry_delay(5) == 0.0
