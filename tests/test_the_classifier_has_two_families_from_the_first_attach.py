"""The action classifier gets its two-family chain in a fresh process.

The classifier picks its models from the endpoint's served roster. The roster
lived only in an in-process cache that nothing had filled when the client
attached the classifier, so every fresh process -- every headless and swarm
session -- ran the classifier on its default model alone with no confirmation
chain (measured 2026-11-11 on real KIT: deepseek alone; with the roster:
mistral first, then deepseek and glm). The attach now fetches the roster once
when the cache is empty, and a failed fetch is retried at most every
_ROSTER_RETRY_S.
"""
from __future__ import annotations

import pytest

from delfin.agent import api_client as A
from delfin.agent import model_capabilities as MC

SERVED = {"kit.mistral-small-4-119b-a8b", "kit.deepseek-v4-flash", "kit.glm-5.3-flash"}


@pytest.fixture
def fresh(monkeypatch, tmp_path):
    monkeypatch.setattr(MC, "_SERVED_CACHE", {})
    monkeypatch.setattr(A, "_ROSTER_TRIED_AT", None)
    monkeypatch.setattr(A, "_classifier_settings",
                        lambda: {"model": "kit.deepseek-v4-flash", "explicit_model": False,
                                 "rules": [], "broad": False})
    calls = []

    def fetch(base_url, api_key):
        calls.append(base_url)
        if fetch.ok:
            MC._remember_served(base_url, [{"id": m} for m in sorted(SERVED)])
    fetch.ok = True
    monkeypatch.setattr(MC, "_fetch_openai_models", fetch)
    client = A.create_client(backend="api", provider="kit", model="kit.glm-5.3-flash",
                             api_key="test-key-not-real", cwd=str(tmp_path))
    perms = A.KitToolPermissions(workspace=tmp_path.resolve(), mode="bypassPermissions")
    return client, perms, calls, fetch


def test_the_first_attach_fetches_the_roster_and_gets_two_families(fresh):
    client, perms, calls, _ = fresh
    client._attach_classifier_runner(perms)
    assert len(calls) == 1
    assert perms.classifier_model == "kit.mistral-small-4-119b-a8b"
    assert perms.classifier_confirm_models == ("kit.deepseek-v4-flash", "kit.glm-5.3-flash")


def test_a_second_attach_uses_the_cache(fresh):
    client, perms, calls, _ = fresh
    client._attach_classifier_runner(perms)
    client._attach_classifier_runner(perms)
    assert len(calls) == 1


def test_a_failed_fetch_is_not_retried_on_every_attach(fresh):
    client, perms, calls, fetch = fresh
    # The client attaches once when it is built; start this case from an
    # empty roster with a failing endpoint.
    MC._SERVED_CACHE.clear()
    A._ROSTER_TRIED_AT = None
    calls.clear()
    fetch.ok = False
    client._attach_classifier_runner(perms)
    client._attach_classifier_runner(perms)
    assert len(calls) == 1
    # without a roster the configured default stays, as before
    assert perms.classifier_model == "kit.deepseek-v4-flash"
