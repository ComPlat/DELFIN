"""Shared test fixtures.

Redirect the subagent state files (live registry, telemetry, finished
sessions) away from the real ``~/.delfin`` so test runs never leave
artifacts that show up in the user's dashboard (live panel, ``/agents``
listing, stats).
"""

from __future__ import annotations

import pytest


@pytest.fixture(autouse=True)
def _isolate_subagent_state(tmp_path, monkeypatch):
    from delfin.agent import subagents as sa
    monkeypatch.setattr(sa, "_RUNNING_DIR",
                        tmp_path / "subagent_running")
    monkeypatch.setattr(sa, "_TELEMETRY_PATH",
                        tmp_path / "subagent_telemetry.jsonl")
    monkeypatch.setattr(sa, "_SESSIONS_DIR",
                        tmp_path / "subagent_sessions")
    yield


# Secret API keys delfin resolves, each from the env OR ~/.delfin/credentials.json.
_SECRET_KEY_VARS = ("KIT_TOOLBOX_API_KEY", "OPENAI_API_KEY", "ANTHROPIC_API_KEY")


@pytest.fixture(autouse=True)
def _ci_parity_no_secrets(tmp_path, monkeypatch):
    """Make the local test environment match CI's: NO API key reachable from
    either the environment or the ``~/.delfin/credentials.json`` store.

    CI runs key-free and fully mocked (see ``.github/workflows/ci.yml`` —
    "no secrets are required"). A test that quietly relies on an ambient key
    therefore passes on a developer box that has one, but fails in CI — the
    exact "lokal grün, CI rot" trap. Stripping the keys here surfaces that
    whole class of failure on the developer's own machine, before the push.

    A test that genuinely needs a key still sets it itself (this autouse
    fixture runs first; the test's own ``monkeypatch.setenv`` then wins).
    """
    for var in _SECRET_KEY_VARS:
        monkeypatch.delenv(var, raising=False)
    try:
        from delfin.agent import credentials as _cred
        monkeypatch.setattr(_cred, "_DEFAULT_PATH",
                            tmp_path / "no_credentials.json")
    except Exception:
        pass
    yield
