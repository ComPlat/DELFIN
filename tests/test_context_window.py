"""The effective context window: ONE source for engine and client.

Control (red before delfin/agent/context_window.py existed): the cap
logic lived in engine.py only; the client derived its in-turn elision
budget from the raw model window and never saw the cap. Both must read
the same function.
"""

from __future__ import annotations

from delfin.agent.context_window import capped_window


def test_env_cap_lowers_the_window(monkeypatch):
    monkeypatch.setenv("DELFIN_CONTEXT_WINDOW_CAP", "40000")
    assert capped_window(131_000) == 40_000


def test_setting_cap_lowers_the_window(monkeypatch):
    monkeypatch.delenv("DELFIN_CONTEXT_WINDOW_CAP", raising=False)
    import delfin.user_settings as us
    monkeypatch.setattr(
        us, "load_settings",
        lambda: {"agent": {"context_window_cap": 50000}})
    assert capped_window(131_000) == 50_000


def test_no_cap_returns_the_raw_window(monkeypatch):
    monkeypatch.delenv("DELFIN_CONTEXT_WINDOW_CAP", raising=False)
    import delfin.user_settings as us
    monkeypatch.setattr(us, "load_settings", lambda: {})
    assert capped_window(131_000) == 131_000


def test_cap_never_raises_the_window(monkeypatch):
    monkeypatch.setenv("DELFIN_CONTEXT_WINDOW_CAP", "990000")
    assert capped_window(131_000) == 131_000


def test_invalid_env_cap_is_ignored(monkeypatch):
    monkeypatch.setenv("DELFIN_CONTEXT_WINDOW_CAP", "not-a-number")
    import delfin.user_settings as us
    monkeypatch.setattr(us, "load_settings", lambda: {})
    assert capped_window(131_000) == 131_000


def test_engine_forwarder_reads_the_same_source(monkeypatch):
    """engine._capped_context_window must stay a thin forwarder so
    existing callers (and tests) keep working with the ONE source."""
    monkeypatch.setenv("DELFIN_CONTEXT_WINDOW_CAP", "40000")
    from delfin.agent.engine import _capped_context_window
    assert _capped_context_window(131_000) == 40_000
