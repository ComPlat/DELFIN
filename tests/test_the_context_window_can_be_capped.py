"""The context window can be capped, so compaction can be made to happen.

Four long supervised sessions on a 131k-token model never reached a
single compaction (2026-09-26); whether a session survives several could
not be observed, and a user on a small local model had no way to make
DELFIN compact earlier.
"""

from __future__ import annotations

from delfin.agent import engine


def test_the_cap_lowers_the_window(monkeypatch):
    monkeypatch.setenv("DELFIN_CONTEXT_WINDOW_CAP", "40000")
    assert engine._capped_context_window(131_072) == 40_000


def test_the_cap_never_raises_it(monkeypatch):
    monkeypatch.setenv("DELFIN_CONTEXT_WINDOW_CAP", "500000")
    assert engine._capped_context_window(131_072) == 131_072


def test_no_cap_and_nonsense_leave_the_window(monkeypatch):
    monkeypatch.delenv("DELFIN_CONTEXT_WINDOW_CAP", raising=False)
    monkeypatch.setattr("delfin.user_settings.load_settings", lambda: {})
    assert engine._capped_context_window(131_072) == 131_072
    monkeypatch.setenv("DELFIN_CONTEXT_WINDOW_CAP", "lots")
    assert engine._capped_context_window(131_072) == 131_072


def test_the_setting_caps_it(monkeypatch):
    monkeypatch.delenv("DELFIN_CONTEXT_WINDOW_CAP", raising=False)
    monkeypatch.setattr("delfin.user_settings.load_settings",
                        lambda: {"agent": {"context_window_cap": 32000}})
    assert engine._capped_context_window(131_072) == 32_000
