"""Per-turn tool-round budget is settings-driven (agent.max_tool_rounds).

Bug report 20260619-172400 (ka_xn0397): "nach 50 turns stoppt der agent
automatisch" — long multi-file tasks hit the old hard-coded 50-round cap
mid-work, forcing manual 'continue' nudges. The cap is now a setting with a
high default (500); the cost circuit-breaker + consecutive-fail abort remain
the real stops.
"""

from __future__ import annotations

from unittest.mock import patch

from delfin.agent import api_client as A


def _with_settings(settings: dict):
    return patch("delfin.user_settings.load_settings", lambda: settings)


def test_default_is_500_when_unset():
    with _with_settings({"agent": {}}):
        assert A._resolve_max_tool_rounds() == 500


def test_reads_configured_value():
    with _with_settings({"agent": {"max_tool_rounds": 1200}}):
        assert A._resolve_max_tool_rounds() == 1200


def test_zero_means_uncapped():
    # 0 → effectively unlimited rounds (cost cap / fail-abort are the stops).
    with _with_settings({"agent": {"max_tool_rounds": 0}}):
        assert A._resolve_max_tool_rounds() == 100_000


def test_negative_also_uncapped():
    with _with_settings({"agent": {"max_tool_rounds": -1}}):
        assert A._resolve_max_tool_rounds() == 100_000


def test_falls_back_to_500_on_error():
    def _boom():
        raise RuntimeError("corrupt settings")
    with patch("delfin.user_settings.load_settings", _boom):
        assert A._resolve_max_tool_rounds() == 500


def test_default_settings_carries_the_key():
    # The dashboard Settings tab reads/writes agent.max_tool_rounds; the
    # baseline default must exist so the widget shows the real fallback.
    from delfin.user_settings import DEFAULT_SETTINGS
    assert DEFAULT_SETTINGS["agent"]["max_tool_rounds"] == 500
