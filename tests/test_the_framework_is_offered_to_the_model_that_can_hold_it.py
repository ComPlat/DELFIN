"""Sixty-one of DELFIN's own tools were never offered to the model.

A user asked "Hallo" on 2026-09-08 and got a warning back:

    MCP tool surface over budget — 61 MCP tools were not advertised: the
    schema budget of 12000 chars is spent.

Their question in the report was the right one: why on a simple hello?

Measured the same day: DELFIN's own servers offer 93 tools costing 54241
characters of schema — about 13560 tokens. The flat 12000-character budget
advertised 32 of them. On the 524288-token window GLM and DeepSeek now
have, the whole surface is 2.6% of the context; the budget was withholding
two thirds of the framework to save that. On a 32k model the same surface
is 42% and the cap is exactly right — so the budget has to follow the
window rather than be one number for every model.
"""

import pytest

from delfin.agent.api_client import (_MCP_SCHEMA_FLOOR_CHARS,
                                     _mcp_schema_budget_chars)

# What DELFIN's own servers cost, measured 2026-09-08.
_MEASURED_SURFACE_CHARS = 54_241


def test_a_large_window_is_offered_the_whole_framework():
    """524288 tokens is what the KIT flagships have."""
    assert _mcp_schema_budget_chars(524_288) >= _MEASURED_SURFACE_CHARS


def test_a_small_window_keeps_the_cap_it_had():
    """The same surface is 42% of a 32k context. Nothing shrinks, and
    nothing grows either."""
    assert _mcp_schema_budget_chars(32_768) == _MCP_SCHEMA_FLOOR_CHARS
    assert _mcp_schema_budget_chars(8_192) == _MCP_SCHEMA_FLOOR_CHARS


def test_an_unknown_window_falls_back_to_the_old_number():
    for value in (0, None, -1, "", "nonsense"):
        assert _mcp_schema_budget_chars(value) == _MCP_SCHEMA_FLOOR_CHARS


def test_the_budget_grows_with_the_window_and_never_shrinks():
    windows = [8_192, 32_768, 128_000, 262_144, 524_288]
    budgets = [_mcp_schema_budget_chars(w) for w in windows]
    assert budgets == sorted(budgets)
    assert min(budgets) == _MCP_SCHEMA_FLOOR_CHARS


def test_an_explicit_setting_still_wins(tmp_path, monkeypatch):
    """Someone who has measured their own surface must not be
    second-guessed by a share of the window."""
    import json

    path = tmp_path / "settings.json"
    path.write_text(json.dumps({"agent": {"mcp_schema_budget_chars": 3000}}))
    monkeypatch.setattr("delfin.user_settings.get_settings_path",
                        lambda *a, **k: path)
    assert _mcp_schema_budget_chars(524_288) == 3000


def test_the_share_stays_a_small_part_of_the_window():
    """A budget is not a licence: the tool surface must not become the
    request."""
    for window in (128_000, 262_144, 524_288):
        budget_tokens = _mcp_schema_budget_chars(window) / 4
        assert budget_tokens <= window * 0.06, window


def test_the_caller_passes_the_resolved_window():
    """A budget that scales with the window does nothing if the call site
    keeps asking without one."""
    import inspect

    from delfin.agent import api_client

    src = inspect.getsource(api_client.OpenAIClient.stream_message)
    assert "_mcp_schema_budget_chars(" in src
    i = src.index("_mcp_schema_budget_chars(")
    assert "context_window" in src[i:i + 160]
