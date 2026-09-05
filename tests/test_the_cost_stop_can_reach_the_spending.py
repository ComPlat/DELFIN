"""The two stops that were meant to bound a runaway turn, and could not.

THE COST BREAKER RAN ONCE, AFTER THE MONEY. It lives in the engine's
`message_delta` handler, and every shipped client emits `message_delta`
exactly once, at the end of a turn. So a turn could run its full round
budget -- up to 500 rounds, each a fresh request carrying the whole
accumulated tool context -- and the ceiling was evaluated afterwards, as
a post-mortem. The same defect disabled the scheduler's per-entry
`budget_usd`, which is implemented by overriding that same cap.

THE PER-MODEL ROUND CAPS WERE UNREACHABLE. `_resolve_max_tool_rounds`
prefers an explicit `agent.max_tool_rounds` setting and falls back to the
model profile -- but the shipped defaults always supply 500, so "unset"
was indistinguishable from "set to 500" and the fallback never ran.
Measured before the fix: profiles declaring 10, 20 and 50 rounds all
resolved to 500. The stated purpose -- a weak model's degenerate loop
dying early instead of burning hundreds of rounds -- had never once
happened.
"""

from __future__ import annotations

import contextlib
from unittest.mock import MagicMock, patch

import pytest

from delfin.agent import api_client as A
from delfin.agent import engine as E


@contextlib.contextmanager
def _settings(payload):
    with patch("delfin.user_settings.load_settings", return_value=payload):
        yield


# ---------------------------------------------------------------------------
# The breaker can be consulted where the spending happens
# ---------------------------------------------------------------------------

def test_the_loop_can_ask_for_the_cap():
    client = A.OpenAIClient.__new__(A.OpenAIClient)
    assert client._turn_cost_cap() == 0.0, "no cap installed means no cap"

    client.turn_cost_cap = lambda: 12.5
    assert client._turn_cost_cap() == 12.5


def test_a_broken_cap_probe_does_not_end_the_turn():
    client = A.OpenAIClient.__new__(A.OpenAIClient)

    def boom():
        raise RuntimeError("engine gone")

    client.turn_cost_cap = boom
    assert client._turn_cost_cap() == 0.0


def test_the_engine_installs_it(tmp_path):
    from delfin.agent.api_client import StreamEvent

    with patch("delfin.agent.engine.create_client", return_value=MagicMock()):
        eng = E.AgentEngine(repo_dir=tmp_path, backend="api", provider="kit",
                            model="kit.qwen3.5-397b-A17b", mode="solo")

    def fine(**kw):
        yield StreamEvent(type="text_delta", text="ok")

    eng.client.stream_message = MagicMock(side_effect=fine)
    eng.stream_response("hi")
    probe = getattr(eng.client, "turn_cost_cap", None)
    assert callable(probe)
    assert probe() == pytest.approx(eng._cost_hard_cap())


def test_the_loop_actually_checks_it():
    """A cap nobody reads inside the loop is the defect this replaces."""
    import pathlib
    src = pathlib.Path(A.__file__).read_text(encoding="utf-8")
    code = "\n".join(l for l in src.splitlines()
                     if not l.lstrip().startswith("#"))
    start = code.index("for _round in range(_MAX_TOOL_ROUNDS + 1):")
    window = code[start:start + 3000]
    assert "_turn_cost_cap" in window, (
        "the cost ceiling is still only evaluated after the turn")


# ---------------------------------------------------------------------------
# The per-model round caps are reachable again
# ---------------------------------------------------------------------------

def test_an_unset_setting_lets_the_model_profile_decide():
    with _settings({"agent": {}}):
        assert A._resolve_max_tool_rounds("llama3:8b") < 100, (
            "a weak model still gets the frontier round budget")


def test_the_shipped_default_does_not_mask_the_profile():
    """The bug: DEFAULT_SETTINGS always supplied a number, so "unset"
    could not be expressed and the profile branch was dead code."""
    from delfin import user_settings
    assert user_settings.DEFAULT_SETTINGS["agent"]["max_tool_rounds"] is None


def test_a_real_setting_still_wins():
    with _settings({"agent": {"max_tool_rounds": 1200}}):
        assert A._resolve_max_tool_rounds("llama3:8b") == 1200


def test_zero_still_means_no_round_cap():
    with _settings({"agent": {"max_tool_rounds": 0}}):
        assert A._resolve_max_tool_rounds("llama3:8b") >= 100_000


def test_an_unknown_model_still_gets_a_bound():
    with _settings({"agent": {}}):
        n = A._resolve_max_tool_rounds("some-model-nobody-profiled")
        assert 0 < n <= 500
