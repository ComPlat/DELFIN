"""One bad moment at startup should not size the context for a session.

``resolve`` asks the endpoint what window a model has, falls back to a
curated static table when the wire is unavailable, and to a name
heuristic when even that has no entry. Then it cached the answer — any
answer — for the life of the process. So a probe that timed out, 401'd
or raced a slow endpoint pinned a GUESS until the process exited, with
nothing to make it try again.

The numbers make it concrete. Both KIT flagships declare 512K. The static
fallback for them is a deliberate, conservative 131072, and the heuristic
below that is 32768. Compaction fires at 95% of whatever this returns, so
a session that started badly compacts at ~124k or ~31k instead of ~498k —
four to sixteen times too early. On GLM every compaction is a cold prompt
head worth about 200 seconds, so the cost of one unlucky probe is paid
for the rest of the day.

Observed 2026-09-08: two engines built seconds apart, one resolving
524288 and the other 32768.

So a live answer is cached as before, and a guess is cached only long
enough to stop a burst of calls hammering an endpoint that is down.
"""

from __future__ import annotations

import time

import pytest

from delfin.agent import model_capabilities as MC


@pytest.fixture(autouse=True)
def _clean_cache():
    MC.clear_cache()
    yield
    MC.clear_cache()


_URL = "https://endpoint.invalid/v1"


def _probe_returning(value, calls):
    def _fake(base_url, model, api_key=""):
        calls.append(model)
        return dict(value) if value else {}
    return _fake


def test_a_live_answer_is_kept(monkeypatch):
    calls = []
    monkeypatch.setattr(MC, "_discover_openai_models",
                        _probe_returning({"context_window": 524_288}, calls))
    first = MC.resolve("kit", "kit.glm-5.3", _URL, api_key="k")
    second = MC.resolve("kit", "kit.glm-5.3", _URL, api_key="k")
    assert first.source == "live" and first.context_window == 524_288
    assert second.context_window == 524_288
    assert len(calls) == 1, "a live answer was re-probed"


def test_a_guess_is_not_kept_for_the_session(monkeypatch):
    calls = []
    monkeypatch.setattr(MC, "_discover_openai_models",
                        _probe_returning(None, calls))
    first = MC.resolve("kit", "kit.glm-5.3", _URL, api_key="k")
    assert first.source != "live"

    # The endpoint comes back.
    monkeypatch.setattr(MC, "_discover_openai_models",
                        _probe_returning({"context_window": 524_288}, calls))
    monkeypatch.setattr(MC, "_PROVISIONAL_TTL_S", 0.0)
    again = MC.resolve("kit", "kit.glm-5.3", _URL, api_key="k")
    assert again.source == "live", "the failed probe was never retried"
    assert again.context_window == 524_288


def test_a_burst_does_not_hammer_an_endpoint_that_is_down(monkeypatch):
    """The cache still has a job: an endpoint that is down must not be
    asked once per call."""
    calls = []
    monkeypatch.setattr(MC, "_discover_openai_models",
                        _probe_returning(None, calls))
    for _ in range(5):
        MC.resolve("kit", "kit.glm-5.3", _URL, api_key="k")
    assert len(calls) == 1, f"probed {len(calls)} times while down"


def test_the_provisional_window_is_still_a_real_number(monkeypatch):
    """Retrying later must not mean answering with nothing now."""
    monkeypatch.setattr(MC, "_discover_openai_models",
                        _probe_returning(None, []))
    caps = MC.resolve("kit", "kit.glm-5.3", _URL, api_key="k")
    assert caps.context_window >= 32_768


def test_only_a_live_answer_reaches_the_disk(monkeypatch, tmp_path):
    """The disk cache outlives the process; a guess written there would
    outlive the outage that produced it."""
    saved = []
    monkeypatch.setattr(MC, "_save_disk_cache", lambda: saved.append(1))
    monkeypatch.setattr(MC, "_discover_openai_models", _probe_returning(None, []))
    MC.resolve("kit", "kit.glm-5.3", _URL, api_key="k")
    assert not saved
    MC.clear_cache()
    monkeypatch.setattr(MC, "_discover_openai_models",
                        _probe_returning({"context_window": 524_288}, []))
    MC.resolve("kit", "kit.glm-5.3", _URL, api_key="k")
    assert saved


# ---------------------------------------------------------------------------
# ...and something has to ask again
# ---------------------------------------------------------------------------

@pytest.fixture
def agent_tree(tmp_path):
    """The minimal pack the engine needs to compose a prompt at all.

    Same shape as tests/test_agent_engine.py: a `pack` for the role and
    shared files, a `pack_lite` holding the manifest and modes.
    """
    import textwrap
    shared = tmp_path / "pack" / "shared"
    shared.mkdir(parents=True)
    agents = tmp_path / "pack" / "agents"
    agents.mkdir()
    for name in ("delfin_context.md", "work_cycle_rules.md",
                 "goal_decomposition_rules.md",
                 "universal_input_template.md", "minimal_final_verdict.md"):
        (shared / name).write_text("#")
    (agents / "solo_agent.md").write_text("# Solo Agent")
    modes = tmp_path / "pack_lite" / "modes"
    modes.mkdir(parents=True)
    (modes / "solo.md").write_text("# solo mode")
    (tmp_path / "pack_lite" / "manifest.yaml").write_text(textwrap.dedent("""\
        pack_name: TEST
        version: 1
        modes:
          - id: solo
            file: modes/solo.md
            route:
              - solo_agent
        """))
    return tmp_path


def test_a_turn_asks_again_while_the_window_is_still_a_guess(agent_tree,
                                                             monkeypatch):
    """The capability layer lets a guess go stale. That heals nothing on
    its own: the engine probes once, at construction, on a daemon thread.
    If that lost a race with a slow endpoint it kept the fallback for the
    session, and compaction fires at 95% of it."""
    from unittest.mock import MagicMock, patch

    from delfin.agent.api_client import StreamEvent
    from delfin.agent.engine import AgentEngine

    def _stream(system, messages, max_tokens=4096, session_id="",
                thinking_budget=0):
        yield StreamEvent(type="text_delta", text="ok")
        yield StreamEvent(type="message_delta", output_tokens=2)

    client = MagicMock()
    client.stream_message = MagicMock(side_effect=_stream)
    with patch("delfin.agent.engine.create_client", return_value=client):
        engine = AgentEngine(repo_dir=agent_tree, backend="cli",
                             mode="solo", pack_dir=agent_tree)

    refreshed: list[bool] = []
    monkeypatch.setattr(
        engine, "_refresh_context_window",
        lambda *, background=False: refreshed.append(background))

    engine._active_capabilities = None            # the probe never landed
    engine.stream_response("hallo")
    assert refreshed == [True], "a guessed window was never re-asked"

    refreshed.clear()
    engine._active_capabilities = MC.ModelCapabilities(
        model="m", provider="kit", context_window=524_288,
        max_output_tokens=8192, supports_tools=True, supports_vision=False,
        is_reasoning=False, thinking_tagged=False,
        recommended_effort="medium", num_ctx_override=None,
        source="live", discovered_at=0.0, note="")
    engine.stream_response("nochmal")
    assert refreshed == [], "a settled window was probed again"
