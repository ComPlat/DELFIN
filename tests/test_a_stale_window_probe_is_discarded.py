"""A background window probe never replaces a newer window.

``AgentEngine.__init__`` (and every turn whose capabilities are still a
guess) resolves the model's context window on a daemon thread. The
thread used to assign ``context_window_tokens`` and
``_active_capabilities`` whenever it finished, so a window set while it
was in flight -- by a model switch and its synchronous refresh, by hand,
or by a caller after ``restore_state`` -- was replaced by the late
result: the previous model's window after a quick switch, or the
resolved window over a test's 4000.

Contract pinned here: a background result is applied only if no window
or capability assignment happened since the probe started and the client
still talks to the model the probe resolved for; otherwise it is
discarded, the window AND the capabilities both.
"""
from __future__ import annotations

import textwrap
import threading
from types import SimpleNamespace
from unittest.mock import MagicMock, patch

import pytest

from delfin.agent import engine as engine_mod
from delfin.agent import model_capabilities as mc
from delfin.agent.engine import AgentEngine

PROBE_WINDOW = 262_144


@pytest.fixture
def agent_tree(tmp_path):
    agent_dir = tmp_path / "pack"
    shared = agent_dir / "shared"
    shared.mkdir(parents=True)
    agents = agent_dir / "agents"
    agents.mkdir()
    (shared / "delfin_context.md").write_text("# DELFIN Context\nTest.")
    (shared / "work_cycle_rules.md").write_text("# Work Cycle Rules\nRule 1.")
    (shared / "universal_input_template.md").write_text("# Input Template")
    (shared / "minimal_final_verdict.md").write_text("# Verdict")
    (agents / "solo_agent.md").write_text("# Solo Agent\nYou are solo.")
    modes = tmp_path / "pack_lite" / "modes"
    modes.mkdir(parents=True)
    (modes / "solo.md").write_text("# Mode: quick\nQuick mode.")
    (tmp_path / "pack_lite" / "manifest.yaml").write_text(textwrap.dedent("""\
        pack_name: DELFIN_AGENT_LITE
        version: 1
        recommended_default_mode: solo
        modes:
          - id: solo
            file: modes/solo.md
            route:
              - solo_agent
    """))
    return tmp_path


def _caps(model: str, window: int):
    return SimpleNamespace(model=model, context_window=window,
                           max_output_tokens=8192, supports_tools=True,
                           supports_vision=False, source="live")


@pytest.fixture
def gated(monkeypatch):
    """A resolver whose answer for ``model-a`` waits for ``release``;
    every other model answers at once with a window of its own."""
    release = threading.Event()
    answers = {"model-a": _caps("model-a", PROBE_WINDOW),
               "model-b": _caps("model-b", 65_536)}

    def resolve(provider, model, base_url="", **kw):
        answer = answers[model]          # what the backend said at the ask
        if model == "model-a":
            release.wait(timeout=5)
        return answer

    monkeypatch.setattr(mc, "resolve", resolve)
    # The learned floor and the settings cap read the user's stores; the
    # contract under test is the ordering, not those two numbers.
    from delfin.agent import observed_window as ow
    monkeypatch.setattr(ow, "at_least", lambda model, window: window)
    monkeypatch.setattr(engine_mod, "_capped_context_window", lambda w: w)
    return SimpleNamespace(release=release, answers=answers)


def _engine(agent_tree, model="model-a"):
    client = MagicMock()
    client.session_id = "s"
    client.model = model
    with patch("delfin.agent.engine.create_client", return_value=client):
        engine = AgentEngine(repo_dir=agent_tree, backend="cli",
                             mode="quick", pack_dir=agent_tree)
    return engine, client


def _land(gated):
    """Let every probe in flight finish."""
    gated.release.set()
    for t in threading.enumerate():
        if t.name == "caps-prewarm":
            t.join(timeout=5)


def test_the_probe_applies_when_nothing_changed(agent_tree, gated):
    engine, _ = _engine(agent_tree)
    assert engine.context_window_tokens == 100_000
    _land(gated)
    assert engine.context_window_tokens == PROBE_WINDOW
    assert engine._active_capabilities is gated.answers["model-a"]


def test_a_direct_assignment_outlives_the_probe(agent_tree, gated):
    engine, _ = _engine(agent_tree)
    engine.context_window_tokens = 4000
    _land(gated)
    assert engine.context_window_tokens == 4000
    # The capabilities belong to the discarded result: they stay unset
    # rather than describing a window the engine does not use.
    assert engine._active_capabilities is None


def test_an_assigned_capability_outlives_the_probe(agent_tree, gated):
    engine, _ = _engine(agent_tree)
    mine = _caps("model-a", 32_768)
    engine._active_capabilities = mine
    _land(gated)
    assert engine._active_capabilities is mine
    assert engine.context_window_tokens == 100_000


def test_a_switch_to_another_model_outlives_the_old_probe(agent_tree, gated):
    engine, client = _engine(agent_tree)
    client.model = "model-b"
    engine._refresh_context_window()          # the switch path's refresh
    assert engine.context_window_tokens == 65_536
    _land(gated)
    assert engine.context_window_tokens == 65_536
    assert engine._active_capabilities is gated.answers["model-b"]


def test_a_switch_and_back_still_discards_the_first_probe(agent_tree, gated):
    """Same model at the end, but a synchronous refresh happened since
    the probe started: its answer is the newer one."""
    engine, client = _engine(agent_tree)
    client.model = "model-b"
    engine._refresh_context_window()
    gated.answers["model-a"] = _caps("model-a", 131_072)
    client.model = "model-a"
    gated.release.set()
    engine._refresh_context_window()
    _land(gated)                 # the first probe answers PROBE_WINDOW
    assert engine.context_window_tokens == 131_072


def test_a_model_changed_without_a_refresh_discards_the_probe(
        agent_tree, gated):
    """The model check alone: nothing assigned a window, but the client
    now talks to another model, so the answer describes the wrong one."""
    engine, client = _engine(agent_tree)
    client.model = "model-b"
    _land(gated)
    assert engine.context_window_tokens == 100_000
    assert engine._active_capabilities is None


def test_restore_then_a_late_probe_keeps_the_restored_window(
        agent_tree, gated):
    engine, _ = _engine(agent_tree)
    saved = engine.export_state()
    engine.context_window_tokens = 8192       # the caller's window
    engine.restore_state(saved)
    _land(gated)
    assert engine.context_window_tokens == 8192


def test_restore_leaves_a_valid_probe_valid(agent_tree, gated):
    """A restore changes the conversation, not the model: the probe for
    the same model still lands, or a resumed long session would run its
    first turns at the 100k fallback and compact at 95k."""
    engine, _ = _engine(agent_tree)
    engine.restore_state(engine.export_state())
    _land(gated)
    assert engine.context_window_tokens == PROBE_WINDOW


def test_a_turn_reprobe_does_not_replace_an_assigned_window(
        agent_tree, gated):
    """The per-turn re-ask (capabilities still a guess) starts AFTER the
    assignment, so the generation alone would let it through."""
    engine, _ = _engine(agent_tree)
    _land(gated)
    engine._active_capabilities = None
    engine.context_window_tokens = 4000
    gated.release.clear()
    engine._refresh_context_window(background=True)
    _land(gated)
    assert engine.context_window_tokens == 4000


def test_an_engine_built_without_init_reads_the_default():
    engine = object.__new__(AgentEngine)
    assert engine.context_window_tokens == 100_000
    assert engine._active_capabilities is None
    engine.context_window_tokens = 1234
    assert engine.context_window_tokens == 1234
    del engine.context_window_tokens
    assert engine.context_window_tokens == 100_000
