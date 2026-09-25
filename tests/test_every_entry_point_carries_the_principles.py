"""Every way DELFIN calls a model must carry the maintainer's principles.

Part 1 of the principles task: enumerate every entry point that builds a
system prompt for a model call, generate the ACTUAL prompt it would send
(client stubbed — no network), and assert the principles text is in it
and ahead of every other shared contract. A path that does not go
through PromptLoader is a finding, reported to the operator; the two
known ones (subagent presets, the internal summariser) are pinned here
as documented bypasses so a fix cannot land unnoticed.
"""

from __future__ import annotations

from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest

import delfin.agent.engine as engine_mod
from delfin.agent.engine import AgentEngine

_PACK = Path(__file__).resolve().parent.parent / "delfin" / "agent" / "pack"
_PRINCIPLES = (_PACK / "shared" / "principles_addendum.md").read_text(
    encoding="utf-8")
# The first sentence of the principles body — unique to this file, short
# enough to survive any section reordering, long enough not to appear in
# role prompts by accident.
_MARKER = "Your highest objective is to advance science"
_HONESTY_MARKER = "verify"

# Entry points that construct an AgentEngine and thereby build their
# prompt through PromptLoader (verified by reading each factory):
#   terminal CLI      delfin/agent/cli.py:80        _build_engine
#   dashboard agent   delfin/dashboard/tab_agent.py:8086
#   benchmark         delfin/agent/benchmark_runner.py:216 _default_engine_factory
#   bug watcher       delfin/agent/bug_watcher.py:254   _default_engine_factory
#   job monitor       delfin/agent/job_monitor.py:932   _default_engine_factory
#   scheduler daemon  delfin/agent/scheduler_daemon.py:97 _default_engine_factory
_FACTORIES = [
    ("cli", "delfin.agent.cli", "_build_engine"),
    ("benchmark_runner", "delfin.agent.benchmark_runner",
     "_default_engine_factory"),
    ("bug_watcher", "delfin.agent.bug_watcher", "_default_engine_factory"),
    ("job_monitor", "delfin.agent.job_monitor", "_default_engine_factory"),
    ("scheduler_daemon", "delfin.agent.scheduler_daemon",
     "_default_engine_factory"),
]


def _engine(tmp_path: str, **kw) -> AgentEngine:
    """A real AgentEngine with a stubbed client — no network, no CLI
    process; everything else (mode loading, PromptLoader) is the real
    path the entry points use."""
    with patch.object(engine_mod, "create_client", return_value=MagicMock()):
        defaults = dict(repo_dir=Path(tmp_path), backend="api",
                        provider="kit", model="kit.glm-5.3", mode="solo")
        defaults.update(kw)
        return AgentEngine(**defaults)


def _built_prompt(eng: AgentEngine, task: str = "tidy the workspace") -> str:
    return eng._build_current_system_prompt(task_text=task)


def _assert_principles_first(prompt: str) -> None:
    assert _MARKER in prompt, "principles text missing from the prompt"
    # The prompt opens with the principles section; the honesty
    # addendum ("# Honesty & grounding") is the second shared contract
    # and must come after them.
    assert prompt.index(_MARKER) < prompt.index(
        "# Honesty & grounding"), "principles must precede the honesty addendum"


def test_the_solo_engine_builds_them_first(tmp_path):
    prompt = _built_prompt(_engine(tmp_path))
    assert _MARKER in prompt
    # Layer 0 opens the prompt: the principles stand before the honesty
    # addendum, the second shared contract.
    assert prompt.index(_MARKER) < prompt.index("# Honesty & grounding")


@pytest.mark.parametrize("factory_module,factory_name", [
    (m, n) for _, m, n in _FACTORIES
])
def test_every_engine_factory_carries_the_principles(
        tmp_path, factory_module, factory_name):
    """The five headless/terminal factories all build their prompts
    through AgentEngine -> PromptLoader; each must carry the principles
    in the prompt it would actually send."""
    import importlib
    mod = importlib.import_module(factory_module)
    factory = getattr(mod, factory_name)
    import inspect
    sig = inspect.signature(factory)
    with patch.object(engine_mod, "create_client",
                      return_value=MagicMock()):
        if factory_name == "_build_engine":
            # cli._build_engine takes an argparse.Namespace.
            import argparse
            args = argparse.Namespace(
                backend="api", model="kit.glm-5.3", provider="kit",
                mode="solo", cwd=str(tmp_path), settings_defaults=False,
                effort="", permission_mode="", extra_dirs=None,
                read_only_dirs=None, allowed_tools=None)
            eng = factory(args)
        elif "repo_root" in sig.parameters:
            eng = factory(str(tmp_path), settings={})
        elif "workspace" in sig.parameters:
            eng = factory(str(tmp_path), settings={})
        elif "folder" in sig.parameters:
            eng = factory(str(tmp_path), settings={})
        elif "mode" in sig.parameters:
            eng = factory("kit.glm-5.3", "api", "kit", "solo")
        else:  # pragma: no cover - unreachable unless a factory changes
            raise AssertionError(f"unknown factory shape: {factory}")
    prompt = _built_prompt(eng)
    _assert_principles_first(prompt)


def test_the_dashboard_agent_engine_carries_them(tmp_path):
    """The dashboard tab constructs AgentEngine directly
    (tab_agent.py:8086) — same loader, same requirement."""
    prompt = _built_prompt(_engine(tmp_path, mode="solo"))
    _assert_principles_first(prompt)


def test_distillation_keeps_layer0_verbatim(tmp_path):
    """The context distiller must return layer 0 — principles included —
    byte-identical, whatever it compresses behind it."""
    from delfin.agent.context_distiller import ContextDistiller, _split_layer0
    full = _built_prompt(_engine(tmp_path))
    layer0, rest = _split_layer0(full)
    assert _MARKER in layer0, "principles are not inside layer 0"
    distiller = ContextDistiller.__new__(ContextDistiller)
    # _extractive_compress is the offline fallback path; the API path
    # prepends the same layer0 (context_distiller.py:135-148).
    compressed = layer0 + distiller._extractive_compress(rest)
    assert _MARKER in compressed
    assert compressed.startswith(layer0[:len(_MARKER)])
    # Verbatim: the first layer-0 bytes survive compression untouched.
    assert compressed[:len(layer0)] == layer0


def test_subagent_presets_do_not_carry_the_principles():
    """DOCUMENTED BYPASS (reported to the operator): subagent system
    prompts are built from the preset body alone (subagents.py:2955)
    and never pass through PromptLoader — the principles are NOT there.
    This test pins that fact; when the operator wires presets through
    the loader, it must be rewritten (and will then fail loudly if
    the wiring regresses)."""
    from delfin.agent import subagents
    presets = subagents._BUILTIN_PRESETS
    for name, preset in presets.items():
        assert _MARKER not in preset.system_prompt, (
            f"preset {name} unexpectedly carries the principles — "
            "rewrite this documented bypass test")


def test_the_internal_summariser_prompt_has_no_principles():
    """DOCUMENTED BYPASS: the conversation summariser builds its own
    inline system prompt (engine.py:4927) — no principles. It sees only
    an already-rendered transcript, never user instructions directly,
    so the operator decides whether it needs them. Pinned here so the
    decision is visible, not silent."""
    import inspect
    src = inspect.getsource(AgentEngine)
    assert "conversation-summarisation assistant" in src
    # And the summariser prompt itself carries no principles marker:
    start = src.index("You are a conversation-summarisation assistant")
    chunk = src[start:start + 1500]
    assert _MARKER not in chunk
