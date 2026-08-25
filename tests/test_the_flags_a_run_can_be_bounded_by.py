"""The bounds a run can be given on the command line, and the ones it cannot.

Every flag here was named in the documentation and declared nowhere, which
is the quieter half of the same defect `--color` and `--effort` had: a flag
that parses and reaches nothing. So each one is asserted twice — the parser
takes it, AND the value arrives at the object that enforces it — because a
namespace assertion alone is exactly what those two flags would have passed.

Where a flag cannot deliver what its name promises, the run has to SAY so
at startup, the way `--isolate` does when bubblewrap is missing. Those
paths are asserted too: silent non-delivery is the specific failure this
whole area exists to close.
"""

from __future__ import annotations

import argparse

import pytest

from delfin.agent import cli as agent_cli


class _Recorder:
    """Stands in for AgentEngine and keeps what it was constructed with."""

    last: dict = {}

    def __init__(self, **kwargs):
        type(self).last = dict(kwargs)
        self.client = type("C", (), {"model": "kit.qwen3-coder-480b"})()
        self.provider = kwargs.get("provider", "") or "kit"


@pytest.fixture
def recorded(monkeypatch):
    import delfin.agent.engine as engine_mod
    monkeypatch.setattr(engine_mod, "AgentEngine", _Recorder)
    _Recorder.last = {}
    return _Recorder


def _args(**over) -> argparse.Namespace:
    base = dict(backend="", provider="", model="", mode="", cwd="",
                effort="", permission_mode="", settings_defaults=False)
    base.update(over)
    return argparse.Namespace(**base)


def _parse(*argv):
    return agent_cli.build_parser().parse_args(["chat", *argv])


# ---------------------------------------------------------------------------
# --max-budget-usd / --max-run-seconds: the ceiling for the whole session
# ---------------------------------------------------------------------------

def test_the_budget_flags_are_offered_on_the_command_line():
    args = _parse("--max-budget-usd", "2.50", "--max-run-seconds", "900")
    assert args.max_budget_usd == pytest.approx(2.50)
    assert args.max_run_seconds == pytest.approx(900.0)


def test_the_budget_reaches_the_attribute_the_engine_enforces(recorded, tmp_path):
    """`_run_budget` reads these two off the instance, above the settings
    file — the precedence the scheduler daemon already uses per entry."""
    engine = agent_cli._build_engine(
        _args(cwd=str(tmp_path), max_budget_usd=3.0, max_run_seconds=120.0))
    assert engine.run_budget_usd == pytest.approx(3.0)
    assert engine.run_budget_s == pytest.approx(120.0)


def test_the_real_engine_reads_what_the_flag_wrote(tmp_path):
    """Against `AgentEngine._run_budget` itself, not a stand-in.

    The attribute name is the whole contract between this file and the
    enforcement, so it is asserted through the method that reads it.
    """
    from delfin.agent.engine import AgentEngine

    engine = AgentEngine.__new__(AgentEngine)      # no client, no network
    agent_cli._apply_run_budget(
        engine, _args(max_budget_usd=7.5, max_run_seconds=60.0))
    assert AgentEngine._run_budget(engine) == (pytest.approx(7.5),
                                               pytest.approx(60.0))


def test_no_budget_flag_leaves_the_configured_one_in_charge(recorded, tmp_path):
    engine = agent_cli._build_engine(_args(cwd=str(tmp_path)))
    assert not hasattr(engine, "run_budget_usd")
    assert not hasattr(engine, "run_budget_s")


def test_a_namespace_that_predates_the_flags_still_builds(recorded, tmp_path):
    ns = argparse.Namespace(backend="", provider="", model="", mode="",
                            cwd=str(tmp_path))
    agent_cli._build_engine(ns)
    assert _Recorder.last.get("mode") == "solo"


def test_a_usd_ceiling_on_an_unpriced_model_says_it_cannot_fire(monkeypatch):
    """The can't-deliver path.

    `cost_usd` sums only the turns whose price could be looked up, so on a
    model with no published rate the fraction spent stays at 0% for a run
    of any size and the ceiling never fires.
    """
    monkeypatch.setattr(agent_cli, "_usd_ceiling_measurable", lambda e: False)
    engine = type("E", (), {"client": type("C", (), {"model": "mystery-7b"})(),
                            "provider": "ollama"})()
    notes = agent_cli._bounding_notices(_args(max_budget_usd=5.0), engine)
    assert any("not enforceable" in n and "mystery-7b" in n for n in notes)
    assert any("--max-run-seconds" in n for n in notes)


def test_a_usd_ceiling_on_a_priced_model_says_nothing(monkeypatch):
    monkeypatch.setattr(agent_cli, "_usd_ceiling_measurable", lambda e: True)
    engine = type("E", (), {"client": type("C", (), {"model": "m"})(),
                            "provider": "claude"})()
    assert agent_cli._bounding_notices(_args(max_budget_usd=5.0), engine) == []


def test_the_price_question_is_asked_of_the_pricing_table(monkeypatch):
    """`_usd_ceiling_measurable` must consult the real rate lookup —
    None from `price_for` means the price is unknown, never that it is
    free."""
    from delfin.agent import pricing

    engine = type("E", (), {"client": type("C", (), {"model": "x"})(),
                            "provider": "kit"})()
    monkeypatch.setattr(pricing, "price_for", lambda m, p: None)
    assert agent_cli._usd_ceiling_measurable(engine) is False
    monkeypatch.setattr(pricing, "price_for", lambda m, p: (3.0, 15.0))
    assert agent_cli._usd_ceiling_measurable(engine) is True
