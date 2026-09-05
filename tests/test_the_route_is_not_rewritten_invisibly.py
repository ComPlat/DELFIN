"""A recorded success rate rewrote the route, invisibly and from nowhere.

`recommend_task_route` escalates quick -> reviewed whenever a stored
task-success rate is below 0.75. No setting gates it, and the reason it
appends to the result is written by the engine and read by nothing.

Measured on this machine:

  * ~/.delfin/outcome_history.jsonl  55 records, 55 PASS, modes solo and
    dashboard only -- not one `quick` record exists.
  * ~/.delfin/provider_profile_state.json  kit: solo 0.005, coding 0.007;
    claude: quick 0.021.

Two accumulators for the same fact, with no reconciliation:
`update_from_outcome` never reads the history, and the two writes sit in
separate try/except blocks so the profile mutates even when the audit
write fails. The EMA seeds at 0.5, which is below the threshold, so a
task class recovering from a single FAIL escalates until several
consecutive passes.

And the escalation reads the wrong provider. The lookup falls back to
"claude" when no class-level hint was set, so a session running on KIT is
routed by numbers recorded for Anthropic -- which is why the live reason
string says 17% while KIT's own recorded rate for the same class is 0.7%.

It has already corrupted two tests, both of which carry written-up
workarounds saying the assertion measures the box rather than the code.

An adaptive mechanism has to be visible, resettable, and reading its own
provider. Until it is, it is off.
"""

from __future__ import annotations

import pytest

from delfin.agent.engine import AgentEngine as E
from delfin.agent import provider_profile as P


@pytest.fixture(autouse=True)
def neutral_profile(monkeypatch):
    """Isolate from whatever this installation has recorded."""
    monkeypatch.setattr(
        "delfin.agent.provider_profile.load_provider_profile",
        lambda *a, **kw: {"task_performance": {"coding": {"success_rate": 0.1}}})


def _route(**kw):
    # A plain coding task with no chemistry terms: those escalate on
    # their own, which would leave mode != "quick" before the adaptive
    # block runs and make this test prove nothing.
    return E.recommend_task_route(
        "fix the off-by-one in agent_metrics.py",
        current_mode="dashboard", is_delfin_workspace=True, **kw)


# ---------------------------------------------------------------------------
# Off by default
# ---------------------------------------------------------------------------

def test_a_low_recorded_rate_does_not_silently_escalate(monkeypatch):
    monkeypatch.setattr("delfin.user_settings.load_settings",
                        lambda *a, **kw: {})
    out = _route()
    assert not any("adaptive" in r for r in out.get("reasons", [])), out


def test_it_can_be_switched_on(monkeypatch):
    monkeypatch.setattr(
        "delfin.user_settings.load_settings",
        lambda *a, **kw: {"agent": {"routing": {"adaptive_escalation": True}}})
    out = _route()
    assert any("adaptive" in r for r in out.get("reasons", [])), out


def test_switching_it_on_still_reads_its_own_provider(monkeypatch):
    """The lookup fell back to "claude", so a KIT session was routed by
    numbers recorded for Anthropic."""
    seen = {}

    def spy(provider="", *a, **kw):
        seen["provider"] = provider
        return {"task_performance": {"coding": {"success_rate": 0.1}}}

    monkeypatch.setattr(
        "delfin.agent.provider_profile.load_provider_profile", spy)
    monkeypatch.setattr(
        "delfin.user_settings.load_settings",
        lambda *a, **kw: {"agent": {"routing": {"adaptive_escalation": True}}})
    monkeypatch.setattr(E, "_active_provider", "kit", raising=False)
    _route()
    assert seen.get("provider") == "kit", seen


# ---------------------------------------------------------------------------
# The learned state can be seen and undone
# ---------------------------------------------------------------------------

def test_reset_profile_exists_and_is_reachable():
    """The module docstring promises "fully reversible: reset_profile()".
    It had zero callers repo-wide -- no command, no flag, no UI."""
    assert callable(getattr(P, "reset_profile", None))
    import pathlib
    tab = (pathlib.Path(__file__).resolve().parents[1] / "delfin"
           / "dashboard" / "tab_agent.py").read_text(encoding="utf-8")
    assert '"/profile"' in tab, "the learned routing state is still unviewable"
    assert "reset_profile" in tab, "the learned routing state is still unresettable"


def test_the_command_can_show_and_reset():
    import pathlib
    tab = (pathlib.Path(__file__).resolve().parents[1] / "delfin"
           / "dashboard" / "tab_agent.py").read_text(encoding="utf-8")
    i = tab.index('if cmd.startswith("/profile")')
    body = tab[i:i + 2500]
    assert "reset" in body
    assert "load_provider_profile" in body
