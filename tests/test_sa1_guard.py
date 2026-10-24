"""SA1 phase 1: the explore-pilot guard.

The pilot lets ONE narrow use of the existing ``explore`` sub-agent through
its measurement harness and refuses everything else exactly as today:

- only the ``explore`` type is allowed,
- at most one sub-agent running per session,
- no nesting (a delegate may not start another delegate),
- every refusal names the reason (wave-15b: a refusal names a ROUTE).

These tests drive the PUBLIC entry, ``explore_pilot.run_pilot_subagent``,
the function the harness (phase 2) calls — never an inner helper.
"""
import pytest

from delfin.agent import explore_pilot


def _kwargs(**over):
    kw = dict(subagent_type="explore", description="locate X",
              prompt="find where DELFIN decides Y, file:line", model="",
              inherit_context=True)
    kw.update(over)
    return kw


class TestOnlyExploreType:
    def test_explore_is_allowed(self):
        ok, why = explore_pilot.check_delegate("explore")
        assert ok, why

    def test_plan_is_refused(self):
        ok, why = explore_pilot.check_delegate("plan")
        assert not ok
        assert "explore" in why

    def test_general_purpose_is_refused(self):
        ok, why = explore_pilot.check_delegate("general-purpose")
        assert not ok

    def test_unknown_type_is_refused(self):
        ok, why = explore_pilot.check_delegate("nonsense")
        assert not ok


class TestOneAtATime:
    def test_first_call_allowed(self):
        with explore_pilot.pilot_slot() as slot:
            ok, why = explore_pilot.check_delegate("explore", slot=slot)
            assert ok, why

    def test_second_concurrent_call_refused(self):
        with explore_pilot.pilot_slot() as slot:
            explore_pilot.claim_slot(slot)
            ok, why = explore_pilot.check_delegate("explore", slot=slot)
            assert not ok
            assert "running" in why

    def test_slot_released_after_run(self):
        with explore_pilot.pilot_slot() as slot:
            explore_pilot.claim_slot(slot)
        with explore_pilot.pilot_slot() as slot:
            ok, why = explore_pilot.check_delegate("explore", slot=slot)
            assert ok, why


class TestNoNesting:
    def test_nested_call_refused(self):
        ok, why = explore_pilot.check_delegate("explore", depth=1)
        assert not ok
        assert "nest" in why

    def test_top_level_allowed(self):
        ok, why = explore_pilot.check_delegate("explore", depth=0)
        assert ok, why


class TestRefusalsNeverCallTheRunner:
    def test_refused_type_never_reaches_run_subagent(self, monkeypatch):
        called = []
        monkeypatch.setattr(explore_pilot, "_real_run_subagent",
                            lambda **kw: called.append(kw) or "RAN")
        result = explore_pilot.run_pilot_subagent(
            parent_client=object(), parent_perms=object(),
            **_kwargs(subagent_type="plan"))
        assert not called
        assert "explore" in result.error

    def test_allowed_call_reaches_run_subagent(self, monkeypatch):
        monkeypatch.setattr(explore_pilot, "_real_run_subagent",
                            lambda **kw: "RAN")
        result = explore_pilot.run_pilot_subagent(
            parent_client=object(), parent_perms=object(), **_kwargs())
        assert result == "RAN"

    def test_refused_because_busy_never_reaches_run_subagent(self, monkeypatch):
        with explore_pilot.pilot_slot() as slot:
            explore_pilot.claim_slot(slot)
            called = []
            monkeypatch.setattr(explore_pilot, "_real_run_subagent",
                                lambda **kw: called.append(kw) or "RAN")
            result = explore_pilot.run_pilot_subagent(
                parent_client=object(), parent_perms=object(), **_kwargs(),
                _pilot_slot=slot)
            assert not called
            assert "running" in result.error
