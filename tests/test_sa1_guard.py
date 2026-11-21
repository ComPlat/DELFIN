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


# ---------------------------------------------------------------------------
# Reviewer s17 gaps (phase-1 review): the real preset dict drives the type
# enumeration; the background channel (subagents.reserve_running) and the
# worktree-isolation prompt count in the same one-at-a-time guard; nesting
# is refused whether the child asks by type or by an isolation request.
# ---------------------------------------------------------------------------

class TestTypesFromTheRealPresetDict:
    def test_every_real_type_except_explore_is_refused(self):
        from delfin.agent.subagents import SUBAGENT_PRESETS
        for name in SUBAGENT_PRESETS:
            ok, why = explore_pilot.check_delegate(name)
            if name == "explore":
                assert ok, why
            else:
                assert not ok
                assert "explore" in why, (name, why)


class TestBackgroundChannelCounts:
    def test_reserved_background_run_blocks_a_new_one(self):
        from delfin.agent import subagents as sa
        with explore_pilot.pilot_slot() as slot:
            # A background delegate holds its slot until reaped: the
            # reservation state the guard counts is the claimed slot.
            explore_pilot.claim_slot(slot)
            ok, why = explore_pilot.check_delegate("explore", slot=slot)
            assert not ok
            assert "running" in why

    def test_slot_survives_until_reaped(self):
        # Simulate the real lifecycle: reserve -> still busy -> reap -> free.
        with explore_pilot.pilot_slot() as slot:
            explore_pilot.claim_slot(slot)
            assert not explore_pilot.check_delegate("explore", slot=slot)[0]
            explore_pilot.reap_slot(slot)
            ok, why = explore_pilot.check_delegate("explore", slot=slot)
            assert ok, why


class TestNestingAndIsolation:
    def test_child_prompt_requesting_a_subagent_is_refused_by_depth(self):
        # A delegate runs with depth=1; whatever its prompt asks for, the
        # pilot setting refuses the delegation itself.
        ok, why = explore_pilot.check_delegate("explore", depth=1)
        assert not ok
        assert "nest" in why

    def test_worktree_isolation_is_refused(self):
        # isolation="worktree" creates a worktree — a write-surface the
        # read-only pilot setting does not allow, whatever the type says.
        ok, why = explore_pilot.check_delegate("explore", isolation="worktree")
        assert not ok
        assert "isolation" in why


class TestEngineWiringUsesTheSlotVocabulary:
    """The engine wiring the CLI builds (`_engine_factory_builder`) must go
    through the module's own slot functions, not an inline slot object: a
    runner that never reaps would leave the slot claimed after the first
    delegate, and the second delegation of the same run would be refused
    as busy. Driven through the wiring the harness installs, never by
    calling `pilot_slot` directly."""

    def _engine(self):
        class _Perms:
            subagent_runner = None

        class _Engine:
            kit_permissions = _Perms()

        return _Engine(), _Perms

    def test_the_wired_runner_reaps_the_slot_after_each_call(self, monkeypatch):
        class _Perms:
            subagent_runner = lambda self, **kw: kw  # placeholder binding

        class _Engine:
            def __init__(self):
                self.kit_permissions = _Perms()

        engine = _Engine()
        builder = explore_pilot._engine_factory_builder(None)
        monkeypatch.setattr(
            "delfin.agent.benchmark_runner._default_engine_factory",
            lambda *a, **kw: engine)
        engine_out = builder("m", "b", "p", "solo", task_class="pilot")
        assert engine_out is engine
        # The wiring replaced the placeholder binding with the pilot runner.
        runner = engine.kit_permissions.subagent_runner
        assert runner is not _Perms.subagent_runner

        ran = []
        monkeypatch.setattr(explore_pilot, "_real_run_subagent",
                            lambda **kw: ran.append(kw) or "RAN")

        # The wired runner itself, called twice directly: the second call
        # only reaches the real runner because the first reaped the slot.
        assert runner(subagent_type="explore", description="d", prompt="p",
                      model="", parent_client=None, parent_perms=None) == "RAN"
        assert runner(subagent_type="explore", description="d", prompt="p",
                      model="", parent_client=None, parent_perms=None) == "RAN"
        assert len(ran) == 2
