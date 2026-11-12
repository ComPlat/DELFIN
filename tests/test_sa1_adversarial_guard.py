"""SA1 adversarial review (nacht-s17): the explore-pilot guard.

Attacks the guard contract from the package brief, through the PUBLIC
entry the harness will call (``run_pilot_subagent``), never only the
inner helpers the builder's own test pins:

- only ``explore`` exactly -- spelling tricks are refused
- at most one running per session -- the slot must survive an
  exception in the delegate (no leak) and be re-claimable afterwards
- no nesting -- detected through the REAL run path: a delegate whose
  runner itself tries to delegate is refused
- every refusal names the route (non-empty reason)
"""
import pytest

from delfin.agent import explore_pilot


def _kw(**over):
    kw = dict(subagent_type="explore", description="locate X",
              prompt="find where DELFIN decides Y, file:line",
              model="", inherit_context=True)
    kw.update(over)
    return kw


class TestExactTypeOnly:
    @pytest.mark.parametrize("bad", ["Explore", "EXPLORE", " explore",
                                     "explore ", "explore\t", "ｅxplore",
                                     "explore-pilot", "ex;plore", ""])
    def test_lookalike_types_refused(self, bad):
        ok, why = explore_pilot.check_delegate(bad)
        assert not ok, f"lookalike {bad!r} slipped through"
        assert why, "refusal must name the reason"

    def test_exact_explore_allowed(self):
        ok, why = explore_pilot.check_delegate("explore")
        assert ok, why


class TestSlotSurvivesDelegateFailure:
    def test_exception_in_delegate_releases_slot(self, monkeypatch):
        def boom(**kw):
            raise RuntimeError("delegate crashed mid-run")
        monkeypatch.setattr(explore_pilot, "_real_run_subagent", boom)
        with explore_pilot.pilot_slot() as slot:
            result = explore_pilot.run_pilot_subagent(
                parent_client=object(), parent_perms=object(),
                _pilot_slot=slot, **_kw())
            assert result is not None
            # The slot is NOT stuck busy: a second delegate may run.
            with explore_pilot.pilot_slot() as slot2:
                monkeypatch.setattr(
                    explore_pilot, "_real_run_subagent", lambda **kw: "RAN")
                result2 = explore_pilot.run_pilot_subagent(
                    parent_client=object(), parent_perms=object(),
                    _pilot_slot=slot2, **_kw())
                assert result2 == "RAN", (
                    "a crashed delegate must not keep the slot busy")


class TestNestingThroughTheRealPath:
    def test_delegate_that_delegates_is_refused(self, monkeypatch):
        """The guard must catch nesting on the RUN path, not only via a
        depth kwarg: the fake delegate's runner calls the pilot entry
        again, and THAT inner call must be refused."""
        inner = {}
        monkeypatch.setattr(
            explore_pilot, "_real_run_subagent",
            lambda **kw: explore_pilot.run_pilot_subagent(
                parent_client=object(), parent_perms=object(),
                _pilot_slot=kw.get("_pilot_slot"), **_kw()) or "INNER-DONE")
        with explore_pilot.pilot_slot() as slot:
            result = explore_pilot.run_pilot_subagent(
                parent_client=object(), parent_perms=object(),
                _pilot_slot=slot, **_kw(prompt="outer, will try to nest"))
            inner["result"] = result
        # Either the outer call reports the nested refusal, or the nested
        # call itself returned an error result -- but NEVER "INNER-DONE"
        # reaching a second real delegate.
        text = str(inner["result"])
        assert "INNER-DONE" not in text, (
            "a nested delegate ran: nesting was not refused on the run path")


class TestRefusalAlwaysNamesTheRoute:
    @pytest.mark.parametrize("mode", ["type", "busy", "nest"])
    def test_error_nonempty(self, monkeypatch, mode):
        monkeypatch.setattr(explore_pilot, "_real_run_subagent",
                            lambda **kw: "RAN")
        if mode == "type":
            result = explore_pilot.run_pilot_subagent(
                parent_client=object(), parent_perms=object(),
                **_kw(subagent_type="plan"))
        else:
            with explore_pilot.pilot_slot() as slot:
                explore_pilot.claim_slot(slot)
                result = explore_pilot.run_pilot_subagent(
                    parent_client=object(), parent_perms=object(),
                    _pilot_slot=slot, **_kw())
                if mode == "nest":
                    result = explore_pilot.run_pilot_subagent(
                        parent_client=object(), parent_perms=object(),
                        _pilot_slot=slot, **_kw(), _pilot_depth=1)
        assert getattr(result, "error", "") or "RAN" != str(result), (
            f"refusal mode {mode} produced no named reason")
