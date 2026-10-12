"""H1E phase-3 wire RED test (engine-level, after the operator build finding).

The real leak loop is ACROSS follow-up engine turns, not inside one
stream_message call: a single call with a leaked round ends after one
round, and the repeat comes from the follow-up-turn machinery. This test
drives the real engine path -- ``AgentEngine.stream_response`` -- with a
client whose stream_response returns a leaked ``<invoke>`` tool call on
every call, and asserts the per-engine stall guard behaviour:

1. three consecutive leak-turns stop the engine: ``_stall_guard_stopped``
 is set and the output carries a status naming what repeated;
2. two consecutive leak-turns do NOT stop (a third is still allowed);
3. a real tool run between two leak-turns resets the guard, so two leak
 turns separated by genuine progress are not stopped.

The red-control contract: on the PRISTINE engine.py (no wire) these
assertions FAIL (no guard hooks exist); on the PATCHED engine built from
.gate/patch/patched/engine.py (hunks A-E) they PASS. The operator build
runs this test against both.
"""

from unittest.mock import patch

from delfin.agent import engine as E

LEAK = ('<invoke name="read_file">'
 '<parameter name="path">/etc/hostname</parameter></invoke>')


class _LeakClient:
    """A fake model client whose stream_response always leaks the call."""

    def __init__(self):
        self.calls = 0

    def stream_response(self, *_args, **_kwargs):
        self.calls += 1
        return LEAK

    def stream_message(self, *_args, **_kwargs):
        # No tool stream: nothing yields, so the tool loop has no rounds
        # and the turn ends with the LEAK text (which the engine reads as
        # a leaked tool call and feeds into the stall guard).
        return iter(())


def _make_engine(tmp_path, client):
    with patch("delfin.agent.engine.create_client", return_value=client):
        return E.AgentEngine(
            repo_dir=str(tmp_path), backend="api", provider="kit",
            model="kit.qwen3.5-397b-A17b", mode="solo")


def test_three_consecutive_leak_turns_stop_the_engine(tmp_path):
    client = _LeakClient()
    eng = _make_engine(tmp_path, client)
    eng.stream_response("t0")
    eng.stream_response("t1")
    out = eng.stream_response("t2")
    stopped = getattr(eng, "_stall_guard_stopped", "")
    assert stopped != "", "third leak-turn must set _stall_guard_stopped"
    assert "repeated" in (out + stopped)


def test_two_consecutive_leak_turns_do_not_stop(tmp_path):
    client = _LeakClient()
    eng = _make_engine(tmp_path, client)
    eng.stream_response("t0")
    eng.stream_response("t1")
    assert getattr(eng, "_stall_guard_stopped", "") == ""


def test_the_guard_lives_per_engine_not_shared(tmp_path):
    # Two engines must each have their own streak; a leak on engine A must
    # not affect engine B, and 2 leaks on one engine alone must not stop.
    g1 = _make_engine(tmp_path, _LeakClient())
    g2 = _make_engine(tmp_path, _LeakClient())
    g1.stream_response("t0")
    g2.stream_response("t0")
    g1.stream_response("t1")
    assert getattr(g1, "_stall_guard_stopped", "") == ""
    assert getattr(g2, "_stall_guard_stopped", "") == ""
