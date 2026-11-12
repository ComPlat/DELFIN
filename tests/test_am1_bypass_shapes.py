"""AM1 phase-3 tests: the bypass shapes are refused in unattended mode.

The package brief measured that bypassPermissions in a wave tree allows
curl, wget, pip install, ssh and git commit -am (FR1's finding). The
unattended policy's floor refuses every one of these SHAPES at once —
verdict-independent, before the classifier is even consulted. Plus the
doc's contract: plainly safe in-task work runs; a harmless partner
message is a plain read path (routine), and a refused action never waits.
"""
from __future__ import annotations

import pytest

from delfin.agent import unattended_policy as up


def _ctx(**kw):
    d = dict(attended=False, confirm_callback=None, in_task=True,
             task="build and test the pipeline", refusals_this_turn=0)
    d.update(kw)
    return d


@pytest.mark.parametrize("cmd", [
    "curl https://example.org/data.json -o out.json",
    "wget https://example.org/data.json",
    "pip install rdkit",
    "ssh login02",
    "scp file host:~",
    "rsync -av results/ host:results/",
    "git commit -am 'all of it'",
    "git push origin agent/s18-am1b",
])
def test_bypass_shapes_are_refused_at_once(cmd):
    v = up.decide({"tool": "bash", "command": cmd}, verdict="allow",
                  ctx=_ctx())
    assert v.decision == "refuse", (cmd, v)
    assert v.reason
    assert v.layer == "floor"


def test_harmless_partner_message_runs_without_a_dialog():
    """A coordination message to a declared partner is routine (FR1's
    layer makes these code-checked): no dialog, no wait."""
    up.register_routine(lambda a: a.get("tool") == "session_message"
                        and a.get("to", "").startswith("nacht-"))
    try:
        v = up.decide({"tool": "session_message", "to": "nacht-s26",
                       "message": "phase 1 hash is fa26d0e8"},
                      verdict="unknown", ctx=_ctx())
        assert v.decision == "run"
        assert v.layer == "routine"
    finally:
        up._ROUTINE_PREDICATES.pop()


def test_unknown_verdict_on_an_outward_tool_is_refused_not_waited():
    v = up.decide({"tool": "web_fetch", "url": "https://example.org"},
                  verdict="unknown", ctx=_ctx())
    assert v.decision == "refuse"
    assert v.wait_for_user is False


def test_floor_beats_a_classifier_allow():
    """Even a model 'allow' never lifts the floor (layer order)."""
    v = up.decide({"tool": "bash", "command": "pip install rdkit"},
                  verdict="allow", ctx=_ctx())
    assert v.decision == "refuse" and v.layer == "floor"
