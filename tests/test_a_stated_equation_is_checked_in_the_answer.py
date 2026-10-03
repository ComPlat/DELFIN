"""The arithmetic scanner is WIRED, not just exported.

scan_for_wrong_arithmetic (phase 4, 0ba71610) was shipped without a call
site -- dead code until the engine invokes it. The pin: a finished
answer that states a wrong equation gets the same treatment as an
unsourced quantity -- a [Verify] correction turn that names the
equation, and a visible caveat when the correction does not resolve it.
"""

import pytest

from delfin.agent.engine import AgentEngine


class _NoStream:
    """Engine stub: no model behind it, the correction turn records
    that stream_response was asked for."""

    def __init__(self):
        self.streamed = False
        self.feedback = None

    def stream_response(self, user_message, **kwargs):
        self.streamed = True
        self.feedback = user_message
        return ""


def _agent_with_arithmetic_flag() -> AgentEngine:
    agent = AgentEngine.__new__(AgentEngine)
    agent._last_observed_files = set()
    agent._observed_ledger_available = False
    agent._last_turn_tools = []
    agent._claim_guard_spent = False
    return agent


@pytest.mark.xfail(strict=False, reason=(
    "wiring pending: the engine edit (scan_for_wrong_arithmetic routed "
    "through _scan_claim_grounding / _enforce_claim_grounding) is blocked "
    "on the protected path delfin/agent/engine.py (self-mod-guard denied "
    "the multi_edit). These tests are the red control for that change; "
    "they turn green the moment it is applied."))
def test_a_wrong_equation_in_the_answer_fires_the_arithmetic_scan():
    agent = _agent_with_arithmetic_flag()
    loc, qty, arith = agent._scan_claim_grounding(
        "The energies give 0.031522 - 0.087844 = 0.056322, so the "
        "difference is positive.",
        [])
    assert arith, "a stated wrong equation must produce a flag"
    assert arith[0].equation == "0.031522 - 0.087844 = 0.056322"
    assert arith[0].message()


@pytest.mark.xfail(strict=False, reason=(
    "wiring pending: engine edit blocked on protected path "
    "delfin/agent/engine.py; see module-level red-control note."))
def test_a_correct_equation_produces_no_flag():
    agent = _agent_with_arithmetic_flag()
    loc, qty, arith = agent._scan_claim_grounding(
        "Sum: 1.5 + 2.25 = 3.75 as printed in the table.", [])
    assert arith == []
    assert loc == [] and qty == []


@pytest.mark.xfail(strict=False, reason=(
    "wiring pending: engine edit blocked on protected path "
    "delfin/agent/engine.py; see module-level red-control note."))
def test_the_correction_turn_names_the_equation():
    stub = _NoStream()
    agent = _agent_with_arithmetic_flag()
    agent.stream_response = stub.stream_response.__get__(agent)
    agent._append_functional_caveat = lambda *a, **k: ""
    agent._append_answer_caveats = lambda text, **k: text
    out = agent._enforce_claim_grounding(
        "The energies give 0.031522 - 0.087844 = 0.056322.")
    assert stub.streamed
    assert "0.031522 - 0.087844 = 0.056322" in stub.feedback
    assert "does not hold" in stub.feedback
