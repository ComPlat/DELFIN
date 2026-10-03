"""Package B / phase-4 wiring control, operator ruling 2026-10-04.

Drives AgentEngine.stream_response (the public path that funnels
dashboard/solo/benchmark alike) and asserts the arithmetic scanner is
ROUTED, not merely implemented:

- a wrong stated equation ("A - B = C" that does not hold) must enter
  the same forced-correction path as an ungrounded code location --
  the feedback text is the scanner's own arithmetic_feedback wording,
  which no other scanner emits;
- a correct stated equation must produce no correction at all.

RED without .gate/engine_wiring.patch (the wiring commit the operator
applies on integrate/welle-11): _scan_claim_grounding still returns a
2-tuple, so the unpack in _enforce_claim_grounding raises TypeError and
the engine path fails. GREEN after it. The engine.py edit itself is the
operator's -- this file must not depend on local edits to the protected
path."""

from unittest.mock import MagicMock, patch

import pytest

from delfin.agent.task_evidence import check_completion_claim


def test_grounding_caveat_takes_arithmetic_flags_as_third_arg():
    """Operator part (1): the caveat is backward compatible -- old
    two-argument callers keep working, the third argument defaults to
    None and adds the equation to the named items when given."""
    from delfin.agent import verify_guard as vg

    loc = [vg.LocationClaimFlag(claim="Zeile 26: class AgentEngine",
                                path="", line=26, kind="bare_line")]
    qty = [vg.QuantityClaimFlag(quantity="2.31 eV", unit="eV")]
    arith = [vg.ArithmeticFlag(equation="1 + 1 = 3", claimed="3",
                               actual="2")]

    two_arg = vg.grounding_caveat(loc, qty)
    three_arg = vg.grounding_caveat(loc, qty, arith)
    assert two_arg                      # old signature intact
    assert "1 + 1 = 3" in three_arg     # equation joins the named items
    assert "1 + 1 = 3" not in two_arg
    assert vg.grounding_caveat(loc, qty, None) == two_arg


def test_arithmetic_feedback_names_the_equation():
    from delfin.agent import verify_guard as vg

    flag = vg.ArithmeticFlag(equation="0.031522 - 0.087844 = 0.056322",
                             claimed="0.056322", actual="-0.056322")
    text = vg.arithmetic_feedback([flag])
    assert "0.031522 - 0.087844 = 0.056322" in text
    assert "recompute" in text          # the fix it demands is recomputing


# --- public-path wiring ----------------------------------------------------

def _engine(agent_tree, client):
    from delfin.agent.engine import AgentEngine
    with patch("delfin.agent.engine.create_client", return_value=client):
        return AgentEngine(repo_dir=agent_tree, backend="cli",
                           mode="quick", pack_dir=agent_tree)


def _claims_client(replies):
    """Fake backend client: each stream_message call yields the next
    reply as text events (same shape as test_grounding_and_budget)."""
    from delfin.agent.api_client import StreamEvent
    fake = MagicMock()
    fake._observed_files_session = set()
    calls = {"n": 0}

    def _stream(*a, **k):
        i = calls["n"]
        calls["n"] += 1
        text = replies[min(i, len(replies) - 1)]
        yield StreamEvent(type="text_delta", text=text)
        yield StreamEvent(type="message_delta", output_tokens=5,
                          cost_usd=0.0)

    fake.stream_message = MagicMock(side_effect=_stream)
    return fake


@pytest.fixture
def agent_tree(tmp_path):
    lite_dir = tmp_path / "pack_lite"
    modes_dir = lite_dir / "modes"
    modes_dir.mkdir(parents=True)
    (modes_dir / "solo.md").write_text("# solo mode")
    manifest = """\
pack_name: DELFIN_AGENT_LITE
version: 1
modes:
  - id: solo
    file: modes/solo.md
    route:
      - session_manager
"""
    (lite_dir / "manifest.yaml").write_text(manifest)
    return tmp_path


def test_a_wrong_equation_in_the_answer_fires_the_correction_turn(
        agent_tree):
    """RED without the operator's wiring patch: the unpack of the third
    scan result raises TypeError inside _enforce_claim_grounding before
    any correction can run. GREEN after: the engine path routes the
    arithmetic flag through the SAME forced-correction loop as an
    ungrounded code location, and the feedback is the scanner's own
    arithmetic_feedback wording."""
    fake = _claims_client([
        "The energies give 0.031522 - 0.087844 = 0.056322, so the "
        "difference is positive.",
        "Recomputed from the cited outputs: 0.031522 - 0.087844 = "
        "-0.056322. The difference is negative, not positive."])
    engine = _engine(agent_tree, fake)
    out = engine.stream_response("what is the energy difference?")

    assert fake.stream_message.call_count == 2
    feedback = [m for m in engine.messages
                if m.get("role") == "user"
                and "[Verify]" in str(m.get("content", ""))]
    assert feedback, "a wrong stated equation must force a correction"
    assert "0.031522 - 0.087844 = 0.056322" in str(feedback[0]["content"])
    # arithmetic_feedback's own wording, no other scanner emits it
    assert "do not hold" in str(feedback[0]["content"])


def test_a_correct_equation_in_the_answer_runs_no_correction(agent_tree):
    fake = _claims_client([
        "Sum: 1.5 + 2.25 = 3.75 as printed in the table."])
    engine = _engine(agent_tree, fake)
    out = engine.stream_response("what is the sum?")

    assert fake.stream_message.call_count == 1
    assert out == "Sum: 1.5 + 2.25 = 3.75 as printed in the table."
