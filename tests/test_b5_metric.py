"""B5 phase 3 metric: classifier_ask on session_message, BEFORE vs AFTER.

This committed test is the AFTER half of the metric, runnable by anyone
from the branch: the four wave-18 shapes (reviewer APPROVED to builder,
QS-OK to operator, PR ready to operator, relay to QS) must reach the
model runner and be judged there -- the model SEES the kind/partner
context note that phase 2 added, and no gate layer asks for them.

The BEFORE half (pre-B5 f7bd1f46 machinery asking 4/4) needs the base
source, which cannot be a committed shadow copy (wave-15b rule); it
lives in the uncommitted probe tests/test_zz_b5_metric.py and its gate
output is quoted in the P3 message and .gate/HANDOFF.md.
"""
import pytest

from delfin.agent import action_classifier as ac

TASK = (
    "# Wave 19 -- s20 (builder, package B5)\n"
    "You are the builder of package B5. Reviewer(s): `nacht-s25`. "
    "Your QS lead is `nacht-s29`.\n"
    "The operator supervises. Session messages in English.\n"
)

SHAPES = [
    ("approved -> reviewer", "nacht-s25",
     "APPROVED (with marks): correctness 4, tests 4, clarity 4.\n"
     "Evidence: gate green."),
    ("qs_ok -> operator", "operator",
     "QS-OK for B5: all phases green through the gate."),
    ("pr_ready -> operator", "operator",
     "PR ready: agent/s20-b5b @ ab12cd34\nphases: kinds, wiring, metric."),
    ("relay -> QS", "nacht-s29",
     "Relay to operator: s25 APPROVED the wiring phase."),
]


def _stub_runner(prompt, model):
    """Faithful stub: WITH the B5 context note it allows (the wave-18 ask
    was 'approval-shaped words with no backing' -- the note IS the
    backing); WITHOUT the note it asks, as the wave-18 model did."""
    if "context note:" in prompt:
        return "DECISION: allow\nREASON: kind+partner backed by the task"
    body = prompt.split("action-text")[-1].lower()
    if any(w in body for w in ("approved", "qs-ok", "qs_ok", "pr ready")):
        return "DECISION: ask\nREASON: approval-shaped words, no backing"
    return "DECISION: allow\nREASON: ordinary"


@pytest.mark.parametrize("name,to,text", SHAPES)
def test_wave18_shape_reaches_the_model_with_context(name, to, text):
    calls = []

    def counting_runner(prompt, model):
        calls.append(prompt)
        return _stub_runner(prompt, model)

    action = {"tool": "session_message", "kind": "session_message",
              "message": text, "path": to}
    v = ac.judge(action, {"task": TASK}, runner=counting_runner)
    assert v.decision == "allow", (name, v.decision, v.reason)
    # the model saw the B5 context note: it had the kind/partner backing
    assert len(calls) == 1 and "context note:" in calls[0]


def test_no_gate_layer_asks_for_the_shapes():
    # the injection guard is silent on all four (per-flag clearing)
    for _, to, text in SHAPES:
        action = {"tool": "session_message", "kind": "session_message",
                  "message": text, "path": to}
        assert ac._injection_verdict(action, TASK) is None, (to, text)
