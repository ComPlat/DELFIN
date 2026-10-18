"""s25 adversarial tests for package B1 (task presence).

Reviewer layer on top of the builder's controls
(tests/test_b1_task_presence.py). Attacking what the builder's suite
does not:
A. HarnessText wake prompt reaching run_turn — the mark must survive the
   hook (the sent user message must stay harness-marked, and the wake
   text must never count as a task).
B. A person's message arriving mid-session (after a missing-history
   turn) flips the status to present; the notice stops.
C. First turn of a real session: EMPTY history + a plain typed prompt
   must be present, not missing (a real first task is not scolded).
D. None history and None workspace must not raise (assess is called on
   every turn; robustness of the pure module).
E. The hook must not change what the model receives on present (the
   prompt passes through byte-identical).

Wave-15b rule: real call sites (run_turn / prepare_fresh_turn) with the
fake engine/sink; the engine records what it was actually sent.
"""

from __future__ import annotations

import types

import pytest

from delfin.agent import repl as repl_mod
from delfin.agent.harness_rows import HARNESS_KEY, HarnessText

pytestmark = pytest.mark.timeout(60)


class _FakeEngine:
    def __init__(self, history=None):
        self.messages = list(history or [])
        self.session_id = "test-session"
        self.kit_permissions = types.SimpleNamespace(workspace=None)
        self.start_fresh = None
        self.last_turn_stop_reason = ""
        self._refusal_entries = lambda: None
        self._last_input_tokens = 0

    def stream_response(self, user_message=None, **kwargs):
        self.sent = user_message
        return "ack"

    def get_status(self):
        return {"input_tokens": 0, "output_tokens": 0, "cost_usd": 0.0}

    def estimate_context_tokens(self):
        return 0


class _Sink:
    def __init__(self):
        self.items = []

    def __call__(self, item):
        self.items.append(item)


def test_a_wake_prompt_stays_marked_after_the_hook():
    """Fused HarnessText wake text keeps the harness mark after the hook.

    If the hook re-wraps the prompt as a plain str, the wake text would
    look like a person message to every LATER turn — the exact bug the
    wave-17 finding describes.
    """
    history = [{"role": "user", "content": "wake", HARNESS_KEY: True}]
    engine = _FakeEngine(history)
    out = repl_mod.run_turn(engine, HarnessText("[operator] wake"), sink=_Sink())
    assert "no task from the person" in str(engine.sent)
    # The mark must survive whatever the hook did to the prompt.
    assert isinstance(engine.sent, HarnessText), type(engine.sent).__name__


def test_b_person_speaks_mid_session_notice_stops():
    """After a missing turn, a typed person message flips to present."""
    from delfin.agent import task_presence
    harness_only = [{"role": "user", "content": "[wake]", HARNESS_KEY: True}]
    status, _ = task_presence.assess(harness_only, workspace=None)
    assert status == "missing"
    with_person = harness_only + [{"role": "user",
                                   "content": "now the real task: fix X"}]
    status, _ = task_presence.assess(with_person, workspace=None)
    assert status == "present"


def test_c_first_turn_empty_history_plain_prompt_is_present():
    """A real first task must not be told it has no task."""
    from delfin.agent import task_presence
    status, reason = task_presence.assess([], workspace=None)
    # Empty history alone is missing...
    assert status == "missing", reason
    # ...but assess with a plain typed prompt counts it (prompt kwarg).
    status, _ = task_presence.assess([], workspace=None, prompt="Fix bug 1234")
    assert status == "present"


def test_d_assess_survives_none_history_and_none_workspace():
    from delfin.agent import task_presence
    status, reason = task_presence.assess(None, workspace=None)
    assert status == "missing"
    assert isinstance(reason, str)


def test_e_present_passes_prompt_through_unchanged():
    """On present the hook must not touch the prompt at all."""
    history = [{"role": "user", "content": "the task: run the suite"}]
    engine = _FakeEngine(history)
    sink = _Sink()
    out = repl_mod.prepare_fresh_turn(engine, "continue phase 2")
    assert out == "continue phase 2", out
    repl_mod.run_turn(engine, "continue phase 2", sink=sink)
    assert engine.sent == "continue phase 2", engine.sent
