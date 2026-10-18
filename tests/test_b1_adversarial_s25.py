"""Reviewer (s25) adversarial cases for B1, folded into the builder's file.

Source: git show 28fba2c4:tests/test_b1_adversarial_s25.py (the reviewer's
branch; the cases, not the branch, are merged here). 4 were red on the
pre-fix tree; they pin:

A. a HarnessText wake stays HarnessText after the hook (the mark must
   survive the fusion — a re-wrapped plain str would look like the person
   to every later turn, the wave-17 bug again);
B. a person speaking mid-session flips missing -> present and the notice
   stops;
C. empty history + plain typed prompt = present (a real first task is
   not scolded); empty history alone = missing;
D. assess(None history) must not raise;
E. on present the prompt passes through byte-identical.
"""

from __future__ import annotations

import pytest

from delfin.agent import repl as repl_mod
from delfin.agent.harness_rows import HARNESS_KEY, HarnessText

pytestmark = pytest.mark.timeout(60)


def _engine(history):
    from tests.test_b1_task_presence import _FakeEngine
    return _FakeEngine(history)


def _sink():
    from tests.test_b1_task_presence import _Sink
    return _Sink()


def test_a_wake_prompt_stays_marked_after_the_hook():
    history = [{"role": "user", "content": "wake", HARNESS_KEY: True}]
    engine = _engine(history)
    repl_mod.run_turn(engine, HarnessText("[operator] wake"), sink=_sink())
    assert "no task from the person" in str(engine.sent)
    assert isinstance(engine.sent, HarnessText), type(engine.sent).__name__


def test_b_person_speaks_mid_session_notice_stops():
    from delfin.agent import task_presence
    harness_only = [{"role": "user", "content": "[wake]", HARNESS_KEY: True}]
    status, _ = task_presence.assess(harness_only, workspace=None)
    assert status == "missing"
    with_person = harness_only + [{"role": "user",
                                   "content": "now the real task: fix X"}]
    status, _ = task_presence.assess(with_person, workspace=None)
    assert status == "present"


def test_c_first_turn_empty_history_plain_prompt_is_present():
    from delfin.agent import task_presence
    status, reason = task_presence.assess([], workspace=None)
    assert status == "missing", reason
    status, _ = task_presence.assess([], workspace=None, prompt="Fix bug 1234")
    assert status == "present"


def test_d_assess_survives_none_history_and_none_workspace():
    from delfin.agent import task_presence
    status, reason = task_presence.assess(None, workspace=None)
    assert status == "missing"
    assert isinstance(reason, str)


def test_e_present_passes_prompt_through_unchanged():
    history = [{"role": "user", "content": "the task: run the suite"}]
    engine = _engine(history)
    sink = _sink()
    out = repl_mod.prepare_fresh_turn(engine, "continue phase 2")
    assert out == "continue phase 2", out
    repl_mod.run_turn(engine, "continue phase 2", sink=sink)
    assert engine.sent == "continue phase 2", engine.sent
