"""Package T4, reviewer (nacht-s19) — control for the failure-budget DISPATCH.

The pure FailureBudget class is in-tree and green (see
test_t4_adversarial_budget.py). Phase 3 closes only when the budget is WIRED
into the tool dispatch. This test drives the REAL public entry point
(DocToolExecutor.execute -> _execute_pinned -> failure hook) with a
deterministic missing-required-argument failure on get_calc_info (no calc_id),
so it proves the dispatch behaviour, not a too-weak probe against the pure
class.

Contract under test (package-T4 spec verbatim): "the same failing call (tool +
normalized args) a THIRD time in one task returns a hint to stop and ask."
So the 3rd identical failing call carries the stop-and-ask hint; the first two
are bare errors. The wired hint carries the literal message "... STOP and ask
the user for help." (GroundingHint in action_grounding.py); we pin on that
stable text.

NOTE ON THE REWRITE (was a red control for the 4th call): the pre-wiring tree
returned no hint at all, so this test as originally written expected the hint
on the 4th call (first three bare). With the operator's correct wiring the
budget feeds BEFORE the repeat-check, so the THIRD call (streak 3 >=
identical_fail_limit 3) already fires the hint. The test is rewritten along
the true spec (third call fires), which the original had misread as the fourth.

Determinism: a tool with a missing REQUIRED argument (get_calc_info without
calc_id) fails structurally at _missing_required_argument, BEFORE any dispatch,
machine- and filesystem-independently.

Run via the gate:  gate tests/test_t4_adversarial_budget_wiring.py -q | tail -12
"""
from __future__ import annotations

import pytest

import json

from delfin.agent import api_client as A


def _run(architecture: object, name: str, arguments: dict, workspace: object):
    perms = A.KitToolPermissions(workspace=workspace, mode="acceptEdits")
    return architecture.execute(name, arguments, perms)


@pytest.fixture(autouse=True)
def _fresh_budget(monkeypatch):
    """The failure budget is process-wide; a test that counts calls must
    start from an empty one, whatever ran before it."""
    monkeypatch.setattr(A._DocToolExecutor, "_T4_FAIL_BUDGET", None,
                        raising=False)


def test_identical_failing_call_third_time_yields_stop_and_ask(tmp_path):
    """The dispatch must stop a stuck loop: 3 identical missing-arg failures
    on get_calc_info (requires calc_id), and the THIRD returns the grounded
    'stop and ask' instruction (spec: third time) instead of a bare error."""
    ex = A._DocToolExecutor()
    name, arguments = "get_calc_info", {}

    # First two identical failing calls: bare errors, NOT yet a stop-and-ask.
    for _ in range(2):
        out = json.loads(_run(ex, name, arguments, tmp_path))
        assert out.get("error"), out
        assert "STOP and ask" not in json.dumps(out)

    # The THIRD identical call is answered with a stop-and-ask hint.
    third = _run(ex, name, arguments, tmp_path)
    assert "STOP and ask" in third, third
