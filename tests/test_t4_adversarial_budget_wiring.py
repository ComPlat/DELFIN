"""Package T4, reviewer (nacht-s19) — RED control for the failure-budget DISPATCH.

The pure FailureBudget class is in-tree and green (see
test_t4_adversarial_budget.py). But phase 3 is not done until the budget is
WIRED into the tool dispatch: d9bf770b defines the class, nothing in
api_client.py calls it (verified: grep for FailureBudget/repeat_hint in
api_client.py is empty). This test is the reviewer's independent red control
for that wiring, driven through the PUBLIC entry point
(DocToolExecutor.execute -> _execute_pinned -> failure hook ~api_client.py:11529)
so it proves the real dispatch behaviour, not a too-weak probe against the pure
class.

Contract under test (phase 3 spec): the same failing (tool, normalized args)
call a third time in one task must yield a "stop and ask" hint instead of yet
another bare error. The wired hint carries the literal message "... STOP and ask
the user for help." (GroundingHint in action_grounding.py), so the assertion
pins on that stable text rather than on a specific JSON field name the wiring
might choose.

Determinism: a tool with a missing REQUIRED argument (get_calc_info without
calc_id) fails structurally at _missing_required_argument, BEFORE any dispatch,
machine- and filesystem-independently. Three such identical calls are exactly
the stuck loop phase 3 targets.

RED today: the wiring is absent, so all calls come back as bare errors and the
"STOP and ask" text never appears -> the 4th-call assertion fails (the control).
The success-reset hunk of the wiring is already pinned at the PURE level
(test_t4_adversarial_budget.py test_failure_signature...); it cannot be driven
to a guaranteed SUCCESS through the public path without a filesystem-mocked
tool, so it is re-checked against the landed patch instead.

Run via the gate:  gate tests/test_t4_adversarial_budget_wiring.py -q | tail -12
"""
from __future__ import annotations

import json

from delfin.agent import api_client as A


def _run(architecture: object, name: str, arguments: dict, workspace: object):
    perms = A.KitToolPermissions(workspace=workspace, mode="acceptEdits")
    return architecture.execute(name, arguments, perms)


def test_identical_failing_call_third_time_yields_stop_and_ask(tmp_path):
    """The dispatch must stop a stuck loop: 4 identical missing-arg failures
    on get_calc_info (requires calc_id), and the 4th returns the grounded
    'stop and ask' instruction instead of a fourth bare error."""
    ex = A._DocToolExecutor()
    name, arguments = "get_calc_info", {}

    # First three identical failing calls: bare errors, NOT yet a stop-and-ask.
    for _ in range(3):
        out = json.loads(_run(ex, name, arguments, tmp_path))
        assert out.get("error"), out
        assert "STOP and ask" not in json.dumps(out)

    # The fourth identical call must be answered with a stop-and-ask hint,
    # not a fourth bare error. RED on the unwired tree.
    fourth = _run(ex, name, arguments, tmp_path)
    assert "STOP and ask" in fourth, fourth
